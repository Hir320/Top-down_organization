rm(list =ls())
library(tidyverse)
library(ggfortify)
library(ggpubr)
library(gt)


# Define the directory where you want to search for the files
setwd('/home/rstudio/')

# 関数 ----------------------------------------------------------------------

rawdata_to_dataframe <- function(rawdata) {
  df <- rawdata %>% 
    filter(Color == "blue") %>% 
    group_by(filename_plus10) %>%
    slice_min(filename_minus55, with_ties = FALSE) %>%
    ungroup() %>% 
    group_by(filename_minus55) %>%
    slice_min(filename_plus10, with_ties = FALSE) %>%
    ungroup()
  return(df)
}

plot_ei_balance <- function(data) {
  comparisons <- list(c("ACC", "RSC"))
  # Create the plot
  p <- ggplot(data = data, aes(x = Region, y = EIbalance, colour = Region)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
    stat_compare_means(comparisons = comparisons, 
                       method = "wilcox.test", 
                       label = "p.signif", 
                       hide.ns = TRUE) +
    theme_minimal() +
    #labs(title = "EI balance when EPSC ≈ 40 pA", y = "IPSC/EPSC", x = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_classic() +
    theme(legend.position = "None",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 20)) +
    scale_color_manual(values = c('#03af7a', '#ff4b00'))
  
  # Return the plot object
  plot(p)
}

fit_and_plot_gamma_model <- function(data) {
  # Fit the Gamma GLM
  glm_model <- glm(
    EIbalance ~ Region + RoughAP,
    data = data,
    family = Gamma(link = "log"),
    control = glm.control(maxit = 100)
  )
  
  # Create an autoplot of the GLM model
  model_autoplot <- autoplot(glm_model)
  print(model_autoplot)
  
  # Print the model summary
  model_summary <- summary(glm_model)
  print(model_summary)
  
  # Return a list with the model, its summary, and the autoplot object
  return(list(model = glm_model, summary = model_summary, autoplot = model_autoplot))
}
make_gt_table <- function(glm_large_results) {
  # Tidy the model output into a data frame
  tidy_results <- tidy(glm_large_results$model) %>%
    mutate(term = ifelse(term == "RoughAP", "AP coord", term))
  
  # Create a formatted table using gt
  gt_table <- tidy_results %>%
    gt() %>%
    tab_header(
      title = "Gamma GLM Summary",
      subtitle = "EIbalance ~ Region + AP coord"
    ) %>%
    fmt_number(
      columns = vars(estimate, std.error, statistic, p.value),
      decimals = 3
    ) %>%
    cols_label(
      term = "Parameter",
      std.error = "Std. Error",
      statistic = "t Value",
      p.value = "p Value"
    )
  
  return(gt_table)
}


# データ読み込み -----------------------------------------------------------------

rawdata <- read.csv("./Data/PatchClamp_escape//EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_GtCCR4.csv")

data <- rawdata_to_dataframe(rawdata) %>% 
  filter(!is.na(PeakAmp_minus55_mean), !is.na(PeakAmp_plus10_mean)) %>% 
  mutate(EIbalance = PeakAmp_plus10_mean / -PeakAmp_minus55_mean,
         log_EIbalance = log(EIbalance),
         uniqueCellID = paste(BrainID, SliceID, CellID, sep = "_"))


# anterior ACC vs posterior RSC -------------------------------------------

small <- data %>% 
  filter(!is.na(RoughAP),
         PeakAmp_minus55_mean < -10,
         PeakAmp_minus55_mean > -50) %>% 
  group_by(uniqueCellID) %>%
  slice(which.min(abs(PeakAmp_minus55_mean - (-40))))

large <- data %>% 
  filter(!is.na(RoughAP),
         PeakAmp_minus55_mean < -70,
         PeakAmp_minus55_mean > -120) %>%  # Filter for rows where Im_scaled_E is below -80
  group_by(uniqueCellID) %>%
  slice_min(PeakAmp_minus55_mean) %>%  # dplyr::select the row with the smallest Im_scaled_E for each Cell_ID
  ungroup()

ggplot(small, aes(x = RoughAP, y = EIbalance, color = Region)) +
  geom_point() +
  geom_smooth(method="lm")
glm_small_results <- fit_and_plot_gamma_model(small)

ggplot(large, aes(x = RoughAP, y = EIbalance, color = Region)) +
  geom_point() +
  geom_smooth(method="lm")
glm_large_results <- fit_and_plot_gamma_model(large)

plot_ei_balance(small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/GtCCR4_small_EIbalance.png", 
       width = 65, height = 50, units = "mm")
plot_ei_balance(large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/GtCCR4_large_EIbalance.png", 
       width = 65, height = 50, units = "mm")

## Gt tableにGLMまとめる --------------------------------------------------------
gt_table_small <- make_gt_table(glm_small_results)
gtsave(gt_table_small, filename = "./Data/PatchClamp_escape/EIkinetics_APaxis/plots/glm_table_GtCCR4_small.png", vwidth = 600, vheight = 200, zoom = 2)
gt_table_large <- make_gt_table(glm_large_results)
gtsave(gt_table_large, filename = "./Data/PatchClamp_escape/EIkinetics_APaxis/plots/glm_table_GtCCR4_large.png", vwidth = 600, vheight = 200, zoom = 2)


# Wilcox rank sum test ----------------------------------------------------
# Wilcoxon rank sum test for EIbalance between Regions in aACCpRSC_small
wilcox_test_result <- wilcox.test(EIbalance ~ Region, data = small)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(EIbalance ~ Region, data = large)
print(wilcox_test_result)
