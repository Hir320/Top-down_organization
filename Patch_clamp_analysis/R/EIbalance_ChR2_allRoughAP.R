rm(list =ls())
library(tidyverse)
library(ggfortify)
library(ggpubr)
library(broom)
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
  # Create the plot
  p <- ggplot(data = data, aes(x = Region, y = EIbalance, colour = Region)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
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

raw_onlyChR2 <- read.csv("./Data/PatchClamp_escape//EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_ChR2.csv")
raw_dualChR2 <- read.csv("./Data/PatchClamp_escape//EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_Dual_Injection.csv")

data_onlyChR2 <- rawdata_to_dataframe(raw_onlyChR2)
data_dualChR2 <- rawdata_to_dataframe(raw_dualChR2) %>% 
  mutate(Opsin = "ChR2",
         Region = str_sub(Region, 1, 3))

data_allChR2 <- bind_rows(data_onlyChR2, data_dualChR2) %>% 
  filter(!is.na(PeakAmp_minus55_mean), !is.na(PeakAmp_plus10_mean)) %>% 
  mutate(EIbalance = PeakAmp_plus10_mean / -PeakAmp_minus55_mean,
         log_EIbalance = log(EIbalance),
         uniqueCellID = paste(BrainID, SliceID, CellID, sep = "_"))


# ACC vs RSC in all regions -------------------------------------------

allChR2_small <- data_allChR2 %>% 
  filter(!is.na(RoughAP),
         PeakAmp_minus55_mean < -10,
         PeakAmp_minus55_mean > -50) %>% 
  group_by(uniqueCellID) %>%
  slice(which.min(abs(PeakAmp_minus55_mean - (-40))))

allChR2_large <- data_allChR2 %>% 
  filter(!is.na(RoughAP),
         PeakAmp_minus55_mean < -70,
         PeakAmp_minus55_mean > -120) %>%  # Filter for rows where Im_scaled_E is below -80
  group_by(uniqueCellID) %>%
  slice_min(PeakAmp_minus55_mean) %>%  # dplyr::select the row with the smallest Im_scaled_E for each Cell_ID
  ungroup()

ggplot(allChR2_small, aes(x = RoughAP, y = EIbalance, color = Region)) +
  geom_point() +
  geom_smooth(method="lm")
glm_small_results <- fit_and_plot_gamma_model(allChR2_small)

ggplot(allChR2_large, aes(x = RoughAP, y = EIbalance, color = Region)) +
  geom_point() +
  geom_smooth(method="lm")
glm_large_results <- fit_and_plot_gamma_model(allChR2_large)


## Gt tableにGLMまとめる --------------------------------------------------------
gt_table_small <- make_gt_table(glm_small_results)
gtsave(gt_table_small, filename = "./Data/PatchClamp_escape/EIkinetics_APaxis/plots/glm_table_small.png", vwidth = 600, vheight = 200, zoom = 2)
gt_table_large <- make_gt_table(glm_large_results)
gtsave(gt_table_large, filename = "./Data/PatchClamp_escape/EIkinetics_APaxis/plots/glm_table_large.png", vwidth = 600, vheight = 200, zoom = 2)

# GLM model fit調べる -----------------------------------------------------------------

library(DHARMa)

# 3. Run DHARMa residual diagnostics
dharma_res <- simulateResiduals(fittedModel = glm_large_results$model, n = 250) 
# n = number of simulations, adjust if needed

# 4. Plot DHARMa diagnostics
# This opens a series of diagnostic plots in the viewer:
plot(dharma_res)  

# 5. Additional DHARMa tests:
#    - Check if residuals are uniformly distributed
testUniformity(dharma_res)

#    - Check for over/under-dispersion
testDispersion(dharma_res)

# 多群間多重比較 -----------------------------------------------------------------
# 
# library(ggplot2)
# library(dplyr)
# library(rstatix)
# library(ggpubr)
# 
# # 1. Region と APregion の組み合わせを表す因子変数を作成
# allChR2_small <- allChR2_small %>%
#   mutate(Region_AP = interaction(Region, APregion, sep = " - "))
# 
# # 2. 箱ひげ図の作成
# p <- ggplot(allChR2_small, aes(x = Region_AP, y = EIbalance, fill = Region)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "EIbalanceの比較 (RegionとAPregionの組み合わせ)",
#        x = "Region - APregion", y = "EIbalance")
# 
# # 3. pairwise Wilcox検定を実施（rstatix::compare_means を利用）
# # 1. グループ情報を解除して data.frame に変換
# df <- allChR2_small %>% 
#   ungroup() %>% 
#   as.data.frame()
# 
# # 2. compare_means() を data 引数として df を渡して実行
# pairwise_results <- compare_means(EIbalance ~ Region_AP, 
#                                   data = df, 
#                                   method = "wilcox.test", 
#                                   p.adjust.method = "bonferroni")
# 
# # 4. 検定結果の表示位置(y.position)を設定（例：全データの最大値の 1.05 倍）
# pairwise_results <- pairwise_results %>%
#   mutate(y.position = max(allChR2_small$EIbalance) * 1.05)
# 
# # 5. 箱ひげ図にペアごとのp値を表示
# p + stat_pvalue_manual(
#   pairwise_results, 
#   label = "p.adj", 
#   tip.length = 0.01, 
#   hide.ns = TRUE,
#   inherit.aes = FALSE
# )
