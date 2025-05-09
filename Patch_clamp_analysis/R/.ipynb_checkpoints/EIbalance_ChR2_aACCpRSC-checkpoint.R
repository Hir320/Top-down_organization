rm(list =ls())
library(tidyverse)
library(ggfortify)
library(ggpubr)
library(rlang)


# Define the directory where you want to search for the files
setwd('/home/rstudio/')

fig_family <- "Arial"

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

plot_ei <- function(data, y) {
  comparisons <- list(c("ACC", "RSC"))
  y_quo <- enquo(y)
  
  p <- ggplot(data, aes(x = Region, y = !!y_quo, colour = Region)) +
    geom_boxplot(outliers = FALSE, width = 0.5, linewidth = 0.2) +
    geom_jitter(width = 0.1, height=0, size = 0.25, alpha = 0.7) +
    # stat_compare_means(comparisons = comparisons, 
    #                    method = "wilcox.test", 
    #                    label = "p.signif", 
    #                    amily='mono',
    #                    hide.ns = TRUE) +
    #theme_minimal() +
    #labs(title = "EI balance when EPSC ≈ 40 pA", y = "IPSC/EPSC", x = "") +
    #scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    # theme_void() +
    # theme(legend.position = "None",
    #       axis.title = element_blank(),
    #       axis.text.x = element_blank(),
    #       axis.text.y = element_text(size = 6),  # Add y-axis text back
    #       axis.ticks.length.y.left =  unit(1, "mm"),
    #       axis.ticks.y.left =  element_line(linewidth = 0.25)
    # ) +
    theme_classic() +
    theme(legend.position = "None",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8, family = fig_family),  # Add y-axis text back
          axis.line = element_line(linewidth = 1/2.13), # 全体の軸線に適用
          axis.ticks.length.x.bottom = unit(0, "mm"),
          axis.ticks.length.y.left =  unit(1, "mm"),
          # Y軸左側の目盛りの太さを設定
          axis.ticks.y.left =  element_line(linewidth = 1/2.13), # y軸左側の目盛りの線に適用
          # プロット全体を透過。color = NAで枠線も消す
          plot.background = element_rect(fill = "transparent", color = NA),
          # パネル部分を透過
          panel.background = element_rect(fill = "transparent")
    ) +
    scale_color_manual(values = c('#ff4b00', '#03af7a'))
  
  # Return the plot object
  return(p)
}

fit_and_plot_gamma_model <- function(data) {
  # Fit the Gamma GLM
  glm_model <- glm(
    EIbalance ~ Region,
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



# データ読み込み -----------------------------------------------------------------

raw_onlyChR2 <- read.csv("./Data/PatchClamp_escape//EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_ChR2.csv")
raw_dualChR2 <- read.csv("./Data/PatchClamp_escape//EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_Dual_Injection.csv")

data_onlyChR2 <- rawdata_to_dataframe(raw_onlyChR2)
data_dualChR2 <- rawdata_to_dataframe(raw_dualChR2) %>% 
  mutate(Opsin = "ChR2",
         Region = str_sub(Region, 1, 3))

data_allChR2 <- bind_rows(data_onlyChR2, data_dualChR2) %>% 
  filter(!is.na(PeakAmp_minus55_mean), 
         !is.na(PeakAmp_plus10_mean),
         n_minus55 > 5,
         Ra_minus55_mean < 45) %>% 
  mutate(EIbalance = PeakAmp_plus10_mean / -PeakAmp_minus55_mean,
         DecayTau_minus55_mean = 1000 * DecayTau_minus55_mean,
         DecayTau_minus55_std = 1000 * DecayTau_minus55_std,
         DecayTau_plus10_mean = 1000 * DecayTau_plus10_mean,
         DecayTau_plus10_std = 1000 * DecayTau_plus10_std,
         log_EIbalance = log(EIbalance),
         uniqueCellID = paste(BrainID, SliceID, CellID, sep = "_"))


# anterior ACC vs posterior RSC -------------------------------------------

RSC_posterior <- data_allChR2 %>% 
  filter(Region == "RSC", APregion == "posterior")
ACC_anterior <- data_allChR2 %>% 
  filter(Region == "ACC", APregion == "anterior")
aACCpRSC <- bind_rows(RSC_posterior, ACC_anterior) %>% 
  mutate(Region = fct_rev(Region))
aACCpRSC_small <- aACCpRSC %>% 
  filter(PeakAmp_minus55_mean < -10,
         PeakAmp_minus55_mean > -50) %>% 
  group_by(uniqueCellID) %>%
  slice(which.min(abs(PeakAmp_minus55_mean - (-40))))
aACCpRSC_large <- aACCpRSC %>% 
  filter(PeakAmp_minus55_mean < -70,
         PeakAmp_minus55_mean > -120) %>%  # Filter for rows where Im_scaled_E is below -80
  group_by(uniqueCellID) %>%
  slice_min(PeakAmp_minus55_mean) %>%  # dplyr::select the row with the smallest Im_scaled_E for each Cell_ID
  ungroup()

ggplot(aACCpRSC_small, aes(x = Latency_ms_plus10_mean)) +
  geom_histogram() +
  facet_wrap(~Region, ncol = 1)
ggplot(aACCpRSC_large, aes(x = EIbalance)) +
  geom_histogram() +
  facet_wrap(~Region, ncol = 1)

EIbalance_small <- plot_ei(aACCpRSC_small, EIbalance) + 
  scale_y_continuous(limits = c(0, 20.5), breaks = seq(0, 20, by = 5))
plot(EIbalance_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_EIbalance.svg", 
       plot = EIbalance_small,
       width = 24, height = 36, units = "mm")
EIbalance_large <- plot_ei(aACCpRSC_large, EIbalance) + 
  scale_y_continuous(limits = c(0, 20.5), breaks = seq(0, 20, by = 5))
plot(EIbalance_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_EIbalance.svg", 
       plot = EIbalance_large,
       width = 24, height = 36, units = "mm")
IPSCCV_small <- plot_ei(aACCpRSC_small, PeakAmp_plus10_CV) + 
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2, by = 0.5))
plot(IPSCCV_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_IPSCampCV.svg", 
       plot = IPSCCV_small,
       width = 24, height = 35, units = "mm")
IPSCCV_large <- plot_ei(aACCpRSC_large, PeakAmp_plus10_CV) + 
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2, by = 0.5))
plot(IPSCCV_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_IPSCampCV.svg", 
       plot = IPSCCV_large,
       width = 24, height = 35, units = "mm")

EPSCCV_small <- plot_ei(aACCpRSC_small, abs(PeakAmp_minus55_CV)) + 
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
plot(EPSCCV_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_EPSCampCV.svg", 
       plot = EPSCCV_small,
       width = 25, height = 30, units = "mm")
EPSCCV_large <- plot_ei(aACCpRSC_large, abs(PeakAmp_minus55_CV))
plot(EPSCCV_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_EPSCampCV.svg", 
       plot = EPSCCV_large,
       width = 25, height = 30, units = "mm")


aACCpRSC_small_IPSCkinetics <- aACCpRSC_small %>% 
  filter(R2_plus10_mean > 0.5,
         DecayTau_plus10_mean > 0)
aACCpRSC_large_IPSCkinetics <- aACCpRSC_large %>% 
  filter(R2_plus10_mean > 0.5,
         DecayTau_plus10_mean > 0)
aACCpRSC_small_EPSCkinetics <- aACCpRSC_small %>% 
  filter(R2_minus55_mean > 0.5,
         DecayTau_minus55_mean > 0)
aACCpRSC_large_EPSCkinetics <- aACCpRSC_large %>% 
  filter(R2_minus55_mean > 0.5,
         DecayTau_minus55_mean > 0)



IPSCDecay_small <- plot_ei(aACCpRSC_small_IPSCkinetics, DecayTau_plus10_mean) + 
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 200))
plot(IPSCDecay_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_IPSCDecay.svg", 
       plot = IPSCDecay_small,
       width = 24, height = 35, units = "mm")
IPSCDecay_large <- plot_ei(aACCpRSC_large_IPSCkinetics, DecayTau_plus10_mean) + 
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 200))
plot(IPSCDecay_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_IPSCDecay.svg", 
       plot = IPSCDecay_large,
       width = 24, height = 35, units = "mm")
IPSCLatency_small <- plot_ei(aACCpRSC_small_IPSCkinetics, Latency_ms_plus10_mean) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 15, by = 5))
plot(IPSCLatency_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_IPSCLatency.svg", 
       plot = IPSCLatency_small,
       width = 24, height = 35, units = "mm")
IPSCLatency_large <- plot_ei(aACCpRSC_large_IPSCkinetics, Latency_ms_plus10_mean) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 15, by = 5))
plot(IPSCLatency_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_IPSCLatency.svg", 
       plot = IPSCLatency_large,
       width = 24, height = 35, units = "mm")
IPSCRiseTime_small <- plot_ei(aACCpRSC_small_IPSCkinetics, RiseTime_ms_plus10_mean)+
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 50, by = 10))
plot(IPSCRiseTime_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_IPSCRiseTime.svg", 
       plot = IPSCRiseTime_small,
       width = 24, height = 35, units = "mm")
IPSCRiseTime_large <- plot_ei(aACCpRSC_large_IPSCkinetics, RiseTime_ms_plus10_mean) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 50, by = 10))
plot(IPSCRiseTime_large) 
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_IPSCRiseTime.svg", 
       plot = IPSCRiseTime_large,
       width = 24, height = 35, units = "mm")

aACCpRSC_small_EPSCkinetics_Decay <- aACCpRSC_small_EPSCkinetics %>% 
  filter(DecayTau_minus55_std < 100)
EPSCDecay_small <- plot_ei(aACCpRSC_small_EPSCkinetics_Decay, DecayTau_minus55_mean) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20))
plot(EPSCDecay_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_EPSCDecay.svg", 
       plot = EPSCDecay_small,
       width = 25, height = 30, units = "mm")
EPSCDecay_large <- plot_ei(aACCpRSC_large_EPSCkinetics, DecayTau_minus55_mean)
plot(EPSCDecay_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_EPSCDecay.svg", 
       plot = EPSCDecay_large,
       width = 25, height = 30, units = "mm")
EPSCLatency_small <- plot_ei(aACCpRSC_small_EPSCkinetics, Latency_ms_minus55_mean) +
  scale_y_continuous(limits = c(0, 0.01))
plot(EPSCLatency_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_EPSCLatency.svg", 
       plot = EPSCLatency_small,
       width = 25, height = 30, units = "mm")
EPSCLatency_large <- plot_ei(aACCpRSC_large_EPSCkinetics, Latency_ms_minus55_mean)
ggplot(aACCpRSC_large_EPSCkinetics, aes(x = Region, y = Latency_ms_minus55_mean)) +
  geom_point()
plot(EPSCLatency_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_EPSCLatency.svg", 
       plot = EPSCLatency_large,
       width = 25, height = 30, units = "mm")
EPSCRiseTime_small <- plot_ei(aACCpRSC_small_EPSCkinetics, RiseTime_ms_minus55_mean)
plot(EPSCRiseTime_small)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_small_EPSCRiseTime.svg", 
       plot = EPSCRiseTime_small,
       width = 25, height = 30, units = "mm")
EPSCRiseTime_large <- plot_ei(aACCpRSC_large_EPSCkinetics, RiseTime_ms_minus55_mean)
plot(EPSCLatency_large)
ggsave("./Data/PatchClamp_escape/EIkinetics_APaxis/plots/aACCpRSC_large_EPSCRiseTime.svg", 
       plot = EPSCLatency_large,
       width = 25, height = 30, units = "mm")

# GLM ---------------------------------------------------------------------

glm_aACCpRSC_small_results <- fit_and_plot_gamma_model(aACCpRSC_small)
glm_aACCpRSC_large_results <-fit_and_plot_gamma_model(aACCpRSC_large)

# Wilcox rank sum test ----------------------------------------------------
# Wilcoxon rank sum test for EIbalance between Regions in aACCpRSC_small
wilcox_test_result <- wilcox.test(EIbalance ~ Region, data = aACCpRSC_small)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(EIbalance ~ Region, data = aACCpRSC_large)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(PeakAmp_plus10_CV ~ Region, data = aACCpRSC_small)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(PeakAmp_plus10_CV ~ Region, data = aACCpRSC_large)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(DecayTau_plus10_mean ~ Region, data = aACCpRSC_small_IPSCkinetics)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(DecayTau_plus10_mean ~ Region, data = aACCpRSC_large_IPSCkinetics)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(Latency_ms_plus10_mean ~ Region, data = aACCpRSC_small_IPSCkinetics)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(Latency_ms_plus10_mean ~ Region, data = aACCpRSC_large_IPSCkinetics)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(RiseTime_ms_plus10_mean ~ Region, data = aACCpRSC_small_IPSCkinetics)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(RiseTime_ms_plus10_mean ~ Region, data = aACCpRSC_large_IPSCkinetics)
print(wilcox_test_result)

*wilcox_test_result <- wilcox.test(RiseTime_ms_minus55_mean ~ Region, data = aACCpRSC_small_EPSCkinetics)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(PeakAmp_minus55_CV ~ Region, data = aACCpRSC_small)
print(wilcox_test_result)
wilcox_test_result <- wilcox.test(DecayTau_minus55_mean ~ Region, data = aACCpRSC_small_EPSCkinetics_Decay)
print(wilcox_test_result)


# Pythonにabf file nameわたす -------------------------------------------------
# Data/PatchClamp_escape/EIkinetics_APaxis/abfピーク確認.ipynbにterminalに出るfile名をコピペする
filenames <- unique(aACCpRSC_large$filename_plus10)

cat(sprintf("filenames = [\n%s\n]\n",
            paste0('    "', filenames, '"', collapse = ",\n")))

