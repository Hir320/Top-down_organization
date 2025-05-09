rm(list =ls())
library(tidyverse)


# Functions ---------------------------------------------------------------

plot_intensity_ratio <- function(Intensity_rawdata) {
  # Ensure required libraries are loaded
  
  # Split the data for ACC and RSC by patch and interpatch
  ACC_patch <- Intensity_rawdata %>% 
    filter(Region == "ACC", PatchInterpatch == "patch")
  ACC_Interpatch <- Intensity_rawdata %>% 
    filter(Region == "ACC", PatchInterpatch == "interpatch")
  RSC_patch <- Intensity_rawdata %>% 
    filter(Region == "RSC", PatchInterpatch == "patch")
  RSC_Interpatch <- Intensity_rawdata %>% 
    filter(Region == "RSC", PatchInterpatch == "interpatch")
  
  # Join patch and interpatch data for ACC and RSC
  ACC <- left_join(ACC_patch, ACC_Interpatch, by = c("MouseID", "Region"), 
                   suffix = c("_patch", "_interpatch"))
  RSC <- left_join(RSC_patch, RSC_Interpatch, by = c("MouseID", "Region"), 
                   suffix = c("_patch", "_interpatch"))
  
  # Combine the data and calculate the intensity ratio
  Data <- bind_rows(ACC, RSC) %>% 
    mutate(Intensity_ratio = Mean_patch / Mean_interpatch)
  
  # Select only the necessary columns
  Data_select <- Data %>% 
    select(MouseID, Region, Intensity_ratio)
  
  # Create the ggplot object
  p <- ggplot(Data_select, aes(x = Region, y = Intensity_ratio, color = Region)) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1) +
    stat_summary(fun = mean, geom = "point", size = 2, shape = 21, fill = "white") +
    geom_jitter(width = 0.05, size = 1) +
    scale_color_manual(values = c("RSC" = "#ff4b00", "ACC" = "#03af7a")) +
    theme_classic() +
    labs(x = "Region", y = "Intensity Ratio") +
    theme(axis.title = element_blank(),
          legend.position = "none",
          axis.ticks.length.x.bottom = unit(0, "mm"),
          axis.ticks.length.y.left = unit(0.5, "mm"),
          axis.text.y              = element_text(size = 8, family = "Arial"),
          axis.text.x              = element_blank(),
          axis.line                = element_line(linewidth = 0.5/2.13),
          axis.ticks               = element_line(linewidth = 0.5/2.13),
          plot.background          = element_rect(fill = "transparent", color = NA),
          panel.background         = element_rect(fill = "transparent"))
  
  # Return the plot 
  return(p)
}


# Ai34 --------------------------------------------------------------------

Intensity_rawdata_Ai34 <- read.csv("./Data/Histology/Retinotopic_quantification/Patch_Interpatch_intensity/Ai34/Intensity.csv") %>% 
  separate(
    col = Label, 
    into = c("MouseID", "Region", "PatchInterpatch"),
    sep = "_",
    remove = FALSE) %>% 
  select(-Label) %>% 
  mutate(Region = fct_rev(Region))

Ai34_plot <- plot_intensity_ratio(Intensity_rawdata_Ai34)
plot(Ai34_plot)
ggsave("./Data/Histology/Retinotopic_quantification/Patch_Interpatch_intensity/Ai34/IntensityRatio.svg", 
       width = 30, height = 32, units = "mm")

# PV-Cre --------------------------------------------------------------------
Intensity_rawdata_PV <- read.csv("./Data/Histology/Retinotopic_quantification/Patch_Interpatch_intensity/PV-Cre/Intensity.csv") %>% 
  separate(
    col = Label, 
    into = c("MouseID", "Region", "PatchInterpatch"),
    sep = "_",
    remove = FALSE) %>% 
  select(-Label) %>% 
  mutate(Region = fct_rev(Region))

PV_plot <- plot_intensity_ratio(Intensity_rawdata_PV) +
  scale_y_continuous(limits = c(0.4, 1.6), breaks = seq(0, 2, by = 0.5))
plot(PV_plot)
ggsave("./Data/Histology/Retinotopic_quantification/Patch_Interpatch_intensity/PV-Cre/IntensityRatio.svg", 
       width = 45, height = 34, units = "mm")
   