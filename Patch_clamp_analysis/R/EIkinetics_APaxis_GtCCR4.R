rm(list =ls())
library(tidyverse)
library(broom)
library(lme4)
library(lmerTest)
library(ggpubr)

# Define the directory where you want to search for the files
setwd('/home/rstudio/')

df_GtCCR4 <- read.csv("./Data/PatchClamp/EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_GtCCR4.csv") %>% 
  filter(Color == "blue") %>% 
  group_by(filename_plus10) %>%
  slice_min(filename_minus55, with_ties = FALSE) %>%
  ungroup() %>% 
  group_by(filename_minus55) %>%
  slice_min(filename_plus10, with_ties = FALSE) %>%
  ungroup()

write.csv(df_GtCCR4, "./Data/PatchClamp/EIkinetics_APaxis/sorted_directory/Result_of_EIkinetics_GtCCR4_filtered.csv")

df_GtCCR4_EPSCfiltered <- df_GtCCR4 %>% 
  filter(PeakAmp_minus55_mean < -50) %>% 
  group_by(Opsin, Region, BrainID, SliceID, CellID) %>%
  slice(which.min(abs(PeakAmp_minus55_mean - (-100)))) %>% 
  mutate(EIbalance = abs(PeakAmp_plus10_mean / PeakAmp_minus55_mean))

df_GtCCR4_EPSCfiltered %>%
  group_by(Region) %>%
  do(
    broom::tidy(
      cor.test(~ EIbalance + RoughAP, data = .)
    )
  )

ggplot(df_GtCCR4_EPSCfiltered, aes(x = RoughAP, y = EIbalance, color = BrainID)) +
  # Points
  geom_point() +
  # Linear fit
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  # Add correlation coefficients and p-values
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  # Separate plots by region
  facet_wrap(~ Region, scales = "free") +
  theme_classic()

df_rep <- df_GtCCR4 %>%
  group_by(Opsin, Region, BrainID, SliceID, CellID) %>%
  filter(StimPower == 100) %>% 
  mutate(dist_duration = abs(StimDuration - median(StimDuration, na.rm = TRUE))) %>%
  # Keep the row with the smallest distance
  slice_min(dist_duration, n = 1) %>%
  # If there’s a tie, slice_min() keeps all ties; slice(1) would keep only the first
  ungroup()

# Example of random intercept for BrainID & SliceID
model_mixed <- lmer(
  PeakAmp_minus55_mean ~ Region + RoughAP + StimDuration + 
    (1 | BrainID) + (1 | SliceID),
  data = df_rep
)

summary(model_mixed)
anova(model_mixed, type = "III")  # type III if you'd like to see main effects

ggplot(df_rep, aes(x = RoughAP, y = PeakAmp_plus10_mean, color = Region)) +
  geom_point(aes(size = StimDuration)) +
  facet_wrap(~ Region) +
  theme_classic()


df_mode_stim <- df_GtCCR4 %>%
  group_by(Opsin, Region, BrainID, StimDuration) %>%
  summarize(count = n()) %>%
  # For each BrainID, pick StimDuration with the highest count
  slice_max(count, with_ties = FALSE) %>%
  rename(mode_StimDuration = StimDuration) %>%
  ungroup()

df_consistent <- df_GtCCR4 %>%
  filter(StimPower == 100) %>% 
  left_join(df_mode_stim, by = "BrainID") %>%
  filter(StimDuration == mode_StimDuration) %>%
  select(-count, -mode_StimDuration)  # clean up extra columns if you like

df_norm_z <- df_GtCCR4_EPSCfiltered %>%
  group_by(BrainID) %>%
  mutate(
    EIbalance_z = (EIbalance - mean(EIbalance, na.rm = TRUE)) / sd(EIbalance, na.rm = TRUE)
  ) %>%
  ungroup() 

# Then do correlation on EIbalance_z vs RoughAP
df_norm_z %>%
  group_by(Region) %>%
  do(
    broom::tidy(
      cor.test(~ EIbalance_z + RoughAP, data = .)
    )
  )
ggplot(df_norm_z, aes(x = RoughAP, y = EIbalance_z)) +
  # Points
  #geom_point(aes(color = BrainID)) +
  geom_point() +
  # Linear fit
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  # Add correlation coefficients and p-values
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  # Separate plots by region
  facet_wrap(~ Region, scales = "free") +
  theme_classic() +
  theme(legend.position = "None")
