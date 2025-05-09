# ======================================================================
#  Enhancement Index  – ACSF  vs.  25 µM PTX
#    • ACC-R → RSC-2
#    • RSC-R → ACC-2
#  (plot creation separated from saving)
# ======================================================================
rm(list =ls())
# ---- 1. packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)   # readr, dplyr, ggplot2, stringr, purrr, tidyr
  library(broom)
})

# ---- 2. project root --------------------------------------------------
setwd("/home/rstudio/")            # ← change if necessary

# ---- 3. helpers -------------------------------------------------------
prefer_edited <- function(paths){
  edited <- sub("merged_mean_peak.csv", "merged_mean_peak_edited.csv", paths)
  ifelse(file.exists(edited), edited, paths)
}

wrangle_integration <- function(dir_path, extra_paths = character(0)){
  # 3-A read ------------------------------------------------------------
  files <- list.files(dir_path, "merged_mean_peak.csv",
                      recursive = TRUE, full.names = TRUE) |>
    prefer_edited()
  files <- c(files, extra_paths[file.exists(extra_paths)])
  
  raw <- map_dfr(files,
                 read_csv,
                 col_types = cols(.default = col_character()),
                 show_col_types = FALSE) |>
    select(-matches("^\\.\\.\\."))              # drop unnamed cols
  
  # 3-B clean -----------------------------------------------------------
  stim_regex <- "\\b[A-Z]_\\d+ A (\\d+\\.?\\d*%?|\\d*\\.?\\d+%?) \\d+\\.?\\d* ms\\b"
  
  raw |>
    mutate(
      across(c(Interval,
               `Actual Peak Value`,
               `Estimated Peak Value`), as.numeric),
      
      PTX = str_detect(log, regex("PTX", ignore_case = TRUE)),
      PTX_conc_num = parse_number(
        str_extract(log,
                    regex("PTX\\s+[0-9]+(?:\\.[0-9]+)?\\s*(?:uM|mM)", ignore_case = TRUE))
      ),
      
      Cell_ID = paste(str_sub(filename, 1, 5),
                      str_sub(SliceID, -1, -1),
                      str_sub(CellID , -1, -1), sep = "_"),
      
      Matches = str_extract_all(log, stim_regex),
      Stim1   = map_chr(Matches, 1, .default = NA),
      Stim2   = map_chr(Matches, 2, .default = NA)
    ) |>
    filter(Interval > 0,
           str_detect(log, "-75 mV"),
           !is.na(Stim1), !is.na(Stim2)) |>
    select(Cell_ID, PTX, PTX_conc_num, Interval,
           Stim1, Stim2,
           `Actual Peak Value`, `Estimated Peak Value`) |>
    group_by(Cell_ID, PTX, PTX_conc_num, Interval, Stim1, Stim2) |>
    summarise(.groups = "drop",
              Actual_Peak    = mean(`Actual Peak Value`,    na.rm = TRUE),
              Estimated_Peak = mean(`Estimated Peak Value`, na.rm = TRUE)) |>
    mutate(Enhancement_index = Actual_Peak / Estimated_Peak)
}

pair_acsf_ptx25 <- function(tbl) {
  
  acsf <- tbl %>%
    filter(PTX == FALSE) %>%
    select(Cell_ID, Interval, Stim1, Stim2, Enhancement_index) %>%
    rename(Enh_acsf = Enhancement_index)          # ← ここで改名
  
  ptx  <- tbl %>%
    filter(PTX == TRUE,
           between(PTX_conc_num, 23, 27)) %>%     # 25 µM ± 2
    select(Cell_ID, Interval, Stim1, Stim2, Enhancement_index) %>%
    rename(Enh_ptx = Enhancement_index)           # ← ここで改名
  
  inner_join(acsf, ptx,
             by = c("Cell_ID", "Interval", "Stim1", "Stim2")) %>%
    # ---- drop rows where either value is NA ---------------------------
  filter(!is.na(Enh_acsf), !is.na(Enh_ptx))
}
make_paired_panel <- function(paired_tbl, y_lim = c(0, 1.5)) {
  
  paired_tbl |>
    pivot_longer(c(Enh_acsf, Enh_ptx),
                 names_to  = "Cond",
                 values_to = "Enh") |>
    mutate(Cond = recode(Cond,
                         Enh_acsf = "ACSF",
                         Enh_ptx  = "PTX 25 µM")) |>
    ggplot(aes(Cond, Enh, group = Cell_ID)) +
    geom_line(colour = "grey60", linewidth = 0.3) +
    geom_point(size = 1) +
    facet_wrap(
      ~ Interval, nrow = 1,
      labeller = labeller(
        Interval = function(x) paste0(as.numeric(x) * 1000, " ms")  # ← 修正点
      )
    ) +
    ylim(y_lim) +
    theme_classic() +
    theme(
      axis.ticks.length.x.bottom = unit(0, "mm"),
      axis.ticks.length.y.left   = unit(1, "mm"),
      axis.line                = element_line(linewidth = 1/2.13),
      axis.ticks               = element_line(linewidth = 1/2.13),
      axis.title  = element_blank(),
      #axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      strip.background = element_blank(),
      strip.text  = element_text(size = 8),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent")
    )
}
calc_paired_t <- function(paired_tbl) {
  paired_tbl |>
    group_by(Interval) |>
    summarise(
      n_pairs   = n(),
      mean_acsf = mean(Enh_acsf, na.rm = TRUE),
      mean_ptx  = mean(Enh_ptx,  na.rm = TRUE),
      .test     = list(t.test(Enh_acsf, Enh_ptx, paired = TRUE)),
      .groups   = "drop"
    ) |>
    mutate(stats = map(.test, broom::tidy)) |>
    unnest(stats) |>
    ## ── 自由度 (parameter) を df に改名 ─────────────────────────
    rename(df = parameter) |>
    select(Interval, n_pairs,
           mean_acsf, mean_ptx,
           estimate, statistic, df, p.value)
}

# ---- 4. load & organise ----------------------------------------------
# 4-A  ACC-R → RSC-2
acc_dir <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-2_ACC-R/"
ACC_R_cell <- wrangle_integration(acc_dir)

# 4-B  RSC-R → ACC-2
rsc_dir <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/"
manual_path <- file.path(rsc_dir,
                         "240817(ACC ChR2+ RSC ChrimsonR)",
                         "merged_mean_peak_edited.csv")
RSC_R_cell <- wrangle_integration(rsc_dir, extra_paths = manual_path)

paired_ACC_R <- pair_acsf_ptx25(ACC_R_cell)
paired_RSC_R <- pair_acsf_ptx25(RSC_R_cell)

# ---- 5. build ggplot objects (no saving yet) -------------------------
p_acc_ptx25 <- make_paired_panel(paired_ACC_R)
p_rsc_ptx25 <- make_paired_panel(paired_RSC_R) +
  scale_y_continuous(limits = c(0, 1.5),
                     breaks  = seq(0, 2, 0.5),
                     expand  = expansion(mult = c(0, 0.02)))

# ---- 6. display (interactive) ----------------------------------------
print(p_acc_ptx25)
print(p_rsc_ptx25)

# ---- 7. save explicitly ----------------------------------------------
ggsave("./Data/PatchClamp/graphs/TwoColour/ACC-R_RSC-2_PTX25.svg",
       p_acc_ptx25, width = 70, height = 35, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_PTX25.svg",
       p_rsc_ptx25, width = 100, height = 45, units = "mm")
# RSC-R → ACC-2
t_stat_rsc <- calc_paired_t(paired_RSC_R)
print(t_stat_rsc)