# ======================================================================
#  Enhancement Index  – ACSF  vs.  4-AP (≈3 mM)
#    • ACC-R → RSC-2
#    • RSC-R → ACC-2
#  (plot creation separated from saving)
# ======================================================================
rm(list = ls())

# ---- 1. packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# ---- 2. project root --------------------------------------------------
setwd("/home/rstudio/")            # ← 必要なら変更

# ---- 3. helpers -------------------------------------------------------
wrangle_integration <- function(files, extra_paths = character(0)){
  raw <- map_dfr(c(files, extra_paths),
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
      
      ## ---- condition flags -------------------------------------------
      PTX  = str_detect(log, regex("PTX", ignore_case = TRUE)),
      AP4  = str_detect(log, regex("4\\s*-?\\s*AP", ignore_case = TRUE)),
      
      ## ---- concentration numbers -------------------------------------
      PTX_conc_num = parse_number(
        str_extract(log,
                    regex("PTX\\s+[0-9]+(?:\\.[0-9]+)?\\s*(?:uM|mM)", ignore_case = TRUE))
      ),
      AP4_conc_num = parse_number(
        str_extract(log,
                    regex("4\\s*-?\\s*AP\\s+[0-9]+(?:\\.[0-9]+)?\\s*(?:uM|mM)", ignore_case = TRUE))
      ),
      
      ## ---- cell & stim descriptors -----------------------------------
      Cell_ID = paste(str_sub(filename, 1, 5),
                      str_sub(SliceID, -1, -1),
                      str_sub(CellID , -1, -1), sep = "_"),
      
      Matches = str_extract_all(log, stim_regex),
      Stim1   = map_chr(Matches, 1, .default = NA_character_),
      Stim2   = map_chr(Matches, 2, .default = NA_character_)
    ) |>
    filter(Interval > 0,
           str_detect(log, "-75 mV"),         # 保持電位フィルタ
           !is.na(Stim1), !is.na(Stim2)) |>
    select(Cell_ID, PTX, PTX_conc_num,
           AP4, AP4_conc_num,
           Interval, Stim1, Stim2,
           `Actual Peak Value`, `Estimated Peak Value`) |>
    group_by(Cell_ID, PTX, PTX_conc_num,
             AP4, AP4_conc_num,
             Interval, Stim1, Stim2) |>
    summarise(.groups = "drop",
              Actual_Peak    = mean(`Actual Peak Value`,    na.rm = TRUE),
              Estimated_Peak = mean(`Estimated Peak Value`, na.rm = TRUE)) |>
    mutate(Enhancement_index = Actual_Peak / Estimated_Peak)
}

## ---- ACSF vs 4-AP (≈3 mM) --------------------------------------------
pair_acsf_4ap <- function(tbl) {
  
  acsf <- tbl |>
    filter(AP4 == FALSE) |>
    select(Cell_ID, Interval, Stim1, Stim2, Enhancement_index) |>
    rename(Enh_acsf = Enhancement_index)
  
  ap4 <- tbl |>
    filter(AP4 == TRUE,
           between(AP4_conc_num, 2, 4)) |>  # 3 mM ±1 mM を許容
    select(Cell_ID, Interval, Stim1, Stim2, Enhancement_index) |>
    rename(Enh_ap4 = Enhancement_index)
  
  inner_join(acsf, ap4,
             by = c("Cell_ID", "Interval", "Stim1", "Stim2")) |>
    filter(!is.na(Enh_acsf), !is.na(Enh_ap4))
}

## ---- plotting helper (変更ナシ) ---------------------------------------
make_paired_panel <- function(paired_tbl, y_lim = c(0, 2)) {
  
  paired_tbl |>
    pivot_longer(c(Enh_acsf, Enh_ap4),
                 names_to  = "Cond",
                 values_to = "Enh") |>
    mutate(Cond = fct_rev(recode(Cond,
                         Enh_acsf = "ACSF",
                         Enh_ap4  = "4-AP 3 mM"))) |>
    ggplot(aes(Cond, Enh, group = Cell_ID)) +
    geom_line(colour = "grey60", linewidth = 0.3) +
    geom_point(size = 1) +
    facet_wrap(
      ~ Interval, nrow = 1,
      labeller = labeller(
        Interval = function(x) paste0(as.numeric(x) * 1000, " ms")
      )
    ) +
    ylim(y_lim) +
    theme_classic() +
    theme(
      axis.ticks.length.x.bottom = unit(0, "mm"),
      axis.ticks.length.y.left   = unit(1, "mm"),
      axis.line   = element_line(linewidth = 1/2.13),
      axis.ticks  = element_line(linewidth = 1/2.13),
      axis.title  = element_blank(),
      axis.text.x = element_blank(),
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
      mean_ap4  = mean(Enh_ap4,  na.rm = TRUE),
      .test     = list(t.test(Enh_acsf, Enh_ap4, paired = TRUE)),
      .groups   = "drop"
    ) |>
    mutate(stats = map(.test, broom::tidy)) |>
    unnest(stats) |>
    ## ── 自由度 (parameter) を df に改名 ─────────────────────────
    rename(df = parameter) |>
    select(Interval, n_pairs,
           mean_acsf, mean_ap4,
           estimate, statistic, df, p.value)
}

# ---- 4. load & organise ----------------------------------------------
file_paths <- list(
  "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/250222(RSC ChrimsonR ACC-ChR2)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/250224(RSC ChrimsonR ACC-ChR2)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/250226(RSC ChrimsonR ACC-ChR2)/merged_mean_peak_edited.csv"
)

## extra_paths を使う場合はベクトルで渡す。無ければ character(0) で OK
#ACC_R_cell <- wrangle_integration(file_paths)   # 例として読み直し
RSC_R_cell <- wrangle_integration(file_paths)   # 本来は別セット

# ---- 5. pairing -------------------------------------------------------
#paired_ACC_R <- pair_acsf_4ap(ACC_R_cell)
paired_RSC_R <- pair_acsf_4ap(RSC_R_cell)

# ---- 6. build ggplot objects -----------------------------------------
#p_acc_4ap <- make_paired_panel(paired_ACC_R)
p_rsc_4ap <- make_paired_panel(paired_RSC_R)

#print(p_acc_4ap)
print(p_rsc_4ap)

paired_RSC_R_no25 <- paired_RSC_R |> 
  filter(Interval != 0.025)

p_rsc_4ap_no25 <- make_paired_panel(paired_RSC_R_no25)+
  scale_y_continuous(limits = c(0, 1.5),
                     breaks  = seq(0, 2, 0.5),
                     expand  = expansion(mult = c(0, 0.02)))

print(p_rsc_4ap_no25)




# ---- 7. save ----------------------------------------------------------
#ggsave("./Data/PatchClamp/graphs/TwoColour/ACC-R_RSC-2_4AP.svg", p_acc_4ap, width = 70, height = 35, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_4AP.svg",
       p_rsc_4ap, width = 50, height = 40, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_4AP_no25ms.svg",
       p_rsc_4ap_no25, width = 44, height = 45, units = "mm") 

# RSC-R → ACC-2
t_stat_rsc <- calc_paired_t(paired_RSC_R)
print(t_stat_rsc)