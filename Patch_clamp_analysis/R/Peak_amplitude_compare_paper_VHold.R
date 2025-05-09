rm(list =ls())
# ────────────────────────────────────────────────────────────────────────
#  0.  Packages
# ────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)     # readr, dplyr, ggplot2, stringr, purrr, tibble …
  library(janitor)       # clean_names()
  library(plotrix)       # std.error()
  library(ggbreak)       # axis break
})

# ────────────────────────────────────────────────────────────────────────
#  1.  Global STYLE object  (edit here ☞ plots update everywhere)
# ────────────────────────────────────────────────────────────────────────
sty <- list(
  # axis breaks (ggbreak)
  x_break        = c(0.26, 0.75),
  x_limits       = c(0, 0.85),
  x_breaks       = c(0.025, 0.05, 0.1, 0.25, 0.8),
  x_labels       = c("25", "50", "100", "250", "800"),
  
  # y-axis
  y_limits       = c(0.5, 1.3),
  
  # points & lines
  cell_pt_size   = 0.5,
  cell_pt_alpha  = 0.10,
  mean_pt_size   = 1.5,
  mean_line_lwd  = 1/2.13,
  err_bar_width  = 0.01,
  err_bar_lwd  = 1/2.13,
  
  axis_line_lwd  = 1/2.13,
  tick_line_lwd  = 1/2.13,
  
  col_vhold      = c(`-75 mV` = "black",`-55 mV` = "red"),
  
  # fonts
  base_font_sz   = 8,
  base_font_fam  = "Arial",
  
  # export size
  fig_w_mm       = 60,
  fig_h_mm       = 46.5
)

# ────────────────────────────────────────────────────────────────────────
#  2.  Helper – read ONE file per leaf directory
#      (“*_edited.csv” wins over plain “*.csv”)
# ────────────────────────────────────────────────────────────────────────
load_peak_files <- function(dir_path) {
  
  # 1. Collect both edited / un-edited ------------------------------------------------
  list.files(
    path       = dir_path,
    pattern    = "merged_mean_peak(_edited)?\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  ) |>
    # 2. Keep ONLY ONE per leaf directory --------------------------------------------
  tibble(path = _) |>
    mutate(
      leaf_dir = dirname(path),                     # e.g. ".../240419(...)"
      is_edit  = str_detect(path, "_edited\\.csv$")
    ) |>
    group_by(leaf_dir) |>
    # slice_max on is_edit = TRUE wins; if none TRUE, takes the plain file
    slice_max(is_edit, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(path) |>
    
    # 3. Read & clean ---------------------------------------------------------------
  map_dfr(read_csv, show_col_types = FALSE) |>
    clean_names()                                # spaces → snake_case
}

# ────────────────────────────────────────────────────────────────────────
#  3.  Wrangle one directory  →  “cell-level” tibble
# ────────────────────────────────────────────────────────────────────────
wrangle_dir <- function(dir_path, vhold_label) {
  
  regex_stim <- "\\b[A-Z]_\\d+ A (\\d+\\.?\\d*%?|\\d*\\.?\\d+%?) \\d+\\.?\\d* ms\\b"
  
  load_peak_files(dir_path) |>
    mutate(
      ptx       = str_detect(log,   regex("PTX", ignore_case = TRUE)),
      cell_id   = paste(
        str_sub(filename,  1, 5),
        str_sub(slice_id, -1, -1),
        str_sub(cell_id,  -1, -1),
        sep = "_"
      ),
      matches   = str_extract_all(log, regex_stim),
      stim1     = map_chr(matches, 1, .default = NA),
      stim2     = map_chr(matches, 2, .default = NA),
      vhold     = vhold_label,
      
      # keep numeric columns numeric
      across(
        c(interval,
          first_peak_value_combined,
          third_peak_value_combined,
          actual_peak_value,
          estimated_peak_value), as.numeric
      )
    ) |>
    ## ── NEW: keep only stim1 that begins with "R_1" ───────────────
    filter(interval > 0,
           !is.na(stim1), !is.na(stim2),
           str_detect(stim1, "^R_1")) |>      #  ← 追加ここだけ
    ## ───────────────────────────────────────────────────────────────
    select(cell_id, ptx, vhold, interval, stim1, stim2,
           first_peak_value_combined, third_peak_value_combined,
           actual_peak_value, estimated_peak_value) |>
    group_by(ptx, interval, cell_id, stim1, stim2, vhold) |>
    summarise(across(
      c(first_peak_value_combined,
        third_peak_value_combined,
        actual_peak_value,
        estimated_peak_value),
      mean, na.rm = TRUE),
      .groups = "drop") |>
    mutate(enhancement_index = actual_peak_value / estimated_peak_value)
}

# ────────────────────────────────────────────────────────────────────────
#  4.  Summary to “mean table” (per interval, per Vhold, ±SEM)
# ────────────────────────────────────────────────────────────────────────
make_mean_tbl <- function(cell_tbl) {
  cell_tbl |>
    group_by(vhold, ptx, interval) |>
    summarise(
      actual_peak             = mean(actual_peak_value),
      estimated_peak          = mean(estimated_peak_value),
      enhancement_index_mean  = mean(enhancement_index),
      enhancement_index_sem   = std.error(enhancement_index),
      n = n(),
      .groups = "drop"
    )
}

# ────────────────────────────────────────────────────────────────────────
#  5-A.  Plot helper: cell-level spaghetti
# ────────────────────────────────────────────────────────────────────────
plot_cells <- function(cell_tbl, S = sty) {
  
  ggplot(cell_tbl,
         aes(ptx, enhancement_index, group = cell_id)) +
    geom_line(size = S$mean_line_lwd, alpha = 0.5) +
    facet_wrap(~ interval) +
    ylim(S$y_limits) +
    theme_minimal(base_size = S$base_font_sz, base_family = S$base_font_fam)
}

# ────────────────────────────────────────────────────────────────────────
#  5-B.  Plot helper: mean ± SEM  (optionally ACSF only, with ggbreak)
# ────────────────────────────────────────────────────────────────────────
plot_mean <- function(mean_tbl, S = sty,
                      acsf_only = FALSE,
                      line_cut  = 0.25) {
  
  df <- if (acsf_only) filter(mean_tbl, ptx == FALSE) else mean_tbl
  df_line  <- filter(df, interval <= line_cut)
  df_point <- df
  
  ggplot() +
    # 点
    geom_point(data = df_point,
               aes(interval, enhancement_index_mean,
                   colour = vhold, group = vhold),
               size = S$mean_pt_size) +
    # 線
    geom_line(data = df_line,
              aes(interval, enhancement_index_mean,
                  colour = vhold, group = vhold),
              linewidth = S$mean_line_lwd) +
    # 誤差バー
    geom_errorbar(data = df_point,
                  aes(interval,
                      ymin = enhancement_index_mean - enhancement_index_sem,
                      ymax = enhancement_index_mean + enhancement_index_sem,
                      colour = vhold, group = vhold),
                  width = S$err_bar_width,
                  linewidth = S$err_bar_lwd) +
    
    ## —— 手動色指定 & 凡例非表示 —— ##
    scale_color_manual(values = S$col_vhold, guide = "none") +
    
    # 軸などは今まで通り
    scale_x_break(breaks = S$x_break, scales = 0.4, space = 0.02) +
    scale_x_continuous(limits = S$x_limits,
                       breaks  = S$x_breaks,
                       labels  = S$x_labels) +
    ylim(S$y_limits) +
    labs(x = "Interval (s)", y = "Enhancement index") +
    theme_classic(base_size = S$base_font_sz, base_family = S$base_font_fam) +
    theme(
      axis.text  = element_text(size = S$axis_font_sz,
                                family = S$axis_font_fam),
      axis.title = element_blank(),
      axis.line  = element_line(linewidth = S$axis_line_lwd),
      axis.ticks = element_line(linewidth = S$tick_line_lwd),
      axis.ticks.length.x.bottom = unit(1, "mm"),
      axis.ticks.length.y.left   = unit(1, "mm"),
      axis.line.x.top   = element_blank(),
      axis.ticks.x.top  = element_blank(),
      axis.text.x.top   = element_blank(),
      axis.title.x.top  = element_blank(),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent")
    )
}

# ────────────────────────────────────────────────────────────────────────
#  6.  Import → process
# ────────────────────────────────────────────────────────────────────────
dir55 <- "./Data/PatchClamp/Integration/Integration_VHold-55mV/RSC-R_ACC-2/"
dir75 <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/"

cell55 <- wrangle_dir(dir55, "-55 mV")
cell75 <- wrangle_dir(dir75, "-75 mV")

cell_all  <- bind_rows(cell55, cell75)
mean_all  <- make_mean_tbl(cell_all)

# ────────────────────────────────────────────────────────────────────────
#  7.  Make plots
# ────────────────────────────────────────────────────────────────────────
p_spaghetti <- plot_cells(cell_all)
p_mean_all  <- plot_mean(mean_all)
p_mean_acsf <- plot_mean(mean_all, acsf_only = TRUE) 
print(p_mean_acsf)

# ────────────────────────────────────────────────────────────────────────
#  8.  Save
# ────────────────────────────────────────────────────────────────────────
ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_spaghetti_Vhold.svg",
       p_spaghetti,
       width = sty$fig_w_mm, height = sty$fig_h_mm, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_mean_Vhold.svg",
       p_mean_all,
       width = sty$fig_w_mm, height = sty$fig_h_mm, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_mean_ACSF_Vhold.svg",
       p_mean_acsf,
       width = sty$fig_w_mm, height = sty$fig_h_mm, units = "mm")
