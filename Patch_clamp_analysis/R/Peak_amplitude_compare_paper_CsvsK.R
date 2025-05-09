rm(list =ls())
# ────────────────────────────────────────────────────────────────────────
#  0.  Packages
# ────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)      # readr, dplyr, ggplot2, stringr, purrr, tibble …
  library(janitor)        # clean_names()
  library(plotrix)        # std.error()
  library(ggbreak)        # axis break
})

# ────────────────────────────────────────────────────────────────────────
#  1.  GLOBAL STYLE (edit just this section to restyle everything)
# ────────────────────────────────────────────────────────────────────────
sty <- list(
  # axis breaks
  x_break        = c(0.26, 0.75),
  x_limits       = c(0, 0.85),
  x_breaks       = c(0.025, 0.05, 0.1, 0.25, 0.8),
  x_labels       = c("25", "50", "100", "250", "800"),
  
  # y-axis
  y_limits       = c(0.5, 1.2),
  
  # glyphs & lines
  cell_pt_size   = 0.5,
  cell_pt_alpha  = 0.10,
  mean_pt_size   = 1.5,
  mean_line_lwd  = 1 / 2.13,
  err_bar_width  = 0.01,
  err_bar_lwd    = 1/2.13,
  
  axis_line_lwd  = 1 / 2.13,          # Illustrator 0.5 pt
  tick_line_lwd  = 1 / 2.13,
  
  # COLOURS — Cs = red, K = black
  col_intra      = c(Cs = "red", K = "black"),
  
  # fonts
  base_font_sz   = 8,
  base_font_fam  = "Arial",
  axis_font_sz   = 8,
  axis_font_fam  = "Arial",
  
  # export size
  fig_w_mm       = 60,
  fig_h_mm       = 46.5
)


# ────────────────────────────────────────────────────────────────────────
#  2-A.  CSV loader for a *directory* (keeps _edited if present)
# ────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────
# 2-A. CSV loader for a *directory*
#      ── keeps **one** file per leaf folder:
#         • merged_mean_peak_edited.csv  (if present)
#         • otherwise merged_mean_peak.csv
# ───────────────────────────────────────────────────────────────────────
load_peak_dir <- function(dir_path) {
  
  list.files(
    path       = dir_path,
    pattern    = "merged_mean_peak(_edited)?\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  ) |>
    # ── group by *leaf directory* ────────────────────────────────
    tibble(path = _) |>
    mutate(leaf_dir = dirname(path),
           is_edit  = str_detect(path, "_edited\\.csv$")) |>
    group_by(leaf_dir) |>
    # if an edited file exists in that folder, keep it; else keep the plain one
    slice_max(is_edit, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(path)
}

# ────────────────────────────────────────────────────────────────────────
#  2-B.  CSV loader for an *explicit vector* of paths
# ────────────────────────────────────────────────────────────────────────
load_peak_vec <- function(vec) vec[file.exists(vec)]

# ────────────────────────────────────────────────────────────────────────
#  3.  Read & wrangle a *set of files* into a “cell-level” tibble
#      (adds a column `intra` = "Cs" or "K")
# ────────────────────────────────────────────────────────────────────────
wrangle_files <- function(file_paths, intra_label) {
  
  stim_regex <- "\\b[A-Z]_\\d+ A (\\d+\\.?\\d*%?|\\d*\\.?\\d+%?) \\d+\\.?\\d* ms\\b"
  
  map_dfr(file_paths, read_csv, show_col_types = FALSE) |>
    clean_names() |>
    mutate(
      ptx       = str_detect(log, regex("PTX", ignore_case = TRUE)),
      cell_id   = paste(
        str_sub(filename,  1, 5),
        str_sub(slice_id, -1, -1),
        str_sub(cell_id,  -1, -1),
        sep = "_"
      ),
      matches   = str_extract_all(log, stim_regex),
      stim1     = map_chr(matches, 1, .default = NA),
      stim2     = map_chr(matches, 2, .default = NA),
      intra     = intra_label,
      
      # keep numeric columns numeric
      across(
        c(interval,
          first_peak_value_combined,
          third_peak_value_combined,
          actual_peak_value,
          estimated_peak_value),
        as.numeric
      )
    ) |>
    filter( interval > 0,
            !is.na(stim1), !is.na(stim2),
            str_detect(stim1, "^R_1") ) |>   #  ← NEW  行
    select(cell_id, ptx, intra, interval, stim1, stim2,
           first_peak_value_combined, third_peak_value_combined,
           actual_peak_value, estimated_peak_value) |>
    group_by(ptx, interval, cell_id, stim1, stim2, intra) |>
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
#  4.  Summarise to “mean table” (± SEM & CI)
# ────────────────────────────────────────────────────────────────────────
make_mean_tbl <- function(cell_tbl) {
  cell_tbl |>
    group_by(intra, ptx, interval) |>
    summarise(
      enhancement_index_mean = mean(enhancement_index),
      enhancement_index_sem  = std.error(enhancement_index),
      n  = n(),
      ci = qt(0.975, df = n - 1) * enhancement_index_sem,
      .groups = "drop"
    )
}

# ────────────────────────────────────────────────────────────────────────
#  5-A.  Plot helper: mean ± SEM (optionally ACSF only)
# ────────────────────────────────────────────────────────────────────────
plot_mean <- function(mean_tbl, S = sty,
                      acsf_only = FALSE,
                      line_cut  = 0.25) {
  
  df <- if (acsf_only) filter(mean_tbl, ptx == FALSE) else mean_tbl
  df_line  <- filter(df, interval <= line_cut)
  df_point <- df
  
  ggplot() +
    geom_point(data = df_point,
               aes(interval, enhancement_index_mean,
                   colour = intra, group = intra),
               size = S$mean_pt_size) +
    geom_line(data = df_line,
              aes(interval, enhancement_index_mean,
                  colour = intra, group = intra),
              linewidth = S$mean_line_lwd) +
    geom_errorbar(data = df_point,
                  aes(interval,
                      ymin = enhancement_index_mean - enhancement_index_sem,
                      ymax = enhancement_index_mean + enhancement_index_sem,
                      colour = intra, group = intra),
                  width = S$err_bar_width,
                  linewidth = S$err_bar_lwd) +
    
    scale_color_manual(values = S$col_intra, guide = "none") +
    
    scale_x_break(breaks = S$x_break, scales = 0.4, space = 0.02) +
    scale_x_continuous(limits = S$x_limits,
                       breaks  = S$x_breaks,
                       labels  = S$x_labels) +
    ylim(S$y_limits) +
    labs(x = "Interval (s)", y = "Enhancement index") +
    theme_classic(base_size = S$base_font_sz, base_family = S$base_font_fam) +
    theme(
      axis.text  = element_text(size   = S$axis_font_sz,
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
#  6.  FILE PATHS
# ────────────────────────────────────────────────────────────────────────
# --- Cs: explicit list ---------------------------------------------------
flist_cs <- c(
  "./Data/PatchClamp/240712(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak.csv",
  "./Data/PatchClamp/240718(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/240719(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/240807(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/240809(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv"
)

# --- K: any file inside this dir (preferring *_edited.csv) --------------
dir_k <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/"

# ────────────────────────────────────────────────────────────────────────
#  7.  IMPORT → wrangle
# ────────────────────────────────────────────────────────────────────────
cell_cs <- wrangle_files(load_peak_vec(flist_cs), "Cs")
cell_k  <- wrangle_files(load_peak_dir(dir_k),  "K")

cell_all <- bind_rows(cell_cs, cell_k)
mean_all <- make_mean_tbl(cell_all)

# ────────────────────────────────────────────────────────────────────────
#  8.  PLOTS
# ────────────────────────────────────────────────────────────────────────
p_mean_all  <- plot_mean(mean_all)                   # PTX+PTX-両方
p_mean_acsf <- plot_mean(mean_all, acsf_only = TRUE) # PTX なしのみ

print(p_mean_acsf)

# ────────────────────────────────────────────────────────────────────────
#  9.  SAVE
# ────────────────────────────────────────────────────────────────────────
ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_Cs_vs_K_all.svg",
       p_mean_all,
       width = sty$fig_w_mm, height = sty$fig_h_mm, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_Cs_vs_K_ACSF.svg",
       p_mean_acsf,
       width = sty$fig_w_mm, height = sty$fig_h_mm, units = "mm")
