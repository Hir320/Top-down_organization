rm(list = ls())

# ───────────────────────────────────────────────────────────────
# 0.  Packages
# ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(plotrix)   # std.error()
})

# ───────────────────────────────────────────────────────────────
# 1.  GLOBAL STYLE
# ───────────────────────────────────────────────────────────────
sty <- list(
  # axis breaks
  x_break        = c(0.26, 0.75),
  x_limits       = c(0, 0.85),
  x_breaks       = c(0.025, 0.05, 0.1, 0.25, 0.8),
  x_labels       = c("25", "50", "100", "250", "800"),
  
  # y-axis
  y_limits       = c(0, 2.1),
  
  # glyphs & lines
  cell_pt_size   = 0.5,
  cell_pt_alpha  = 0.10,
  mean_pt_size   = 1.5,
  mean_line_lwd  = 1 / 2.13,
  err_bar_width  = 0.01,
  err_bar_lwd    = 1/2.13,
  
  axis_line_lwd  = 1 / 2.13,          # Illustrator 0.5 pt
  tick_line_lwd  = 1 / 2.13,
  
  # intervals
  x_order      = c("-75 mV K", "-55 mV K", "-75 mV Cs"),
  x_labels     = c("-75 mV\nK-glu", "-55 mV\nK-glu", "-75 mV\nCs"),
  point_size   = 1.8,
  err_bar_w    = 0.15,
  err_bar_lwd  = 0.4,
  
  # colours
  col_cond     = c(`-75 mV K`  = "black",
                   `-55 mV K`  = "grey40",
                   `-75 mV Cs` = "red"),
  
  # fonts & export
  base_sz      = 8,
  base_fam     = "Arial",
  fig_w_mm     = 100,
  fig_h_mm     = 45
)

theme_manual <- theme(
  strip.background = element_blank(),
  strip.text = element_text(size = sty$base_sz, family = sty$base_fam),
  axis.text.x  = element_blank(),
  axis.text.y  = element_text(size = sty$base_sz, family = sty$base_fam),
  axis.title = element_blank(),
  axis.line  = element_line(linewidth = sty$axis_line_lwd),
  axis.ticks = element_line(linewidth = sty$tick_line_lwd),
  axis.ticks.length.x.bottom = unit(0, "mm"),
  axis.ticks.length.y.left   = unit(1, "mm"),
  plot.background  = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent")
)

# ───────────────────────────────────────────────────────────────
# 2.  Small helpers  (edited)
#     — keep **one** file per leaf folder:
#         • merged_mean_peak_edited.csv  (if it exists)
#         • otherwise merged_mean_peak.csv
# ───────────────────────────────────────────────────────────────

## 2-A. return a character vector of paths  -------------------------------
load_peak_dir <- function(dir_path) {
  list.files(
    path       = dir_path,
    pattern    = "merged_mean_peak(_edited)?\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  ) |>
    tibble(path = _) |>
    mutate(
      leaf_dir = dirname(path),                  # “…/240419(… )”
      is_edit  = str_detect(path, "_edited\\.csv$")
    ) |>
    group_by(leaf_dir) |>
    slice_max(is_edit, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(path)
}

## 2-B. read & clean *all* files inside one directory ---------------------
load_peak_files <- function(dir_path) {
  load_peak_dir(dir_path) |>
    map_dfr(read_csv, show_col_types = FALSE) |>
    clean_names()        # spaces → snake_case
}

## 2-C. explicit-vector loader (unchanged) --------------------------------
load_peak_vec <- function(v) v[file.exists(v)]

wrangle_dir <- function(dir_path, vhold_label) {
  
  regex_stim <- "\\b[A-Z]_\\d+ A (\\d+\\.?\\d*%?|\\d*\\.?\\d+%?) \\d+\\.?\\d* ms\\b"
  
  load_peak_files(dir_path) |>
    mutate(
      ptx       = str_detect(log,   regex("PTX", ignore_case = TRUE)),
      ap4       = str_detect(log, regex("4\\s*-?\\s*AP", ignore_case = TRUE)),
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
           !ptx,
           !ap4,
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

wrangle_files <- function(file_paths, intra_label) {
  
  stim_regex <- "\\b[A-Z]_\\d+ A (\\d+\\.?\\d*%?|\\d*\\.?\\d+%?) \\d+\\.?\\d* ms\\b"
  
  map_dfr(file_paths, read_csv, show_col_types = FALSE) |>
    clean_names() |>
    mutate(
      ptx       = str_detect(log, regex("PTX", ignore_case = TRUE)),
      ap4       = str_detect(log, regex("4\\s*-?\\s*AP", ignore_case = TRUE)),
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
            !ptx,
            !ap4,
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




summary_tbl <- function(cell_tbl) {
  cell_tbl |>
    group_by(interval, cond) |>
    summarise(
      enhancement_index_mean = mean(enhancement_index, na.rm = TRUE),
      enhancement_index_sem  = std.error(enhancement_index),
      n                      = n(),
      ci                     = qt(0.975, df = n - 1) * enhancement_index_sem,
      .groups = "drop"
    ) |>
    mutate(cond = factor(cond, levels = sty$x_order))
}

# ───────────────────────────────────────────────────────────────
# 3.  FILE LISTS
# ───────────────────────────────────────────────────────────────
## 3-A.  Cs (assumed -75 mV)
files_cs <- load_peak_vec(c(
  "./Data/PatchClamp/240712(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak.csv",
  "./Data/PatchClamp/240718(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/240719(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/240807(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv",
  "./Data/PatchClamp/240809(ACC ChR2+ RSC ChrimsonR)/merged_mean_peak_edited.csv"
))

## 3-B.  K, two holding potentials
dir_k75 <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/"
dir_k55 <- "./Data/PatchClamp/Integration/Integration_VHold-55mV/RSC-R_ACC-2/"


# ───────────────────────────────────────────────────────────────
# 4.  IMPORT → COMBINE → SUMMARISE   (REPLACED)
# ───────────────────────────────────────────────────────────────

## 4-A.  -75 mV  Cs
cell_cs75 <- wrangle_files(files_cs, "Cs") |>
  filter(ptx == FALSE) |>
  mutate(vhold = "-75 mV",
         cond  = paste(vhold, intra))

## 4-B.  -75 mV  K
cell_k75  <- wrangle_dir(dir_k75, "-75 mV") |>
  filter(ptx == FALSE) |>
  mutate(intra = "K",
         cond  = paste(vhold, intra))

## 4-C.  -55 mV  K
cell_k55  <- wrangle_dir(dir_k55, "-55 mV") |>
  filter(ptx == FALSE) |>
  mutate(intra = "K",
         cond  = paste(vhold, intra))

## 4-D.  combine
cell_all  <- bind_rows(cell_cs75, cell_k75, cell_k55)

## 4-E.  summarise  (mean ± SEM)
mean_all  <- summary_tbl(cell_all)

# ───────────────────────────────────────────────────────────────
# 5.  PLOT  (mean  ±  SEM   +  every cell)
# ───────────────────────────────────────────────────────────────

plot_mean_cells <- function(cell_tbl, mean_tbl, S = sty) {
  
  cells <- cell_tbl |>                       # every sweep
    filter(ptx == FALSE) |>           # ACSF only
    mutate(cond = factor(cond, levels = S$x_order))
  
  means <- mean_tbl                          # already filtered via summary_tbl()
  
  ggplot() +
    ## per-cell points (faint) -------------------------------------------
  geom_jitter(data = cells,
              aes(cond, enhancement_index,
                  colour = cond),
              size   = 0.5,
              alpha  = 0.12,
              width  = 0.15, height = 0) +
    
    ## mean ± SEM ---------------------------------------------------------
  geom_point(data = means,
             aes(cond, enhancement_index_mean,
                 colour = cond),
             size = S$point_size) +
    geom_errorbar(data = means,
                  aes(cond,
                      ymin = enhancement_index_mean - enhancement_index_sem,
                      ymax = enhancement_index_mean + enhancement_index_sem,
                      colour = cond),
                  width     = S$err_bar_w,
                  linewidth = S$err_bar_lwd) +
    
    ## colour & guides ----------------------------------------------------
  scale_colour_manual(values = S$col_cond, guide = "none") +
    
    ## panels -------------------------------------------------------------
  facet_wrap(~ interval,
             labeller = labeller(interval = function(x)
               paste0(as.numeric(as.character(x)) * 1000, " ms")),
             nrow = 1) +
    
    ## axes, theme --------------------------------------------------------
  ylim(S$y_limits) +
    labs(x = NULL, y = "Enhancement index") +
    theme_classic(base_size = S$base_sz, base_family = S$base_fam) +
    theme_manual
}

plot_mean_cells_ci <- function(cell_tbl, mean_tbl, S = sty) {
  
  cells <- filter(cell_tbl, ptx == FALSE) |>         # ACSF only
    mutate(cond = factor(cond, levels = S$x_order))
  
  means <- mean_tbl                                  # already has ci column
  
  ggplot() +
    geom_jitter(data = cells,
                aes(cond, enhancement_index, colour = cond),
                size  = 0.5, alpha = 0.12,
                width = 0.15, height = 0) +
    
    geom_point(data = means,
               aes(cond, enhancement_index_mean, colour = cond),
               size = S$point_size) +
    geom_errorbar(data = means,
                  aes(cond,
                      ymin = enhancement_index_mean - ci,
                      ymax = enhancement_index_mean + ci,
                      colour = cond),
                  width     = S$err_bar_w,
                  linewidth = S$err_bar_lwd) +
    
    scale_colour_manual(values = S$col_cond, guide = "none") +
    facet_wrap(~ interval,
               labeller = labeller(interval = function(x)
                 paste0(as.numeric(as.character(x)) * 1000, " ms")),
               nrow = 1) +
    ylim(S$y_limits) +
    labs(x = NULL, y = "Enhancement index") +
    theme_classic(base_size = S$base_sz, base_family = S$base_fam) +
    theme_manual
}


# build & show ------------------------------------------------------------
p <- plot_mean_cells(cell_all, mean_all)
print(p)

# --- build & show CI version --------------------------------------------
p_ci <- plot_mean_cells_ci(cell_all, mean_all)+
  scale_y_continuous(limits = c(0, 2.1),
                     breaks  = seq(0, 2, 0.5),
                     expand  = expansion(mult = c(0, 0.02))
                     )
print(p_ci) 



# ───────────────────────────────────────────────────────────────
# 6.  SAVE
# ───────────────────────────────────────────────────────────────
ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_K_vs_Cs_Hold.svg",
       p,
       width  = sty$fig_w_mm,
       height = sty$fig_h_mm,
       units  = "mm")
# --- save ---------------------------------------------------------------
ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_K_vs_Cs_Hold_CI.svg",
       p_ci,
       width  = sty$fig_w_mm,
       height = sty$fig_h_mm,
       units  = "mm")
