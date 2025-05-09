rm(list =ls())
# =======================================================================
#  ACC-R → RSC-2   and   RSC-R → ACC-2
#  Interval-series panels (ACSF only)
#  – plotting separated from saving –
# =======================================================================

# ──────────────────────────────────────────────────────────────────────
#  0.  Centralised plot-style parameters
# ──────────────────────────────────────────────────────────────────────
plt <- list(
  # geometry
  x_break       = c(0.26, 0.75),   # ggbreak() lower/upper bound
  x_limits      = c(0, 0.85),
  x_breaks      = c(0.025, 0.05, 0.1, 0.25, 0.8),
  x_labels      = c("25", "50", "100", "250", "800"),
  y_limits      = c(0.2,1.8),
  
  # point layer
  cell_pt_size  = 0.5,
  cell_pt_alpha = 0.10,
  mean_pt_size  = 1.5,
  
  # fonts & lines
  axis_font_sz  = 8,               # pt
  axis_font_fam = "Arial",
  axis_line_lwd = 1/2.13,        # Illustrator 0.5 pt ≈ 0.235 mm
  tick_line_lwd = 1/2.13,
  
  # output
  fig_w_mm      = 70,
  fig_h_mm      = 55
)


# ---- 1.  Packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)     # readr, dplyr, ggplot2, etc.
  library(plotrix)       # std.error()
  library(ggbreak)
  library(broom)
  library(emmeans)
  library(effectsize)   # Cohen’s d, Hedges g, r, Kendall’s W, etc.
  library(lmerTest)
})

# ---- 2.  Project root --------------------------------------------------
setwd("/home/rstudio/")            # <- change if necessary

# ---- 3.  Helper functions ----------------------------------------------
prefer_edited <- function(paths) {
  edited <- sub("merged_mean_peak.csv", "merged_mean_peak_edited.csv", paths)
  ifelse(file.exists(edited), edited, paths)
}

# ------------------------------------------------------------------
# wrangle_integration(): 先ほどの read_and_clean() をもう一度修正
# ------------------------------------------------------------------
wrangle_integration <- function(dir_path, extra_paths = character(0)) {
  
  files <- list.files(dir_path,
                      pattern   = "merged_mean_peak.csv",
                      recursive = TRUE,
                      full.names = TRUE) |>
    prefer_edited()
  print(files)
  files <- c(files, extra_paths[file.exists(extra_paths)])
  
  stim_regex <- "\\b[A-Z]_\\d+ A (\\d+\\.?\\d*%?|\\d*\\.?\\d+%?) \\d+\\.?\\d* ms\\b"
  
  # ---- 読み込み -------------------------------------------------------
  read_and_clean <- function(fp) {
    # ① すべて文字列で読み込む（型衝突を根絶）
    df <- readr::read_csv(
      fp,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    )
    
    # ② 無名列 “…1” “…2”… を捨て、
    #    解析に必要な最小限の列だけ残す
    df |>
      select(
        -matches("^\\.\\.\\."),                # 無名列全部除去
        any_of(c("filename", "SliceID", "CellID",
                 "log", "Interval",
                 "Actual Peak Value", "Estimated Peak Value"))
      )
  }
  
  raw <- purrr::map_dfr(files, read_and_clean)
  
  # ---- 型をここで整える ----------------------------------------------
  raw <- raw |>
    mutate(
      Interval            = as.numeric(Interval),
      `Actual Peak Value` = as.numeric(`Actual Peak Value`),
      `Estimated Peak Value` = as.numeric(`Estimated Peak Value`)
    )
  
  # ---- 以降は元の処理フロー -------------------------------------------
  raw |>
    mutate(
      PTX = str_detect(log, regex("PTX", ignore_case = TRUE)),
      ## --- 4-AP の有無を検出 ------------------------------------------
      AP4 = str_detect(log, regex("4\\s*-?\\s*AP", ignore_case = TRUE)),
      ## ----------------------------------------------------------------
      Cell_ID = paste(str_sub(filename, 1, 5),
                      str_sub(SliceID, -1, -1),
                      str_sub(CellID,  -1, -1),
                      sep = "_"),
      Matches = str_extract_all(log, stim_regex),
      Stim1   = map_chr(Matches, 1, .default = NA),
      Stim2   = map_chr(Matches, 2, .default = NA)
    ) |>
    filter(
      Interval > 0,
      str_detect(log, "-75 mV"),
      !PTX,          # ← PTX を除外
      !AP4,          # ← 4-AP も除外   ★ 追加 ★
      !is.na(Stim1), !is.na(Stim2)
    ) |>
    select(
      Cell_ID,
      PTX,            # ← 残す
      AP4,            # ← 残す（後で使うなら）
      Interval,
      `Actual Peak Value`,
      `Estimated Peak Value`
    ) |>
    ##  group_by() に PTX を追加 ---------------
    group_by(PTX, Interval, Cell_ID) |>
    summarise(
      .groups = "drop",
      Actual_Peak    = mean(`Actual Peak Value`,    na.rm = TRUE),
      Estimated_Peak = mean(`Estimated Peak Value`, na.rm = TRUE)
    ) |>
    mutate(Enhancement_index = Actual_Peak / Estimated_Peak)
}



make_mean_table <- function(cell_tbl) {
  cell_tbl |>
    group_by(PTX, Interval) |>        # ← PTX を戻す
    summarise(.groups = "drop",
              Enhancement_index_mean = mean(Enhancement_index),
              Enhancement_index_sem  = std.error(Enhancement_index),
              n = n()) |>
    mutate(CI = qt(0.975, df = n - 1) * Enhancement_index_sem)
}


# ---- 3-B.  Plot-only function ------------------------------------------
# ─────────────────────────────────────────────────────────────────────────
#  NEW helper: mean + per-cell points  (only PTX == FALSE, intervals wanted)
# ─────────────────────────────────────────────────────────────────────────
make_interval_plot_cells <- function(cell_tbl, mean_tbl, S = plt, jitter = FALSE) {
  
  cells <- cell_tbl |>
    filter(PTX == FALSE,
           Interval %in% S$x_breaks)
  
  means <- mean_tbl |> filter(PTX == FALSE)
  
  ggplot() +
    {if (jitter)
      geom_jitter(data = cells,
                  aes(Interval, Enhancement_index),
                  size  = S$cell_pt_size,
                  alpha = S$cell_pt_alpha,
                  width = 0, height = 0)   # zero width keeps dots aligned
      else
        geom_point (data = cells,
                    aes(Interval, Enhancement_index),
                    size  = S$cell_pt_size,
                    alpha = S$cell_pt_alpha)} +
    geom_point(data = means,
               aes(Interval, Enhancement_index_mean),
               size = S$mean_pt_size) +
    geom_errorbar(data = means,
                  aes(Interval,
                      ymin = Enhancement_index_mean - CI,
                      ymax = Enhancement_index_mean + CI),
                  width = 0.01,
                  linewidth = 1/2.13) +
    ylim(S$y_limits) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = S$axis_font_sz,
                                family = S$axis_font_fam),
      axis.title = element_blank(),
      axis.line  = element_line(linewidth = S$axis_line_lwd),
      axis.ticks = element_line(linewidth = S$tick_line_lwd),
      axis.line.x.top   = element_blank(),
      axis.ticks.x.top  = element_blank(),
      axis.text.x.top   = element_blank(),
      axis.title.x.top  = element_blank(),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent")
    )
}



# ---- 4.  Data import & processing --------------------------------------
## ACC-R → RSC-2
acc_dir <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-2_ACC-R/"
ACC_R_cell <- wrangle_integration(acc_dir)
ACC_R_mean <- make_mean_table(ACC_R_cell)

## RSC-R → ACC-2
rsc_dir <- "./Data/PatchClamp/Integration/Integration_VHold-75mV/RSC-R_ACC-2/"
manual_path <- file.path(rsc_dir,
                         "240817(ACC ChR2+ RSC ChrimsonR)",
                         "merged_mean_peak_edited.csv")
RSC_R_cell <- wrangle_integration(rsc_dir, extra_paths = manual_path)
RSC_R_mean <- make_mean_table(RSC_R_cell)

# ---- 5.  Create plots (objects only) -----------------------------------
p_acc <- make_interval_plot_cells(ACC_R_cell, ACC_R_mean) +
  scale_x_break(breaks = plt$x_break, scales = 0.4, space = 0.02) +
  scale_x_continuous(limits = plt$x_limits,
                     breaks  = plt$x_breaks,
                     labels  = plt$x_labels)

p_rsc <- make_interval_plot_cells(RSC_R_cell, RSC_R_mean) +
  scale_x_break(breaks = plt$x_break, scales = 0.4, space = 0.02) +
  scale_x_continuous(limits = plt$x_limits,
                     breaks  = plt$x_breaks,
                     labels  = plt$x_labels)

print(p_acc)
print(p_rsc)

# ---- 6.  SAVE to disk (separate step) ----------------------------------
ggsave("./Data/PatchClamp/graphs/TwoColour/ACC-R_RSC-2_interval_panel_cells.svg",
       p_acc, width = plt$fig_w_mm, height = plt$fig_h_mm, units = "mm")

ggsave("./Data/PatchClamp/graphs/TwoColour/RSC-R_ACC-2_interval_panel_cells.svg",
       p_rsc, width = plt$fig_w_mm, height = plt$fig_h_mm, units = "mm")


# ──────────────────────────────────────────────────────────────
#  interval_summary(): mean ± CI  +  effect size ± CI
# ──────────────────────────────────────────────────────────────
interval_summary <- function(cell_tbl, nonparam = FALSE, mu = 1) {
  cell_tbl %>% 
    filter(PTX == FALSE) %>% 
    group_by(Interval) %>% 
    summarise(
      vec = list(Enhancement_index),   # ← 生ベクトルを保持
      n   = n(),
      .groups = "drop"
    ) %>% 
    mutate(
      ## ----- 平均と信頼区間 & p 値  -----------------------------------
      test_obj = map(vec, \(x)
                     if (nonparam)
                       wilcox.test(x, mu = mu, conf.int = TRUE, exact = FALSE)
                     else
                       t.test(x, mu = mu, conf.int = TRUE)
      ),
      mean          = map_dbl(vec, mean),
      mean_CI_low   = map_dbl(test_obj, \(tt) tt$conf.int[1]),
      mean_CI_high  = map_dbl(test_obj, \(tt) tt$conf.int[2]),
      p.value       = map_dbl(test_obj, \(tt) tt$p.value),
      ## ----- 効果量  ---------------------------------------------------
      es_tbl = map(vec, \(x) {
        if (length(x) < 2) {
          tibble(effect = NA_real_, CI_low = NA_real_, CI_high = NA_real_)
        } else if (nonparam) {
          rb <- effectsize::rank_biserial(x, mu = mu, ci = 0.95)
          tibble(effect = rb$r,
                 CI_low = rb$CI_low,
                 CI_high = rb$CI_high)
        } else {
          hg <- effectsize::hedges_g(x, mu = mu, ci = 0.95)
          tibble(effect = hg$Hedges_g,
                 CI_low = hg$CI_low,
                 CI_high = hg$CI_high)
        }
      })
    ) %>% 
    select(-vec, -test_obj) %>% 
    tidyr::unnest(es_tbl, names_sep = "_") %>%  # effect, CI_low/high が展開される
    rename(
      eff_CI_low  = es_tbl_CI_low,
      eff_CI_high = es_tbl_CI_high,
      effect      = es_tbl_effect
    )
}


## ---- VS 1 (parametric = Hedges g) -------------------------------------
acc_summary   <- interval_summary(ACC_R_cell, nonparam = FALSE)
rsc_summary   <- interval_summary(RSC_R_cell, nonparam = FALSE)

## ---- VS 1 (non-param = rank-biserial r) -------------------------------
acc_summary_np <- interval_summary(ACC_R_cell, nonparam = TRUE)
rsc_summary_np <- interval_summary(RSC_R_cell, nonparam = TRUE)

print(acc_summary)
print(rsc_summary)

print(acc_summary_np)
print(rsc_summary_np)

# ──────────────────────────────────────────────────────────────
#  across_interval_summary(): Interval main-effect + pairwise
#                             p & effect size (param/Friedman)
# ──────────────────────────────────────────────────────────────
across_interval_summary <- function(cell_tbl,
                                    method = c("parametric", "friedman")) {
  
  method <- match.arg(method)
  
  df <- cell_tbl %>% 
    filter(PTX == FALSE) %>% 
    mutate(Interval_f = factor(Interval))
  
  if (method == "parametric") {
    
    ## ---------- LMM main effect ---------------------------------------
    fit  <- lmerTest::lmer(Enhancement_index ~ Interval_f + (1|Cell_ID), data = df)
    main <- anova(fit)
    
    ## ---------- LS-means pairwise ------------------------------------
    pw0 <- emmeans(fit, "Interval_f") %>%
      contrast("pairwise") %>%
      summary(infer = TRUE, adjust = "none") %>%     # raw p
      as_tibble()
    
    ## ---------- effect size (Hedges g_av) ----------------------------
    pw  <- pw0 %>%
      mutate(
        # convert t & df to Hedges g + CI
        g_tbl = map2(t.ratio, df, \(tval, dfr)
                     effectsize::t_to_d(tval, df = dfr,
                                        paired = TRUE, hedges = TRUE,
                                        ci = 0.95)
        ),
        effect      = map_dbl(g_tbl, \(x) x$d),
        eff_CI_low  = map_dbl(g_tbl, \(x) x$CI_low),
        eff_CI_high = map_dbl(g_tbl, \(x) x$CI_high),
        p.adj       = p.adjust(p.value, "BH")
      ) %>%
      select(contrast, estimate, p.raw = p.value, p.adj,
             effect, eff_CI_low, eff_CI_high)
    
    return(list(main = main, pairwise = pw))
    
  } else {       # ---------------------- Friedman / Wilcoxon ----------
    df_wide <- df %>%
      select(Cell_ID, Interval_f, Enhancement_index) %>%
      tidyr::pivot_wider(names_from = Interval_f,
                         values_from = Enhancement_index) %>%
      drop_na()
    
    long <- df_wide %>%
      tidyr::pivot_longer(-Cell_ID,
                          names_to = "Interval_f",
                          values_to = "Enhancement_index") %>%
      mutate(Interval_f = factor(Interval_f,
                                 levels = levels(df$Interval_f)))
    
    fri <- with(long,
                friedman.test(Enhancement_index ~ Interval_f | Cell_ID))
    W   <- effectsize::kendalls_w(long$Enhancement_index,
                                 long$Interval_f, long$Cell_ID)
    
    ## pairwise Wilcoxon
    pw_mat <- pairwise.wilcox.test(long$Enhancement_index,
                                   long$Interval_f,
                                   paired = TRUE,
                                   p.adjust.method = "none")$p.value
    
    pw <- as.data.frame(pw_mat) %>%
      rownames_to_column("i1") %>%
      pivot_longer(-i1, names_to = "i2", values_to = "p.raw") %>%
      filter(!is.na(p.raw)) %>%
      mutate(
        contrast = paste(i1, i2, sep = " – "),
        # rank-biserial r for each pair
        rb_tbl = pmap(list(i1, i2), \(a, b) {
          x <- long %>% filter(Interval_f == a) %>% pull(Enhancement_index)
          y <- long %>% filter(Interval_f == b) %>% pull(Enhancement_index)
          effectsize::rank_biserial(x, y, paired = TRUE, ci = 0.95)
        }),
        effect      = map_dbl(rb_tbl, \(x) x$r),
        eff_CI_low  = map_dbl(rb_tbl, \(x) x$CI_low),
        eff_CI_high = map_dbl(rb_tbl, \(x) x$CI_high),
        p.adj       = p.adjust(p.raw, "BH")
      ) %>%
      select(contrast, p.raw, p.adj,
             effect, eff_CI_low, eff_CI_high)
    
    return(list(main = fri, W = W, pairwise = pw))
  }
}


acc_int_res <- across_interval_summary(ACC_R_cell, method = "parametric")
rsc_int_res <- across_interval_summary(RSC_R_cell, method = "parametric")

print(acc_int_res$main)      # F, df, p
print(acc_int_res$pairwise)  # Δ, p, g_av ± CI

print(rsc_int_res$main)      # F, df, p
print(rsc_int_res$pairwise)  # Δ, p, g_av ± CI

# ノンパラ版
acc_int_np <- across_interval_summary(ACC_R_cell, method = "friedman")
print(acc_int_np$main)       # Friedman χ²
print(acc_int_np$W)          # Kendall's W (effect size)
print(acc_int_np$pairwise)   # Wilcoxon p + r ± CI
