rm(list = ls())
# ライブラリを読み込む
library(tidyverse)
library(plotrix)



# ディレクトリ設定 ----------------------------------------------------------------

# Define the directory where you want to search for the files
setwd('/home/rstudio/')

# ワーキングディレクトリを設定
Path = "./Data/Histology/Retinotopic_quantification/APgradation/"
setwd(Path)


# 関数 ----------------------------------------------------------------------

scale_intensity <- function(csv_data) {
  result <- csv_data %>% 
    group_by(Region, MouseID) %>%
    mutate(
      Mean_Scaled = (Mean - min(Mean, na.rm = TRUE)) / (max(Mean, na.rm = TRUE) - min(Mean, na.rm = TRUE))
    ) %>%
    ungroup()
  
  return(result)
}

# 1) Summarize scaled intensity
summarize_scaled_intensity <- function(data) {
  data_scaled <- scale_intensity(data)  # assumes you already have scale_intensity()
  
  summary_df <- data_scaled %>%
    group_by(Region, degree) %>%
    summarise(
      mean_intensity = mean(Mean_Scaled, na.rm = TRUE),
      SEM            = std.error(Mean_Scaled, na.rm = TRUE),
      n              = n(),
      .groups        = "drop"
    )
  
  return(summary_df)
}


# 2) Plot from that summary
plot_summary_intensity_ACCRSC <- function(summary_df,
                                          regions        = c("RSC", "ACC"),
                                          region_colors  = c(RSC = "#ff4b00", ACC = "#03af7a"),
                                          errorbar_size  = 0.1,
                                          x_step         = NULL,
                                          x_keep         = TRUE,
                                          y_step         = NULL,
                                          y_keep         = FALSE) {
  # 1) filter
  df_plot <- summary_df %>%
    filter(Region %in% regions)
  
  # 2) base plot
  p <- ggplot(df_plot, aes(x      = degree,
                           y      = mean_intensity,
                           fill   = Region,
                           colour = Region)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.3) +
    geom_errorbar(aes(ymax = mean_intensity + SEM,
                      ymin = mean_intensity - SEM),
                  position = position_dodge(0.8),
                  width     = 0.3,
                  size      = errorbar_size) +
    theme_classic()
  
  # 3) optional x-axis label thinning
  if (!is.null(x_step)) {
    lbl_x <- if (x_keep) {
      # keep only every x_step-th label
      function(x) {
        labs <- as.character(x)
        labs[seq_along(labs) %% x_step != 0] <- ""
        labs
      }
    } else {
      # blank every x_step-th label
      function(x) {
        labs <- as.character(x)
        labs[seq_along(labs) %% x_step == 0] <- ""
        labs
      }
    }
    p <- p + scale_x_discrete(labels = lbl_x)
  }
  
  # 4) optional y-axis label thinning
  y_breaks <- seq(0, 1, by = 0.2)
  if (!is.null(y_step)) {
    lbl_y <- if (y_keep) {
      # keep only every y_step-th label
      function(y) {
        labs <- as.character(y)
        labs[seq_along(labs) %% y_step != 0] <- ""
        labs
      }
    } else {
      # blank every y_step-th label
      function(y) {
        labs <- as.character(y)
        labs[seq_along(labs) %% y_step == 0] <- ""
        labs
      }
    }
    p <- p + scale_y_continuous(breaks = y_breaks,
                                labels = lbl_y,
                                expand = expansion(mult = c(0, 0.1)))
  } else {
    p <- p + scale_y_continuous(breaks = y_breaks,
                                expand = expansion(mult = c(0, 0.1)))
  }
  
  # 5) finish styling
  p + scale_fill_manual(values  = region_colors) +
    scale_color_manual(values = region_colors) +
    theme(axis.title               = element_blank(),
          legend.position          = "none",
          axis.ticks.length.x.bottom = unit(0, "mm"),
          axis.ticks.length.y.left   = unit(0.5, "mm"),
          axis.text.y              = element_text(size = 8, family = "Arial"),
          axis.text.x              = element_text(size = 8, family = "Arial"),
          axis.line                = element_line(linewidth = 0.5/2.13),
          axis.ticks               = element_line(linewidth = 0.5/2.13),
          plot.background          = element_rect(fill = "transparent", color = NA),
          panel.background         = element_rect(fill = "transparent"))
}
plot_intensity_LGN <- function(data,
                               x_step     = 5,
                               x_keep     = TRUE,
                               y_step     = 2,
                               y_keep     = FALSE,
                               errorbar_size = 0.1) {
  
  # base ggplot
  p <- ggplot(data, aes(x = degree, y = mean_intensity)) +
    geom_col(position = position_dodge(width = 0.8),
             width    = 0.7,
             fill     = 'black') +
    geom_errorbar(aes(ymax = mean_intensity + SEM,
                      ymin = mean_intensity - SEM),
                  color    = 'black',
                  position = position_dodge(0.8),
                  width    = 0.3,
                  size     = errorbar_size) +
    theme_classic()
  
  # X-axis thinning
  if (!is.null(x_step)) {
    lbl_x <- if (x_keep) {
      # keep only every x_step-th
      function(x) {
        labs <- as.character(x)
        labs[seq_along(labs) %% x_step != 0] <- ""
        labs
      }
    } else {
      # blank every x_step-th
      function(x) {
        labs <- as.character(x)
        labs[seq_along(labs) %% x_step == 0] <- ""
        labs
      }
    }
    p <- p + scale_x_discrete(labels = lbl_x)
  }
  
  # Y-axis thinning (with your original breaks/expand)
  y_breaks <- seq(0, 1, by = 0.2)
  if (!is.null(y_step)) {
    lbl_y <- if (y_keep) {
      # keep only every y_step-th
      function(y) {
        labs <- as.character(y)
        labs[seq_along(labs) %% y_step != 0] <- ""
        labs
      }
    } else {
      # blank every y_step-th
      function(y) {
        labs <- as.character(y)
        labs[seq_along(labs) %% y_step == 0] <- ""
        labs
      }
    }
    p <- p + scale_y_continuous(breaks = y_breaks,
                                labels = lbl_y,
                                expand = expansion(mult = c(0, 0.1)))
  } else {
    p <- p + scale_y_continuous(breaks = y_breaks,
                                expand = expansion(mult = c(0, 0.1)))
  }
  
  # final theme tweaks
  p + theme(
    axis.title               = element_blank(),
    axis.ticks.length.x.bottom = unit(0, "mm"),
    axis.ticks.length.y.left = unit(0.5, "mm"),
    axis.text.y              = element_text(size = 8, family = "Arial"),
    axis.text.x              = element_text(size = 8, family = "Arial"),
    axis.line                = element_line(linewidth = 0.5/2.13),
    axis.ticks               = element_line(linewidth = 0.5/2.13),
    plot.background          = element_rect(fill = "transparent", color = NA),
    panel.background         = element_rect(fill = "transparent")
  )
}


# csv読み込み -----------------------------------------------------------------

# List all CSV files in the directory
Layer1_intensity <- read.csv("./APaxonintensity/Layer1_intensity.csv") %>% 
  select(-Min, -Max, -MinThr, -MaxThr) %>% 
  mutate(ROIID = as.numeric(str_match(ROIID, "grid[^:]*\\.tif:(\\d+):\\d")[, 2]))

ROIid_to_degree_altitude <- read.csv("./unique_roi_ids_altitude.csv")
ROIid_to_degree_azimuth <- read.csv("./unique_roi_ids_azimuth.csv")

Altitude_data <- Layer1_intensity %>% 
  filter(Axis == "altitude") %>% 
  left_join(
    ROIid_to_degree_altitude,
    by = c("ROIID" = "Unique_ROIID") # Specify joining columns and their mapping
  ) %>% 
  mutate(degree = as.factor(degree))
  
Azimuth_data <- Layer1_intensity %>% 
  filter(Axis == "azimuth") %>% 
  left_join(
    ROIid_to_degree_azimuth,
    by = c("ROIID" = "Unique_ROIID") # Specify joining columns and their mapping
  ) %>% 
  mutate(degree = as.factor(degree))



# Min-Max normalizationとplotting ------------------------------------------

# 1) サマリー作成
sum_alt <- summarize_scaled_intensity(Altitude_data)
sum_alt_LGN <- filter(sum_alt, Region == "LGN")
sum_azm <- summarize_scaled_intensity(Azimuth_data)
sum_azm_LGN <- filter(sum_azm, Region == "LGN")

# 2) プロット作成

# ■ 毎5つ目だけ表示（5以外は空白）
p_azm_ACCRSC <- plot_summary_intensity_ACCRSC(sum_azm,
                                    x_step = 5,
                                    x_keep = TRUE,
                                    y_step = 2,
                                    y_keep = FALSE)

# ■ 2の倍数だけ空白にする（1,3,5…は表示）
p_alt_ACCRSC <- plot_summary_intensity_ACCRSC(sum_alt,
                                    x_step = 2,
                                    x_keep = FALSE,
                                    y_step = 2,
                                    y_keep = FALSE)

p_azm_LGN <- # show *only* every 3rd x-label, blank all others:
  plot_intensity_LGN(sum_azm_LGN, x_step = 5, x_keep = TRUE)
p_alt_LGN <- plot_intensity_LGN(sum_alt_LGN, x_step = 2, x_keep = FALSE)


plot(p_alt_ACCRSC)
plot(p_alt_LGN)
plot(p_azm_ACCRSC)
plot(p_azm_LGN)


# 3) プロット表示＆保存
ggsave("./plots/ACCvsRSC_altitude.svg", p_alt_ACCRSC,width = 35, height = 32, units = "mm")
ggsave("./plots/LGNlayer1_altitude.svg", p_alt_LGN,width = 35, height = 32, units = "mm")
ggsave("./plots/ACCvsRSC_azimuth.svg", p_azm_ACCRSC, width = 35, height = 32, units = "mm")
ggsave("./plots/LGNlayer1_azimuth.svg", p_azm_LGN, width = 35, height = 32, units = "mm")

