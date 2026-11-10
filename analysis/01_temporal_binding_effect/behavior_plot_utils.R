suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gghalves)
  library(readr)
  library(dplyr)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_plot_pdf <- function(p, stem, dir_pdf, w = 8, h = 6) {
  ensure_dir(dir_pdf)

  ggsave(
    filename = file.path(dir_pdf, paste0(stem, ".pdf")),
    plot = p,
    width = w,
    height = h,
    units = "in",
    device = cairo_pdf
  )
}

to_long <- function(df, cols, levels, labels) {

  df %>%
    pivot_longer(
      cols = all_of(cols),
      names_to = "type",
      values_to = "value"
    ) %>%
    mutate(
      type = factor(type, levels = levels, labels = labels)
    )
}

plot_halfviolin <- function(df_long, title,
                            y_breaks = NULL,
                            y_label = "Perceptual Shift (ms)",
                            palette = NULL, fill_palette = palette,
                            point_alpha = .60, fill_alpha = .60, point_size = .65,
                            font_size = 6, font_family = "Arial") {
  df_plot <- df_long %>%
    mutate(
      type_id = as.numeric(type),
      point_x = type_id - 0.16
    )

  p <- ggplot(df_plot, aes(x = type_id, y = value, fill = type, color = type)) +
    geom_half_violin(
      side = "R",
      alpha = fill_alpha,
      linewidth = .22,
      width = 0.6,
      position = position_nudge(x = 0.1, y = 0)
    ) +
    geom_half_boxplot(
      alpha = fill_alpha,
      linewidth = .22,
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.25
    ) +
    geom_point(
      aes(x = point_x),
      alpha = point_alpha,
      position = position_jitter(width = 0.04, height = 0),
      size = point_size
    ) +
    scale_x_continuous(
      breaks = seq_along(levels(df_plot$type)),
      labels = levels(df_plot$type)
    ) +
    scale_fill_manual(values = fill_palette) +
    scale_color_manual(values = palette) +
    labs(
      x = NULL,
      y = y_label,
      title = title
    ) +
    theme_classic(base_size=font_size, base_family=font_family) +
    theme(
      plot.title = element_text(hjust = 0.5, size = font_size+1, face="bold"),
      axis.title.y = element_text(size = font_size+1, face="bold"),
      axis.text.x = element_text(size = font_size, color="black"),
      axis.text.y = element_text(size = font_size, color="black"),
      axis.line = element_line(linewidth=.25),
      axis.ticks = element_line(linewidth=.20),
      legend.position = "none"
    )

  if (!is.null(y_breaks)) {
    p <- p + scale_y_continuous(breaks = y_breaks)
  }

  return(p)
}
