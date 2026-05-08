# analysis/01_temporal_binding_effect/temporal_binding_effect_utils.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gghalves)
  library(readr)
  library(dplyr)
})

# ---------- I/O ----------
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

# ---------- data reshape ----------
to_long <- function(df, cols, levels, labels) {
  stopifnot(
    length(cols) == length(levels),
    length(levels) == length(labels)
  )
  
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

# ---------- plotting ----------
plot_halfviolin <- function(df_long, title,
                            y_breaks = NULL,
                            y_label = "Perceptual Shift (ms)",
                            palette = NULL) {
  
  p <- ggplot(df_long, aes(x = type, y = value, fill = type, color = type)) +
    geom_half_violin(
      side = "R",
      alpha = 0.6,
      width = 0.6,
      position = position_nudge(x = 0.1, y = 0)
    ) +
    geom_half_boxplot(
      alpha = 0.6,
      width = 0.6
    ) +
    geom_half_point(
      alpha = 0.6,
      side = "L",
      position = position_nudge(x = -0.1, y = 0)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(
      x = NULL,
      y = y_label,
      title = title
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(2)),
      axis.title.y = element_text(size = rel(2)),
      axis.text.x = element_text(size = rel(2)),
      axis.text.y = element_text(size = rel(2)),
      legend.position = "none"
    )
  
  if (!is.null(y_breaks)) {
    p <- p + scale_y_continuous(breaks = y_breaks)
  }
  
  return(p)
}

# ---------- stats ----------
one_sample_tests <- function(x, label, mu = 0) {
  tt <- t.test(x, mu = mu)
  wx <- wilcox.test(x, mu = mu, exact = FALSE)
  
  tibble(
    test = label,
    n = length(x),
    t_stat = unname(tt$statistic),
    t_p = tt$p.value,
    w_stat = unname(wx$statistic),
    w_p = wx$p.value
  )
}