# pearson_utils.R
# Plotting utilities for iEEG behavioural response-time panels

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
label_pretty <- function(x) {
  gsub("chan", "", x)
}

# ------------------------------------------------------------
# Plot function: response-time distribution by subject/contact
# - jittered trial points
# - median horizontal line
# - IQR vertical line
# - different shapes for different contacts
# ------------------------------------------------------------
plot_half_violin_subjects <- function(
    behavior_df,
    feature_df = NULL,
    condition_pick,
    map_df,
    subjects_order = NULL,
    fill_vals = NULL,
    point_shapes = NULL,
    y_limits = c(-200, 600),
    y_breaks = c(-200, 0, 200, 400, 600),
    
    # Font sizes
    base_size = 8,
    x_text_size = 10,
    y_text_size = 9,
    y_title_size = 10.5,
    
    # Points
    point_size = 2.1,
    point_alpha = 0.72,
    point_nudge = -0.06,
    jitter_width = 0.075,
    point_stroke = 0.20,
    
    # Summary: median + IQR
    summary_nudge = 0.13,
    summary_linewidth = 0.75,
    median_linewidth = 1.35,
    median_width = 0.30,
    
    # Zero line
    zero_line_color = "grey72",
    zero_line_linetype = "dashed",
    zero_line_linewidth = 0.45
) {
  
  stopifnot(all(c("subject", "condition", "behavior_value") %in% names(behavior_df)))
  stopifnot(all(c("subject", "label") %in% names(map_df)))
  
  map_df2 <- map_df %>%
    mutate(
      subject = as.character(subject),
      label = as.character(label)
    )
  
  if (is.null(subjects_order)) {
    subjects_order <- unique(map_df2$subject)
  } else {
    subjects_order <- as.character(subjects_order)
  }
  
  label_df <- map_df2 %>%
    filter(subject %in% subjects_order) %>%
    mutate(
      subject = factor(subject, levels = subjects_order),
      x_label = label_pretty(label)
    ) %>%
    arrange(subject)
  
  x_levels <- label_df$x_label
  
  plot_data <- behavior_df %>%
    mutate(
      subject = as.character(subject),
      condition = as.character(condition)
    ) %>%
    filter(condition == condition_pick) %>%
    inner_join(map_df2, by = "subject") %>%
    mutate(
      x_label = label_pretty(label),
      x_label = factor(x_label, levels = x_levels),
      x_num = as.numeric(x_label),
      x_pt = x_num + point_nudge
    ) %>%
    filter(!is.na(x_label), is.finite(behavior_value))
  
  if (nrow(plot_data) == 0) {
    stop("No behavior data after filtering. Check condition_pick, map_df, and subjects_order.")
  }
  
  if (is.null(fill_vals)) {
    fill_vals <- c("#02c39a", "#ff758f")
    fill_vals <- rep(fill_vals, length.out = length(x_levels))
    names(fill_vals) <- x_levels
  } else {
    if (is.null(names(fill_vals))) {
      names(fill_vals) <- x_levels
    }
    
    missing_cols <- setdiff(x_levels, names(fill_vals))
    if (length(missing_cols) > 0) {
      stop("fill_vals missing names for: ", paste(missing_cols, collapse = ", "))
    }
    
    fill_vals <- fill_vals[x_levels]
  }
  
  # Use filled shapes by default, so a white edge can be added.
  if (is.null(point_shapes)) {
    point_shapes <- c(21, 24)
    point_shapes <- rep(point_shapes, length.out = length(x_levels))
    names(point_shapes) <- x_levels
  } else {
    if (is.null(names(point_shapes))) {
      names(point_shapes) <- x_levels
    }
    
    missing_shapes <- setdiff(x_levels, names(point_shapes))
    if (length(missing_shapes) > 0) {
      stop("point_shapes missing names for: ", paste(missing_shapes, collapse = ", "))
    }
    
    point_shapes <- point_shapes[x_levels]
  }
  
  sum_df <- plot_data %>%
    group_by(x_label, x_num) %>%
    summarise(
      y_med = median(behavior_value, na.rm = TRUE),
      y_q1  = quantile(behavior_value, 0.25, na.rm = TRUE),
      y_q3  = quantile(behavior_value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_sum = x_num + summary_nudge
    )
  
  ggplot() +
    geom_segment(
      aes(
        x = 0.58,
        xend = length(x_levels) + 0.42,
        y = 0,
        yend = 0
      ),
      color = zero_line_color,
      linetype = zero_line_linetype,
      linewidth = zero_line_linewidth
    ) +
    geom_point(
      data = plot_data,
      aes(
        x = x_pt,
        y = behavior_value,
        fill = x_label,
        shape = x_label
      ),
      color = "white",
      stroke = point_stroke,
      alpha = point_alpha,
      size = point_size,
      position = position_jitter(width = jitter_width, height = 0)
    ) +
    geom_linerange(
      data = sum_df,
      aes(
        x = x_sum,
        ymin = y_q1,
        ymax = y_q3,
        color = x_label
      ),
      linewidth = summary_linewidth
    ) +
    geom_segment(
      data = sum_df,
      aes(
        x = x_sum - median_width / 2,
        xend = x_sum + median_width / 2,
        y = y_med,
        yend = y_med,
        color = x_label
      ),
      linewidth = median_linewidth,
      lineend = "round"
    ) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_color_manual(values = fill_vals, guide = "none") +
    scale_shape_manual(values = point_shapes, guide = "none") +
    scale_x_continuous(
      breaks = seq_along(x_levels),
      labels = x_levels,
      limits = c(0.55, length(x_levels) + 0.45),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x = NULL,
      y = "Sound response time (ms)"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      axis.line.y = element_line(linewidth = 0.65, color = "black"),
      axis.line.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 0.55, color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(
        size = x_text_size,
        face = "bold",
        color = "black",
        margin = margin(t = 4)
      ),
      axis.text.y = element_text(
        size = y_text_size,
        color = "black"
      ),
      axis.title.y = element_text(
        size = y_title_size,
        face = "bold",
        margin = margin(r = 6)
      ),
      plot.margin = margin(4, 4, 4, 4),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}