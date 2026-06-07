# pearson_utils.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

label_pretty <- function(x) {
  gsub("chan", "", x)
}

build_itpc_xlabels <- function(
    feature_df,
    condition_pick,
    map_df,
    itpc_cols = c(Delta = "Delta_bandITPC", Theta = "Theta_bandITPC"),
    digits = 3
) {
  map_df2 <- map_df %>%
    mutate(
      subject = as.character(subject),
      label = as.character(label)
    )
  
  feature_df2 <- feature_df %>%
    mutate(
      subject = as.character(subject),
      condition = as.character(condition),
      label = as.character(label)
    )
  
  need_cols <- unname(itpc_cols)
  miss <- setdiff(need_cols, names(feature_df2))
  
  if (length(miss) > 0) {
    stop("feature_df missing ITPC columns: ", paste(miss, collapse = ", "))
  }
  
  sum_df <- feature_df2 %>%
    filter(condition == condition_pick) %>%
    inner_join(map_df2, by = c("subject", "label")) %>%
    group_by(subject, label) %>%
    summarise(
      across(all_of(need_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  label_df <- map_df2 %>%
    left_join(sum_df, by = c("subject", "label"))
  
  if (any(!complete.cases(label_df[, need_cols, drop = FALSE]))) {
    stop("Missing ITPC values after matching feature_df with map_df.")
  }
  
  x_labels <- character(nrow(label_df))
  
  for (i in seq_len(nrow(label_df))) {
    base <- label_pretty(label_df$label[i])
    
    for (nm in names(itpc_cols)) {
      col_name <- itpc_cols[[nm]]
      val <- label_df[[col_name]][i]
      
      base <- paste0(
        base,
        "\n", nm, "-ITPC = ",
        sprintf(paste0("%.", digits, "f"), val)
      )
    }
    
    x_labels[i] <- base
  }
  
  setNames(x_labels, label_df$subject)
}

resolve_plot_values <- function(
    values,
    x_levels,
    x_base_levels,
    default_values,
    value_name
) {
  if (is.null(values)) {
    values <- rep(default_values, length.out = length(x_levels))
    names(values) <- x_levels
    return(values)
  }
  
  if (is.null(names(values))) {
    values <- rep(values, length.out = length(x_levels))
    names(values) <- x_levels
    return(values)
  }
  
  if (all(x_levels %in% names(values))) {
    return(values[x_levels])
  }
  
  if (all(x_base_levels %in% names(values))) {
    return(setNames(values[x_base_levels], x_levels))
  }
  
  stop(
    value_name,
    " names must match either full x-axis labels or base channel labels: ",
    paste(x_base_levels, collapse = ", ")
  )
}

plot_half_violin_subjects <- function(
    behavior_df,
    feature_df,
    condition_pick,
    map_df,
    subjects_order = NULL,
    fill_vals = NULL,
    point_shapes = NULL,
    itpc_cols = c(Delta = "Delta_bandITPC", Theta = "Theta_bandITPC"),
    itpc_digits = 3,
    y_limits = c(-200, 600),
    y_breaks = c(-200, 0, 200, 400, 600),
    base_size = 8,
    x_text_size = 10,
    y_text_size = 9,
    y_title_size = 10.5,
    point_size = 2.1,
    point_alpha = 0.72,
    point_nudge = -0.06,
    jitter_width = 0.075,
    point_stroke = 0.20,
    summary_nudge = 0.13,
    summary_linewidth = 0.75,
    median_linewidth = 1.35,
    median_width = 0.30,
    range_cap_width = 0.18,
    zero_line_color = "grey72",
    zero_line_linetype = "dashed",
    zero_line_linewidth = 0.45
) {
  stopifnot(all(c("subject", "condition", "trial", "behavior_value") %in% names(behavior_df)))
  stopifnot(all(c("subject", "label", "condition", "trial") %in% names(feature_df)))
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
  
  map_df2 <- map_df2 %>%
    filter(subject %in% subjects_order)
  
  x_lab_map <- build_itpc_xlabels(
    feature_df = feature_df,
    condition_pick = condition_pick,
    map_df = map_df2,
    itpc_cols = itpc_cols,
    digits = itpc_digits
  )
  
  label_df <- map_df2 %>%
    mutate(
      subject = factor(subject, levels = subjects_order),
      x_base = label_pretty(label),
      x_label = unname(x_lab_map[as.character(subject)])
    ) %>%
    arrange(subject)
  
  x_levels <- label_df$x_label
  x_base_levels <- label_df$x_base

  feature_trial_keys <- feature_df %>%
    mutate(
      subject = as.character(subject),
      label = as.character(label),
      condition = as.character(condition),
      trial = as.character(trial)
    ) %>%
    filter(condition == condition_pick) %>%
    inner_join(map_df2, by = c("subject", "label")) %>%
    distinct(subject, label, condition, trial)
  
  plot_data <- behavior_df %>%
    mutate(
      subject = as.character(subject),
      condition = as.character(condition),
      trial = as.character(trial)
    ) %>%
    filter(condition == condition_pick) %>%
    inner_join(map_df2, by = "subject") %>%
    inner_join(
      feature_trial_keys,
      by = c("subject", "label", "condition", "trial")
    ) %>%
    mutate(
      x_label = unname(x_lab_map[subject]),
      x_label = factor(x_label, levels = x_levels),
      x_num = as.numeric(x_label),
      x_pt = x_num + point_nudge
    ) %>%
    filter(!is.na(x_label), is.finite(behavior_value))
  
  if (nrow(plot_data) == 0) {
    stop("No behavior data after filtering. Check condition_pick, map_df, and subjects_order.")
  }
  
  fill_vals <- resolve_plot_values(
    values = fill_vals,
    x_levels = x_levels,
    x_base_levels = x_base_levels,
    default_values = c("#02c39a", "#ff758f"),
    value_name = "fill_vals"
  )
  
  point_shapes <- resolve_plot_values(
    values = point_shapes,
    x_levels = x_levels,
    x_base_levels = x_base_levels,
    default_values = c(21, 24),
    value_name = "point_shapes"
  )
  
  sum_df <- plot_data %>%
    group_by(x_label, x_num) %>%
    summarise(
      y_med = median(behavior_value, na.rm = TRUE),
      y_q1 = quantile(behavior_value, 0.25, na.rm = TRUE),
      y_q3 = quantile(behavior_value, 0.75, na.rm = TRUE),
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
        x = x_sum - range_cap_width / 2,
        xend = x_sum + range_cap_width / 2,
        y = y_q1,
        yend = y_q1,
        color = x_label
      ),
      linewidth = summary_linewidth,
      lineend = "round"
    ) +
    geom_segment(
      data = sum_df,
      aes(
        x = x_sum - range_cap_width / 2,
        xend = x_sum + range_cap_width / 2,
        y = y_q3,
        yend = y_q3,
        color = x_label
      ),
      linewidth = summary_linewidth,
      lineend = "round"
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
