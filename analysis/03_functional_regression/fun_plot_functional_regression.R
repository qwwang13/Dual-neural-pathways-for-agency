# fun_plot_binding_by_columns.R
# Horizontal layout:
# Col 1: original amplitude
# Col 2: zoomed amplitude
# Col 3: estimated coefficients
# Each row = one ROI

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
  library(scales)
})

build_zoom_x_scales <- function(border_df) {
  lst <- vector("list", nrow(border_df))
  
  for (i in seq_len(nrow(border_df))) {
    rg <- as.character(border_df$region[i])
    x1 <- border_df$xmin[i]
    x2 <- border_df$xmax[i]
    
    brks <- pretty(c(x1, x2), n = 4)
    brks <- brks[brks >= x1 & brks <= x2]
    
    lst[[i]] <- eval(parse(text =
                             paste0(
                               "region == '", rg, "' ~ scale_x_continuous(limits = c(",
                               x1, ", ", x2, "), breaks = c(",
                               paste(brks, collapse = ", "),
                               "))"
                             )
    ))
  }
  
  lst
}

build_overlap_df <- function(coef_df,
                             regions_name,
                             time_min = 0,
                             time_max = 1000,
                             overlap_step_ms = 10) {
  time_points <- seq(time_min, time_max, by = overlap_step_ms)
  
  map_dfr(regions_name, function(rg) {
    df_rg <- coef_df %>%
      filter(region == rg)
    
    if (nrow(df_rg) == 0) return(tibble())
    
    overlap_count <- sapply(time_points, function(t) {
      sum(t >= df_rg$lower & t < df_rg$upper)
    })
    
    tibble(
      lower   = time_points,
      upper   = time_points + overlap_step_ms,
      Overlap = as.character(overlap_count),
      region  = rg
    ) %>%
      filter(lower < time_max, Overlap != "0")
  }) %>%
    mutate(region = factor(region, levels = regions_name))
}

build_coef_points_and_smooth <- function(coef_df,
                                         border_df,
                                         regions_name,
                                         smooth_span) {
  step_ms <- 10
  
  roi <- border_df %>%
    filter(region %in% regions_name) %>%
    mutate(region = as.character(region))
  
  coef_df2 <- coef_df %>%
    mutate(
      region = as.character(region),
      time_window = as.character(time_window),
      Time = map2_dbl(lower, upper, ~ (.x + .y) / 2)
    )
  
  coef_grid <- tidyr::crossing(
    time_window = unique(coef_df2$time_window),
    region      = unique(coef_df2$region)
  ) %>%
    left_join(
      roi %>% transmute(region, xmin, xmax),
      by = "region"
    ) %>%
    mutate(
      lower = map2(xmin - step_ms, xmax, ~ seq(.x, .y, by = step_ms))
    ) %>%
    unnest(lower) %>%
    mutate(
      upper = lower + step_ms
    ) %>%
    select(time_window, region, lower, upper)
  
  region_ranges <- roi %>%
    transmute(
      region,
      min_time = xmin - step_ms,
      max_time = xmax + step_ms
    )
  
  # 注意：
  # action_binding.csv / outcome_binding.csv 保存前已经做过：
  # coef = coef / (as.numeric(time_window) / 10)
  # 所以这里不能再除一次，直接使用 csv 里的 coef。
  coef_points <- coef_df2 %>%
    left_join(region_ranges, by = "region") %>%
    filter(lower >= min_time, upper <= max_time) %>%
    select(time_window, region, lower, upper, coef) %>%
    rename(y_value_raw = coef) %>%
    mutate(
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      new_lower = map2(lower, upper, ~ seq(.x, .y - step_ms, by = step_ms)),
      y_value = y_value_raw
    ) %>%
    select(-lower, -upper, -y_value_raw) %>%
    unnest(new_lower) %>%
    mutate(
      lower = new_lower,
      upper = lower + step_ms
    ) %>%
    select(-new_lower) %>%
    right_join(
      coef_grid,
      by = c("time_window", "region", "lower", "upper")
    ) %>%
    mutate(
      y_value = if_else(is.na(y_value), 0, y_value)
    ) %>%
    left_join(region_ranges, by = "region") %>%
    mutate(
      Time = pmap_dbl(
        list(lower, upper, min_time, max_time),
        function(x, y, min_t, max_t) {
          if (x == min_t) {
            min_t + step_ms
          } else if (y == max_t) {
            max_t - step_ms
          } else {
            (x + y) / 2
          }
        }
      )
    ) %>%
    select(-min_time, -max_time) %>%
    mutate(
      region = factor(region, levels = regions_name),
      time_window = factor(time_window, levels = c("10", "20", "40"))
    )
  
  smooth_lines <- map_dfr(regions_name, function(rg) {
    df_rg <- coef_points %>%
      filter(region == rg)
    
    sp <- smooth_span[[rg]]
    if (is.null(sp) || is.na(sp)) sp <- 0.6
    
    fit <- loess(y_value ~ Time, data = df_rg, span = sp)
    
    x0 <- roi$xmin[roi$region == rg]
    x1 <- roi$xmax[roi$region == rg]
    
    grid_x <- seq(x0, x1, length.out = 200)
    pred_y <- predict(fit, newdata = data.frame(Time = grid_x))
    
    tibble(
      region = factor(rg, levels = regions_name),
      Time = grid_x,
      y_value = pred_y
    )
  })
  
  list(
    coef_points = coef_points,
    smooth_lines = smooth_lines
  )
}

theme_white_compact <- function(base_size = 8.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = base_size + 0.6,
        color = "black"
      ),
      
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.40
      ),
      panel.spacing = unit(5, "pt"),
      
      # 避免 axis.line 和 panel.border 叠加导致边框粗细不一致
      axis.line = element_blank(),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.3
      ),
      axis.ticks.length = unit(1.1, "mm"),
      
      axis.title = element_text(
        size = base_size + 0.2,
        color = "black"
      ),
      axis.text = element_text(
        size = base_size - 0.1,
        color = "black"
      ),
      
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(
        face = "bold",
        size = base_size
      ),
      legend.text = element_text(
        size = base_size - 0.1
      ),
      
      plot.title = element_text(
        face = "bold",
        size = base_size + 1,
        hjust = 0.5
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
}

plot_binding_by_columns <- function(df_meg_mean_raw,
                                    coef_df_raw,
                                    regions_name,
                                    border_df,
                                    condition_order,
                                    style,
                                    smooth_span,
                                    time_min = 0,
                                    time_max = 1000,
                                    overlap_step_ms = 10,
                                    reverse_amp_y = TRUE,
                                    reverse_coef_y = TRUE,
                                    plot_widths = c(1.42, 0.7, 0.7),
                                    base_size = 8.5,
                                    base_family = "Arial") {
  
  df_meg_mean <- df_meg_mean_raw %>%
    rename(
      time = Time,
      amp  = Value
    ) %>%
    mutate(
      region = factor(region, levels = regions_name),
      condition = factor(condition, levels = condition_order)
    )
  
  coef_df <- coef_df_raw %>%
    mutate(
      region = factor(region, levels = regions_name),
      time_window = as.character(time_window)
    )
  
  border_df <- border_df %>%
    mutate(
      region = factor(region, levels = regions_name)
    )
  
  df_full <- df_meg_mean %>%
    filter(
      region %in% regions_name,
      time >= time_min,
      time <= time_max,
      condition %in% condition_order
    ) %>%
    mutate(
      region = factor(region, levels = regions_name)
    )
  
  df_zoom <- df_full %>%
    left_join(border_df, by = "region") %>%
    filter(
      time >= xmin,
      time <= xmax
    ) %>%
    mutate(
      region = factor(region, levels = regions_name)
    )
  
  coef_use <- coef_df %>%
    filter(region %in% regions_name) %>%
    mutate(
      region = factor(region, levels = regions_name)
    )
  
  coef_build <- build_coef_points_and_smooth(
    coef_df = coef_use,
    border_df = border_df,
    regions_name = regions_name,
    smooth_span = smooth_span
  )
  
  coef_points  <- coef_build$coef_points
  smooth_lines <- coef_build$smooth_lines
  
  overlap_df <- build_overlap_df(
    coef_df = coef_use,
    regions_name = regions_name,
    time_min = time_min,
    time_max = time_max,
    overlap_step_ms = overlap_step_ms
  )
  
  zoom_x_scales <- build_zoom_x_scales(border_df)
  
  # Col 1: original amplitude
  y1_upper  <- ceiling(max(df_full$amp, na.rm = TRUE))
  y1_lower  <- floor(min(df_full$amp, na.rm = TRUE))
  y1_limits <- c(y1_upper, y1_lower)
  y1_breaks <- seq(y1_lower, y1_upper, by = 1)
  
  # Col 2: zoomed amplitude
  y2_upper  <- ceiling(max(df_zoom$amp, na.rm = TRUE))
  y2_lower  <- floor(min(df_zoom$amp, na.rm = TRUE))
  y2_limits <- c(y2_upper, y2_lower)
  y2_breaks <- pretty(c(y2_lower, y2_upper), n = 5)
  
  p_col1 <- ggplot() +
    geom_rect(
      data = overlap_df,
      aes(xmin = lower, xmax = upper, fill = Overlap),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.55
    ) +
    geom_line(
      data = df_full,
      aes(
        x = time,
        y = amp,
        color = condition,
        linetype = condition
      ),
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.25
    ) +
    geom_vline(
      data = border_df,
      aes(xintercept = xmin),
      inherit.aes = FALSE,
      linewidth = 0.28,
      linetype = "dotted",
      colour = "grey65"
    ) +
    geom_vline(
      data = border_df,
      aes(xintercept = xmax),
      inherit.aes = FALSE,
      linewidth = 0.28,
      linetype = "dotted",
      colour = "grey65"
    ) +
    ggh4x::facet_wrap2(
      ~ region,
      ncol = 1,
      scales = "fixed",
      drop = FALSE,
      axes = "x",
      remove_labels = "none"
    ) +
    scale_fill_manual(
      name = "Overlap",
      values = style$overlap_fill,
      guide = guide_legend(order = 1)
    ) +
    scale_color_manual(
      name = "Condition",
      values = style$line_colors,
      labels = style$line_labels,
      guide = guide_legend(order = 2)
    ) +
    scale_linetype_manual(
      name = "Condition",
      values = style$line_lty,
      labels = style$line_labels,
      guide = guide_legend(order = 2)
    ) +
    scale_x_continuous(
      breaks = seq(time_min, time_max, by = 200),
      limits = c(time_min, time_max)
    ) +
    labs(
      title = "Original amplitude",
      x = "Time (ms)",
      y = "Amplitude (z-score)"
    ) +
    theme_white_compact(
      base_size = base_size,
      base_family = base_family
    )
  
  if (reverse_amp_y) {
    p_col1 <- p_col1 +
      scale_y_reverse(
        limits = y1_limits,
        breaks = y1_breaks
      )
  } else {
    p_col1 <- p_col1 +
      scale_y_continuous(
        limits = rev(y1_limits),
        breaks = y1_breaks
      )
  }
  
  p_col2 <- ggplot(
    df_zoom,
    aes(
      x = time,
      y = amp,
      color = condition,
      linetype = condition
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.25
    ) +
    facet_wrap(
      ~ region,
      ncol = 1,
      scales = "free_x",
      drop = FALSE
    ) +
    ggh4x::facetted_pos_scales(
      x = zoom_x_scales
    ) +
    scale_color_manual(
      name = "Condition",
      values = style$line_colors,
      labels = style$line_labels,
      guide = "none"
    ) +
    scale_linetype_manual(
      name = "Condition",
      values = style$line_lty,
      labels = style$line_labels,
      guide = "none"
    ) +
    labs(
      title = "Selected time window",
      x = "Time (ms)",
      y = "Amplitude (z-score)"
    ) +
    theme_white_compact(
      base_size = base_size,
      base_family = base_family
    )
  
  if (reverse_amp_y) {
    p_col2 <- p_col2 +
      scale_y_reverse(
        limits = y2_limits,
        breaks = y2_breaks
      )
  } else {
    p_col2 <- p_col2 +
      scale_y_continuous(
        limits = rev(y2_limits),
        breaks = y2_breaks
      )
  }
  
  p_col3 <- ggplot() +
    geom_point(
      data = coef_points,
      aes(
        x = Time,
        y = y_value,
        shape = time_window
      ),
      color = style$coef_point_color,
      size = 0.95,
      alpha = 0.75
    ) +
    geom_line(
      data = smooth_lines,
      aes(
        x = Time,
        y = y_value,
        group = region
      ),
      color = style$smooth_color,
      linewidth = 1.0
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.25
    ) +
    facet_wrap(
      ~ region,
      ncol = 1,
      scales = "free",
      drop = FALSE
    ) +
    ggh4x::facetted_pos_scales(
      x = zoom_x_scales
    ) +
    scale_shape_manual(
      name = "Window",
      values = style$coef_shapes,
      breaks = c("10", "20", "40"),
      labels = c("10 ms", "20 ms", "40 ms"),
      guide = guide_legend(order = 3)
    ) +
    labs(
      title = "Estimated coefficients",
      x = "Time (ms)",
      y = "Estimated coefficient"
    ) +
    theme_white_compact(
      base_size = base_size,
      base_family = base_family
    )
  
  if (reverse_coef_y) {
    p_col3 <- p_col3 +
      scale_y_reverse(
        breaks = function(x) pretty(x, n = 4)
      )
  } else {
    p_col3 <- p_col3 +
      scale_y_continuous(
        breaks = function(x) pretty(x, n = 4)
      )
  }
  
  final_plot <- (p_col1 | p_col2 | p_col3) +
    patchwork::plot_layout(
      widths = plot_widths,
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  return(final_plot)
}