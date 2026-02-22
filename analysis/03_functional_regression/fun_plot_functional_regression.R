# 03_functional_regression/plot_binding_panels.R
# Minimal, reusable: top panel + bottom panel + combine/save

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
  library(scales)
})

# ============================================
# 1) Top panel: overlap background + mean lines
# ============================================
plot_top_overlap_panel <- function(df_meg_mean,
                                   coef_df,
                                   regions_name,
                                   time_min = 0,
                                   time_max = 1000,
                                   overlap_step_ms = 10,
                                   desired_order_condition = c("vaa_ba","iaa_bi","vaa_iaa"),
                                   style = list(
                                     overlap_fill = c("1"="grey90","2"="grey80","3"="grey40"),
                                     line_colors  = c('iaa_bi'="#42a5f5",'vaa_ba'="#344fac",'vaa_iaa'="#ff4d6d"),
                                     line_lty     = c('iaa_bi'=5,'vaa_ba'=5,'vaa_iaa'=1),
                                     line_labels  = c('iaa_bi'="Involuntary",'vaa_ba'="Voluntary",'vaa_iaa'="Difference")
                                   )) {
  
  desired_order <- regions_name
  
  df_amp <- df_meg_mean %>%
    filter(region %in% regions_name, Time >= time_min, Time <= time_max) %>%
    mutate(
      region = factor(region, levels = desired_order),
      condition = factor(condition, levels = desired_order_condition)
    )
  
  y_upper <- ceiling(max(df_amp$Value, na.rm = TRUE))
  y_lower <- floor(min(df_amp$Value, na.rm = TRUE))
  y_limits <- c(y_upper, y_lower)
  y_breaks <- seq(y_lower, y_upper, by = 1)
  
  coef_use <- coef_df %>%
    mutate(region = factor(region, levels = desired_order)) %>%
    filter(region %in% regions_name)
  
  time_points <- seq(time_min, time_max, by = overlap_step_ms)
  
  overlap_df <- map_dfr(regions_name, function(rg) {
    df_rg <- coef_use %>% filter(region == rg)
    if (nrow(df_rg) == 0) return(tibble())
    
    overlap_count <- sapply(time_points, function(t) {
      sum(t >= df_rg$lower & t < df_rg$upper)
    })
    
    tibble(
      lower = time_points,
      upper = time_points + overlap_step_ms,
      Overlap = as.character(overlap_count),
      region = rg
    ) %>%
      filter(lower < time_max, Overlap != "0")
  }) %>%
    mutate(region = factor(region, levels = desired_order))
  
  ggplot() +
    geom_rect(
      data = overlap_df,
      aes(xmin = lower, xmax = upper, fill = Overlap),
      ymin = -Inf, ymax = Inf, alpha = 0.6
    ) +
    geom_line(
      data = df_amp,
      aes(x = Time, y = Value, group = condition, color = condition, linetype = condition),
      linewidth = 0.8, alpha = 1
    ) +
    facet_wrap(~ region, nrow = 1, scales = "fixed") +
    scale_fill_manual(name = "Overlap", values = style$overlap_fill) +
    scale_y_reverse(name = "Amplitude (z-score)", limits = y_limits, breaks = y_breaks) +
    scale_color_manual(name = NULL, values = style$line_colors, labels = style$line_labels) +
    scale_linetype_manual(name = NULL, values = style$line_lty, labels = style$line_labels) +
    scale_x_continuous(name = NULL, breaks = seq(time_min, time_max, by = 200), limits = c(time_min, time_max)) +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5, linewidth = 0.8) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(1.5), color = "black"),
      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
      panel.spacing = unit(10, "pt"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.box = "horizontal",
      legend.spacing = unit(0.5, "cm"),
      legend.title = element_text(face = "bold", size = rel(1.2)),
      legend.text = element_text(size = rel(1.2)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title.y = element_text(size = rel(1.5)),
      axis.text.x  = element_text(size = rel(1.5)),
      axis.text.y  = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5))
    )
}

# ==========================================================
# 2) Bottom panel: nested facets (Amplitude + Coefficients)
# ==========================================================

plot_bottom_nested <- function(df_meg_mean,
                               coef_df,      # filtered_models_final: region,time_window,lower,upper,coef
                               border_df,    # tibble: region,xmin,xmax
                               regions_name = border_df$region,
                               smooth_span = c("SMA-R"=0.75,"PreSMA-R"=0.40,"PreSMA-L"=0.75),
                               x_break = seq(180, 280, by = 20),
                               strip_fill   = "grey90",
                               smooth_color = "#9b5de5",
                               style = list(
                                 line_colors  = c('iaa_bi'="#42a5f5",'vaa_ba'="#344fac",'vaa_iaa'="#ff4d6d"),
                                 coef_shapes  = c("vaa_ba"=NA,"iaa_bi"=NA,"vaa_iaa"=NA),
                                 linetype_map = c("vaa_ba"="dashed","iaa_bi"="dashed","vaa_iaa"="solid")
                               )) {
  
  # --- ROI table + display label ---
  roi <- border_df %>%
    filter(region %in% regions_name) %>%
    mutate(
      region = as.character(region),
      display_group = paste0(xmin, "-", xmax, "ms")
    )
  
  region_levels <- roi$display_group
  tw_levels <- sort(unique(as.character(coef_df$time_window)))
  
  # --- Amplitude (clip to ROI) ---
  amp_df <- df_meg_mean %>%
    filter(region %in% regions_name) %>%
    mutate(region = as.character(region)) %>%
    left_join(roi, by = "region") %>%
    filter(Time >= xmin, Time <= xmax) %>%
    transmute(
      display_group,
      plot_type = "Amplitude",
      Time,
      y_value = Value,
      condition
    )
  
  # --- Coef points: expand to 10ms bins, fill missing with 0, and clip by ROI (per region) ---
  step_ms <- 10
  
  coef_df2 <- coef_df %>%
    mutate(Time = pmap_dbl(list(lower, upper), \(x, y) (x + y) / 2))
  
  coef_grid <- crossing(
    time_window = unique(coef_df2$time_window),
    region      = unique(coef_df2$region)
  ) %>%
    left_join(
      border_df %>%
        transmute(
          region = as.character(region),
          xmin   = as.numeric(xmin),
          xmax   = as.numeric(xmax)
        ),
      by = "region"
    ) %>%
    mutate(
      lower = map2(xmin - step_ms, xmax, \(a, b) seq(a, b, by = step_ms))
    ) %>%
    unnest(lower) %>%
    mutate(upper = lower + step_ms) %>%
    select(time_window, region, lower, upper)
  
  region_ranges <- border_df %>%
    transmute(
      region   = as.character(region),
      min_time = as.numeric(xmin) - step_ms,
      max_time = as.numeric(xmax) + step_ms
    )
  
  coef_points <- coef_df2 %>%
    merge(region_ranges, by="region") %>%
    filter(lower >= min_time, upper <= max_time) %>%
    select(time_window, region, lower, upper, coef) %>%
    rename(
      y_value_raw = coef
    ) %>%
    mutate(
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      new_lower = pmap(
        list(lower, upper),
        \(lwr, upr) seq(lwr, upr - step_ms, by = step_ms)
      ),
      y_value = pmap_dbl(
        list(y_value_raw, time_window),
        \(x, tw) x / (as.numeric(tw) / step_ms)
      )
    ) %>%
    select(-c(lower, upper, y_value_raw)) %>%
    unnest(new_lower) %>%
    mutate(
      lower = new_lower,
      upper = lower + step_ms
    ) %>%
    select(-new_lower) %>%
    merge(
      coef_grid,
      by = c("time_window", "region", "lower", "upper"),
      all.y = TRUE
    ) %>%
    mutate(
      y_value = if_else(is.na(y_value), 0, y_value)
    ) %>%
    merge(region_ranges, by = "region", all.x = TRUE) %>%
    mutate(
      Time = pmap_dbl(
        list(lower, upper, min_time, max_time),
        \(x, y, min_t, max_t) {
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
    select(-c(min_time, max_time)) %>%
    mutate(region = factor(region, levels = regions_name))%>%
    left_join(
      roi %>% select(region, display_group),
      by = "region"
    ) %>% 
    mutate(plot_type = "Coefficients",
           condition = as.character(time_window))
  
  
  
  # --- Smooth lines (loess on coef_points, per ROI) ---
  if (is.null(smooth_span)) {
    smooth_span <- stats::setNames(rep(0.6, length(regions_name)), regions_name)
  }
  
  smooth_lines <- purrr::map_dfr(regions_name, function(rg) {
    df_rg <- coef_points %>% filter(region == rg)
    sp <- smooth_span[[rg]]
    if (is.null(sp) || is.na(sp)) sp <- 0.6
    
    fit <- loess(y_value ~ Time, data = df_rg, span = sp)
    
    x0 <- roi$xmin[roi$region == rg]
    x1 <- roi$xmax[roi$region == rg]
    grid_x <- seq(x0, x1, length.out = 200)
    pred_y <- predict(fit, newdata = data.frame(Time = grid_x))
    
    tibble(
      display_group = roi$display_group[roi$region == rg],
      plot_type = "Coefficients",
      Time = grid_x,
      y_value = pred_y
    )
  })
  
  # --- Combine for plotting ---
  combined <- bind_rows(
    amp_df,
    coef_points %>% select(display_group, plot_type, Time, y_value, condition)
  )
  
  combined$display_group <- factor(combined$display_group, levels = region_levels)
  combined$plot_type <- factor(combined$plot_type, levels = c("Amplitude", "Coefficients"))
  
  ggplot(combined, aes(x = Time, y = y_value)) +
    geom_line(
      data = combined %>% filter(plot_type == "Amplitude"),
      aes(group = interaction(condition, display_group),
          color = condition,
          linetype = condition),
      linewidth = 0.8,
      alpha = 1
    ) +
    geom_point(
      data = combined %>% filter(plot_type == "Coefficients"),
      aes(color = condition, shape = condition),
      size = 3,
      alpha = 0.9
    ) +
    geom_line(
      data = smooth_lines,
      aes(group = display_group),
      color = smooth_color,
      linewidth = 1.2,
      alpha = 1
    ) +
    ggh4x::facet_nested(
      ~ display_group + plot_type,
      scales = "free_x",
      space = "fixed",
      nest_line = element_line(color = "gray60", linewidth = 0.5)
    ) +
    scale_y_reverse(
      name = NULL,
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = 0.05)
    ) +
    scale_x_continuous(
      name = "Time (ms)",
      breaks = x_break,
      expand = expansion(mult = 0.02)
    ) +
    scale_color_manual(
      name = NULL,
      values = c(style$line_colors, "10"="#FF8C00","20"="#3399FF","40"="#00CC99"),
      breaks = names(c("10"="10 ms","20"="20 ms","40"="40 ms")),
      labels = unname(c("10"="10 ms","20"="20 ms","40"="40 ms"))
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(style$coef_shapes, "10"=15,"20"=16,"40"=17),
      breaks = names(c("10"="10 ms","20"="20 ms","40"="40 ms")),
      labels = unname(c("10"="10 ms","20"="20 ms","40"="40 ms"))
    ) +
    scale_linetype_manual(name = NULL, values = c(style$linetype_map,"10"="blank","20"="blank","40"="blank"), guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = strip_fill, color = strip_fill, linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.5), color = "black"),
      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
      panel.spacing = unit(15, "pt"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = rel(1.2)),
      legend.text  = element_text(size = rel(1.2)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title.y = element_text(size = rel(1.5)),
      axis.text.x  = element_text(size = rel(1.5)),
      axis.text.y  = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5))
    )
}


# =============================================================
# Plot STGR
# =============================================================

# ============================================
# 1) Left panel: overlap background + mean lines
# ============================================
plot_left_overlap_panel_50 <- function(df_meg_mean,
                                   coef_df,
                                   regions_name,
                                   time_min = 0,
                                   time_max = 1000,
                                   overlap_step_ms = 10,
                                   desired_order_condition = c("vaa_ba","iaa_bi","vaa_iaa"),
                                   style = list(
                                     overlap_fill = c("1"="grey90","2"="grey80","3"="grey40"),
                                     line_colors  = c('iaa_bi'="#42a5f5",'vaa_ba'="#344fac",'vaa_iaa'="#ff4d6d"),
                                     line_lty     = c('iaa_bi'=5,'vaa_ba'=5,'vaa_iaa'=1),
                                     line_labels  = c('iaa_bi'="Involuntary",'vaa_ba'="Voluntary",'vaa_iaa'="Difference")
                                   )) {
  
  desired_order <- regions_name
  
  df_amp <- df_meg_mean %>%
    filter(region %in% regions_name, Time >= time_min, Time <= time_max) %>%
    mutate(
      region = factor(region, levels = desired_order),
      condition = factor(condition, levels = desired_order_condition)
    )
  
  y_upper <- ceiling(max(df_amp$Value, na.rm = TRUE))
  y_lower <- floor(min(df_amp$Value, na.rm = TRUE))
  y_limits <- c(y_upper, y_lower)
  y_breaks <- seq(y_lower, y_upper, by = 2)
  
  coef_use <- coef_df %>%
    mutate(region = factor(region, levels = desired_order)) %>%
    filter(region %in% regions_name)
  
  time_points <- seq(time_min, time_max, by = overlap_step_ms)
  
  overlap_df <- map_dfr(regions_name, function(rg) {
    df_rg <- coef_use %>% filter(region == rg)
    if (nrow(df_rg) == 0) return(tibble())
    
    overlap_count <- sapply(time_points, function(t) {
      sum(t >= df_rg$lower & t < df_rg$upper)
    })
    
    tibble(
      lower = time_points,
      upper = time_points + overlap_step_ms,
      Overlap = as.character(overlap_count),
      region = rg
    ) %>%
      filter(lower < time_max, Overlap != "0")
  }) %>%
    mutate(region = factor(region, levels = desired_order))
  
  ggplot() +
    geom_rect(
      data = overlap_df,
      aes(xmin = lower, xmax = upper, fill = Overlap),
      ymin = -Inf, ymax = Inf, alpha = 0.6
    ) +
    geom_line(
      data = df_amp,
      aes(x = Time, y = Value, group = condition, color = condition, linetype = condition),
      linewidth = 0.8, alpha = 1
    ) +
    facet_wrap(~ region, nrow = 1, scales = "fixed") +
    scale_fill_manual(name = "Overlap", values = style$overlap_fill,
                      labels = function(x) paste0(x),
                      guide = guide_legend(
                        override.aes = list(alpha = 0.6),
                        order = 1  # Overlap图例放在第一个
                      )) +
    scale_y_reverse(name = "Amplitude (z-score)", limits = y_limits, breaks = y_breaks) +
    scale_color_manual(name = NULL, values = style$line_colors, labels = style$line_labels,
                       guide = guide_legend(
                         nrow = 1,  # 一行显示
                         order = 2,  # 线条图例放在第二个
                         title.position = "top"
                       )) +
    scale_linetype_manual(name = NULL, values = style$line_lty, labels = style$line_labels,
                          guide = guide_legend(
                            nrow = 1,  # 一行显示
                            order = 2  # 线条图例放在第二个
                          )) +
    scale_x_continuous(name = NULL, breaks = seq(time_min, time_max, by = 200), limits = c(time_min, time_max)) +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5, linewidth = 0.8) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(1.5), color = "black"),
      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
      panel.spacing = unit(10, "pt"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.justification = c(0.5, 1),  # 居中，顶部对齐
      legend.background = element_blank(),
      legend.box = "vertical",  
      legend.box.just = "top",  
      legend.spacing = unit(0.2, "cm"),
      legend.spacing.x = unit(0.5, "cm"), 
      legend.spacing.y = unit(0.3, "cm"), 
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.8, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.title = element_text(face = "bold", size = rel(1.2)),
      legend.text = element_text(size = rel(1.2)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title.y = element_text(size = rel(1.5)),
      axis.text.x  = element_text(size = rel(1.5)),
      axis.text.y  = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5))
    )+
    # 12. 通过guides进一步控制图例
    guides(
      fill = guide_legend(
        title = "Overlap",
        direction = "vertical",  
        ncol = 1,  
        title.position = "top",
        title.hjust = 0.5,
        order = 1,
        position = "left"
      ),
      color = guide_legend(
        title = NULL,
        direction = "horizontal",  
        nrow = 1,  
        title.position = "top",
        title.hjust = 0.5,
        order = 2,
        position = "top"
      ),
      linetype = guide_legend(
        title = NULL,
        direction = "horizontal",  
        nrow = 1,  
        title.position = "top",
        title.hjust = 0.5,
        order = 2,
        position = "top"
      )
    ) 
}

# ==========================================================
# 2) Right panel: nested facets (Amplitude + Coefficients)
# ==========================================================

plot_right_nested_50 <- function(df_meg_mean,
                               coef_df,      # filtered_models_final: region,time_window,lower,upper,coef
                               border_df,    # tibble: region,xmin,xmax
                               regions_name = border_df$region,
                               smooth_span = c("SMA-R"=0.75,"PreSMA-R"=0.40,"PreSMA-L"=0.75),
                               x_break = seq(180, 280, by = 20),
                               strip_fill   = "grey90") {
  
  # --- ROI table + display label ---
  roi <- border_df %>%
    filter(region %in% regions_name) %>%
    mutate(
      region = as.character(region),
      display_group = paste0(xmin, "-", xmax, "ms")
    )
  
  region_levels <- roi$display_group
  tw_levels <- sort(unique(as.character(coef_df$time_window)))
  
  # --- Amplitude (clip to ROI) ---
  amp_df <- df_meg_mean %>%
    filter(region %in% regions_name) %>%
    mutate(region = as.character(region)) %>%
    left_join(roi, by = "region") %>%
    filter(Time >= xmin, Time <= xmax) %>%
    transmute(
      display_group,
      plot_type = "Amplitude",
      Time,
      y_value = Value,
      condition
    )
  
  # --- Coef points: expand to 10ms bins, fill missing with 0, and clip by ROI (per region) ---
  step_ms <- 10
  
  coef_df2 <- coef_df %>%
    mutate(Time = pmap_dbl(list(lower, upper), \(x, y) (x + y) / 2))
  
  coef_grid <- crossing(
    time_window = unique(coef_df2$time_window),
    region      = unique(coef_df2$region)
  ) %>%
    left_join(
      border_df %>%
        transmute(
          region = as.character(region),
          xmin   = as.numeric(xmin),
          xmax   = as.numeric(xmax)
        ),
      by = "region"
    ) %>%
    mutate(
      lower = map2(xmin - step_ms, xmax, \(a, b) seq(a, b, by = step_ms))
    ) %>%
    unnest(lower) %>%
    mutate(upper = lower + step_ms) %>%
    select(time_window, region, lower, upper)
  
  region_ranges <- border_df %>%
    transmute(
      region   = as.character(region),
      min_time = as.numeric(xmin) - step_ms,
      max_time = as.numeric(xmax) + step_ms
    )
  
  coef_points <- coef_df2 %>%
    merge(region_ranges, by="region") %>%
    filter(lower >= min_time, upper <= max_time) %>%
    select(time_window, region, lower, upper, coef) %>%
    rename(
      y_value_raw = coef
    ) %>%
    mutate(
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      new_lower = pmap(
        list(lower, upper),
        \(lwr, upr) seq(lwr, upr - step_ms, by = step_ms)
      ),
      y_value = pmap_dbl(
        list(y_value_raw, time_window),
        \(x, tw) x / (as.numeric(tw) / step_ms)
      )
    ) %>%
    select(-c(lower, upper, y_value_raw)) %>%
    unnest(new_lower) %>%
    mutate(
      lower = new_lower,
      upper = lower + step_ms
    ) %>%
    select(-new_lower) %>%
    merge(
      coef_grid,
      by = c("time_window", "region", "lower", "upper"),
      all.y = TRUE
    ) %>%
    mutate(
      y_value = if_else(is.na(y_value), 0, y_value)
    ) %>%
    merge(region_ranges, by = "region", all.x = TRUE) %>%
    mutate(
      Time = pmap_dbl(
        list(lower, upper, min_time, max_time),
        \(x, y, min_t, max_t) {
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
    select(-c(min_time, max_time)) %>%
    mutate(region = factor(region, levels = regions_name))%>%
    left_join(
      roi %>% select(region, display_group),
      by = "region"
    ) %>% 
    mutate(plot_type = "Coefficients",
           condition = as.character(time_window))
  
  

  # --- Smooth lines (loess on coef_points, per ROI) ---
  if (is.null(smooth_span)) {
    smooth_span <- stats::setNames(rep(0.6, length(regions_name)), regions_name)
  }
  
  smooth_lines <- purrr::map_dfr(regions_name, function(rg) {
    df_rg <- coef_points %>% filter(region == rg)
    sp <- smooth_span[[rg]]
    if (is.null(sp) || is.na(sp)) sp <- 0.6
    
    fit <- loess(y_value ~ Time, data = df_rg, span = sp)
    
    x0 <- roi$xmin[roi$region == rg]
    x1 <- roi$xmax[roi$region == rg]
    grid_x <- seq(x0, x1, length.out = 200)
    pred_y <- predict(fit, newdata = data.frame(Time = grid_x))
    
    tibble(
      display_group = roi$display_group[roi$region == rg],
      plot_type = "Coefficients",
      Time = grid_x,
      y_value = pred_y
    )
  })
  
  # --- Combine for plotting ---
  combined <- bind_rows(
    amp_df,
    coef_points %>% select(display_group, plot_type, Time, y_value, condition)
  )
  
  combined$display_group <- factor(combined$display_group, levels = region_levels)
  combined$plot_type <- factor(combined$plot_type, levels = c("Amplitude", "Coefficients"))
  
  amplitude_plot <- ggplot(amp_df, aes(x = Time, y = y_value)) +
    geom_line(
      aes(group = interaction(condition, display_group), 
          color = condition, 
          linetype = condition),
      linewidth = 0.8,
      alpha = 1
    ) +
    facet_wrap(~ display_group, scales = "fixed") +
    scale_y_reverse(
      name = "Amplitude",
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = 0.05)
    ) +
    scale_x_continuous(
      name = "Time (ms)",
      breaks = seq(0, 1500, by = 10),
      expand = expansion(mult = 0.02)
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        'a50_an_ba' = "#42a5f5",
        'a50_ay_ba' = "#344fac",
        'a50_ay_an' = "#ff4d6d"
      ),
      labels = c("a50_an_ba", "a50_ay_ba", "a50_ay_an")
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "a50_ay_ba" = "dashed",
        "a50_an_ba" = "dashed",
        "a50_ay_an" = "solid"
      ),
      labels = c("a50_an_ba", "a50_ay_ba", "a50_ay_an")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = strip_fill, color = strip_fill, linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.5), color = "black"),
      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
      panel.spacing = unit(15, "pt"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "none",  
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title.y = element_text(size = rel(1.5)),
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5))
    )
  

  coefficient_plot <- ggplot(combined %>% filter(plot_type == "Coefficients"), aes(x = Time, y = y_value)) +
    geom_point(
      aes(color = condition, shape = condition),
      size = 3,
      alpha = 0.9
    ) +
    geom_line(
      data = smooth_lines,
      aes(group = display_group),
      color = "#9b5de5",
      linewidth = 1.2,
      alpha = 1
    ) +
    facet_wrap(~ display_group, scales = "fixed") +
    scale_y_continuous(
      name = "Coefficients",
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = 0.05)
    ) +
    scale_x_continuous(
      name = "Time (ms)",
      breaks = seq(0, 1500, by = 10),
      expand = expansion(mult = 0.02)
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "10" = "#FF8C00", 
        "20" = "#3399FF", 
        "40" = "#00CC99"
      ),
      breaks = c("10", "20", "40"),
      labels = c("10 ms", "20 ms", "40 ms")
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "10" = 15,
        "20" = 16,
        "40" = 17
      ),
      breaks = c("10", "20", "40"),
      labels = c("10 ms", "20 ms", "40 ms")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = strip_fill, color = strip_fill, linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(1.5), color = "black"),
      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
      panel.spacing = unit(15, "pt"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.spacing = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = rel(1.2)),
      legend.text = element_text(size = rel(1.2)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.title.y = element_text(size = rel(1.5)),
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5))
    )
  

  final_match_image_plot <- amplitude_plot + coefficient_plot +
    plot_layout(nrow = 1, widths = c(1, 1)+ theme(text = element_text(family = "Arial")))
  
  final_match_image_plot
  
}
