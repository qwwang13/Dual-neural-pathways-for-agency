# =========================================
# fun_plot_lmm_split.R
# Plotting utilities for significant LMM results
# =========================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(ggeffects)
  library(ggnewscale)
})

pretty_term_label <- function(term) {
  term %>%
    gsub("_bandITPC$", "-band ITPC", .) %>%
    gsub("_bandPower$", "-band power", .) %>%
    gsub("^meanAmplitude$", "Mean amplitude", .)
}

# =========================================
# 1) Forest plot (vertical x=term, y=estimate)
# =========================================
plot_forest_vertical_split <- function(
    split_res,
    p_cut = 0.05,
    conditions = c("va_a","ia_a"),
    cond_labels = c(va_a = "Voluntary", ia_a = "Involuntary"),
    terms_select = NULL,
    
    ci_level = 0.95,
    point_size = 4.8,
    line_width  = 1.8,
    dodge_width = 0.50,
    
    point_colors = c(Voluntary="#bc3908", Involuntary="#1F77B4"),
    point_shapes = c(Voluntary=16, Involuntary=17),
    line_types   = c(Voluntary="solid", Involuntary="dashed"),
    
    add_p_labels = TRUE,
    p_digits = 3,
    h_value = 0.20,
    
    zero_linetype = "dashed",
    zero_linewidth = 0.9,
    
    ylim = NULL
) {
  models <- split_res$models
  conditions <- intersect(conditions, names(models))
  stopifnot(length(conditions) > 0)
  
  # ---- collect fixef ----
  get_fixef_df <- function(m, cond) {
    sm <- summary(m)
    coefs <- as.data.frame(sm$coefficients)
    coefs$term <- rownames(coefs)
    rownames(coefs) <- NULL
    
    alpha <- 1 - ci_level
    z <- qnorm(1 - alpha/2)
    
    coefs %>%
      transmute(
        condition_raw = cond,
        term = term,
        estimate = .data[["Estimate"]],
        se = .data[["Std. Error"]],
        p = .data[["Pr(>|t|)"]],
        conf.low  = estimate - z * se,
        conf.high = estimate + z * se
      )
  }
  
  df <- map_dfr(conditions, ~ get_fixef_df(models[[.x]], .x)) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      condition = factor(condition_raw, levels = conditions, labels = cond_labels[conditions]),
      sig = !is.na(p) & p < p_cut,
      est_plot  = ifelse(sig, estimate, NA_real_),
      low_plot  = ifelse(sig, conf.low,  NA_real_),
      high_plot = ifelse(sig, conf.high, NA_real_)
    )
  
  if (!is.null(terms_select)) {
    terms_select <- intersect(terms_select, unique(df$term))
    stopifnot(length(terms_select) > 0)
    df <- df %>% filter(term %in% terms_select)
  }
  
  term_levels <- unique(df$term)
  term_labels <- setNames(sapply(term_levels, pretty_term_label), term_levels)
  df <- df %>% mutate(term_f = factor(term, levels = term_levels))
  
  pd <- position_dodge(width = dodge_width)
  
  p <- ggplot(df, aes(x = term_f, y = est_plot,
                      color = condition, shape = condition, linetype = condition)) +
    geom_hline(yintercept = 0, linewidth = zero_linewidth, linetype = zero_linetype) +
    geom_errorbar(aes(ymin = low_plot, ymax = high_plot),
                  width = 0, linewidth = line_width, position = pd, na.rm = TRUE) +
    geom_point(size = point_size, position = pd, na.rm = TRUE) +
    scale_color_manual(values = point_colors) +
    scale_shape_manual(values = point_shapes) +
    scale_linetype_manual(values = line_types) +
    scale_x_discrete(labels = function(x) str_wrap(term_labels[x], width = 12)) +
    labs(x = NULL, y = "Fixed effects") +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 28, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
      axis.text.y = element_text(size = 30, color = "black", face = "bold"),
      axis.title.y = element_text(size = 30, face = "bold", color = "black"),
      axis.line = element_line(linewidth = 1.1),
      axis.ticks = element_line(linewidth = 1.0)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          linetype = c("solid", "dashed"),
          shape = c(16, 17)
        )
      ),
      shape = "none",
      linetype = "none"
    )
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  if (isTRUE(add_p_labels)) {
    df_annot <- df %>%
      filter(sig) %>%
      mutate(
        lab = case_when(
          is.na(p) ~ "p=NA",
          p < 0.001 ~ "p<0.001",
          TRUE ~ paste0("p=", formatC(p, format="f", digits=p_digits))
        )
      )
    
    p <- p +
      geom_text(
        data = df_annot,
        aes(label = lab),
        position = pd,
        vjust = -1,
        hjust = h_value,
        size = 8,
        show.legend = FALSE,
        angle = 90,
        fontface = "bold"
      )
  }
  
  p
}

# =========================================
# 2) Partial residual overlay (two conditions)
# =========================================
plot_partial_overlay_split <- function(
    split_res,
    term,
    conditions = c("va_a","ia_a"),
    xlim = NULL,
    
    subj_mean_points = TRUE,
    n_grid = 120,
    ci_level = 0.95,
    
    point_size = 5,
    point_alpha = 1,
    line_width = 2,
    ribbon_alpha = 0.30,
    
    point_colors  = c(va_a = "#bc3908", ia_a = "#1F77B4"),
    line_colors   = c(va_a = "#fb8500", ia_a = "#08519C"),
    ribbon_colors = c(va_a = "#ffd29d", ia_a = "#BDD7E7"),
    point_shapes  = c(va_a = 16, ia_a = 17),
    line_types    = c(va_a = "solid", ia_a = "dashed"),
    point_label   = c(va_a = "Voluntary", ia_a = "Involuntary"),
    
    legend_position = c(0.02, 0.98),
    add_annot = FALSE,
    annotate_x = NULL,
    annotate_y = NULL,
    annot_text = NULL
) {
  df <- split_res$data
  models <- split_res$models
  conditions <- intersect(conditions, names(models))
  stopifnot(length(conditions) > 0)
  
  # ---- points ----
  pts_all <- map_dfr(conditions, function(cond) {
    m <- models[[cond]]
    d <- df %>% filter(condition == cond)
    stopifnot(term %in% names(d))
    stopifnot(term %in% names(fixef(m)))
    
    b <- unname(fixef(m)[term])
    
    pts <- d %>%
      transmute(
        condition = cond,
        subject = subject,
        x = .data[[term]],
        y = residuals(m) + b * .data[[term]]
      )
    
    if (isTRUE(subj_mean_points)) {
      pts <- pts %>%
        group_by(condition, subject) %>%
        summarise(x = mean(x, na.rm = TRUE),
                  y = mean(y, na.rm = TRUE),
                  .groups = "drop")
    }
    pts
  })
  
  if (is.null(xlim)) xlim <- range(pts_all$x, na.rm = TRUE)
  
  # ---- effect lines/bands ----
  eff_all <- map_dfr(conditions, function(cond) {
    m <- models[[cond]]
    stopifnot(term %in% names(fixef(m)))
    
    x_seq <- seq(xlim[1], xlim[2], length.out = n_grid)
    newdata <- setNames(data.frame(x_seq), term)
    
    pred <- ggpredict(
      m,
      terms = newdata,
      type = "fixed",
      typical = "zero",
      ci_level = ci_level
    ) %>% as.data.frame()
    
    b0 <- unname(fixef(m)["(Intercept)"])
    
    pred %>%
      transmute(
        condition = cond,
        x = x,
        y = predicted - b0,
        ymin = conf.low - b0,
        ymax = conf.high - b0
      )
  })
  
  # remap display labels
  pts_all$condition <- factor(pts_all$condition, levels = conditions, labels = point_label[conditions])
  eff_all$condition <- factor(eff_all$condition, levels = conditions, labels = point_label[conditions])
  
  # remap style vectors to display labels
  remap <- function(v) { out <- v[conditions]; names(out) <- point_label[conditions]; out }
  pc <- remap(point_colors)
  lc <- remap(line_colors)
  rc <- remap(ribbon_colors)
  ps <- remap(point_shapes)
  lt <- remap(line_types)
  
  p <- ggplot() +
    geom_ribbon(
      data = eff_all,
      aes(x = x, ymin = ymin, ymax = ymax, fill = condition),
      alpha = ribbon_alpha,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = rc, guide = "none") +
    
    ggnewscale::new_scale_color() +
    geom_line(
      data = eff_all,
      aes(x = x, y = y, color = condition, linetype = condition),
      linewidth = line_width,
      show.legend = FALSE
    ) +
    scale_color_manual(values = lc, guide = "none") +
    scale_linetype_manual(values = lt, guide = "none") +
    
    ggnewscale::new_scale_color() +
    geom_point(
      data = pts_all,
      aes(x = x, y = y, color = condition, shape = condition),
      size = point_size,
      alpha = point_alpha
    ) +
    scale_color_manual(values = pc, name = NULL) +
    scale_shape_manual(values = ps, name = NULL) +
    coord_cartesian(xlim = xlim) +
    labs(x = paste0(pretty_term_label(term), " (z-score)"), y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1.1),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.title.x = element_text(size = 30, face = "bold"),
      axis.text.x  = element_text(size = 30, color = "black", face = "bold"),
      axis.text.y  = element_text(size = 30, face = "bold", color = "black"),
      axis.title.y = element_text(size = 30, face = "bold"),
      legend.position = legend_position,
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_blank(),
      legend.text = element_text(size = 28)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(linetype = 0, size = point_size, alpha = 1, shape = unname(ps))
      ),
      shape = "none"
    )
  
  if (isTRUE(add_annot)) {
    ax <- if (!is.null(annotate_x)) annotate_x else (xlim[1] + 0.05 * diff(xlim))
    ay <- if (!is.null(annotate_y)) annotate_y else max(pts_all$y, na.rm = TRUE)
    lab <- if (!is.null(annot_text)) annot_text else ""
    p <- p + annotate("text", x = ax, y = ay, label = lab, hjust = 0, size = 8, fontface = "bold")
  }
  
  p
}

# =========================================
# 3) Partial residual (significant term)
# =========================================
plot_partial_residual <- function(
    split_res,
    p_cut = 0.05,
    terms_select = NULL,
    cond_select = NULL,
    xlim_list = NULL,
    
    point_colors,
    line_colors,
    ribbon_colors,
    point_shapes = NULL,
    line_types   = NULL,
    
    point_size = 2.2,
    point_alpha = 1,
    line_width = 1.2,
    ribbon_alpha = 0.3,
    
    n_grid = 100,
    ci_level = 0.95,
    subj_mean_points = TRUE,
    
    point_label = NULL,
    legend_position_list = NULL,
    annotate_x = NULL,
    annotate_y = NULL
) {
  
  df <- split_res$data
  models <- split_res$models
  fixed_table <- split_res$fixed_table
  
  sig_tbl <- fixed_table %>%
    filter(term != "(Intercept)", !is.na(p_value), p_value < p_cut)
  
  if (!is.null(terms_select)) {
    sig_tbl <- sig_tbl %>% filter(term %in% terms_select)
  }
  
  if (!is.null(cond_select)) {
    stopifnot(length(cond_select) == 1)
    sig_tbl <- sig_tbl %>% filter(condition == cond_select)
  }
  
  if (nrow(sig_tbl) == 0) {
    warning("No significant predictors found under current filters.")
    return(NULL)
  }
  
  plot_one <- function(cond, term) {
    if (!cond %in% names(models)) return(NULL)
    
    m <- models[[cond]]
    d <- df %>% filter(condition == cond)
    
    if (!term %in% names(d)) return(NULL)
    if (!term %in% names(fixef(m))) return(NULL)
    
    b0 <- unname(fixef(m)["(Intercept)"])
    b  <- unname(fixef(m)[term])
    
    pt_col <- point_colors[[cond]]
    ln_col <- line_colors[[cond]]
    rb_col <- ribbon_colors[[cond]]
    
    shp <- if (!is.null(point_shapes) && cond %in% names(point_shapes)) point_shapes[[cond]] else 16
    lty <- if (!is.null(line_types)   && cond %in% names(line_types))   line_types[[cond]]   else "solid"
    
    pts <- d %>%
      transmute(
        subject = subject,
        condition = cond,
        x = .data[[term]],
        y = residuals(m) + b * .data[[term]]
      )
    
    if (isTRUE(subj_mean_points)) {
      pts <- pts %>%
        group_by(subject, condition) %>%
        summarise(
          x = mean(x, na.rm = TRUE),
          y = mean(y, na.rm = TRUE),
          .groups = "drop"
        )
    }
    
    x_rng <- if (!is.null(xlim_list) && term %in% names(xlim_list)) {
      xlim_list[[term]]
    } else {
      range(pts$x, na.rm = TRUE)
    }
    
    x_seq <- seq(x_rng[1], x_rng[2], length.out = n_grid)
    newdata <- setNames(data.frame(x_seq), term)
    
    pred <- ggpredict(
      m,
      terms = newdata,
      type = "fixed",
      typical = "zero",
      ci_level = ci_level
    ) %>% as.data.frame()
    
    if (nrow(pred) == 0) return(NULL)
    
    eff <- pred %>%
      transmute(
        x = x,
        y = predicted - b0,
        ymin = conf.low - b0,
        ymax = conf.high - b0
      )
    
    cond_label <- if (!is.null(point_label) && cond %in% names(point_label)) point_label[[cond]] else cond
    leg_pos <- if (!is.null(legend_position_list) && term %in% names(legend_position_list)) {
      legend_position_list[[term]]
    } else c(0.02, 0.98)
    
    ax <- if (!is.null(annotate_x) && term %in% names(annotate_x)) annotate_x[[term]] else (x_rng[1] + 0.65 * diff(x_rng))
    ay <- if (!is.null(annotate_y) && term %in% names(annotate_y)) annotate_y[[term]] else (min(pts$y, na.rm = TRUE) + 0.2 * diff(range(pts$y, na.rm = TRUE)))
    
    ggplot() +
      geom_ribbon(
        data = eff,
        aes(x = x, ymin = ymin, ymax = ymax),
        fill = rb_col,
        alpha = ribbon_alpha,
        show.legend = FALSE
      ) +
      geom_line(
        data = eff,
        aes(x = x, y = y),
        color = ln_col,
        linetype = lty,
        linewidth = line_width,
        show.legend = FALSE
      ) +
      geom_point(
        data = pts,
        aes(x = x, y = y, color = condition),
        shape = shp,
        size = point_size,
        alpha = point_alpha
      ) +
      coord_cartesian(xlim = x_rng) +
      scale_color_manual(
        values = setNames(pt_col, cond),
        breaks = cond,
        labels = cond_label,
        name = NULL
      ) +
      guides(
        color = guide_legend(
          override.aes = list(shape = shp, size = point_size, alpha = 1),
          title = NULL
        )
      ) +
      labs(
        x = paste0(pretty_term_label(term), " (z-score)"),
        y = NULL
      ) +
      annotate(
        "text",
        x = ax,
        y = ay,
        hjust = 0,
        label = paste0("Slope = ", round(b, 3)),
        size = 8,
        fontface = "bold"
      ) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1.1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.text.x  = element_text(size = 30, color = "black", face = "bold"),
        axis.text.y  = element_text(size = 30, face = "bold", color = "black"),
        axis.title.y = element_text(size = 30, face = "bold"),
        legend.position = leg_pos,
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 28)
      )
  }
  
  sig_pairs <- sig_tbl %>%
    distinct(condition, term) %>%
    arrange(term, condition)
  
  plots <- sig_pairs %>%
    mutate(
      key  = paste0(condition, "_", term),
      plot = map2(condition, term, plot_one)
    ) %>%
    filter(!map_lgl(plot, is.null)) %>%
    { setNames(.$plot, .$key) }
  
  plots
}


plot_split_term_overlay <- function(
    split_res,
    term,
    conditions = NULL,
    xlim = NULL,
    n_grid = 120,
    ci_level = 0.95,
    
    subj_mean_points = TRUE,
    point_size = 2.5,
    point_alpha = 1,
    line_width = 1.2,
    
    point_colors  = NULL,
    line_colors   = NULL,
    ribbon_colors = NULL,
    ribbon_alpha  = 0.30,
    
    point_shapes = NULL,
    line_types   = NULL,
    
    legend_title = NULL,
    legend_position = c(0.02, 0.98),
    point_label = NULL,
    
    add_annot = TRUE,
    annotate_x = NULL,
    annotate_y = NULL,
    annot_digits = 3,
    
    h_value = 0.2
) {

  df <- split_res$data
  models <- split_res$models

  
  if (is.null(conditions)) conditions <- names(models)
  conditions <- intersect(conditions, names(models))

  
  # ---- build partial-residual points (per condition) ----
  pts_all <- map_dfr(conditions, function(cond) {
    m <- models[[cond]]
    d <- df %>% filter(condition == cond)
    
    if (!term %in% names(d)) return(NULL)
    if (!term %in% names(fixef(m))) return(NULL)
    
    b <- unname(fixef(m)[term])
    
    pts <- d %>%
      transmute(
        condition = cond,
        subject = subject,
        x = .data[[term]],
        y = residuals(m) + b * .data[[term]]
      )
    
    if (isTRUE(subj_mean_points)) {
      pts <- pts %>%
        group_by(condition, subject) %>%
        summarise(
          x = mean(x, na.rm = TRUE),
          y = mean(y, na.rm = TRUE),
          .groups = "drop"
        )
    }
    pts
  })
  
  if (nrow(pts_all) == 0) return(NULL)
  
  # ---- x-range ----
  if (is.null(xlim)) xlim <- range(pts_all$x, na.rm = TRUE)
  
  # ---- fitted effect lines + CI ribbons (per condition) ----
  eff_all <- map_dfr(conditions, function(cond) {
    m <- models[[cond]]
    d <- df %>% filter(condition == cond)
    
    if (!term %in% names(d)) return(NULL)
    if (!term %in% names(fixef(m))) return(NULL)
    
    x_seq <- seq(xlim[1], xlim[2], length.out = n_grid)
    newdata <- setNames(data.frame(x_seq), term)
    
    pred <- ggpredict(
      m,
      terms = newdata,
      type = "fixed",
      typical = "zero",
      ci_level = ci_level
    ) %>% as.data.frame()
    
    if (nrow(pred) == 0) return(NULL)
    
    b0 <- unname(fixef(m)["(Intercept)"])
    
    pred %>%
      transmute(
        condition = cond,
        x = x,
        y = predicted - b0,
        ymin = conf.low - b0,
        ymax = conf.high - b0
      )
  })
  
  if (nrow(eff_all) == 0) return(NULL)
  
  # ---- label mapping (optional) ----
  cond_show_levels <- conditions
  if (!is.null(point_label)) {
    pts_all$condition <- factor(pts_all$condition, levels = conditions, labels = point_label[conditions])
    eff_all$condition <- factor(eff_all$condition, levels = conditions, labels = point_label[conditions])
    cond_show_levels <- point_label[conditions]
  } else {
    pts_all$condition <- factor(pts_all$condition, levels = conditions)
    eff_all$condition <- factor(eff_all$condition, levels = conditions)
  }
  
  # ---- remap aesthetics to displayed levels ----
  remap_vec <- function(v) {
    if (is.null(v)) return(NULL)
    out <- v[conditions]
    if (!is.null(point_label)) names(out) <- point_label[conditions]
    out
  }
  
  pc <- remap_vec(point_colors)
  lc <- remap_vec(line_colors)
  rc <- remap_vec(ribbon_colors)
  ps <- remap_vec(point_shapes)
  lt <- remap_vec(line_types)
  
  # ---- plot ----
  p <- ggplot() +
    geom_ribbon(
      data = eff_all,
      aes(x = x, ymin = ymin, ymax = ymax, fill = condition),
      alpha = ribbon_alpha,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = rc, guide = "none") +
    
    ggnewscale::new_scale_color() +
    geom_line(
      data = eff_all,
      aes(x = x, y = y, color = condition, linetype = condition),
      linewidth = line_width,
      show.legend = FALSE
    ) +
    scale_color_manual(values = lc, guide = "none") +
    scale_linetype_manual(values = lt, guide = "none") +
    
    ggnewscale::new_scale_color() +
    geom_point(
      data = pts_all,
      aes(x = x, y = y, color = condition, shape = condition),
      size = point_size,
      alpha = point_alpha
    ) +
    scale_color_manual(values = pc, name = legend_title) +
    scale_shape_manual(values = ps, name = legend_title) +
    coord_cartesian(xlim = xlim) +
    labs(
      title = NULL,
      x = paste0(pretty_term_label(term), " (z-score)"),
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1.1),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.title.x = element_text(size = 30, face = "bold"),
      axis.text.x  = element_text(size = 30, face = "bold", color = "black"),
      axis.text.y  = element_text(size = 30, face = "bold", color = "black"),
      axis.title.y = element_text(size = 30, face = "bold"),
      legend.position = legend_position,
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.title = element_blank(),
      legend.text = element_text(size = 28)
    ) +
    guides(
      color = guide_legend(
        title = legend_title,
        override.aes = list(
          linetype = 0,
          shape = unname(ps),
          size = point_size,
          alpha = 1
        )
      ),
      shape = "none"
    )
  
  # ---- slope annotation (optional) ----
  if (isTRUE(add_annot)) {
    slope_df <- map_dfr(conditions, function(cond) {
      m <- models[[cond]]
      if (!term %in% names(fixef(m))) return(NULL)
      data.frame(
        cond = cond,
        slope = unname(fixef(m)[term]),
        stringsAsFactors = FALSE
      )
    })
    
    if (nrow(slope_df) > 0) {
      if (!is.null(point_label)) slope_df$cond <- point_label[slope_df$cond]
      slope_df$cond <- factor(slope_df$cond, levels = cond_show_levels)
      
      slope_df <- slope_df %>%
        arrange(cond) %>%
        mutate(line = paste0(as.character(cond), ": slope = ", format(round(slope, annot_digits), nsmall = annot_digits))) 
      
      label_txt <- paste(slope_df$line, collapse = "\n")
      
      ax <- if (!is.null(annotate_x)) annotate_x else (xlim[1] + 0.05 * diff(xlim))
      ay <- if (!is.null(annotate_y)) annotate_y else max(pts_all$y, na.rm = TRUE)
      
      p <- p + annotate(
        "text",
        x = ax, y = ay,
        label = label_txt,
        hjust = h_value,
        size = 8,
        fontface = "bold"
      )
    }
  }
  
  p
}