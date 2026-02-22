# pearson_utils.R
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggExtra)
  library(stringr)
  library(readr)
  library(ggplotify)
  library(gghalves)
})

# ----------------------------
# Helpers
# ----------------------------
label_pretty <- function(x) gsub("chan", "", x)

pretty_feature_name <- function(feature_name) {
  x <- feature_name
  x <- gsub("_bandPower$", "-band power", x)
  x <- gsub("^meanAmplitude$", "Mean amplitude", x)
  x <- gsub("_", " ", x)
  paste0(x, " (z-score)")
}

format_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) "p < 0.001" else paste0("p = ", formatC(p, format = "f", digits = 3))
}

# ----------------------------
# Data IO + Join
# ----------------------------
read_behavior_feature <- function(
    behavior_csv,
    feature_csv,
    feature_encoding = "UTF-8"
) {
  be_df <- read_csv(behavior_csv,
                    locale = locale(encoding = feature_encoding),
                    show_col_types = FALSE)
  
  feature_df <- read_csv(
    feature_csv,
    locale = locale(encoding = feature_encoding),
    show_col_types = FALSE
  )
  
  # 统一类型，避免 join 出问题
  be_df <- be_df %>%
    mutate(
      subject = as.character(subject),
      trial   = as.character(trial),
      condition = as.character(condition)
    )
  
  feature_df <- feature_df %>%
    mutate(
      subject = as.character(subject),
      trial   = as.character(trial),
      condition = as.character(condition),
      label   = as.character(label)
    )
  
  list(be_df = be_df, feature_df = feature_df)
}

make_long_joined_df <- function(
    be_df,
    feature_df,
    feature_cols,
    condition_keep = NULL
) {
  if (!all(feature_cols %in% names(feature_df))) {
    miss <- setdiff(feature_cols, names(feature_df))
    stop("feature_df missing feature columns: ", paste(miss, collapse = ", "))
  }
  
  if (!is.null(condition_keep)) {
    condition_keep <- as.character(condition_keep)
    be_df <- be_df %>% filter(condition %in% condition_keep)
    feature_df <- feature_df %>% filter(condition %in% condition_keep)
  }
  
  beh_tbl <- be_df %>% select(subject, trial, condition, y = behavior_value)
  
  feat_long <- feature_df %>%
    select(subject, trial, condition, label, all_of(feature_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(feature_cols),
      names_to = "feature",
      values_to = "x"
    )
  
  dat <- feat_long %>%
    inner_join(beh_tbl, by = c("subject", "trial", "condition")) %>%
    filter(is.finite(x), is.finite(y))
  
  dat
}

# ----------------------------
# Pearson results
# ----------------------------
compute_pearson_results <- function(
    joined_long_df,
    min_n = 5
) {
  grp <- c("condition", "feature")

  
  joined_long_df %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      n = sum(is.finite(x) & is.finite(y)),
      r = {
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) >= min_n && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
          unname(cor(x[ok], y[ok], method = "pearson"))
        } else {
          NA_real_
        }
      },
      p = {
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) >= min_n && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
          suppressWarnings(cor.test(x[ok], y[ok], method = "pearson")$p.value)
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(p_text = vapply(p, format_p, character(1)))
}

# ----------------------------
# Plot: ggMarginal scatter for ONE condition + optional label subset
# ----------------------------
make_scatter_marginal_plot <- function(
    joined_long_df,
    condition_pick,
    feature_pick,
    show_y_title = TRUE,
    show_right_margin = TRUE,
    min_n = 5,
    y_limits = c(-200, 600),
    y_breaks = seq(-200, 600, by = 200),
    point_size = 5,
    axis_size = 30,
    annotate_size = 8,
    label_cols = NULL,
    label_shapes = NULL,
    show_legend = NULL
    
) {
  dat <- joined_long_df %>%
    filter(condition == condition_pick, feature == feature_pick) %>%
    filter(is.finite(x), is.finite(y))
  
  # 1) Pearson on RAW x
    ct <- cor.test(dat$x, dat$y, method = "pearson")
    r_value <- unname(ct$estimate)
    p_value <- ct$p.value
    p_show  <- format_p(p_value)

  # 2) Z-score x for plotting + fitting
  mx <- mean(dat$x, na.rm = TRUE)
  sx <- sd(dat$x, na.rm = TRUE)
  dat <- dat %>% mutate(x_std = (x - mx) / sx)
  
  # X ticks (step = 1) on x_std
  x_min <- floor(min(dat$x_std) / 1) * 1
  x_max <- ceiling(max(dat$x_std) / 1) * 1
  x_breaks <- seq(x_min, x_max, by = 1)
  if (length(x_breaks) < 3) {
    x_breaks <- seq(x_min - 1, x_max + 1, by = 1)
    x_min <- min(x_breaks); x_max <- max(x_breaks)
  }
  
  # 3) Fit line on z-scored x
  can_fit <- nrow(dat) >= 2 &&
    is.finite(sd(dat$x_std)) && sd(dat$x_std) > 0 &&
    is.finite(sd(dat$y)) && sd(dat$y) > 0
  
  if (can_fit) {
    lm_fit <- lm(y ~ x_std, data = dat)
    a <- unname(coef(lm_fit)[1])
    b <- unname(coef(lm_fit)[2])
  }
  
  # 4) Plot
  ann_label <-paste0("r = ", sprintf("%.2f", r_value), "\n", p_show)

  
  p <- ggplot(dat, aes(x = x_std, y = y, color = label, shape = label)) +
    geom_point(size = point_size, alpha = 0.8) +
    {if (can_fit) geom_abline(intercept = a, slope = b, color = "black", linewidth = 1)} +
    scale_x_continuous(breaks = x_breaks, limits = c(x_min, x_max)) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1.1),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.title.x = element_text(size = axis_size, face = "bold", margin = margin(t = 10)),
      axis.text.x  = element_text(size = axis_size, color = "black", face = "bold"),
      axis.text.y  = element_text(size = axis_size, face = "bold", color = "black"),
      axis.title.y = if (show_y_title) element_text(size = axis_size, face = "bold", margin = margin(r = 10)) else element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      x = pretty_feature_name(feature_pick),  # keeps "(z-score)" label
      y = if (show_y_title) "Sound response time (ms)" else NULL
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = ann_label,
      hjust = 1.1, vjust = 1.2,
      size = annotate_size,
      color = "black",
      fontface = "bold"
    )
  
  if (!is.null(label_cols)) {
    p <- p + scale_color_manual(values = label_cols, labels = label_pretty(names(label_cols)), drop = FALSE)
  }
  if (!is.null(label_shapes)) {
    p <- p + scale_shape_manual(values = label_shapes, labels = label_pretty(names(label_shapes)), drop = FALSE)
  }
  
  pm <- ggMarginal(
    p,
    type = "density",
    groupColour = TRUE,
    groupFill = TRUE,
    size = 5,
    margins = if (isTRUE(show_right_margin)) "both" else "x"
  )
  
  pm
}



# ---------------------------------------------------
# ----------plot_half_violin_subjects----------------
# ---------------------------------------------------
# Helper: "chan81-chan82" -> "81-82"
label_pretty <- function(x) gsub("chan", "", x)

# Build x-axis labels that include ITPC summaries
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
      label   = as.character(label)
    )
  
  df2 <- feature_df %>%
    filter(condition == condition_pick) %>%
    mutate(
      subject = as.character(subject),
      label   = as.character(label)
    ) %>%
    inner_join(map_df2, by = c("subject", "label"))
  
  need_cols <- unname(itpc_cols)
  miss <- setdiff(need_cols, names(df2))
  
  sum_df <- df2 %>%
    group_by(subject, label) %>%
    summarise(
      across(all_of(need_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Build multiline labels safely
  x_labels <- character(nrow(sum_df))
  
  for (i in seq_len(nrow(sum_df))) {
    
    base <- label_pretty(sum_df$label[i])
    
    for (nm in names(itpc_cols)) {
      col_name <- itpc_cols[[nm]]
      val <- sum_df[[col_name]][i]
      
      base <- paste0(
        base,
        "\n", nm, "-band ITPC=",
        sprintf(paste0("%.", digits, "f"), val)
      )
    }
    
    x_labels[i] <- base
  }
  
  setNames(x_labels, sum_df$subject)
}

# Main plot function (half violin + half boxplot + left jitter points)
plot_half_violin_subjects <- function(
    behavior_df,
    feature_df,
    condition_pick,
    map_df,                       # tibble/data.frame with columns: subject, label
    subjects_order = NULL,        # e.g. c("4","2")
    fill_vals = NULL,             # named vector; names must match x_label levels
    itpc_cols = c(Delta = "Delta_bandITPC", Theta = "Theta_bandITPC"),
    y_limits = c(-200, 600),
    y_breaks = c(-200, 0, 200, 400, 600),
    base_size = 16,
    x_text_size = 22,
    y_text_size = 22,
    y_title_size = 24,
    violin_width = 0.5,
    violin_alpha = 0.6,
    box_width = 0.5,
    box_alpha = 0.75,
    point_size = 2,
    point_alpha = 0.6
) {
  stopifnot(all(c("subject", "condition", "behavior_value") %in% names(behavior_df)))
  
  map_df2 <- map_df %>%
    mutate(subject = as.character(subject), label = as.character(label))
  
  # Build x-axis labels with ITPC summaries
  x_lab_map <- build_itpc_xlabels(
    feature_df = feature_df,
    condition_pick = condition_pick,
    map_df = map_df2,
    itpc_cols = itpc_cols
  )
  
  # Plot data from behavior_df
  plot_data <- behavior_df %>%
    filter(condition == condition_pick) %>%
    mutate(subject = as.character(subject)) %>%
    filter(subject %in% map_df2$subject) %>%
    mutate(
      x_label = x_lab_map[subject]
    ) %>%
    filter(!is.na(x_label)) %>%
    mutate(x_label = as.character(x_label))
  
  if (nrow(plot_data) == 0) stop("No behavior data after filtering. Check condition_pick and map_df subjects.")
  
  # Ordering on x-axis
  if (is.null(subjects_order)) {
    subjects_order <- unique(map_df2$subject)
  } else {
    subjects_order <- as.character(subjects_order)
  }
  
  x_levels <- unname(x_lab_map[subjects_order])
  plot_data <- plot_data %>%
    mutate(
      x_label = factor(x_label, levels = x_levels),
      x_num = as.numeric(x_label)
    )
  
  # Default colors (if not provided): follow x_label order
  if (is.null(fill_vals)) {
    fill_vals <- setNames(rep("#999999", length(levels(plot_data$x_label))), levels(plot_data$x_label))
  } else {
    # Ensure names match x_label levels
    if (is.null(names(fill_vals))) stop("fill_vals must be a named vector with names matching x_label levels.")
    missc <- setdiff(levels(plot_data$x_label), names(fill_vals))
    if (length(missc) > 0) stop("fill_vals missing names for x_label levels: ", paste(missc, collapse = ", "))
    fill_vals <- fill_vals[levels(plot_data$x_label)]
  }
  
  ggplot(plot_data, aes(x = x_num, y = behavior_value, fill = x_label, color = x_label)) +
    geom_half_violin(
      side = "R",
      alpha = violin_alpha,
      width = violin_width,
      position = position_nudge(x = 0.05, y = 0),
      trim = TRUE
    ) +
    geom_half_boxplot(
      side = "R",
      alpha = box_alpha,
      width = box_width,
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    geom_point(
      aes(x = x_num - 0.12),
      alpha = point_alpha,
      size = point_size,
      position = position_jitter(width = 0.03, height = 0)
    ) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_color_manual(values = fill_vals, guide = "none") +
    scale_x_continuous(
      breaks = seq_along(levels(plot_data$x_label)),
      labels = levels(plot_data$x_label)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks
    ) +
    labs(x = NULL, y = "Sound response time (ms)") +
    theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 1.1),
      axis.ticks = element_line(linewidth = 1),
      axis.text.x  = element_text(size = x_text_size, face = "bold", color = "black"),
      axis.text.y  = element_text(size = y_text_size, face = "bold", color = "black"),
      axis.title.y = element_text(size = y_title_size, face = "bold", margin = margin(r = 10)),
      plot.margin = margin(8, 8, 8, 8),
      plot.background = element_rect(fill = "white", color = NA)
    )
}