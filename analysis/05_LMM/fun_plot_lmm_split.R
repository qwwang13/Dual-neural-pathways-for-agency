# =========================================
# fun_plot_lmm_split.R
# ITPC-only coefficient forest plots for LMM results
# =========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# -------------------------
# Default colors
# -------------------------
col_ns   <- "#B8B8B8"
col_zero <- "#B0B0B0"
col_text <- "#222222"

# -------------------------
# ITPC terms
# -------------------------
itpc_terms_all <- c(
  "Delta_bandITPC",
  "Theta_bandITPC",
  "Alpha_bandITPC",
  "Beta_bandITPC"
)

itpc_term_labels <- c(
  Delta_bandITPC = "Delta",
  Theta_bandITPC = "Theta",
  Alpha_bandITPC = "Alpha",
  Beta_bandITPC  = "Beta"
)

# -------------------------
# String helpers
# -------------------------
clean_name <- function(x) {
  x %>%
    as.character() %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "")
}

label_for_condition <- function(cond_labels, cond) {
  if (!is.null(names(cond_labels)) && cond %in% names(cond_labels)) {
    return(unname(cond_labels[[cond]]))
  }
  cond
}

stat_match <- function(x_clean, stat) {
  if (stat == "estimate") {
    return(str_detect(x_clean, "estimate|coefficient|coef|beta|^b$"))
  }
  if (stat == "std.error") {
    return(str_detect(x_clean, "stderr|stderror|standarderror|se"))
  }
  if (stat == "p.value") {
    return(str_detect(x_clean, "pvalue|pval|^p$|^pr"))
  }
  if (stat == "df") {
    return(str_detect(x_clean, "^df|dof|degreesfreedom"))
  }
  if (stat == "conf.low") {
    return(str_detect(x_clean, "conflow|cilow|lower|lowerci|lwr"))
  }
  if (stat == "conf.high") {
    return(str_detect(x_clean, "confhigh|cihigh|upper|upperci|upr"))
  }
  
  rep(FALSE, length(x_clean))
}

find_first_col <- function(nms, patterns) {
  x_clean <- clean_name(nms)
  
  idx <- which(x_clean %in% clean_name(patterns))
  if (length(idx) > 0) return(nms[idx[1]])
  
  idx <- which(str_detect(x_clean, paste(clean_name(patterns), collapse = "|")))
  if (length(idx) > 0) return(nms[idx[1]])
  
  NA_character_
}

find_stat_cond_col <- function(nms, stat, cond, cond_label = NULL) {
  x_clean <- clean_name(nms)
  cond_keys <- unique(clean_name(c(cond, cond_label)))
  cond_keys <- cond_keys[!is.na(cond_keys) & cond_keys != ""]
  
  stat_idx <- which(stat_match(x_clean, stat))
  cond_idx <- which(map_lgl(x_clean, function(z) {
    any(str_detect(z, fixed(cond_keys)))
  }))
  
  idx <- intersect(stat_idx, cond_idx)
  
  if (length(idx) > 0) return(nms[idx[1]])
  NA_character_
}

# -------------------------
# Standardize fixed-effect tables
# -------------------------
standardize_long_table <- function(tab, roi, conditions, cond_labels) {
  nms <- names(tab)
  
  term_col <- find_first_col(nms, c("term", "effect", "predictor", "variable"))
  cond_col <- find_first_col(nms, c("condition", "cond", "group"))
  
  if (is.na(term_col) || is.na(cond_col)) {
    return(NULL)
  }
  
  estimate_col <- find_first_col(nms, c("estimate", "beta", "coef", "coefficient"))
  se_col       <- find_first_col(nms, c("std.error", "stderr", "standarderror", "se"))
  df_col       <- find_first_col(nms, c("df", "dof", "degrees.freedom", "degreesfreedom"))
  p_col        <- find_first_col(nms, c("p.value", "pvalue", "pval", "p", "Pr(>|t|)"))
  low_col      <- find_first_col(nms, c("conf.low", "conflow", "ci.low", "cilow", "lower", "lowerci"))
  high_col     <- find_first_col(nms, c("conf.high", "confhigh", "ci.high", "cihigh", "upper", "upperci"))
  
  out <- tibble(
    ROI = roi,
    condition_raw = as.character(tab[[cond_col]]),
    term = as.character(tab[[term_col]]),
    estimate = if (!is.na(estimate_col)) suppressWarnings(as.numeric(tab[[estimate_col]])) else NA_real_,
    std.error = if (!is.na(se_col)) suppressWarnings(as.numeric(tab[[se_col]])) else NA_real_,
    df = if (!is.na(df_col)) suppressWarnings(as.numeric(tab[[df_col]])) else NA_real_,
    p.value = if (!is.na(p_col)) suppressWarnings(as.numeric(tab[[p_col]])) else NA_real_,
    conf.low = if (!is.na(low_col)) suppressWarnings(as.numeric(tab[[low_col]])) else NA_real_,
    conf.high = if (!is.na(high_col)) suppressWarnings(as.numeric(tab[[high_col]])) else NA_real_
  )
  
  out %>%
    mutate(
      condition = map_chr(condition_raw, function(x) {
        xc <- clean_name(x)
        
        hit <- conditions[clean_name(conditions) == xc]
        if (length(hit) == 1) return(hit)
        
        hit <- conditions[clean_name(unname(cond_labels[conditions])) == xc]
        if (length(hit) == 1) return(hit)
        
        x
      })
    ) %>%
    filter(condition %in% conditions) %>%
    select(ROI, condition, term, estimate, std.error, df, p.value, conf.low, conf.high)
}

standardize_wide_table <- function(tab, roi, conditions, cond_labels) {
  nms <- names(tab)
  
  term_col <- find_first_col(nms, c("term", "effect", "predictor", "variable"))
  
  if (is.na(term_col)) {
    if (!is.null(rownames(tab)) && !all(rownames(tab) == as.character(seq_len(nrow(tab))))) {
      tab <- tab %>% mutate(term = rownames(tab), .before = 1)
      term_col <- "term"
      nms <- names(tab)
    } else {
      stop("Cannot find a term/effect/predictor column in fixed_table_wide.")
    }
  }
  
  map_dfr(conditions, function(cond) {
    cond_label <- label_for_condition(cond_labels, cond)
    
    est_col  <- find_stat_cond_col(nms, "estimate",  cond, cond_label)
    se_col   <- find_stat_cond_col(nms, "std.error", cond, cond_label)
    df_col   <- find_stat_cond_col(nms, "df",        cond, cond_label)
    p_col    <- find_stat_cond_col(nms, "p.value",   cond, cond_label)
    low_col  <- find_stat_cond_col(nms, "conf.low",  cond, cond_label)
    high_col <- find_stat_cond_col(nms, "conf.high", cond, cond_label)
    
    tibble(
      ROI = roi,
      condition = cond,
      term = as.character(tab[[term_col]]),
      estimate = if (!is.na(est_col)) suppressWarnings(as.numeric(tab[[est_col]])) else NA_real_,
      std.error = if (!is.na(se_col)) suppressWarnings(as.numeric(tab[[se_col]])) else NA_real_,
      df = if (!is.na(df_col)) suppressWarnings(as.numeric(tab[[df_col]])) else NA_real_,
      p.value = if (!is.na(p_col)) suppressWarnings(as.numeric(tab[[p_col]])) else NA_real_,
      conf.low = if (!is.na(low_col)) suppressWarnings(as.numeric(tab[[low_col]])) else NA_real_,
      conf.high = if (!is.na(high_col)) suppressWarnings(as.numeric(tab[[high_col]])) else NA_real_
    )
  })
}

extract_fixed_long <- function(split_res, roi, conditions, cond_labels) {
  if (!is.null(split_res$fixed_table_wide)) {
    tab <- split_res$fixed_table_wide
    
    long_try <- standardize_long_table(tab, roi, conditions, cond_labels)
    if (!is.null(long_try) && nrow(long_try) > 0) {
      out <- long_try
    } else {
      out <- standardize_wide_table(tab, roi, conditions, cond_labels)
    }
    
  } else if (!is.null(split_res$fixed_table)) {
    tab <- split_res$fixed_table
    
    long_try <- standardize_long_table(tab, roi, conditions, cond_labels)
    if (!is.null(long_try) && nrow(long_try) > 0) {
      out <- long_try
    } else {
      out <- standardize_wide_table(tab, roi, conditions, cond_labels)
    }
    
  } else {
    stop("split_res does not contain fixed_table_wide or fixed_table.")
  }
  
  out %>%
    mutate(
      ci_crit = if_else(
        !is.na(df) & df > 0,
        qt(0.975, df),
        qnorm(0.975)
      ),
      conf.low = if_else(
        is.na(conf.low) & !is.na(estimate) & !is.na(std.error),
        estimate - ci_crit * std.error,
        conf.low
      ),
      conf.high = if_else(
        is.na(conf.high) & !is.na(estimate) & !is.na(std.error),
        estimate + ci_crit * std.error,
        conf.high
      )
    ) %>%
    select(-ci_crit)
}

match_itpc_term <- function(term_vec, terms_select) {
  term_key <- clean_name(term_vec)
  target_key <- clean_name(terms_select)
  
  map_chr(term_key, function(x) {
    hit <- target_key[str_detect(x, fixed(target_key))]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  })
}

# -------------------------
# Build plotting data
# -------------------------
make_itpc_plot_data <- function(res_list, roi_order, conditions, cond_labels) {
  term_lookup <- tibble(
    term = itpc_terms_all,
    term_key = clean_name(itpc_terms_all),
    term_label = unname(itpc_term_labels[itpc_terms_all]),
    y = rev(seq_along(itpc_terms_all))
  )
  
  raw <- map_dfr(roi_order, function(roi) {
    extract_fixed_long(
      split_res = res_list[[roi]],
      roi = roi,
      conditions = conditions,
      cond_labels = cond_labels
    )
  })
  
  raw2 <- raw %>%
    mutate(term_key = match_itpc_term(term, itpc_terms_all)) %>%
    filter(!is.na(term_key)) %>%
    group_by(ROI, condition, term_key) %>%
    slice(1) %>%
    ungroup()
  
  grid <- tidyr::expand_grid(
    ROI = roi_order,
    condition = conditions,
    term_key = term_lookup$term_key
  )
  
  grid %>%
    left_join(raw2, by = c("ROI", "condition", "term_key")) %>%
    left_join(term_lookup, by = "term_key") %>%
    mutate(
      ROI = factor(ROI, levels = roi_order),
      condition = factor(condition, levels = conditions),
      condition_label = map_chr(as.character(condition), ~ label_for_condition(cond_labels, .x)),
      sig = !is.na(p.value) & p.value < 0.05
    )
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(p < 0.001, "p < .001", sprintf("p = %.3f", p))
  )
}

# -------------------------
# Single condition forest panel
# -------------------------
plot_condition_forest <- function(df,
                                  roi,
                                  cond,
                                  cond_label,
                                  cond_color,
                                  show_y = TRUE,
                                  title_size = 8.5,
                                  cap_height = 0.16) {
  
  d <- df %>%
    filter(as.character(ROI) == roi, as.character(condition) == cond) %>%
    arrange(desc(y)) %>%
    mutate(
      plot_col = if_else(sig, cond_color, col_ns),
      p_lab = if_else(sig, format_p(p.value), "")
    )
  
  xs <- c(d$conf.low, d$conf.high, d$estimate, 0)
  xs <- xs[is.finite(xs)]
  
  if (length(xs) == 0) {
    xlim <- c(-1, 1)
  } else {
    xlim <- range(xs, na.rm = TRUE)
    
    if (diff(xlim) == 0) {
      pad <- max(0.15, abs(xlim[1]) * 0.3)
    } else {
      pad <- diff(xlim) * 0.18
    }
    
    xlim <- c(xlim[1] - pad, xlim[2] + pad)
  }
  
  d <- d %>%
    mutate(
      p_x = case_when(
        sig & !is.na(conf.high) & estimate >= 0 ~ conf.high + diff(xlim) * 0.04,
        sig & !is.na(conf.low)  & estimate <  0 ~ conf.low  - diff(xlim) * 0.04,
        TRUE ~ NA_real_
      ),
      p_hjust = if_else(estimate >= 0, 0, 1)
    )
  
  ggplot(d, aes(y = y)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      color = col_zero
    ) +
    geom_segment(
      data = d %>% filter(is.finite(conf.low), is.finite(conf.high)),
      aes(
        x = conf.low,
        xend = conf.high,
        yend = y,
        color = plot_col
      ),
      linewidth = 0.58,
      lineend = "round"
    ) +
    geom_segment(
      data = d %>% filter(is.finite(conf.low)),
      aes(
        x = conf.low,
        xend = conf.low,
        y = y - cap_height,
        yend = y + cap_height,
        color = plot_col
      ),
      linewidth = 0.58,
      lineend = "round"
    ) +
    geom_segment(
      data = d %>% filter(is.finite(conf.high)),
      aes(
        x = conf.high,
        xend = conf.high,
        y = y - cap_height,
        yend = y + cap_height,
        color = plot_col
      ),
      linewidth = 0.58,
      lineend = "round"
    ) +
    geom_point(
      data = d %>% filter(is.finite(estimate)),
      aes(
        x = estimate,
        fill = plot_col,
        color = plot_col
      ),
      shape = 22,
      size = 2.8,
      stroke = 0.25
    ) +
    geom_text(
      data = d %>% filter(sig, is.finite(p_x)),
      aes(
        x = p_x,
        label = p_lab,
        hjust = p_hjust,
        color = plot_col
      ),
      size = 2.5,
      vjust = 0.45
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_y_continuous(
      breaks = d$y,
      labels = if (show_y) d$term_label else rep("", nrow(d)),
      limits = c(0.5, length(itpc_terms_all) + 0.5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(
      title = cond_label,
      x = "Estimated coefficient",
      y = NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 8.5) +
    theme(
      plot.title = element_text(
        color = cond_color,
        face = "bold",
        hjust = 0,
        size = title_size,
        margin = margin(b = 2)
      ),
      axis.title.x = element_text(
        size = 8,
        margin = margin(t = 4)
      ),
      axis.text.x = element_text(
        size = 7,
        color = col_text
      ),
      axis.text.y = element_text(
        size = 7.5,
        color = col_text
      ),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(
        linewidth = 0.35,
        color = "#444444"
      ),
      axis.ticks.x = element_line(
        linewidth = 0.3,
        color = "#444444"
      ),
      plot.margin = margin(4, 18, 4, 4)
    )
}

make_roi_label_plot <- function(roi) {
  ggplot() +
    annotate(
      "text",
      x = 1,
      y = 0.5,
      label = roi,
      hjust = 1,
      vjust = 0.5,
      family = "Arial",
      fontface = "bold",
      size = 3.6,
      color = col_text
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(4, 4, 4, 0)
    )
}

plot_roi_pair <- function(df, roi, conditions, cond_labels, cond_colors) {
  cond1 <- conditions[1]
  cond2 <- conditions[2]
  
  p1 <- plot_condition_forest(
    df = df,
    roi = roi,
    cond = cond1,
    cond_label = label_for_condition(cond_labels, cond1),
    cond_color = cond_colors[[cond1]],
    show_y = TRUE
  )
  
  p2 <- plot_condition_forest(
    df = df,
    roi = roi,
    cond = cond2,
    cond_label = label_for_condition(cond_labels, cond2),
    cond_color = cond_colors[[cond2]],
    show_y = TRUE
  )
  
  roi_lab <- make_roi_label_plot(roi)
  coef_block <- (p1 / p2) + plot_layout(heights = c(1, 1))
  
  (roi_lab | coef_block) +
    plot_layout(widths = c(0.18, 1))
}

plot_group_itpc <- function(res_list,
                            roi_order,
                            conditions,
                            cond_labels,
                            cond_colors,
                            group_title) {
  
  df <- make_itpc_plot_data(
    res_list = res_list,
    roi_order = roi_order,
    conditions = conditions,
    cond_labels = cond_labels
  )
  
  roi_plots <- map(roi_order, function(roi) {
    plot_roi_pair(
      df = df,
      roi = roi,
      conditions = conditions,
      cond_labels = cond_labels,
      cond_colors = cond_colors
    )
  })
  
  wrap_plots(roi_plots, ncol = 1) +
    plot_annotation(
      title = group_title,
      theme = theme(
        plot.title = element_text(
          family = "Arial",
          face = "bold",
          size = 13,
          hjust = 0,
          margin = margin(b = 8)
        )
      )
    )
}

