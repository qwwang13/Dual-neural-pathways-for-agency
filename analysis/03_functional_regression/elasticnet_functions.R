# R/elasticnet_functions.R
# =========================================
# ElasticNet
# =========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(glmnet)
})

# Convert region label like "PreSMAL" -> "PreSMA-L"
make_region_label <- function(region_raw) {
  paste0(
    substr(region_raw, 1, nchar(region_raw) - 1),
    "-",
    substr(region_raw, nchar(region_raw), nchar(region_raw))
  )
}

# Map each time point to a bin index (p_number) in [1, n_bins]
assign_time_bin <- function(time_ms, n_bins, lower_ms, upper_ms) {
  breaks <- seq(lower_ms, upper_ms, length.out = n_bins + 1)
  k <- findInterval(time_ms, breaks, rightmost.closed = TRUE)
  k[time_ms < lower_ms | time_ms > upper_ms] <- NA_integer_
  k
}

# Given p_number and time_window, compute [lower, upper] bounds (ms)
bin_time_bounds <- function(p_number, time_window_ms, lower0_ms) {
  lower <- (p_number - 1) * time_window_ms + lower0_ms
  upper <- p_number * time_window_ms + lower0_ms
  tibble(lower = lower, upper = upper)
}

# R^2 for each lambda (columns of fitted matrix)
r2_by_lambda <- function(fit_values, y_true) {
  sst <- sum((y_true - mean(y_true))^2)
  apply(fit_values, 2, function(pred) {
    sse <- sum((y_true - pred)^2)
    1 - sse / sst
  })
}

# Strict: require consecutive bins with step=1 (e.g., 5,6,7,...)
has_consecutive_strict <- function(pos_vec, p) {
  if (length(pos_vec) < p) return(FALSE)
  diffs <- diff(sort(pos_vec))
  r <- rle(diffs)
  any(r$values == 1 & r$lengths >= (p - 1))
}

# Relaxed: allow steps in allowed_steps (default 1 or 2),
# and require p points in a "near-consecutive" run
has_consecutive_relaxed <- function(pos_vec, p, allowed_steps = c(1, 2)) {
  if (length(pos_vec) < p) return(FALSE)
  diffs <- diff(sort(pos_vec))
  
  run_len <- 0
  for (d in diffs) {
    if (d %in% allowed_steps) {
      run_len <- run_len + 1
      if (run_len >= (p - 1)) return(TRUE)
    } else {
      run_len <- 0
    }
  }
  FALSE
}


# Build binned brain-wave features:
# Output columns include: region, subject, condition, time_window, p_number, Value
prep_brain_features <- function(brain_long,
                                time_windows_ms,
                                time_lower_ms = 0,
                                time_upper_ms = 1000) {
  expanded <- map_df(time_windows_ms, ~mutate(brain_long, time_window = .x))
  
  expanded %>%
    mutate(
      dimension = (time_upper_ms - time_lower_ms) / time_window,
      p_number = pmap_int(
        list(Time, dimension),
        ~assign_time_bin(..1, n_bins = ..2, lower_ms = time_lower_ms, upper_ms = time_upper_ms)
      )
    ) %>%
    filter(!is.na(p_number)) %>%
    group_by(region, subject, condition, time_window, dimension, p_number) %>%
    summarise(Value = mean(Value), .groups = "drop")
}

# Fit ElasticNet for ONE (region, time_window, alpha, intercept)
# Returns a tibble with one row per lambda, including selected bins and coefficients (list columns)
fit_elasticnet_one <- function(feature_df,
                               binding_df,
                               region_name,
                               time_window_ms,
                               alpha_value,
                               intercept_flag,
                               binding_col,
                               n_lambda = 100,
                               min_n = 5) {
  df <- feature_df %>%
    filter(region == region_name, time_window == time_window_ms) %>%
    select(subject, p_number, Value) %>%
    pivot_wider(names_from = p_number, values_from = Value) %>%
    inner_join(binding_df, by = "subject")
  
  if (nrow(df) < min_n) return(tibble())
  
  x <- df %>%
    select(-subject, -all_of(binding_col)) %>%
    as.matrix()
  
  y <- df[[binding_col]]
  
  fit <- glmnet(
    x = x, y = y,
    alpha = alpha_value,
    family = "gaussian",
    nlambda = n_lambda,
    intercept = intercept_flag
  )
  
  lambdas <- fit$lambda
  fitted_vals <- predict(fit, newx = x, type = "response")  # n x nlambda
  r2 <- r2_by_lambda(fitted_vals, y_true = y)
  
  coef_mat <- as.matrix(coef(fit))  # (p+1) x nlambda
  
  map_dfr(seq_along(lambdas), function(i) {
    lam <- lambdas[i]
    beta <- coef_mat[-1, i]  # drop intercept row
    pos <- which(beta != 0)
    co <- beta[pos]
    
    tibble(
      region = region_name,
      time_window = time_window_ms,
      alpha = alpha_value,
      intercept = intercept_flag,
      lambda = lam,
      r2 = unname(r2[i]),
      n_selected = length(pos),
      pos = list(as.integer(pos)),
      coef = list(as.numeric(co))
    )
  })
}

# Run full ElasticNet grid and return all lambda-level results
run_elasticnet_grid <- function(feature_df,
                                binding_df,
                                regions,
                                time_windows_ms,
                                alpha_grid = seq(0.1, 1, 0.1),
                                intercepts = c(TRUE, FALSE),
                                binding_col,
                                n_lambda = 100) {
  settings <- tidyr::crossing(
    region = regions,
    time_window = time_windows_ms,
    alpha = alpha_grid,
    intercept = intercepts
  )
  
  pmap_dfr(
    list(settings$region, settings$time_window, settings$alpha, settings$intercept),
    ~fit_elasticnet_one(
      feature_df = feature_df,
      binding_df = binding_df,
      region_name = ..1,
      time_window_ms = ..2,
      alpha_value = ..3,
      intercept_flag = ..4,
      binding_col = binding_col,
      n_lambda = n_lambda
    )
  )
}

# Filter results:
# For each (region,time_window), keep the "best" model subject to consecutive constraint.
# "Best" = sparsest (min n_selected), tie-break by highest r2.
# Output is expanded per selected p_number with time bounds.
filter_elasticnet_by_consecutive <- function(all_models_df,
                                             conti_p_list,
                                             time_lower_ms = 0,
                                             time_upper_ms = 1000,
                                             consecutive_mode = c("strict", "relaxed"),
                                             allowed_steps = c(1, 2)) {
  
  keys <- all_models_df %>% distinct(region, time_window)
  
  map_dfr(seq_len(nrow(keys)), function(i) {
    rg <- keys$region[i]
    tw <- keys$time_window[i]
    con_p <- conti_p_list[[rg]][[as.character(tw)]]
    
    cand <- all_models_df %>%
      filter(region == rg, time_window == tw) %>%
      mutate(
        ok = map_lgl(pos, function(pos_vec) {
          if (consecutive_mode == "strict") {
            has_consecutive_strict(pos_vec, p = con_p)
          } else {
            has_consecutive_relaxed(pos_vec, p = con_p, allowed_steps = allowed_steps)
          }
        })
      ) %>%
      filter(ok)
    
    if (nrow(cand) == 0) return(tibble())
    
    # Choose best: sparsest (min n_selected), tie-break by max r2
    # best <- cand %>%
    #   arrange(n_selected, desc(r2)) %>%
    #   slice(1)
    
    best <- cand %>% 
      arrange(n_selected)%>%
      slice(1)
    
    pos_vec  <- best$pos[[1]]
    coef_vec <- best$coef[[1]]
    
    bounds <- bin_time_bounds(pos_vec, time_window_ms = tw, lower0_ms = time_lower_ms) %>%
      mutate(upper = pmin(upper, time_upper_ms))
    
    tibble(
      region = rg,
      time_window = tw,
      alpha = best$alpha,
      intercept = best$intercept,
      lambda = best$lambda,
      r2 = best$r2,
      n_selected = best$n_selected,
      consecutive_mode = consecutive_mode,
      allowed_steps = list(1),
      p_number = pos_vec,
      coef = coef_vec
    ) %>%
      bind_cols(bounds)
  })
}

# ---------------------------
# Post-filtering utilities
# ---------------------------
# Merge overlapping/adjacent intervals and return merged tibble(lower, upper)
merge_intervals <- function(df, gap_ms = 0) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(lower = numeric(), upper = numeric()))
  }
  
  x <- df %>%
    transmute(
      lower = suppressWarnings(as.numeric(lower)),
      upper = suppressWarnings(as.numeric(upper))
    ) %>%
    filter(!is.na(lower), !is.na(upper), upper > lower) %>%
    arrange(lower, upper)
  
  if (nrow(x) == 0) {
    return(tibble(lower = numeric(), upper = numeric()))
  }
  
  # If only one interval, return it directly
  if (nrow(x) == 1) {
    return(x)
  }
  
  out_lower <- numeric(0)
  out_upper <- numeric(0)
  
  cur_l <- x$lower[1]
  cur_u <- x$upper[1]
  
  # Safe index sequence: 2..n only when n>=2
  for (i in seq(2, nrow(x))) {
    l <- x$lower[i]
    u <- x$upper[i]
    
    if (l <= cur_u + gap_ms) {
      cur_u <- max(cur_u, u)
    } else {
      out_lower <- c(out_lower, cur_l)
      out_upper <- c(out_upper, cur_u)
      cur_l <- l
      cur_u <- u
    }
  }
  
  out_lower <- c(out_lower, cur_l)
  out_upper <- c(out_upper, cur_u)
  
  tibble(lower = out_lower, upper = out_upper)
}


# Total covered length after merging (ms)
covered_length_ms <- function(df, gap_ms = 0) {
  merged <- merge_intervals(df, gap_ms = gap_ms)
  if (nrow(merged) == 0) return(0)
  sum(merged$upper - merged$lower)
}


# Intersect two interval sets: A and B are tibbles(lower, upper)
intersect_two_interval_sets <- function(A, B) {
  if (nrow(A) == 0 || nrow(B) == 0) {
    return(tibble(lower = numeric(), upper = numeric()))
  }
  
  A <- A %>% arrange(lower, upper)
  B <- B %>% arrange(lower, upper)
  
  i <- 1
  j <- 1
  out <- vector("list", 0)
  
  while (i <= nrow(A) && j <= nrow(B)) {
    a1 <- A$lower[i]; a2 <- A$upper[i]
    b1 <- B$lower[j]; b2 <- B$upper[j]
    
    lo <- max(a1, b1)
    hi <- min(a2, b2)
    
    if (hi > lo) {
      out[[length(out) + 1]] <- c(lo, hi)
    }
    
    # Move the pointer with smaller endpoint
    if (a2 < b2) {
      i <- i + 1
    } else {
      j <- j + 1
    }
  }
  
  if (length(out) == 0) {
    tibble(lower = numeric(), upper = numeric())
  } else {
    mat <- do.call(rbind, out)
    tibble(lower = mat[, 1], upper = mat[, 2])
  }
}

# Get triple-overlap interval set across multiple time_windows (default 10/20/40)
triple_overlap_intervals <- function(df_region,
                                     time_windows = c(10, 20, 40),
                                     gap_ms_merge_within_window = 0) {
  # df_region should contain: time_window, lower, upper, coef (non-zero selection)
  # Make per-window merged interval sets
  interval_sets <- lapply(time_windows, function(tw) {
    merge_intervals(
      df_region %>%
        filter(time_window == tw) %>%
        filter(!is.na(coef), coef != 0) %>%
        select(lower, upper),
      gap_ms = gap_ms_merge_within_window
    )
  })
  
  # If any window is empty -> no overlap
  if (any(vapply(interval_sets, nrow, integer(1)) == 0)) {
    return(tibble(lower = numeric(), upper = numeric()))
  }
  
  # Intersect iteratively: (((I10 ∩ I20) ∩ I40) ...)
  ov <- interval_sets[[1]]
  for (k in 2:length(interval_sets)) {
    ov <- intersect_two_interval_sets(ov, interval_sets[[k]])
    if (nrow(ov) == 0) break
  }
  
  ov %>% arrange(lower, upper)
}

# Check if there exists a "continuous overlap run" of length >= min_run_ms,
# allowing gaps up to allowed_gap_bins * base_bin_ms (base_bin_ms default=min time_window)
has_overlap_run <- function(overlap_df,
                            min_run_ms = 40,
                            allowed_gap_bins = 0,
                            base_bin_ms = 10) {
  if (nrow(overlap_df) == 0) return(FALSE)
  
  allowed_gap_ms <- allowed_gap_bins * base_bin_ms
  
  # Merge overlap intervals allowing small gaps (treat small gaps as continuous)
  merged_run <- merge_intervals(overlap_df, gap_ms = allowed_gap_ms)
  
  any((merged_run$upper - merged_run$lower) >= min_run_ms)
}


post_filter_sparse_and_overlap_run <- function(filtered_models,
                                               time_windows = c(10, 20, 40),
                                               max_window_ms = 100,
                                               min_overlap_run_ms = 40,
                                               allowed_gap_bins = 0,     # 0/1/2
                                               base_bin_ms = NULL,       # default=min(time_windows)
                                               gap_ms_merge_within_window = 0) {
  if (is.null(base_bin_ms)) base_bin_ms <- min(time_windows)
  
  regions <- sort(unique(filtered_models$region))
  
  out <- lapply(regions, function(rg) {
    df_rg <- filtered_models %>% filter(region == rg)
    
    # 1) per-window duration <= max_window_ms
    # Keep only results whose non-zero total duration (per time_window) <= max_window_ms
    len_ok <- purrr::map_lgl(time_windows, function(tw) {
      df_tw <- df_rg %>%
        filter(time_window == tw, !is.na(coef), coef != 0) %>%
        select(lower, upper)
      
      covered_length_ms(df_tw, gap_ms = 0) <= max_window_ms
    })
    
    if (!all(len_ok)) return(tibble())
    
    # passed -> keep original (unmodified)
    trimmed <- df_rg
    
    
    # 2) triple-overlap intervals among 10/20/40
    ov_df <- triple_overlap_intervals(
      df_region = trimmed,
      time_windows = time_windows,
      gap_ms_merge_within_window = gap_ms_merge_within_window
    )
    
    ok <- has_overlap_run(
      overlap_df = ov_df,
      min_run_ms = min_overlap_run_ms,
      allowed_gap_bins = allowed_gap_bins,
      base_bin_ms = base_bin_ms
    )
    
    if (!ok) return(tibble())
    
    # attach some meta info (optional)
    trimmed %>%
      mutate(
        overlap_run_required_ms = min_overlap_run_ms,
        allowed_gap_bins = allowed_gap_bins,
        base_bin_ms = base_bin_ms
      )
  })
  
  dplyr::bind_rows(out)
}

