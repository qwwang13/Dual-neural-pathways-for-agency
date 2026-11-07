suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(glmnet)
})

make_region_label <- function(region_raw) {
  paste0(
    substr(region_raw, 1, nchar(region_raw) - 1),
    "-",
    substr(region_raw, nchar(region_raw), nchar(region_raw))
  )
}

assign_time_bin <- function(time_ms, n_bins, lower_ms, upper_ms) {
  breaks <- seq(lower_ms, upper_ms, length.out = n_bins + 1)
  k <- findInterval(time_ms, breaks, rightmost.closed = TRUE)
  k[time_ms < lower_ms | time_ms > upper_ms] <- NA_integer_
  k
}

bin_time_bounds <- function(p_number, time_window_ms, lower0_ms) {
  lower <- (p_number - 1) * time_window_ms + lower0_ms
  upper <- p_number * time_window_ms + lower0_ms
  tibble(lower = lower, upper = upper)
}

r2_by_lambda <- function(fit_values, y_true) {
  sst <- sum((y_true - mean(y_true))^2)
  apply(fit_values, 2, function(pred) {
    sse <- sum((y_true - pred)^2)
    1 - sse / sst
  })
}

has_consecutive_strict <- function(pos_vec, p) {
  if (length(pos_vec) < p) return(FALSE)
  diffs <- diff(sort(pos_vec))
  r <- rle(diffs)
  any(r$values == 1 & r$lengths >= (p - 1))
}

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
  fitted_vals <- predict(fit, newx = x, type = "response")
  r2 <- r2_by_lambda(fitted_vals, y_true = y)

  coef_mat <- as.matrix(coef(fit))

  map_dfr(seq_along(lambdas), function(i) {
    lam <- lambdas[i]
    beta <- coef_mat[-1, i]
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
