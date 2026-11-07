suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lme4)
  library(lmerTest)
})

add_signif <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", "")))))
}

.tidy_fixef <- function(m, condition_raw) {
  coefs <- as.data.frame(summary(m)$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL

  coefs %>%
    transmute(
      condition = condition_raw,
      term = term,
      estimate = .data[["Estimate"]],
      std_error = .data[["Std. Error"]],
      df = .data[["df"]],
      t_value = .data[["t value"]],
      p_value = .data[["Pr(>|t|)"]],
      signif = add_signif(.data[["Pr(>|t|)"]])
    ) %>%
    mutate(
      ci_crit = qt(0.975, df),
      conf_low = estimate - ci_crit * std_error,
      conf_high = estimate + ci_crit * std_error
    ) %>%
    select(
      condition, term, estimate, std_error, df, t_value, p_value,
      conf_low, conf_high, signif
    )
}

run_lmm_split <- function(
    feature_csv,
    behavior_csv,
    scale_features = TRUE,
    condition_levels = NULL,
    exclude_subjects = NULL,
    exclude_cols = NULL,
    encoding = "UTF-8"
) {
  feat <- read_csv(feature_csv, locale = locale(encoding = encoding), show_col_types = FALSE)
  beh  <- read_csv(behavior_csv, locale = locale(encoding = encoding), show_col_types = FALSE)

  df <- feat %>%
    inner_join(beh, by = c("subject","condition","trial")) %>%
    filter(!is.na(behavior_value))

  if (!is.null(exclude_subjects)) {
    df <- df %>% filter(!as.character(subject) %in% as.character(exclude_subjects))
  }

  if (!is.null(condition_levels)) {
    df <- df %>%
      filter(condition %in% condition_levels) %>%
      mutate(condition = factor(condition, levels = condition_levels))
  } else {
    df <- df %>% mutate(condition = factor(condition))
    condition_levels <- levels(df$condition)
  }

  df <- df %>% mutate(subject = factor(subject))

  predictors <- setdiff(
    names(df),
    c("behavior_value", "subject", "condition", "trial", exclude_cols)
  )

  df_list <- split(df, df$condition)
  scaling_parameters <- NULL
  if (isTRUE(scale_features) && length(predictors) > 0) {
    scaling_parameters <- bind_rows(lapply(names(df_list), function(cond) {
      d <- df_list[[cond]]
      centers <- vapply(d[predictors], mean, numeric(1))
      scales <- vapply(d[predictors], sd, numeric(1))
      if (any(!is.finite(centers)) || any(!is.finite(scales) | scales <= 0)) {
        stop("Non-finite or zero-variance predictors in condition: ", cond)
      }
      tibble(condition = cond, term = predictors, center = unname(centers),
             scale = unname(scales), n_trials = nrow(d))
    }))
    df_list <- lapply(df_list, function(d) {
      d %>% mutate(across(all_of(predictors), ~ as.numeric(scale(.))))
    })
    df <- bind_rows(df_list)
  }

  fixed_formula <- reformulate(predictors, response = "behavior_value")
  full_formula  <- update(fixed_formula, . ~ . + (1 | subject))

  models <- lapply(df_list, function(d) lmer(full_formula, data = d, REML = FALSE))

  fixed_table <- bind_rows(lapply(names(models), function(cond) {
    .tidy_fixef(models[[cond]], cond)
  })) %>%
    select(condition, term, estimate, std_error, df, t_value, p_value, conf_low, conf_high, signif)

  fixed_table_wide <- fixed_table %>%
    select(condition, term, estimate, std_error, df, t_value, p_value, conf_low, conf_high, signif) %>%
    pivot_wider(
      names_from = condition,
      values_from = c(estimate, std_error, df, t_value, p_value, conf_low, conf_high, signif),
      names_sep = "_"
    )

  list(
    data = df,
    predictors = predictors,
    scaling_parameters = scaling_parameters,
    models = models,
    fixed_table = fixed_table,
    fixed_table_wide = fixed_table_wide
  )
}
