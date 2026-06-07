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

# ---- internal: extract fixed effect table from one model ----
.tidy_fixef <- function(m, condition_raw) {
  coefs <- as.data.frame(summary(m)$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  
  # lmerTest: columns are typically Estimate / Std. Error / df / t value / Pr(>|t|)
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


# =========================================
# Main API: run LMMs
# =========================================
run_lmm_split <- function(
    feature_csv,
    behavior_csv,
    scale_features = TRUE,
    condition_levels = NULL,      # e.g. c("ia_a","va_a")
    exclude_subjects = NULL,      # e.g. c(19)
    exclude_cols = NULL,          # columns to exclude from predictors
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
  
  if (isTRUE(scale_features) && length(predictors) > 0) {
    df <- df %>% mutate(across(all_of(predictors), ~ as.numeric(scale(.))))
  }
  
  # split by condition -> fit LMM
  fixed_formula <- reformulate(predictors, response = "behavior_value")
  full_formula  <- update(fixed_formula, . ~ . + (1 | subject))
  
  df_list <- split(df, df$condition)
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
    models = models,
    fixed_table = fixed_table,
    fixed_table_wide = fixed_table_wide
  )
}
