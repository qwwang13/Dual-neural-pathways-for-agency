# Reproduce the manuscript's JZS Bayes factors and participant-level ITPC slopes.
# Run this script in RStudio with Source, or with Rscript --vanilla analysis/run.R robustness.
# All analysis inputs are included. Only base R is required.

# Read data ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
source_paths <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
script_path <- if (length(source_paths)) tail(source_paths, 1L)[[1L]] else
  if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else "07_robustness.R"
module_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
project_root <- normalizePath(file.path(module_dir, "../.."), mustWork = TRUE)
data_dir <- file.path(project_root, "data", "07_supplementary_robustness")
output_dir <- file.path(project_root, "results", "07_supplementary_robustness")
output_option <- match("--output", args)
if (!is.na(output_option)) {
  if (output_option == length(args)) stop("Missing directory after --output")
  output_dir <- args[output_option + 1L]
}
options(digits = 17)
checksums <- read.csv(file.path(data_dir, "input_checksums.csv"),
                      stringsAsFactors = FALSE)
# Normalize CRLF to LF so Git text checkout settings do not affect checks.
input_md5 <- function(path) {
  text <- rawToChar(readBin(path, "raw", n = file.info(path)$size))
  text <- gsub("\r\n", "\n", text, fixed = TRUE)
  canonical <- tempfile()
  on.exit(unlink(canonical))
  writeBin(charToRaw(text), canonical)
  unname(tools::md5sum(canonical))
}
actual_md5 <- unname(vapply(file.path(data_dir, checksums$file),
                            input_md5, character(1)))
stopifnot(!anyNA(actual_md5), identical(actual_md5, checksums$md5))
behavior <- read.csv(file.path(project_root,
  "data/01_temporal_binding_effect/perceptual_shift_action_outcome.csv"),
  check.names = FALSE)
behavior <- behavior[order(behavior$subject), ]
behavior_columns <- c("va_ba", "ia_bi", "vs_bs", "is_bs")
stopifnot(nrow(behavior) == 19L, !anyDuplicated(behavior$subject),
          identical(as.integer(behavior$subject), 2:20),
          all(is.finite(as.matrix(behavior[behavior_columns]))))
itpc <- read.csv(file.path(data_dir, "participant_itpc.csv"),
                 stringsAsFactors = FALSE)
stopifnot(nrow(itpc) == 152L, all(is.finite(itpc$itpc)),
          all(itpc$itpc >= 0 & itpc$itpc <= 1),
          !anyDuplicated(itpc[c("subject", "binding", "roi", "band", "condition")]))
draws <- as.matrix(read.csv(file.path(data_dir, "bootstrap_indices.csv"),
                            header = FALSE))
stopifnot(identical(dim(draws), c(20000L, 19L)),
          all(draws >= 1L & draws <= 19L), all(draws == as.integer(draws)))
storage.mode(draws) <- "integer"

# Run Bayesian t tests ----------------------------------------------------
# Paired tests are one-sample tests on within-participant differences.
# The alternative uses delta ~ Cauchy(0, sqrt(2)/2); H0 fixes delta = 0.
jzs_test <- function(values) {
  n <- length(values)
  df <- n - 1L
  prior_scale <- sqrt(2) / 2
  statistic <- mean(values) / (sd(values) / sqrt(n))
  null_kernel <- (1 + statistic^2 / df)^(-n / 2)
  mixture <- function(g) exp(
    -0.5 * log1p(n * g * prior_scale^2) -
      n / 2 * log1p(statistic^2 / (df * (1 + n * g * prior_scale^2))) -
      0.5 * log(2 * pi) - 1.5 * log(g) - 1 / (2 * g))
  bf10 <- integrate(mixture, 0, Inf, rel.tol = 1e-10,
                     subdivisions = 2000L)$value / null_kernel
  # A second integration over the standardized effect verifies BF10.
  kernel <- function(delta) suppressWarnings(
    dt(statistic, df = df, ncp = delta * sqrt(n))) *
    dcauchy(delta, scale = prior_scale)
  positive <- integrate(kernel, 0, Inf, rel.tol = 1e-10,
                         subdivisions = 2000L)$value
  negative <- integrate(kernel, -Inf, 0, rel.tol = 1e-10,
                         subdivisions = 2000L)$value
  null_density <- dt(statistic, df = df)
  direct_bf <- (positive + negative) / null_density
  error <- abs(bf10 - direct_bf)
  stopifnot(error < 1e-7 * max(1, bf10))
  data.frame(n = n, mean_ms = mean(values), sd_ms = sd(values),
             t = statistic, df = df, cohen_dz = mean(values) / sd(values),
             prior_cauchy_scale = prior_scale, BF10 = bf10,
             BF_positive0 = 2 * positive / null_density,
             BF_negative0 = 2 * negative / null_density,
             integration_difference = error)
}
comparisons <- c(as.list(behavior[behavior_columns]), list(
  outcome_voluntary_minus_involuntary = behavior$vs_bs - behavior$is_bs))
bayes <- do.call(rbind, lapply(names(comparisons), function(name) {
  cbind(comparison = name, jzs_test(comparisons[[name]]))
}))
rownames(bayes) <- NULL

# Run paired participant bootstrap ---------------------------------------
# Each condition has a separate OLS model with an estimated intercept.
# ITPC uses its original scale; the response is the participant binding in ms.
specs <- data.frame(
  binding = c("action", "outcome", "outcome", "outcome"),
  roi = c("PreSMA-L", "PreCG-R", "MTG-R", "MTG-R"),
  band = c("Theta", "Alpha", "Delta", "Theta"),
  lower_ms = c(220, 180, 780, 780), upper_ms = c(270, 240, 830, 830),
  voluntary_behavior = c("va_ba", "vs_bs", "vs_bs", "vs_bs"),
  involuntary_behavior = c("ia_bi", "is_bs", "is_bs", "is_bs"),
  stringsAsFactors = FALSE)

fit_condition <- function(spec, condition, response_column) {
  selected <- itpc[itpc$binding == spec$binding & itpc$roi == spec$roi &
                     itpc$band == spec$band & itpc$condition == condition, ]
  selected <- selected[order(selected$subject), ]
  stopifnot(nrow(selected) == 19L,
            identical(as.integer(selected$subject), 2:20),
            all(selected$lower_ms == spec$lower_ms),
            all(selected$upper_ms == spec$upper_ms))
  x <- selected$itpc
  y <- behavior[[response_column]][match(selected$subject, behavior$subject)]
  model <- lm(y ~ x)
  # The closed-form OLS slope is checked against lm() for the observed data.
  xc <- x - mean(x)
  yc <- y - mean(y)
  slope <- sum(xc * yc) / sum(xc^2)
  stopifnot(abs(slope - unname(coef(model)["x"])) < 1e-9)
  xx <- matrix(x[as.vector(draws)], nrow = nrow(draws))
  yy <- matrix(y[as.vector(draws)], nrow = nrow(draws))
  xx <- xx - rowMeans(xx)
  yy <- yy - rowMeans(yy)
  denominators <- rowSums(xx^2)
  stopifnot(all(denominators > 0))
  bootstrap_slopes <- rowSums(xx * yy) / denominators
  # Independently verify sampled fits, including the first and last draws.
  for (index in unique(c(1L, 2L, 100L, 1000L, nrow(draws)))) {
    sampled <- draws[index, ]
    independent_fit <- lm(y[sampled] ~ x[sampled])
    stopifnot(abs(bootstrap_slopes[index] - unname(coef(independent_fit)[2L])) < 1e-8)
  }
  list(intercept = unname(coef(model)[1L]), slope = slope,
       bootstrap_slopes = bootstrap_slopes)
}

bootstrap <- do.call(rbind, lapply(seq_len(nrow(specs)), function(index) {
  spec <- specs[index, ]
  voluntary <- fit_condition(spec, "Voluntary", spec$voluntary_behavior)
  involuntary <- fit_condition(spec, "Involuntary", spec$involuntary_behavior)
  differences <- voluntary$bootstrap_slopes - involuntary$bootstrap_slopes
  # R type 7 is the same linear percentile definition used by numpy.quantile.
  interval <- quantile(differences, c(0.025, 0.975), type = 7, names = FALSE)
  p <- min(1, 2 * min((1 + sum(differences <= 0)) / (nrow(draws) + 1),
                      (1 + sum(differences >= 0)) / (nrow(draws) + 1)))
  data.frame(binding = spec$binding, roi = spec$roi, band = spec$band,
    n = 19L, n_bootstrap = nrow(draws),
    voluntary_intercept = voluntary$intercept,
    involuntary_intercept = involuntary$intercept,
    voluntary_slope = voluntary$slope, involuntary_slope = involuntary$slope,
    difference_ms_per_unit_ITPC = voluntary$slope - involuntary$slope,
    ci_low = interval[1L], ci_high = interval[2L], bootstrap_p = p)
}))

# Check manuscript values and save results --------------------------------
# These numbers are verification targets only; they are not model inputs.
assert_display <- function(observed, reported) {
  stopifnot(identical(sprintf("%.3f", observed), sprintf("%.3f", reported)))
}
assert_display(bayes$BF10, c(0.826, 0.264, 258.260, 130.875, 4.883))
assert_display(bayes$BF_positive0[1L], 1.552)
assert_display(bootstrap$difference_ms_per_unit_ITPC,
               c(374.509, -621.004, 54.763, -641.156))
assert_display(bootstrap$ci_low,
               c(-134.163, -1650.805, -279.925, -1505.191))
assert_display(bootstrap$ci_high, c(869.900, 235.740, 280.104, 581.905))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(bayes, file.path(output_dir, "bayesian_t_tests.csv"), row.names = FALSE)
write.csv(bootstrap, file.path(output_dir, "itpc_bootstrap.csv"), row.names = FALSE)

labels <- c("Voluntary action vs zero", "Involuntary action vs zero",
            "Voluntary outcome vs zero", "Involuntary outcome vs zero",
            "Outcome: voluntary minus involuntary")
lines <- c("# Supplementary robustness results", "",
  "All manuscript Bayes factors, slope differences and percentile confidence intervals reproduce at three decimal places.",
  "", "## Bayesian t tests", "",
  "JZS prior: standardized effect ~ Cauchy(0, sqrt(2)/2). BF10 compares the two-sided alternative with zero; BF+0 and BF-0 use the corresponding truncated and renormalized priors.",
  "", "| Comparison | n | t | df | BF10 | BF+0 | BF-0 |",
  "|---|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(bayes))) lines <- c(lines, sprintf(
  "| %s | %d | %.3f | %d | %.3f | %.3f | %.3f |", labels[i],
  bayes$n[i], bayes$t[i], bayes$df[i], bayes$BF10[i],
  bayes$BF_positive0[i], bayes$BF_negative0[i]))
lines <- c(lines, "", "## Participant-level ITPC slope differences", "",
  "Each condition is fitted separately as binding_ms = intercept + slope * ITPC. Differences are voluntary minus involuntary, in ms per unit ITPC. Each of the 20,000 draws resamples participants together with both conditions. Intervals are 95% percentile intervals.",
  "", "| Region | Band | Voluntary slope | Involuntary slope | Difference | CI low | CI high | Bootstrap p |",
  "|---|---|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(bootstrap))) lines <- c(lines, sprintf(
  "| %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",
  bootstrap$roi[i], bootstrap$band[i], bootstrap$voluntary_slope[i],
  bootstrap$involuntary_slope[i], bootstrap$difference_ms_per_unit_ITPC[i],
  bootstrap$ci_low[i], bootstrap$ci_high[i], bootstrap$bootstrap_p[i]))
lines <- c(lines, "",
  "The supplementary bootstrap p value is twice the smaller zero-crossing tail proportion, with a plus-one correction and a cap of 1. The manuscript reports the percentile confidence intervals. All four intervals include zero.",
  "", "## Reproduction details", "",
  paste("R:", R.version.string),
  "R package dependencies: base and stats, supplied with R.",
  "Bootstrap draws: NumPy PCG64, seed 20261109, shape 20000 x 19; included as one-based participant-row positions.",
  sprintf("Largest absolute difference between the two Bayes-factor integrations: %.3g.",
          max(bayes$integration_difference)),
  "Observed and sampled regression slopes were also checked with lm(), including its intercept.",
  "Input MD5 checks (LF-normalized text):",
  paste0("- ", checksums$file, ": ", actual_md5), "")
writeLines(lines, file.path(output_dir, "RESULTS.md"), useBytes = TRUE)
cat("All manuscript numerical checks passed.\n")
print(bayes[, c("comparison", "BF10", "BF_positive0")], row.names = FALSE, digits = 7)
print(bootstrap[, c("roi", "band", "difference_ms_per_unit_ITPC", "ci_low", "ci_high")],
      row.names = FALSE, digits = 8)
cat("Results:", normalizePath(output_dir), "\n")
