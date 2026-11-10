# The JZS alternative uses delta ~ Cauchy(0, sqrt(2)/2); H0 fixes delta=0.
jzs_test <- function(values) {
  n <- length(values)
  df <- n-1L
  prior_scale <- sqrt(2)/2
  statistic <- mean(values)/(sd(values)/sqrt(n))
  null_kernel <- (1+statistic^2/df)^(-n/2)
  mixture <- function(g) exp(
    -0.5*log1p(n*g*prior_scale^2) -
      n/2*log1p(statistic^2/(df*(1+n*g*prior_scale^2))) -
      0.5*log(2*pi)-1.5*log(g)-1/(2*g))
  bf10 <- integrate(mixture, 0, Inf, rel.tol=1e-10,
                     subdivisions=2000L)$value/null_kernel
  kernel <- function(delta) suppressWarnings(
    dt(statistic, df=df, ncp=delta*sqrt(n)))*dcauchy(delta, scale=prior_scale)
  positive <- integrate(kernel, 0, Inf, rel.tol=1e-10, subdivisions=2000L)$value
  negative <- integrate(kernel, -Inf, 0, rel.tol=1e-10, subdivisions=2000L)$value
  null_density <- dt(statistic, df=df)
  data.frame(n=n, mean_ms=mean(values), sd_ms=sd(values), t=statistic, df=df,
    cohen_dz=mean(values)/sd(values), prior_cauchy_scale=prior_scale, BF10=bf10,
    BF_positive0=2*positive/null_density, BF_negative0=2*negative/null_density)
}

run_bayesian_tests <- function(behavior) {
  comparisons <- c(as.list(behavior[c("va_ba", "ia_bi", "vs_bs", "is_bs")]), list(
    outcome_voluntary_minus_involuntary=behavior$vs_bs-behavior$is_bs))
  result <- do.call(rbind, lapply(names(comparisons), function(name) {
    cbind(comparison=name, jzs_test(comparisons[[name]]))
  }))
  rownames(result) <- NULL
  result
}

