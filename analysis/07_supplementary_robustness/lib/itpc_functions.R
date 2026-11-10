# Prepare one ITPC value per participant and condition from the loaded tables.
prepare_itpc_inputs <- function(behavior, features_by_roi) {
  behavior <- behavior[order(behavior$subject), ]
  subjects <- behavior$subject
  specs <- data.frame(
    binding=c("action", "outcome", "outcome", "outcome"),
    roi_key=c("PreSMAL", "PrecentralR", "MTGR", "MTGR"),
    roi=c("PreSMA-L", "PreCG-R", "MTG-R", "MTG-R"),
    band=c("Theta", "Alpha", "Delta", "Theta"),
    voluntary_behavior=c("va_ba", "vs_bs", "vs_bs", "vs_bs"),
    involuntary_behavior=c("ia_bi", "is_bs", "is_bs", "is_bs"),
    stringsAsFactors=FALSE)
  records <- lapply(seq_len(nrow(specs)), function(index) {
    spec <- specs[index, ]
    features <- features_by_roi[[spec$roi_key]]
    condition_keys <- if(spec$binding=="action") c("va_a", "ia_a") else c("vs_s", "is_s")
    do.call(rbind, lapply(seq_along(condition_keys), function(condition_index) {
      selected <- features[features$condition==condition_keys[condition_index], ]
      values <- selected[match(subjects, selected$subject), paste0(spec$band, "_bandITPC")]
      data.frame(subject=subjects, binding=spec$binding, roi=spec$roi, band=spec$band,
        condition=c("Voluntary", "Involuntary")[condition_index], itpc=values,
        stringsAsFactors=FALSE)
    }))
  })
  list(behavior=behavior, itpc=do.call(rbind, records), specs=specs)
}

# Generate the original PCG64 draws at runtime rather than storing an index table.
paired_bootstrap_draws <- function(generator, n, replicates, seed) {
  python <- Sys.getenv("DUAL_AGENCY_PYTHON", unset="")
  if(!nzchar(python)) python <- unname(Sys.which("python"))
  if(!nzchar(python)) stop("Set DUAL_AGENCY_PYTHON to a Python executable with NumPy.")
  output <- tempfile(fileext=".bin")
  on.exit(unlink(output))
  messages <- system2(python, c(shQuote(generator), "--output", shQuote(output),
    "--participants", as.character(n), "--replicates", as.character(replicates),
    "--seed", as.character(seed)), stdout=TRUE, stderr=TRUE)
  status <- attr(messages, "status")
  if(!is.null(status) && status!=0L) stop(paste(
    "Bootstrap draw generation failed. Check Python and the module's NumPy requirement.",
    paste(messages, collapse="\n"), sep="\n"))
  values <- readBin(output, integer(), n=n*replicates, size=4L, endian="little")
  matrix(values, nrow=replicates, ncol=n, byrow=TRUE)
}

fit_itpc_condition <- function(itpc, behavior, spec, condition, response, draws) {
  selected <- itpc[itpc$binding==spec$binding & itpc$roi==spec$roi &
                     itpc$band==spec$band & itpc$condition==condition, ]
  selected <- selected[order(selected$subject), ]
  x <- selected$itpc
  y <- behavior[[response]][match(selected$subject, behavior$subject)]
  model <- lm(y ~ x)
  xx <- matrix(x[as.vector(draws)], nrow=nrow(draws))
  yy <- matrix(y[as.vector(draws)], nrow=nrow(draws))
  xx <- xx-rowMeans(xx)
  yy <- yy-rowMeans(yy)
  denominator <- rowSums(xx^2)
  list(intercept=unname(coef(model)[1L]), slope=unname(coef(model)[2L]),
       bootstrap_slopes=rowSums(xx*yy)/denominator)
}

run_itpc_bootstrap <- function(inputs, draws) {
  do.call(rbind, lapply(seq_len(nrow(inputs$specs)), function(index) {
    spec <- inputs$specs[index, ]
    voluntary <- fit_itpc_condition(inputs$itpc, inputs$behavior, spec,
                                    "Voluntary", spec$voluntary_behavior, draws)
    involuntary <- fit_itpc_condition(inputs$itpc, inputs$behavior, spec,
                                      "Involuntary", spec$involuntary_behavior, draws)
    difference <- voluntary$bootstrap_slopes-involuntary$bootstrap_slopes
    interval <- quantile(difference, c(0.025, 0.975), type=7, names=FALSE)
    p <- min(1, 2*min((1+sum(difference<=0))/(nrow(draws)+1),
                       (1+sum(difference>=0))/(nrow(draws)+1)))
    data.frame(binding=spec$binding, roi=spec$roi, band=spec$band,
      n=ncol(draws), n_bootstrap=nrow(draws),
      voluntary_intercept=voluntary$intercept, involuntary_intercept=involuntary$intercept,
      voluntary_slope=voluntary$slope, involuntary_slope=involuntary$slope,
      difference_ms_per_unit_ITPC=voluntary$slope-involuntary$slope,
      ci_low=interval[1L], ci_high=interval[2L], bootstrap_p=p)
  }))
}
