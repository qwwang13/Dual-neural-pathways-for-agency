binding_read_inputs <- function(binding, cfg) {
  brain_file <- file.path(data_root,"03_functional_regression",paste0("meg_",binding,"_difference.csv"))
  behavior_file <- file.path(data_root,"03_functional_regression/perceptual_shift_difference.csv")
  brain <- read_csv(brain_file,show_col_types=FALSE) %>%
    filter(condition==cfg$condition, Time>=cfg$lower, Time<=cfg$upper) %>%
    mutate(region=make_region_label(region))
  behavior <- read_csv(behavior_file,show_col_types=FALSE) %>% select(subject,all_of(cfg$binding_col))
  check19(behavior$subject); stopifnot(nrow(behavior)==19, !anyDuplicated(behavior$subject),
    all(is.finite(behavior[[cfg$binding_col]])), all(is.finite(brain$Value)),
    !anyDuplicated(brain[c("subject","region","Time")]))
  for(rg in unique(brain$region)) check19(brain$subject[brain$region==rg])
  list(brain=brain,behavior=behavior)
}

binding_load_paths <- function(project_root, binding, cfg, alpha_grid, n_lambda, inputs=NULL) {
  stopifnot(length(alpha_grid)>0,all(is.finite(alpha_grid)),all(alpha_grid>=0),all(alpha_grid<=1),
    !anyDuplicated(alpha_grid),length(n_lambda)==1,is.finite(n_lambda),n_lambda>=2,n_lambda%%1==0)
  brain_file <- file.path(data_root,"03_functional_regression",paste0("meg_",binding,"_difference.csv"))
  behavior_file <- file.path(data_root,"03_functional_regression/perceptual_shift_difference.csv")
  hashes <- tools::md5sum(c(brain_file, behavior_file,
    file.path(project_root,"analysis/03_functional_regression/elasticnet_fit.R"),
    file.path(project_root,"analysis/03_functional_regression/functional_inputs.R")))
  stopifnot(!anyNA(hashes))
  key <- list(input_hash=hashes, glmnet_version=as.character(packageVersion("glmnet")),
              time_windows=cfg$time_windows, alpha=alpha_grid, nlambda=as.integer(n_lambda), intercept=TRUE)
  cache_dir <- file.path(project_root,"analysis/.cache/functional")
  dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  cache_file <- file.path(cache_dir,paste0(binding,"_model_paths.rds"))
  if(file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if(identical(cached$key,key)) return(list(models=cached$models, fit_version=key$glmnet_version,
      source=cache_file, input_hash=hashes, refitted=FALSE))
  }
  if(key$glmnet_version!="4.1.10") stop("Use the project glmnet 4.1.10 library to retain the verified fit version.")
  if(is.null(inputs)) inputs <- binding_read_inputs(binding,cfg)
  brain <- inputs$brain
  behavior <- inputs$behavior
  features <- prep_brain_features(brain,cfg$time_windows,cfg$lower,cfg$upper)
  for(g in split(features,interaction(features$region,features$time_window,drop=TRUE))) {
    stopifnot(nrow(g)==19*(cfg$upper-cfg$lower)/unique(g$time_window), !anyNA(g$Value))
  }
  models <- run_elasticnet_grid(features,behavior,sort(unique(brain$region)),cfg$time_windows,
    alpha_grid=alpha_grid,intercepts=TRUE,binding_col=cfg$binding_col,n_lambda=n_lambda)
  saveRDS(list(key=key,models=models),cache_file)
  list(models=models,fit_version=key$glmnet_version,source=cache_file,input_hash=hashes,refitted=TRUE)
}
