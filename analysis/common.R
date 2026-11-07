project_root <- Sys.getenv("DUAL_AGENCY_ROOT", unset=getwd())
while(!file.exists(file.path(project_root,"dual-agency-pathways.Rproj"))) {
  parent <- dirname(project_root)
  if(identical(parent,project_root)) stop("Project root not found.")
  project_root <- parent
}
extra_library <- Sys.getenv("DUAL_AGENCY_R_LIBRARY")
if(nzchar(extra_library)) .libPaths(c(extra_library,.libPaths()))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr);library(ggplot2)})
if(.Platform$OS.type=="windows") grDevices::windowsFonts(Arial=grDevices::windowsFont("Arial"))
data_root <- file.path(project_root,"data")
result_root <- file.path(project_root,"results")
figure_root <- file.path(project_root,"figures")
cache_root <- file.path(project_root,"analysis/.cache/derived")
report_end_ms <- 850
manuscript_rois <- c("PreSMAL","PrecentralR","MTGR")
# Short names are used only in figure text.
roi_display <- function(x) {
  aliases <- c("Precentral-R"="PreCG-R","PrecentralR"="PreCG-R",
    "Postcentral-R"="PoCG-R","PostcentralR"="PoCG-R",
    "Postcentral-L"="PoCG-L","PostcentralL"="PoCG-L")
  labels <- as.character(x)
  matched <- labels %in% names(aliases)
  labels[matched] <- unname(aliases[labels[matched]])
  labels
}
check19 <- function(x) stopifnot(setequal(unique(x),2:20),!anyNA(x))
out_dir <- function(part,figures=FALSE) {
  p <- file.path(if(figures)figure_root else result_root,part)
  dir.create(p,recursive=TRUE,showWarnings=FALSE);p
}
knitr::opts_knit$set(root.dir=project_root)
knitr::opts_chunk$set(message=FALSE,warning=TRUE,error=FALSE)

# Read the selected analysis windows from the coefficient exports.
read_functional_windows <- function() {
  bind_rows(lapply(c("action","outcome"),function(b) {
    read_csv(file.path(result_root,"03_functional_regression",paste0(b,"_binding.csv")),
      show_col_types=FALSE) %>%
      distinct(region,primary_lower_ms,primary_upper_ms) %>%
      mutate(binding=b,.before=1)
  })) %>%
    mutate(primary_window=paste0(primary_lower_ms,"–",primary_upper_ms)) %>%
    arrange(binding,primary_lower_ms,region)
}

cache_dir <- function(part) {
  p <- file.path(cache_root,part)
  dir.create(p,recursive=TRUE,showWarnings=FALSE);p
}
