args <- commandArgs(trailingOnly=TRUE)
file_arg <- grep("^--file=",commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=","",file_arg[1])),".."),winslash="/",mustWork=TRUE)
Sys.setenv(DUAL_AGENCY_ROOT=root)
source(file.path(root,"analysis/common.R"))
steps <- c(behavior="01_temporal_binding_effect/01_behavior.Rmd",meg="02_meg_results/02_timecourses.Rmd",
 functional="03_functional_regression/03_functional.Rmd",models="05_LMM/03_models.Rmd",
 coefficients="05_LMM/04_coefficients.Rmd",amplitude="05_LMM/05_amplitude.Rmd",ieeg="06_ieeg/02_ieeg.Rmd")
html <- "--html" %in% args
args <- args[args!="--html"]
step <- if(length(args))args[1]else"behavior"
stopifnot(step %in% names(steps))
folders <- c(behavior="01_temporal_binding_effect",meg="02_meg_results",functional="03_functional_regression",
 models="05_LMM/models",coefficients="05_LMM/coefficients",amplitude="05_LMM/mean_amplitude",ieeg="06_ieeg")
input <- file.path(root,"analysis",steps[[step]])
if(html) {
  rmarkdown::render(input,output_file=paste0(step,".html"),output_dir=out_dir(folders[[step]]),
    knit_root_dir=root,envir=new.env(parent=globalenv()),encoding="UTF-8")
} else {
  # Run the same report code without requiring an HTML renderer.
  setwd(root)
  grDevices::cairo_pdf(tempfile(fileext=".pdf"),family="Arial")
  lines <- readLines(input,encoding="UTF-8")
  e <- new.env(parent=globalenv())
  for(a in grep("^```\\{r",lines)) {
    b <- which(seq_along(lines)>a & grepl("^```\\s*$",lines))[1]
    eval(parse(text=lines[seq.int(a+1,b-1)]),envir=e)
  }
  grDevices::dev.off()
  message("Completed ",step,". Outputs: ",out_dir(folders[[step]]))
}
