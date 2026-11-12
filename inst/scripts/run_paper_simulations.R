# Usage: Rscript run_paper_simulations.R <benchmark|full|restart|summary> [output-dir] [workers] [replications]
# Full runs and summaries default to paper-results in the working directory.
# Set thread limits before loading numerical libraries. Parallelism is across
# simulation replications, not inside each estimation routine.
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1",BLIS_NUM_THREADS="1")
args <- commandArgs(trailingOnly=TRUE)
# EDIT THESE SETTINGS before starting a new experiment. They are recorded in
# the manifest and must remain unchanged when resuming with mode="full".
estimator_args <- list(MLE=list(initial_cov="stationary",pilot="HR1"))
# To fix the initial state covariance at its stabilized pilot value:
# estimator_args <- list(MLE=list(initial_cov="fixed",pilot="HR1"))
if(length(args)<1L || length(args)>4L) stop(paste(
  "Usage: Rscript run_paper_simulations.R <benchmark|full|restart|summary> [output-dir] [workers] [replications]",
  "Example: Rscript run_paper_simulations.R benchmark timing-test 4 5",
  "Example: Rscript run_paper_simulations.R full paper-results 4 1000",sep="\n"))
library(varma)
if(!exists("run_varma_paper_simulation",asNamespace("varma"),inherits=FALSE))
  stop("Install the updated varma package supplied with this driver")
mode <- match.arg(args[1L],c("benchmark","full","restart","summary"))
if(mode != "summary" && utils::packageVersion("varma") < "0.1.3")
  stop("Install varma >= 0.1.3 for CML, HR1 starting values and selectable Kalman initialization")
output <- if(length(args)>=2L) args[2L] else if(mode=="benchmark") "timing-test" else "paper-results"
workers <- if(length(args)>=3L) as.integer(args[3L]) else 1L
repetitions <- if(length(args)>=4L) as.integer(args[4L]) else if(mode=="benchmark") 5L else 1000L
if(is.na(workers) || workers<1L || is.na(repetitions) || repetitions<1L)
  stop("workers and replications must be positive integers")
if(mode=="restart") {
  # Explicit restart preserves one previous run rather than mixing checkpoints.
  # Resolve both paths and verify that the move remains in the same parent.
  if(dir.exists(output)) {
    target <- normalizePath(output,winslash="/",mustWork=TRUE)
    parent <- normalizePath(dirname(target),winslash="/",mustWork=TRUE)
    if(identical(target,parent) || !file.exists(file.path(target,"manifest.rds")))
      stop("restart requires an existing simulation directory with manifest.rds")
    backup <- file.path(parent,paste0(basename(target),"-previous-",format(Sys.time(),"%Y%m%d-%H%M%S")))
    if(!identical(normalizePath(dirname(backup),winslash="/",mustWork=TRUE),parent) ||
       file.exists(backup) || dir.exists(backup)) stop("Invalid or existing backup path")
    if(!file.rename(target,backup)) stop("Could not preserve the previous experiment; close files and retry")
    message("Previous experiment preserved at ",backup)
  }
  mode <- "full"
}
if(mode=="benchmark") {
  result <- benchmark_varma_paper_simulation(output,pilot_nsim=repetitions,workers=workers,
    chunk_size=1L,estimator_args=estimator_args)
  print(result$cell_times)
  if(workers==1L) cat(sprintf("Projected serial duration: %.2f hours (pilot MCSE %.2f hours).\n",
    result$serial_seconds/3600,result$serial_mcse_seconds/3600))
  cat(sprintf("Projection from measured %d-worker throughput: %.2f hours.\n",
    result$workers,result$throughput_seconds/3600))
} else if(mode=="full") {
  result <- run_varma_paper_simulation(output,nsim=repetitions,workers=workers,estimator_args=estimator_args)
  cat("Completed. Checkpoints and summary CSVs:",result$output_dir,"\n")
} else {
  result <- summarize_varma_paper_simulation(output)
  cat("Summary CSVs refreshed in",normalizePath(output),"\n")
}
