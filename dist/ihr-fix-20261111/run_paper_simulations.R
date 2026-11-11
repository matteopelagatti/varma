# Usage: Rscript run_paper_simulations.R <benchmark|full|summary> <output-dir> [workers] [replications]
# Set thread limits before loading numerical libraries. Parallelism is across
# simulation replications, not inside each estimation routine.
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1",BLIS_NUM_THREADS="1")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2L) stop(paste(
  "Usage: Rscript run_paper_simulations.R <benchmark|full|summary> <output-dir> [workers] [replications]",
  "Example: Rscript run_paper_simulations.R benchmark timing-test 4 5",
  "Example: Rscript run_paper_simulations.R full paper-results 4 1000",sep="\n"))
library(varma)
if(!exists("run_varma_paper_simulation",asNamespace("varma"),inherits=FALSE))
  stop("Install the updated varma package supplied with this driver")
mode <- match.arg(args[1L],c("benchmark","full","summary"))
output <- args[2L]
workers <- if(length(args)>=3L) as.integer(args[3L]) else 1L
repetitions <- if(length(args)>=4L) as.integer(args[4L]) else if(mode=="benchmark") 5L else 1000L
if(mode=="benchmark") {
  result <- benchmark_varma_paper_simulation(output,pilot_nsim=repetitions,workers=workers,chunk_size=1L)
  print(result$cell_times)
  if(workers==1L) cat(sprintf("Projected serial duration: %.2f hours (pilot MCSE %.2f hours).\n",
    result$serial_seconds/3600,result$serial_mcse_seconds/3600))
  cat(sprintf("Projection from measured %d-worker throughput: %.2f hours.\n",
    result$workers,result$throughput_seconds/3600))
} else if(mode=="full") {
  result <- run_varma_paper_simulation(output,nsim=repetitions,workers=workers)
  cat("Completed. Checkpoints and summary CSVs:",result$output_dir,"\n")
} else {
  result <- summarize_varma_paper_simulation(output)
  cat("Summary CSVs refreshed in",normalizePath(output),"\n")
}
