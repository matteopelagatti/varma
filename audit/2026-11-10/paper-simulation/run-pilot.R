# Run from the repo root. The native estimation code is the last checked build;
# the new simulation functions are sourced while their package build is pending.
.libPaths(c("C:/Users/Matteo/AppData/Local/Temp/varma-kronid-20261110/varma.Rcheck",.libPaths()))
library(varma)
source("R/paper_simulations.R")
for(nm in c("kalman",".varma_process_mean",".varma_integer"))
  assign(nm,get(nm,asNamespace("varma")))
started<-Sys.time()
result<-benchmark_varma_paper_simulation("audit/2026-11-10/paper-simulation/pilot-sequential",
  pilot_nsim=10L,full_nsim=1000L,workers=1L,chunk_size=2L)
print(result)
cat("Complete benchmark wall seconds:",as.numeric(difftime(Sys.time(),started,units="secs")),"\n")
