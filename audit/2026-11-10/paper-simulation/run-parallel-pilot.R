.libPaths(c("C:/Users/Matteo/AppData/Local/Temp/varma-paper-20261110/varma.Rcheck",.libPaths()))
library(varma)
started<-Sys.time()
result<-benchmark_varma_paper_simulation("audit/2026-11-10/paper-simulation/pilot-parallel",
  pilot_nsim=10L,full_nsim=1000L,workers=4L,chunk_size=2L)
print(result)
cat("Complete benchmark wall seconds:",as.numeric(difftime(Sys.time(),started,units="secs")),"\n")
