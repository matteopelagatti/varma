# Uses the corrected, isolated package installation; run from the repo root.
.libPaths(c("C:/Users/Matteo/AppData/Local/Temp/varma-kronid-20261110/varma.Rcheck",.libPaths()))
library(varma)
source("R/paper_simulations.R")
models <- varma_paper_designs()
rows <- list()
for(design in c("model1","model4")) {
  set.seed(445);Y<-sim_varma(models[[design]],200)
  orders<-dim(models[[design]])
  for(backend in c("skf","fkf","kfas")) {
    t0<-proc.time()[["elapsed"]]
    f<-suppressWarnings(.varma_paper_fit_function("MLE",backend)(Y,p=orders[2],q=orders[3],
      intercept=FALSE,control=list(maxit=500,reltol=1e-8)))
    row<-data.frame(design=design,backend=backend,seconds=proc.time()[["elapsed"]]-t0,
      convergence=f$convergence,loglik=as.numeric(logLik(f)),objective_calls=f$optim.out$counts[1])
    rows[[length(rows)+1L]]<-row
    print(row)
    write.csv(do.call(rbind,rows),"audit/2026-11-10/paper-simulation/mle-backends.csv",row.names=FALSE)
  }
}
