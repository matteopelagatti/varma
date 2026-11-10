# Run from the repository root: Rscript benchmark.R <library> <before|after>
args <- commandArgs(TRUE)
.libPaths(c(args[1],"C:/Users/Matteo/AppData/Local/R/win-library/4.5",.libPaths()))
library(varma)
ns <- asNamespace("varma")
phase <- args[2]
out <- "audit/2026-11-10/performance"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
results <- list(); fitted <- list()
bench <- function(name,fn,iterations=1L,repetitions=5L) {
  invisible(suppressWarnings(fn()))
  # Windows elapsed timers are coarse. Amortise short calls over >=0.1 seconds.
  repeat {
    probe <- system.time(for (i in seq_len(iterations)) suppressWarnings(fn()))[["elapsed"]]
    if (probe >= .1 || iterations >= 10000L) break
    iterations <- min(10000L,iterations*2L)
  }
  times <- numeric(repetitions)
  for (j in seq_len(repetitions)) {
    gc(FALSE)
    times[j] <- system.time(for (i in seq_len(iterations)) suppressWarnings(fn()))[["elapsed"]]/iterations
  }
  results[[length(results)+1L]] <<- data.frame(task=name,phase=phase,
    median_seconds=median(times),min_seconds=min(times),max_seconds=max(times),
    iterations=iterations,repetitions=repetitions)
  cat(name,format(median(times),digits=5),"seconds\n")
}
set.seed(1437)
designs <- list(small=c(2,200,1,1),medium=c(3,400,2,2),larger=c(5,800,2,2))
for (name in names(designs)) {
  z<-designs[[name]];m<-z[1];n<-z[2];p<-z[3];q<-z[4]
  a<-array(0,c(m,m,p));b<-array(0,c(m,m,q))
  a[,,1]<-diag(.35,m)+matrix(rnorm(m*m,sd=.025),m)
  b[,,1]<-diag(.2,m)+matrix(rnorm(m*m,sd=.025),m)
  if (p>1) a[,,2]<-diag(.12,m)
  if (q>1) b[,,2]<-diag(-.08,m)
  model<-varma(ar=a,ma=b,cov=diag(m))
  y<-sim_varma(model,n,burn_in=200)
  bench(paste0("IHR_",name),function() fit_varma_ihr(y,p,q,FALSE,maxit=30),iterations=5)
  bench(paste0("DPF_",name),function() fit_varma_dpf(y,p,q,FALSE),iterations=5)
  if (m<=3) {
    for (backend in c("kfas","fkf","skf")) {
      fn<-get(paste0("fit_varma_",backend),ns)
      bench(paste0("ML_",backend,"_",name),function() fn(y,p,q,FALSE,maxit=10,
          control=list(maxit=10)),repetitions=3)
    }
    bench(paste0("RBY_",name),function() fit_varma_rby(y,p,q,FALSE,maxit=10),repetitions=3)
  }
  fitted[[name]]<-list(y=y,model=model,ihr=fit_varma_ihr(y,p,q,FALSE,maxit=30),
                       dpf=fit_varma_dpf(y,p,q,FALSE))
}
for (d in c(2,8,24,48)) {
  T<-matrix(0,d,d);diag(T)<-.95
  if(d>1) T[cbind(seq_len(d-1L),2:d)]<-.01
  C<-diag(d);sp<-Matrix::Matrix(T,sparse=TRUE)
  bench(paste0("Lyapunov_",d),function() get("solve_dlyap_iter",ns)(sp,C),iterations=100)
  if (exists("solve_discrete_lyapunov",ns,inherits=FALSE))
    bench(paste0("Lyapunov_dense_",d),function() solve_discrete_lyapunov(T,C),iterations=100)
}
# Compare existing Kronecker alternatives with the new state-space covariance route.
for (m in c(2,4,8)) {
  a<-array(0,c(m,m,2));a[,,1]<-diag(.4,m);a[,,2]<-diag(.15,m)
  b<-array(0,c(m,m,2));b[,,1]<-diag(.2,m)
  f<-varma(ar=a,ma=b,cov=diag(m))
  bench(paste0("autocov_vk_",m),function() autocov_vk(f,10),iterations=5)
  bench(paste0("autocov_sp_",m),function() autocov_sp(f,10),iterations=5)
  if(exists("autocov_ss",ns,inherits=FALSE))
    bench(paste0("autocov_ss_",m),function() autocov_ss(f,10),iterations=20)
}
write.csv(do.call(rbind,results),file.path(out,paste0(phase,".csv")),row.names=FALSE)
saveRDS(fitted,file.path(out,paste0(phase,"-fits.rds")))
capture.output(sessionInfo(),file=file.path(out,paste0(phase,"-session.txt")))
