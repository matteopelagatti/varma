# Development benchmark; MTS is used only as the external reference.
# Run from the repository root. Times include each function's input processing.
Rcpp::sourceCpp("src/rcpp_functions.cpp",showOutput=FALSE,
  cacheDir=file.path(Sys.getenv("TEMP",tempdir()),"varma-regression-source-cache"))
for(file in c("estimation_helpers.R","functions.R","echelon.R")) source(file.path("R",file))
stopifnot(requireNamespace("MTS",quietly=TRUE))
bench<-function(fn) {
  invisible(fn());batch<-1L
  repeat {
    elapsed<-system.time(for(i in seq_len(batch)) invisible(fn()))[["elapsed"]]
    if(elapsed>=.15 || batch>=10000L) break
    batch<-min(10000L,2L*batch)
  }
  elapsed<-replicate(5,{
    gc(FALSE)
    system.time(for(i in seq_len(batch)) invisible(fn()))[["elapsed"]]/batch
  })
  c(median_seconds=median(elapsed),min_seconds=min(elapsed),max_seconds=max(elapsed),batch=batch)
}
designs<-list(c(2L,300L),c(3L,1000L),c(5L,2000L),c(8L,1500L))
results<-list()
for(j in seq_along(designs)) {
  m<-designs[[j]][1];n<-designs[[j]][2]
  set.seed(900+j)
  model<-varma(ar=array(diag(seq(.5,.7,length.out=m),m),c(m,m,1)),
                ma=array(diag(.15,m),c(m,m,1)),cov=diag(m))
  Y<-sim_varma(model,n)
  invisible(capture.output(reference<-MTS::Kronid(scale(Y))))
  internal<-select_kronecker_indices(Y)
  stopifnot(identical(as.integer(reference$index),as.integer(internal$index)))
  # MTS has no quiet flag; capture its console output, as the former wrapper did.
  old<-bench(function() invisible(capture.output(MTS::Kronid(scale(Y)))))
  new<-bench(function() select_kronecker_indices(Y))
  result<-data.frame(m=m,n=n,plag=5L,method=c("MTS::Kronid","internal"),
    rbind(old,new),indices=paste(internal$index,collapse=","),row.names=NULL)
  results[[j]]<-result
  print(result)
  cat("Speed ratio:",old[1]/new[1],"\n")
}
write.csv(do.call(rbind,results),"audit/2026-11-10/kronid/benchmark.csv",row.names=FALSE)
