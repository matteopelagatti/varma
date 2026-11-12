if(identical(Sys.getenv("VARMA_TEST_SOURCE"),"1")) {
  Rcpp::sourceCpp("src/rcpp_functions.cpp",showOutput=FALSE,
    cacheDir=file.path(Sys.getenv("TEMP",tempdir()),"varma-regression-source-cache"))
  for(file in c("estimation_helpers.R","functions.R","cml.R","echelon.R")) source(file.path("R",file))
} else library(varma)
check <- function(name,expr) {force(expr);cat("PASS:",name,"\n")}
close <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr,pattern) {
  result<-tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(result,"error"),grepl(pattern,conditionMessage(result),ignore.case=TRUE))
}
fixture_path<-if(file.exists("fixtures/kronid-reference.rds")) "fixtures/kronid-reference.rds" else "tests/fixtures/kronid-reference.rds"
fixtures<-readRDS(fixture_path)$cases

# Independent, deliberately simple oracle: rebuild raw matrices and use base R
# cancor and acf for every test. No compiled basis updates or MTS calls.
reference_selector <- function(Y,plag=5,crit=.05) {
  Y<-as.matrix(scale(Y));n<-nrow(Y);m<-ncol(Y)
  indices<-rep(NA_integer_,m);retained<-matrix(integer(),0,2)
  values<-list();meta<-list()
  for(h in 0:(m*plag)) {
    dates<-seq.int(plag+1L,n-h)
    past<-do.call(cbind,lapply(seq_len(plag),function(lag) Y[dates-lag,,drop=FALSE]))
    for(component in which(is.na(indices))) {
      candidate<-rbind(retained,c(h,component))
      future<-vapply(seq_len(nrow(candidate)),function(j)
        Y[dates+candidate[j,1],candidate[j,2]],numeric(length(dates)))
      cc<-stats::cancor(past,future)
      df<-ncol(future);dp<-ncol(past);rho2<-cc$cor[df]^2;correction<-1
      if(h>0) {
        sx<-past%*%cc$xcoef[,df];sy<-future%*%cc$ycoef[,df]
        ax<-as.numeric(stats::acf(sx,lag.max=h,plot=FALSE)$acf)[-1]
        ay<-as.numeric(stats::acf(sy,lag.max=h,plot=FALSE)$acf)[-1]
        correction<-1+2*sum(ax*ay)
      }
      statistic<--(n-1-(dp+df-1)/2)*log1p(-rho2/correction)
      p<-stats::pchisq(statistic,dp-df+1,lower.tail=FALSE)
      values[[length(values)+1L]]<-c(statistic,dp-df+1,p)
      meta[[length(meta)+1L]]<-c(h,component,rho2,correction)
      if(p>crit) indices[component]<-h else retained<-candidate
    }
    if(all(!is.na(indices))) break
  }
  list(index=indices,tests=do.call(cbind,values),details=do.call(rbind,meta))
}

check("internal selection agrees with saved MTS numerical references",{
  for(case in fixtures) {
    z<-select_kronecker_indices(case$Y,case$plag,case$crit)
    stopifnot(identical(as.integer(z$index),as.integer(case$index)))
    close(z$tests,case$tests)
  }
})

check("cached QR and incremental basis agree with independent cancor tests",{
  for(seed in 101:108) {
    set.seed(seed);m<-2L+seed%%2L
    f<-varma(ar=array(diag(seq(.5,.7,length.out=m),m),c(m,m,1)),
              ma=array(diag(.15,m),c(m,m,1)),cov=diag(m))
    Y<-sim_varma(f,500)
    z<-select_kronecker_indices(Y,plag=4)
    r<-reference_selector(Y,plag=4)
    stopifnot(identical(z$index,r$index));close(z$tests,r$tests)
    close(as.matrix(z$details[,c("horizon","component","canonical_sq","correction")]),r$details)
  }
  # An AR(2) needs more than one future horizon.
  set.seed(611)
  f<-varma(ar=array(c(0,.8),c(1,1,2)),cov=matrix(1))
  Y<-sim_varma(f,1200)
  z<-select_kronecker_indices(Y)
  r<-reference_selector(Y)
  stopifnot(z$index==2L,max(z$details$horizon)>=2)
  close(z$tests,r$tests)
})

check("a zero index removes its future column before the next horizon",{
  Y<-fixtures$zero_first$Y[,2:1]
  z<-select_kronecker_indices(Y)
  r<-reference_selector(Y)
  stopifnot(identical(z$index,c(1L,0L)),
            z$details$future_dim[z$details$horizon==1L][1L]==2L)
  close(z$tests,r$tests)
  # This ordering exercises a former dimension-dropping case as well.
  close(select_kronecker_indices(Y[,2:1])$tests,reference_selector(Y[,2:1])$tests)
})

check("univariate inputs, white noise, dimensions, and silent output",{
  Y<-fixtures$dynamic2$Y[,1]
  z<-select_kronecker_indices(Y)
  stopifnot(identical(z$index,1L),is.matrix(z$tests),nrow(z$tests)==3L)
  close(z$tests,reference_selector(matrix(Y))$tests)
  fit<-fit_varma_echelon(Y,control=list(maxit=500))
  stopifnot(fit$indices_estimated,identical(fit$kronecker_indices,1L),fit$convergence==0L)
  stopifnot(identical(dim(predict(fit,2)$cov),c(1L,1L,2L)))
  Y<-fixtures$white2$Y
  stopifnot(all(select_kronecker_indices(Y,crit=.01)$index==0L))
  stopifnot(length(capture.output(invisible(select_kronecker_indices(Y))))==0L)
  stopifnot(length(capture.output(select_kronecker_indices(Y,verbose=TRUE)))>0L)
})

check("selection preserves column names, units, locations, and the random seed",{
  Y<-fixtures$dynamic2$Y
  colnames(Y)<-c("output","inflation")
  set.seed(1);saved_seed<-.Random.seed
  z<-select_kronecker_indices(Y)
  stopifnot(identical(.Random.seed,saved_seed),identical(names(z$index),colnames(Y)))
  transformed<-sweep(sweep(Y,2,c(1e-5,-1e5),"*"),2,c(1,-20),"+")
  close(select_kronecker_indices(transformed)$tests,z$tests,1e-8)
  close(select_kronecker_indices(Y*1e200)$tests,z$tests,1e-8)
  close(select_kronecker_indices(Y*1e-200)$tests,z$tests,1e-8)
})

check("invalid and unsupported tests fail explicitly",{
  Y<-fixtures$dynamic2$Y
  fails(select_kronecker_indices(Y,plag=0),"plag")
  fails(select_kronecker_indices(Y,plag=1.5),"plag")
  fails(select_kronecker_indices(Y,plag=1e100),"plag")
  fails(select_kronecker_indices(Y,crit=0),"crit")
  fails(select_kronecker_indices(Y,crit=NA_real_),"crit")
  fails(select_kronecker_indices(Y,tol=0),"tol")
  fails(select_kronecker_indices(Y,verbose=NA),"verbose")
  fails(select_kronecker_indices(Y+1i),"real numeric")
  fails(select_kronecker_indices(matrix(NA_real_,100,2)),"finite")
  fails(select_kronecker_indices(Y[1:10,]),"Too few")
  fails(select_kronecker_indices(cbind(Y[,1],1)),"constant")
  fails(select_kronecker_indices(cbind(Y[,1],Y[,1])),"rank deficient")
  fails(select_kronecker_indices(Y,plag=1),"Past dimension exhausted")
})
cat("All internal Kronecker-index selection checks passed.\n")
