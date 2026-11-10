# Regression checks for estimation, likelihood, IRFs, and residual bootstrap.
# Runs against the installed package under R CMD check. For a source-only audit,
# set VARMA_TEST_SOURCE=1 and run from the repository root.
if (identical(Sys.getenv("VARMA_TEST_SOURCE"),"1")) {
  Rcpp::sourceCpp("src/rcpp_functions.cpp",showOutput=FALSE,
    cacheDir=file.path(Sys.getenv("TEMP",tempdir()),"varma-regression-source-cache"))
  source("R/estimation_helpers.R")
  source("R/functions.R")
  source("R/covariance.R")
  source("R/echelon.R")
  internal <- function(x) get(x,envir=.GlobalEnv)
} else {
  library(varma)
  internal <- function(x) get(x,envir=asNamespace("varma"))
}
close <- function(a,b,tol=1e-7) stopifnot(isTRUE(all.equal(unname(a),unname(b),tolerance=tol,check.attributes=FALSE)))
check <- function(name,expr) { force(expr); cat("PASS:",name,"\n") }
set.seed(20261109)
A <- array(matrix(c(.5,.15,-.2,.35),2,byrow=TRUE),c(2,2,1))
B <- array(matrix(c(.2,-.12,.08,.25),2,byrow=TRUE),c(2,2,1))
S <- matrix(c(1,.6,.6,2),2)
Y <- sim_varma(varma(ar=A,ma=B,cov=S),n=240,burn_in=200)
Yc <- sweep(Y,2,colMeans(Y),"-"); Ys <- sweep(Yc,2,c(8,-5),"+")

check("DP estimators are translation equivariant",{
  for (fn in list(fit_varma_dpd,fit_varma_dpf)) {
    a <- fn(Yc,1,1,TRUE,r=8); b <- fn(Ys,1,1,TRUE,r=8)
    close(a$ar,b$ar); close(a$ma,b$ma); close(a$cov,b$cov)
  }
})
check("one-series and pure-MA boundary orders",{
  for (fn in list(fit_varma_ihr,fit_varma_rby,fit_varma_dpd,fit_varma_dpf)) {
    a <- fn(Y,0,1,FALSE,r=8)
    stopifnot(identical(as.integer(dim(a)),c(2L,0L,1L)),all(is.finite(a$cov)))
    b <- fn(Y[,1,drop=FALSE],1,1,TRUE,r=8)
    stopifnot(identical(as.integer(dim(b)),c(1L,1L,1L)),is.finite(as.numeric(logLik(b))))
  }
  for (p in c(1,2)) {
    a <- fit_varma_ihr(Y[,1,drop=FALSE],p,0)
    stopifnot(identical(dim(predict(a,2)$cov),c(1L,1L,2L)))
    invisible(capture.output(print(a)))
  }
})
check("IHR stores innovations consistent with final coefficients",{
  for (iterations in c(1,50)) {
    f <- fit_varma_ihr(Y,1,1,FALSE,maxit=iterations,r=8)
    close(f$residuals,residuals(f)); stopifnot(nrow(f$residuals)==nrow(Y))
  }
  f <- fit_varma_ihr(Y,1,0,TRUE)
  E <- f$residuals[-1,,drop=FALSE]
  close(f$cov,crossprod(E)/nrow(E))
  close(f$loglik,internal(".varma_gaussian_ll")(E,f$cov))
  stopifnot(attr(logLik(f),"df")==9,nobs(f)==239,inherits(logLik(f),"logLik"))
  close(AIC(f),-2*f$loglik+2*f$npar)
})
check("NET honours fixed lambda, SUR refreshes residuals, and scalar input works",{
  a <- fit_varma_net(Y,1,1,FALSE,lambda=1e-6,relax=FALSE,sur=FALSE,r=8,maxit=5)
  b <- fit_varma_net(Y,1,1,FALSE,lambda=1e6,relax=FALSE,sur=FALSE,r=8,maxit=5)
  stopifnot(max(abs(c(a$ar,a$ma)))>.01,max(abs(c(b$ar,b$ma)))<1e-8)
  for (relax in c(FALSE,TRUE)) {
    f <- fit_varma_net(Y,1,1,FALSE,lambda=.05,relax=relax,sur=TRUE,r=8,maxit=5)
    close(f$residuals,residuals(f)); close(f$cov,crossprod(f$residuals[-1,,drop=FALSE])/239)
  }
  f <- fit_varma_net(Y[,1,drop=FALSE],1,0,lambda=.01,relax=FALSE)
  stopifnot(all(is.finite(f$ar)))
  # An equation with no selected regressors must not break the SUR routine.
  beta <- internal("sur_cpp")(list(matrix(0,20,0),matrix(rnorm(40),20,2)),
                               list(matrix(rnorm(20)),matrix(rnorm(20))),diag(2))
  stopifnot(length(beta)==2L)
})
check("RBY score agrees with finite differences of concentrated Gaussian likelihood",{
  aa <- A; bb <- B; cc <- c(.1,-.2)
  step <- internal("rby_optimization_step")(Y,aa,bb,1L,1L,TRUE,cc)
  obj <- function(a,b,c) {
    E <- internal("varma_residuals")(Y,a,b,c)[-1,,drop=FALSE]
    internal(".varma_gaussian_ll")(E,crossprod(E)/nrow(E))
  }
  numerical <- matrix(0,5,2); h <- 1e-6
  for (eq in 1:2) {
    for (j in 1:2) {
      ap<-am<-aa; ap[eq,j,1]<-ap[eq,j,1]+h; am[eq,j,1]<-am[eq,j,1]-h
      numerical[j,eq]<-(obj(ap,bb,cc)-obj(am,bb,cc))/(2*h)
      bp<-bm<-bb; bp[eq,j,1]<-bp[eq,j,1]+h; bm[eq,j,1]<-bm[eq,j,1]-h
      numerical[2+j,eq]<-(obj(aa,bp,cc)-obj(aa,bm,cc))/(2*h)
    }
    cp<-cm<-cc;cp[eq]<-cp[eq]+h;cm[eq]<-cm[eq]-h
    numerical[5,eq]<-(obj(aa,bb,cp)-obj(aa,bb,cm))/(2*h)
  }
  close(as.numeric(step$ZE),as.numeric(numerical),tol=1e-6)
  f <- fit_varma_rby(Y,1,1,TRUE,r=8,maxit=50)
  stopifnot(all(diff(f$loglik_history)>0))
  close(f$loglik,tail(f$loglik_history,1))
})
check("Kalman acceleration waits for full Riccati convergence",{
  bb <- array(0,c(2,2,2));bb[,,2]<-diag(c(.8,.6))
  f <- varma(ma=bb,cov=S);ss<-varma_to_ss(f)
  TT<-Matrix::Matrix(ss$T,sparse=TRUE);RR<-Matrix::Matrix(ss$R,sparse=TRUE)
  P<-internal("solve_dlyap_iter")(TT,ss$R%*%S%*%t(ss$R))
  set.seed(392); yy<-t(sim_varma(f,n=300,burn_in=200))
  fast<-internal("kalmanLogLik")(TT,RR,S,matrix(0,nrow(TT),1),P,yy)
  exact<-internal("kalmanLogLik")(TT,RR,S,matrix(0,nrow(TT),1),P,yy,ss_tol=0)
  close(fast,exact,tol=1e-8)
  Pnear <- internal("solve_dlyap_iter")(Matrix::Matrix(matrix(.999,1,1),sparse=TRUE),matrix(1,1,1))
  close(as.numeric(Pnear),1/(1-.999^2),tol=1e-8)
})
check("likelihood backends return matching coefficients, covariances, and stationary initialisation",{
  # Pure VAR has a well-behaved likelihood and provides a shared backend benchmark.
  fits <- lapply(list(fit_varma_kfas,fit_varma_fkf,fit_varma_skf),function(fn)
    fn(Y,1,0,TRUE,control=list(maxit=300,reltol=1e-10)))
  for (f in fits) {
    ss<-varma_to_ss(f)
    close(f$transition_matrix,ss$T)
    close(f$disturbance_cov,ss$R%*%f$cov%*%t(ss$R))
    close(f$initial_state_cov,ss$T%*%f$initial_state_cov%*%t(ss$T)+f$disturbance_cov)
    stopifnot(f$n==240,nobs(f)==240,f$convergence==0L,!is.null(f[["npar",exact=TRUE]]))
    # Joint one-step errors, not sequentially orthogonalised component errors.
    manual <- sweep(Y[-1,,drop=FALSE],2,f$mean,"-")-
      sweep(Y[-nrow(Y),,drop=FALSE],2,f$mean,"-")%*%t(f$ar[,,1])
    cat("  ",f$estimation_method,"joint innovation max error:",max(abs(f$residuals[-1,,drop=FALSE]-manual)),"\n")
    close(as.numeric(f$residuals[-1,,drop=FALSE]),as.numeric(manual),tol=1e-6)
    invisible(capture.output(print(f)))
    stopifnot(identical(dim(predict(f,1)$cov),c(2L,2L,1L)))
  }
  close(fits[[1]]$loglik,fits[[2]]$loglik,1e-6)
  close(fits[[1]]$loglik,fits[[3]]$loglik,1e-6)
  close(predict(fits[[1]],3)$mean,predict(fits[[3]],3)$mean,1e-4)
  for (fn in list(fit_varma_kfas,fit_varma_fkf,fit_varma_skf)) {
    f <- fn(Y[,1,drop=FALSE],0,1,TRUE,control=list(maxit=300))
    stopifnot(identical(dim(predict(f,2)$cov),c(1L,1L,2L)))
  }
  # Hold coefficients at a common pilot to test higher-order packing independently
  # of optimiser paths and local minima.
  higher <- lapply(list(fit_varma_kfas,fit_varma_fkf,fit_varma_skf),function(fn)
    suppressWarnings(fn(Y,2,2,FALSE,control=list(maxit=0))))
  for (f in higher) close(f$transition_matrix,varma_to_ss(f)$T)
  close(higher[[1]]$loglik,higher[[2]]$loglik,1e-7)
  close(higher[[1]]$loglik,higher[[3]]$loglik,1e-7)
})
check("IRFs, prediction covariance, and identification",{
  f<-varma(ar=A,ma=B,cov=S)
  close(irf(f,0,"generalized")[,,1],S%*%diag(1/sqrt(diag(S))))
  psi<-irf(f,4,"none"); expected<-array(0,dim(psi));expected[,,1]<-diag(2)
  expected[,,2]<-A[,,1]+B[,,1]
  for (i in 3:5) expected[,,i]<-A[,,1]%*%expected[,,i-1]
  close(psi,expected)
  # A and B share a scalar root, but Psi1 is invertible: A=Psi2/Psi1 is unique.
  identified<-varma(ar=array(diag(c(.4,.2)),c(2,2,1)),ma=array(diag(c(.3,-.4)),c(2,2,1)),cov=diag(2))
  stopifnot(is_identified(identified))
  cancelled<-varma(ar=A,ma=-A,cov=S)
  stopifnot(!is_identified(cancelled))
  v<-varma(ar=A,cov=S);w<-varma(ar=A,ma=array(0,c(2,2,1)),cov=S)
  close(varma_to_ss(v,"pearlman")$R,varma_to_ss(w,"pearlman")$R)
  a<-array(0,c(2,2,2));a[1,1,1]<-.3;a[2,2,2]<-.4
  stopifnot(identical(as.integer(kronecker_indices(varma(ar=a,cov=diag(2)))),c(1L,2L)))
  delayed<-array(0,c(2,2,4));delayed[,,4]<-diag(c(.2,.3))
  stopifnot(identical(as.integer(kronecker_indices(varma(ma=delayed,cov=diag(2)))),c(4L,4L)))
})
check("bootstrap retains the fitted mean, centres errors, and skip=0 loses no rows",{
  f<-fit_varma_ihr(Ys,1,0,TRUE)
  capture<-function(Y,...) Y
  set.seed(510); yy<-bootstrap_varma(f,nsim=1,n=5000,fit_fn=capture)[[1]]
  stopifnot(max(abs(colMeans(yy)-internal(".varma_process_mean")(f)))<.2)
  # White noise gives an exact check of centring and inclusion of the first row.
  w<-varma(cov=diag(2),n=20,residuals=cbind(c(100,rep(0,19)),rep(0,20)))
  set.seed(93); zz<-bootstrap_varma(w,nsim=1,n=1000,resids_to_skip=0,burn_in=0,fit_fn=capture)[[1]]
  stopifnot(any(zz[,1]==95),any(zz[,1]==-5))
})
check("order selection handles white noise, common samples, and failed candidates",{
  tab<-select_varma_order(Y,p_max=1,q_max=1,ic="AIC",r=8,maxit=5)
  stopifnot(nrow(tab)==8,all(is.finite(tab$AIC)))
  f<-fit_varma_ihr(Y,1,0,FALSE)
  close(tab$AIC[tab$p==1 & tab$q==0 & !tab$intercept],AIC(f))
  broken <- function(Y,p,q,intercept,...) {
    if (q) stop("intentional failed candidate")
    fit_varma_ihr(Y,p,q,intercept,...)
  }
  tab <- suppressWarnings(select_varma_order(Y,1,1,FALSE,fit_fn=broken))
  stopifnot(sum(is.na(tab$AIC))==2L)
})
if (!identical(Sys.getenv("VARMA_TEST_SOURCE"),"1")) {
  check("parallel bootstrap restores the caller's future plan",{
    old <- future::plan()
    tryCatch({
      future::plan(future::multisession,workers=2L)
      f <- fit_varma_ihr(Y,1,0,FALSE)
      bb <- bootstrap_varma(f,nsim=2,parallel=TRUE,n_cores=1L)
      stopifnot(length(bb)==2L,future::nbrOfWorkers()==2L)
      net <- fit_varma_net(Y,1,0,FALSE,lambda=.05,relax=FALSE,
                          parallel=TRUE,n_cores=2L)
      stopifnot(all(is.finite(net$cov)),future::nbrOfWorkers()==2L)
    },finally=future::plan(old))
  })
}
cat("All estimation regression checks passed.\n")
