library(varma)
internal <- function(name) get(name,asNamespace("varma"))
check <- function(name,expr) { force(expr); cat("PASS:",name,"\n") }
close <- function(a,b,tol=1e-8) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),"error"))

check("full CML and unrestricted echelon coincide in the identifiable paper designs",{
  models <- varma_paper_designs()
  for(i in c(2L,4L)) {
    set.seed(610+i); Y <- sim_varma(models[[i]],250)
    k <- dim(models[[i]])[2L]; m <- ncol(Y)
    for(mean in c(FALSE,TRUE)) {
      a <- fit_varma_cml(Y,k,k,mean,control=list(maxit=500,reltol=1e-9))
      b <- fit_varma_echelon(Y,rep(k,m),mean,control=list(maxit=500,reltol=1e-9))
      close(a$ar,b$ar,1e-12); close(a$ma,b$ma,1e-12)
      close(a$cov,b$cov,1e-12); close(a$loglik,b$loglik,1e-12)
      close(predict(a,3,cov=FALSE),predict(b,3,cov=FALSE),1e-12)
      close(a$loglik,-a$optim.out$value)
      stopifnot(a$npar==m*m*(2*k)+as.integer(mean)*m+m*(m+1)/2,
                a$likelihood_type=="conditional",a$residuals_type=="recursive")
    }
  }
})

check("conditional VAR likelihood agrees with OLS, and unequal/zero orders work",{
  set.seed(84)
  f <- varma(ar=array(c(.3,.05,-.1,.2),c(2,2,1)),intercept=c(.4,-.2),cov=diag(2))
  Y <- sim_varma(f,350)
  fit <- fit_varma_cml(Y,1,0,TRUE,control=list(maxit=500,reltol=1e-12))
  ols <- lm.fit(cbind(1,Y[-nrow(Y),]),Y[-1,])
  close(fit$ar,t(ols$coefficients[-1,,drop=FALSE]),1e-5)
  close(fit$intercept,ols$coefficients[1,],1e-5)
  close(fit$cov,crossprod(ols$residuals)/(nrow(Y)-1),1e-7)
  for(orders in list(c(0,1),c(2,1),c(1,2))) {
    a <- suppressWarnings(fit_varma_cml(Y,orders[1],orders[2],TRUE))
    stopifnot(identical(as.integer(dim(a)),as.integer(c(2,orders))),
      all(is.finite(predict(a,2,cov=FALSE))))
    close(a$loglik,-a$optim.out$value)
  }
  a <- fit_varma_cml(Y[,1,drop=FALSE],0,1,FALSE)
  stopifnot(identical(as.integer(dim(a)),c(1L,0L,1L)))
  fails(fit_varma_cml(Y,-1,1)); fails(fit_varma_cml(Y,1,.5))
})

check("full conditional gradient matches finite differences for unequal orders",{
  set.seed(911); Y <- matrix(rnorm(600),300,2)
  for(orders in list(c(0,2),c(2,0),c(2,1),c(1,2))) {
    pat <- internal(".varma_cml_pattern")(2,orders[1],orders[2],TRUE)
    x <- rnorm(nrow(pat),sd=.03); k <- max(orders)
    fun <- function(v) internal("echelon_css_cpp")(Y,v,pat,k)
    z <- fun(x); stopifnot(z$valid)
    numerical <- vapply(seq_along(x),function(j) {
      u<-v<-x;u[j]<-u[j]+1e-6;v[j]<-v[j]-1e-6
      (fun(u)$value-fun(v)$value)/2e-6
    },numeric(1))
    close(z$gradient,numerical,1e-5)
  }
})

check("nonidentified CML returns transfer estimates without information inversion",{
  set.seed(713); Y <- sim_varma(varma_paper_designs()$model1,200)
  a <- fit_varma_cml(Y,1,1,FALSE,control=list(maxit=500))
  stopifnot(all(is.finite(irf(a,5))),all(is.finite(predict(a,5,cov=FALSE))))
  b <- bootstrap_varma(a,nsim=1,n=150,parallel=FALSE)[[1]]
  stopifnot(b$estimation_method=="fit_varma_cml",b$refit_args$r==a$refit_args$r)
})

check("all state-space backends agree under both initial covariance policies",{
  set.seed(122); Y <- sim_varma(varma_paper_designs()$model2,220)
  backends <- list(fit_varma_skf,fit_varma_kfas,fit_varma_fkf)
  for(policy in c("stationary","fixed")) {
    fits <- lapply(backends,function(fn) suppressWarnings(fn(Y,2,2,FALSE,
      initial_cov=policy,control=list(maxit=0))))
    for(f in fits) {
      stopifnot(f$initial_cov==policy,f$pilot_method=="HR1")
      close(f$initial_state_cov,f$pilot_state_cov)
    }
    for(f in fits[-1]) {
      close(f$loglik,fits[[1]]$loglik,1e-7)
      close(predict(f,3,cov=FALSE),predict(fits[[1]],3,cov=FALSE),1e-7)
    }
  }
  # Estimate a pure VAR to compare converged objectives across independent filters.
  for(policy in c("stationary","fixed")) {
    fits <- lapply(backends,function(fn) fn(Y,1,0,TRUE,initial_cov=policy,
      control=list(maxit=400,reltol=1e-10)))
    for(f in fits[-1]) close(f$loglik,fits[[1]]$loglik,1e-6)
    for(f in fits) {
      if(policy=="fixed") close(f$initial_state_cov,f$pilot_state_cov,1e-12)
      else close(f$initial_state_cov,
        f$transition_matrix%*%f$initial_state_cov%*%t(f$transition_matrix)+f$disturbance_cov)
      stopifnot(f$refit_args$initial_cov==policy,f$refit_args$pilot=="HR1")
    }
  }
  # Compatibility option preserves iterated pilot initialization.
  a <- suppressWarnings(fit_varma_skf(Y,1,1,FALSE,pilot="IHR",maxit=3,
                                     initial_cov="fixed",control=list(maxit=0)))
  pilot <- fit_varma_ihr(Y,1,1,FALSE,maxit=3,ret="regression")
  close(a$cov,pilot$res_cov)
  b <- suppressWarnings(fit_varma_skf(Y,1,1,FALSE,control=list(maxit=0)))
  close(b$cov,fit_varma_ihr(Y,1,1,FALSE,maxit=1,ret="regression")$res_cov)
  fails(fit_varma_skf(Y,initial_cov="unknown")); fails(fit_varma_fkf(Y,pilot="unknown"))
  # Frozen and stationary covariances must differ once dynamics change.
  pack <- function(f) {
    ch<-chol(f$cov);diag(ch)<-log(diag(ch))
    c(as.numeric(aperm(f$ar,c(1,3,2))),as.numeric(aperm(f$ma,c(1,3,2))),ch[upper.tri(ch,TRUE)])
  }
  pars <- pack(b); pars[1:4] <- pars[1:4]*.5
  fixed <- internal("varma_candidate_cpp")(pars,2,1,1,FALSE,b$pilot_state_cov)
  varying <- internal("varma_candidate_cpp")(pars,2,1,1,FALSE)
  stopifnot(fixed$valid,varying$valid,max(abs(fixed$P-varying$P))>1e-5)
  close(fixed$P,b$pilot_state_cov,1e-12)
})
cat("All CML and initialization checks passed.\n")
