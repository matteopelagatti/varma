# Shared conditional Gaussian estimator. Both full and echelon models use the
# same pilot regression, initialization, compiled objective/gradient and BFGS.
.varma_conditional_fit <- function(Y,pattern,p,q,intercept,r,control,method) {
  n <- nrow(Y); m <- ncol(Y); order <- max(p,q)
  if(n-order<=m) stop("The orders leave too few observations to estimate covariance")
  pars <- numeric(nrow(pattern)); names(pars) <- rownames(pattern)
  if(order==0L) {
    if(intercept) pars[] <- colMeans(Y)
  } else {
    if(is.null(r)) r <- .varma_long_order(Y,p,q,intercept)
    r <- .varma_integer(r,"r")
    if(r<order || n-r<=m*r+as.integer(intercept))
      stop("Pilot order r must be at least max(p,q) and leave residual degrees of freedom")
    rows <- seq.int(r+1L,n)
    XP <- .varma_lags(Y,r)[rows,,drop=FALSE]
    if(intercept) XP <- cbind(1,XP)
    pilot <- .varma_regression(XP,Y[rows,,drop=FALSE])
    E <- matrix(0,n,m); E[rows,] <- pilot$residuals
    if(n-r-order<=max(tabulate(pattern[,2L]+1L,nbins=m)))
      stop("Too few observations for the conditional ML pilot regression")
    rows <- seq.int(r+order+1L,n)
    for(i in seq_len(m)) {
      columns <- which(pattern[,2L]==i-1L)
      if(!length(columns)) next
      X <- vapply(columns,function(j) {
        kind<-pattern[j,1L];v<-pattern[j,3L]+1L;lag<-pattern[j,4L]
        switch(kind+1L,rep(1,length(rows)),E[rows,v]-Y[rows,v],Y[rows-lag,v],E[rows-lag,v])
      },numeric(length(rows)))
      pars[columns] <- .varma_regression(X,Y[rows,i,drop=FALSE])$coefficients
    }
  }
  evaluate <- function(x,gradient=TRUE) echelon_css_cpp(Y,x,pattern,order,gradient)
  initial <- evaluate(pars,FALSE)
  dynamic <- pattern[,1L] %in% c(2L,3L)
  for(i in seq_len(200L)) {
    if(isTRUE(initial$valid)) break
    pars[dynamic] <- pars[dynamic]*.8
    initial <- evaluate(pars,FALSE)
  }
  if(!isTRUE(initial$valid)) stop("Could not obtain a stable conditional ML starting point with positive definite covariance")
  last_x <- NULL; last <- NULL
  cached <- function(x) {
    if(is.null(last_x) || !identical(x,last_x)) { last <<- evaluate(x); last_x <<- x }
    last
  }
  if(length(pars) && order>0L) {
    opt <- stats::optim(pars,function(x) cached(x)$value,
                       function(x) as.numeric(cached(x)$gradient),
                       method="BFGS",control=control)
    pars <- opt$par
  } else opt <- list(par=pars,value=initial$value,convergence=0L,
                    counts=c("function"=1L,gradient=0L),message="Closed-form white-noise fit")
  z <- evaluate(pars,FALSE)
  if(!isTRUE(z$valid)) stop("Optimisation did not produce a valid conditional ML model")
  if(opt$convergence!=0L) warning("Conditional ML optimiser did not converge; inspect $optim.out")
  f <- varma(intercept=if(intercept) as.numeric(z$intercept) else NULL,
             ar=if(p) z$ar[,,seq_len(p),drop=FALSE] else NULL,
             ma=if(q) z$ma[,,seq_len(q),drop=FALSE] else NULL)
  f <- .varma_finish(f,Y,method,length(pars),convergence=opt$convergence)
  f$coefficients <- pars; f$optim.out <- opt
  list(fit=f,evaluation=z,r=r)
}

.varma_cml_pattern <- function(m,p,q,intercept) {
  rows <- list(); labels <- character()
  add <- function(type,i,j,lag,label) {
    rows[[length(rows)+1L]] <<- c(type,i-1L,j-1L,lag)
    labels <<- c(labels,label)
  }
  for(i in seq_len(m)) {
    if(intercept) add(0L,i,1L,0L,paste0("c[",i,"]"))
    for(lag in seq_len(p)) for(j in seq_len(m))
      add(2L,i,j,lag,sprintf("AR%d[%d,%d]",lag,i,j))
    for(lag in seq_len(q)) for(j in seq_len(m))
      add(3L,i,j,lag,sprintf("MA%d[%d,%d]",lag,i,j))
  }
  pattern <- do.call(rbind,rows); storage.mode(pattern) <- "integer"
  colnames(pattern) <- c("type","equation","variable","lag")
  rownames(pattern) <- labels
  pattern
}

#' Conditional Gaussian maximum likelihood for a full VARMA model
#'
#' Estimates unrestricted AR and MA coefficient matrices using the same
#' conditional likelihood, analytic gradient, HR1-style pilot and BFGS
#' optimization as fit_varma_echelon, without echelon restrictions.
#' @param Y finite numeric n by m matrix; rows are observations.
#' @param p,q non-negative AR and MA orders; at least one must be positive.
#' @param intercept logical; estimate an intercept vector (default TRUE).
#' @param r order of the pilot long VAR, or NULL for the logarithmic default.
#' @param control controls passed to stats::optim; default maxit is 200.
#' @details The criterion concentrates out the full innovation covariance.
#' Errors through max(p,q) are initialized to zero and the likelihood uses
#' observations max(p,q)+1 through n. Stable AR and invertible MA polynomials
#' are required. The compiled recursion evaluates the criterion and analytic
#' gradient jointly. No Kalman filter or stationary initial covariance is used.
#'
#' Identification is not imposed: coefficient estimates can be nonunique and
#' the information matrix singular. Numerical convergence does not establish
#' global maximization. With p=q=k and every echelon index equal to k, the
#' full and echelon routines use identical patterns and numerical algorithms.
#' @return A varma object with coefficients, optim.out, convergence, conditional
#' loglik, covariance, recursive residuals, and refit_args for bootstrap refits.
#' @seealso \code{\link{fit_varma_echelon}}, \code{\link{fit_varma_skf}}
#' @export
#' @examples
#' set.seed(11)
#' y <- matrix(rnorm(400),200,2)
#' fit <- fit_varma_cml(y,p=1,q=1,intercept=FALSE)
fit_varma_cml <- function(Y,p=1,q=1,intercept=TRUE,r=NULL,control=list(maxit=200)) {
  Y <- .varma_data(Y,p,q)
  if(length(intercept)!=1L || is.na(intercept) || !is.logical(intercept))
    stop("intercept must be TRUE or FALSE")
  pattern <- .varma_cml_pattern(ncol(Y),p,q,intercept)
  out <- .varma_conditional_fit(Y,pattern,p,q,intercept,r,control,"fit_varma_cml")
  f <- out$fit
  f$refit_args <- list(p=p,q=q,intercept=intercept,r=out$r,control=control)
  f
}
