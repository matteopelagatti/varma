#' Select Kronecker indices from time-series observations
#'
#' Uses sequential tests of the smallest canonical correlation between past
#' observations and an expanding set of future observations. This is an
#' internal implementation of the procedure described by Tsay (2014), with
#' no dependency on MTS.
#' @param Y finite numeric n by m data matrix; rows are dates and columns are
#' series. A numeric vector is treated as a single series.
#' @param plag positive number of lags representing the past (default 5).
#' @param crit test size for zero canonical correlation (default 0.05).
#' @param tol relative numerical rank tolerance (default 1e-7).
#' @param verbose logical; print the test table and selected indices.
#' @details Series order is preserved. Each candidate future column is retained
#' only while its smallest canonical correlation remains significant. When the
#' null is not rejected, that series' index is set to the current future horizon
#' and the candidate column is removed, including at horizon zero.
#'
#' The compiled implementation factors the past matrix once per horizon,
#' updates an orthonormal future basis by reorthogonalization, and takes an SVD
#' of the small cross-product of the two bases. Canonical scores are obtained
#' directly from those bases, avoiding coefficient inverses. At positive
#' horizons the squared correlation is divided by
#' \eqn{1+2\sum_{j=1}^h\hat\rho_x(j)\hat\rho_y(j)}. The chi-square statistic
#' retains the original sample-size multiplier used by MTS::Kronid, namely
#' \eqn{n-1-(d_p+d_f-1)/2}, with \eqn{d_p-d_f+1} degrees of freedom.
#'
#' At least \eqn{n-plag>2m\,plag+1} observations are required. Rank-deficient
#' data, an invalid serial-correlation correction, and exhaustion of the past
#' dimension produce explicit errors rather than guessed indices. Increase
#' plag (with sufficient data), reduce redundant series, or supply indices to
#' fit_varma_echelon when the data do not support these tests. Data are centered
#' and scaled internally; constant series are not supported. Single-series
#' inputs and all-zero index vectors are supported.
#'
#' This is a finite sequence of approximate tests at the supplied significance
#' level, not a guarantee of consistent model selection or valid post-selection
#' inference. Results agree with MTS on regular cases but can differ when its
#' horizon-zero retained-column or dimension-dropping edge cases occur.
#' @return A list with index (an integer vector in column order), tests (a
#' three-row matrix of statistics, degrees of freedom and p-values), details
#' (a data frame including horizon, component, sample size, canonical
#' correlation and serial-correlation correction), method, plag, crit, and tol.
#' @references Tsay (2014), Multivariate Time Series Analysis with R and
#' Financial Applications, Chapter 4. MTS::Kronid is a compatibility reference;
#' this implementation does not call or include its source code.
#' @seealso \code{\link{fit_varma_echelon}}, \code{\link{kronecker_indices}} for
#' indices calculated from a specified population model rather than data.
#' @export
#' @examples
#' set.seed(17)
#' select_kronecker_indices(matrix(rnorm(600), 300, 2))$index
select_kronecker_indices <- function(Y,plag=5,crit=.05,tol=1e-7,verbose=FALSE) {
  Y <- as.matrix(Y)
  if(!is.numeric(Y) || is.complex(Y) || !length(Y) || any(!is.finite(Y)))
    stop("Y must be a nonempty finite real numeric matrix")
  if(!is.numeric(plag) || length(plag)!=1L || !is.finite(plag) ||
      plag<1 || plag>nrow(Y) || plag!=floor(plag))
    stop("plag must be a positive integer smaller than the sample size")
  plag <- as.integer(plag)
  if(!is.numeric(crit) || length(crit)!=1L || !is.finite(crit) || crit<=0 || crit>=1)
    stop("crit must lie strictly between zero and one")
  if(!is.numeric(tol) || length(tol)!=1L || !is.finite(tol) || tol<=0 || tol>=1)
    stop("tol must lie strictly between zero and one")
  if(!is.logical(verbose) || length(verbose)!=1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE")
  if(nrow(Y)-plag<=2*as.double(ncol(Y))*plag+1)
    stop("Too few observations for plag; reduce plag or supply kidx")
  out <- kronid_cpp(Y,plag,crit,tol)
  names(out$index) <- colnames(Y)
  rownames(out$tests) <- c("statistic","df","p_value")
  out$details <- as.data.frame(out$details)
  names(out$details) <- c("horizon","component","nobs","past_dim","future_dim",
                         "canonical_sq","correction","adjusted_sq",
                         "statistic","df","p_value")
  out$method <- "varma::select_kronecker_indices"
  out$plag <- plag; out$crit <- crit; out$tol <- tol
  if(verbose) { print(out$details); cat("Kronecker indices:\n"); print(out$index) }
  out
}

.varma_echelon_pattern <- function(kidx,intercept) {
  m <- length(kidx)
  rows <- list(); labels <- character()
  add <- function(type,i,j,lag,label) {
    rows[[length(rows)+1L]] <<- c(type,i-1L,j-1L,lag)
    labels <<- c(labels,label)
  }
  # Row degrees are the supplied Kronecker indices, in data-column order.
  # Lag-zero AR and MA factors coincide and have unit lower-triangular form.
  for(i in seq_len(m)) {
    if(intercept) add(0L,i,1L,0L,paste0("c[",i,"]"))
    for(j in seq_len(m)) if(j<i && kidx[j]>kidx[i])
      add(1L,i,j,0L,sprintf("Phi0[%d,%d]",i,j))
    for(lag in seq_len(kidx[i])) for(j in seq_len(m))
      if(lag>kidx[i]-kidx[j]) add(2L,i,j,lag,sprintf("AR%d[%d,%d]",lag,i,j))
    for(lag in seq_len(kidx[i])) for(j in seq_len(m))
      add(3L,i,j,lag,sprintf("MA%d[%d,%d]",lag,i,j))
  }
  pattern <- if(length(rows)) do.call(rbind,rows) else matrix(integer(),0L,4L)
  storage.mode(pattern) <- "integer"
  colnames(pattern) <- c("type","equation","variable","lag")
  rownames(pattern) <- labels
  pattern
}

#' Estimate a VARMA model in echelon form
#'
#' Fits the echelon restrictions implied by Kronecker indices, preserving the
#' order of the data columns. If kidx is NULL, select_kronecker_indices selects
#' the indices using internal canonical-correlation tests before fitting.
#'
#' @param Y numeric n by m matrix without missing values.
#' @param kidx vector of m non-negative integer Kronecker indices, or NULL to
#' estimate them from the data. Indices are not sorted or reordered.
#' @param intercept logical; estimate a vector of intercepts if TRUE.
#' @param plag number of past lags for index selection (only when kidx is NULL).
#' @param crit test size for automatic index selection (default 0.05).
#' @param r order of the initial long VAR, or NULL for the package's logarithmic
#' pilot-order default. Unused when every index is zero.
#' @param control list of controls for stats::optim (BFGS); default maxit is 200.
#' @param verbose logical; display the index-selection output when TRUE.
#' @details The estimator minimises a conditional Gaussian negative log-likelihood,
#' concentrating out the innovation covariance. Presample errors through
#' max(kidx) are zero. A compiled recursion supplies the analytic gradient,
#' including derivatives of the common nonidentity lag-zero factor. The search
#' requires stable AR and invertible MA polynomials. This is conditional
#' likelihood, not the stationary state-space likelihood used by fit_varma_kfas.
#'
#' The returned AR/MA coefficients and intercept are normalised to lag-zero
#' identity for compatibility with irf, predict, residuals, and bootstrap_varma.
#' The echelon slot retains the common lag-zero factor and the original
#' restricted coefficient matrices. Positive MA signs follow this package's
#' convention. A zero index, an all-zero vector, and a single-series model are
#' supported. Automatic selection uses compiled routines within this package;
#' selection failures are reported rather than silently replaced by a guessed
#' specification. Selecting indices adds model-selection
#' uncertainty; no post-selection inference guarantee is implied.
#'
#' Bootstrap refits condition on the fitted indices by default. Pass kidx=NULL
#' to bootstrap_varma to repeat index selection in each replicate.
#' @return A varma object with kronecker_indices, indices_estimated, selection,
#' echelon, coefficients, optim.out, and convergence fields. npar includes only
#' free echelon coefficients plus innovation covariance parameters. Numerical
#' convergence is not a guarantee that a global optimum was found.
#' @references Tsay (2014), Multivariate Time Series Analysis with R and
#' Financial Applications, Chapter 4.
#' @seealso \code{\link{select_kronecker_indices}}
#' @export
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(400),200,2)
#' fit <- fit_varma_echelon(y,kidx=c(0,0))
#' predict(fit,n.ahead=2)
#' \dontrun{
#' fit_varma_echelon(y,kidx=c(1,1))
#' fit_varma_echelon(y) # estimate the indices internally
#' }
fit_varma_echelon <- function(Y,kidx=NULL,intercept=TRUE,plag=5,crit=.05,
                              r=NULL,control=list(maxit=200),verbose=FALSE) {
  Y <- as.matrix(Y)
  if(!is.numeric(Y) || !length(Y) || any(!is.finite(Y)))
    stop("Y must be a nonempty finite numeric matrix")
  n <- nrow(Y); m <- ncol(Y)
  if(n<=m+1L) stop("Too few observations to estimate an innovation covariance")
  if(length(intercept)!=1L || is.na(intercept) || !is.logical(intercept)) stop("intercept must be TRUE or FALSE")
  selected <- is.null(kidx); selection <- NULL
  if(selected) {
    selection <- select_kronecker_indices(Y,plag=plag,crit=crit,verbose=verbose)
    kidx <- selection$index
  }
  if(!is.numeric(kidx) || length(kidx)!=m || any(!is.finite(kidx)) ||
      any(kidx<0) || any(kidx!=as.integer(kidx)))
    stop("kidx must contain one non-negative integer per data column")
  kidx <- as.integer(kidx); names(kidx) <- colnames(Y)
  order <- max(kidx)
  if(n-order<=m) stop("The indices leave too few observations")
  pattern <- .varma_echelon_pattern(kidx,intercept)
  pars <- numeric(nrow(pattern)); names(pars) <- rownames(pattern)
  if(order==0L) {
    if(intercept) pars[] <- colMeans(Y)
  } else {
    if(is.null(r)) r <- .varma_long_order(Y,order,order,intercept)
    r <- .varma_integer(r,"r")
    if(r<order || n-r<=m*r+as.integer(intercept)) stop("Pilot order r must be at least max(kidx) and leave residual degrees of freedom")
    rows <- seq.int(r+1L,n)
    XP <- .varma_lags(Y,r)[rows,,drop=FALSE]
    if(intercept) XP <- cbind(1,XP)
    pilot <- .varma_regression(XP,Y[rows,,drop=FALSE])
    E <- matrix(0,n,m); E[rows,] <- pilot$residuals
    if(n-r-order<=max(tabulate(pattern[,2L]+1L,nbins=m))) stop("Too few observations for the echelon pilot regression")
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
  if(!isTRUE(initial$valid)) stop("Could not obtain a stable echelon starting point with positive definite covariance")
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
  if(!isTRUE(z$valid)) stop("Optimisation did not produce a valid echelon model")
  if(opt$convergence!=0L) warning("Echelon optimiser did not converge; inspect $optim.out")
  f <- varma(intercept=if(intercept) as.numeric(z$intercept) else NULL,
             ar=if(order) z$ar else NULL,ma=if(order) z$ma else NULL)
  f <- .varma_finish(f,Y,"fit_varma_echelon",length(pars),convergence=opt$convergence)
  raw_ar <- raw_ma <- array(0,c(m,m,order)); raw_c <- numeric(m)
  for(j in seq_along(pars)) {
    kind<-pattern[j,1L];i<-pattern[j,2L]+1L;v<-pattern[j,3L]+1L;lag<-pattern[j,4L]
    if(kind==0L) raw_c[i]<-pars[j]
    if(kind==2L) raw_ar[i,v,lag]<-pars[j]
    if(kind==3L) raw_ma[i,v,lag]<-pars[j]
  }
  f$kronecker_indices <- kidx; f$indices_estimated <- selected; f$selection <- selection
  f$echelon <- list(Phi0=z$Phi0,ar=raw_ar,ma=raw_ma,intercept=raw_c,pattern=pattern)
  f$coefficients <- pars; f$optim.out <- opt
  f$refit_args <- list(kidx=kidx,intercept=intercept,plag=plag,crit=crit,r=r,control=control)
  f
}
