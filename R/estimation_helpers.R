# Shared input checks and bookkeeping for the estimation routines.
.varma_data <- function(Y, p, q) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || !length(Y) || any(!is.finite(Y)))
    stop("Y must be a nonempty numeric matrix without missing or infinite values")
  if (length(p) != 1L || length(q) != 1L || any(!is.finite(c(p, q))) ||
      any(c(p, q) < 0) || any(c(p, q) != as.integer(c(p, q))))
    stop("p and q must be non-negative integers")
  if (p + q == 0L) stop("p = q = 0: no VARMA model specified")
  if (nrow(Y) <= max(p, q)) stop("Not enough observations for the requested orders")
  Y
}

.varma_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.finite(x) || x < minimum || x != as.integer(x))
    stop(name, " must be an integer >= ", minimum)
  as.integer(x)
}

.varma_long_order <- function(Y, p, q, intercept = FALSE) {
  # Logarithmic growth avoids fitting a pilot with a fixed fraction of the data.
  cap <- floor((nrow(Y) - as.integer(intercept) - 1L)/(ncol(Y) + 1L))
  min(max(p + q, ceiling(log(nrow(Y))^1.5)), cap)
}

.varma_lags <- function(Y, order) {
  n <- nrow(Y); m <- ncol(Y)
  ans <- matrix(0, n, m * order)
  for (j in seq_len(order))
    ans[(j + 1L):n, ((j - 1L)*m + 1L):(j*m)] <- Y[seq_len(n - j), , drop = FALSE]
  ans
}

.varma_regression <- function(X, Y) {
  reg <- stats::lm.fit(X, Y)
  # Preserve matrix dimensions for m=1; use a particular least-squares solution
  # when a redundant parametrisation produces aliased columns.
  reg$coefficients <- matrix(reg$coefficients, ncol(X), ncol(Y))
  reg$coefficients[is.na(reg$coefficients)] <- 0
  # lm.fit already computed these using its QR factorisation. Repeating both
  # matrix products is costly in every Hannan--Rissanen iteration.
  reg$residuals <- matrix(reg$residuals,nrow(Y),ncol(Y))
  reg$fitted.values <- matrix(reg$fitted.values,nrow(Y),ncol(Y))
  reg
}

.varma_from_coefficients <- function(coeff, p, q, intercept, m) {
  est <- t(coeff)
  cnst <- if (intercept) as.numeric(est[, 1L]) else NULL
  if (intercept) est <- est[, -1L, drop = FALSE]
  varma(intercept = cnst,
        ar = if (p) array(est[, seq_len(m*p), drop = FALSE], c(m, m, p)) else NULL,
        ma = if (q) array(est[, m*p + seq_len(m*q), drop = FALSE], c(m, m, q)) else NULL)
}

.varma_gaussian_ll <- function(E, S) {
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch) || any(!is.finite(E))) return(-Inf)
  -0.5 * (length(E)*log(2*pi) + 2*nrow(E)*sum(log(diag(ch))) +
            sum(backsolve(ch, t(E), transpose = TRUE)^2))
}

.varma_finish <- function(fit, Y, method, dynamic_npar, start = max(dim(fit)[-1L]),
                          rank = NULL, convergence = NULL, iterations = NULL,
                          residuals = NULL) {
  recursive <- is.null(residuals)
  E <- if (recursive) residuals.varma(fit, Y, recompute=TRUE) else as.matrix(residuals)
  if (!identical(dim(E),dim(Y))) stop("Residuals must be aligned with Y")
  used <- seq.int(start + 1L, nrow(Y))
  S <- crossprod(E[used, , drop = FALSE])/length(used)
  fit$cov <- S
  fit$loglik <- .varma_gaussian_ll(E[used, , drop = FALSE], S)
  fit$n <- nrow(Y)
  fit$nobs <- length(used)  # independent vector observations, not scalar entries
  fit$npar <- dynamic_npar + ncol(Y)*(ncol(Y) + 1L)/2L
  fit$y <- Y
  fit$residuals <- E
  fit$residuals_start <- start + 1L
  fit$residuals_type <- if (recursive) "recursive" else "regression"
  fit$estimation_method <- method
  fit$rank <- rank
  fit$convergence <- convergence
  fit$iterations <- iterations
  fit$stability <- is_stable(fit)
  fit$likelihood_type <- if (recursive) "conditional" else "regression_approximation"
  fit
}

.varma_process_mean <- function(fit) {
  m <- dim(fit)[1L]
  if (!is.null(fit$mean)) return(as.numeric(fit$mean))
  if (is.null(fit$intercept)) return(numeric(m))
  aa <- diag(m)
  for (j in seq_len(dim(fit)[2L])) aa <- aa - matrix(fit$ar[, , j], m, m)
  as.numeric(solve(aa, fit$intercept))
}

.varma_ic <- function(loglik, npar, n, ic) {
  switch(ic, AIC = -2*loglik + 2*npar,
         AICc = if (n <= npar + 1) Inf else -2*loglik + 2*npar + 2*npar*(npar+1)/(n-npar-1),
         BIC = -2*loglik + npar*log(n), HQC = -2*loglik + 2*npar*log(log(n)))
}

# All three likelihood backends use the same parametrisation and a stationary
# covariance recomputed at EVERY candidate, not only at the pilot estimate.
.varma_ml <- function(Y, p, q, intercept, maxit, control, backend) {
  Y <- .varma_data(Y,p,q)
  m <- ncol(Y); n <- nrow(Y); d <- m*max(p,q+1L)
  pilot <- fit_varma_ihr(Y,p,q,intercept,maxit=maxit,ret="regression")
  initial <- .varma_from_coefficients(pilot$coefficients,p,q,intercept,m)
  initial$cov <- pilot$res_cov
  # A feasible starting point is needed even when an unrestricted pilot has
  # roots outside the stable/invertible region.
  for (part in c("ar","ma")) {
    if (length(initial[[part]])) {
      for (i in seq_len(1000L)) {
        roots <- inv_roots(initial,part)[[part]]
        if (all(Mod(roots) < .98)) break
        initial[[part]] <- initial[[part]]*.9
      }
    }
  }
  ss0 <- varma_to_ss(initial)
  ch <- chol(initial$cov)
  diag(ch) <- log(diag(ch))
  pack <- function(x) if (length(x)) as.numeric(aperm(x,c(1,3,2))) else numeric(0)
  inits <- c(if (intercept) colMeans(Y),pack(initial$ar),pack(initial$ma),ch[upper.tri(ch,TRUE)])
  Yt <- t(Y)
  make_candidate <- function(pars) {
    z <- varma_candidate_cpp(pars,m,p,q,intercept)
    if (!isTRUE(z$valid)) return(NULL)
    z$mean <- as.numeric(z$mean)
    z$y <- if (intercept) sweep(Y,2,z$mean,"-") else Y
    z
  }
  first <- make_candidate(inits)
  if (is.null(first)) stop("Could not construct a stable, invertible likelihood starting point")
  if (backend == "kfas") {
    SSMcustom <- KFAS::SSMcustom # KFAS recognises the unqualified formula special.
    model0 <- KFAS::SSModel(first$y ~ 0 + SSMcustom(
      Z=ss0$Z,T=first$T,R=first$R,Q=first$Q,a1=numeric(d),
      P1=first$P,P1inf=matrix(0,d,d),index=seq_len(m)),H=matrix(0,m,m))
  }
  kfas_model <- function(z) {
    mod <- model0
    mod$y[] <- z$y; mod$T[,,1] <- z$T; mod$R[,,1] <- z$R
    mod$Q[,,1] <- z$Q; mod$P1 <- z$P
    mod
  }
  fkf_run <- function(z,verbose=FALSE) FKF.SP::fkf.SP(
    a0=numeric(d),P0=z$P,dt=matrix(0,d,1),ct=matrix(0,m,1),
    Tt=z$T,Zt=ss0$Z,HHt=z$W,GGt=numeric(m),yt=t(z$y),verbose=verbose)
  objective <- function(pars) {
    val <- tryCatch({
      if (backend == "skf") return({
        ll <- varma_loglik_cpp(pars,Yt,p,q,intercept)
        if (is.finite(ll)) -ll else 1e100
      })
      z <- make_candidate(pars)
      if (is.null(z)) return(1e100)
      switch(backend,
        kfas = as.numeric(stats::logLik(kfas_model(z),check.model=FALSE)),
        fkf = as.numeric(fkf_run(z)),
        skf = kalmanLogLik(Matrix::Matrix(z$T,sparse=TRUE),Matrix::Matrix(z$R,sparse=TRUE),
                           z$Q,matrix(0,d,1),z$P,t(z$y)))
    },error=function(e) -Inf)
    if (is.finite(val)) -val else 1e100
  }
  opt <- stats::optim(inits,objective,method="BFGS",control=control)
  z <- make_candidate(opt$par)
  if (is.null(z) || opt$value >= 1e99) stop("Likelihood optimisation did not produce a valid model")
  if (opt$convergence != 0L) warning("Likelihood optimiser did not converge; inspect $optim.out")
  fit <- varma(mean=if (intercept) z$mean else NULL,ar=if (p) z$ar else NULL,
               ma=if (q) z$ma else NULL,cov=z$Q)
  if (backend == "kfas") {
    mod <- kfas_model(z)
    kf <- KFAS::KFS(mod,filtering="state",smoothing="none")
    fit$residuals <- KFAS::mvInnovations(kf)$v
    fit$state_pred_mean <- kf$a[n+1L,]
    fit$state_pred_cov <- matrix(kf$P[,,n+1L],d,d)
    fit$loglik <- as.numeric(stats::logLik(mod))
    raw <- list(model=mod,optim.out=opt)
  } else if (backend == "fkf") {
    raw <- fkf_run(z,TRUE)
    # FKF.SP also processes components sequentially; its v is not the joint
    # vector innovation. at contains the state prediction before observing y_t.
    fit$residuals <- z$y-t(raw$at[seq_len(m),seq_len(n),drop=FALSE])
    fit$state_pred_mean <- raw$at[,n+1L]
    fit$state_pred_cov <- matrix(raw$Pt[,,n+1L],d,d)
    fit$loglik <- -objective(opt$par)
    raw$R <- z$R; raw$Q <- z$Q; raw$HH <- z$W; raw$optim.out <- opt
  } else {
    aa <- matrix(0,d,n+1L)
    PP <- array(0,c(d,d,n+1L)); PP[,,1L] <- z$P
    vv <- matrix(0,m,n); FF <- array(0,c(m,m,n))
    fit$loglik <- kalman(Matrix::Matrix(z$T,sparse=TRUE),Matrix::Matrix(z$R,sparse=TRUE),
                         z$Q,aa,PP,t(z$y),vv,FF)
    fit$residuals <- t(vv)
    fit$state_pred_mean <- aa[,n+1L]; fit$state_pred_cov <- matrix(PP[,,n+1L],d,d)
    raw <- NULL
  }
  fit$transition_matrix <- z$T; fit$disturbance_cov <- z$W
  fit$initial_state_cov <- z$P
  fit$n <- n; fit$nobs <- n; fit$npar <- length(opt$par); fit$y <- Y
  fit$residuals_start <- 1L
  fit$estimation_method <- paste0("fit_varma_",backend)
  fit$convergence <- opt$convergence; fit$optim.out <- opt
  fit$stability <- is_stable(fit)
  fit$likelihood_type <- "stationary"
  list(fit=fit,raw=raw)
}

.varma_dp <- function(Y,p,q,intercept,r,type) {
  Y <- .varma_data(Y,p,q); n <- nrow(Y); m <- ncol(Y)
  name <- if (type == "diagonal") "fit_varma_dpd" else "fit_varma_dpf"
  if (!q) {
    out <- fit_varma_ihr(Y,p,0,intercept); out$estimation_method <- name
    return(out)
  }
  r <- .varma_integer(r,"r")
  if (r < max(p,q) || n-r <= m*r || n-r-max(p,q) <= m*p+q)
    stop("Pilot order r leaves too few observations")
  mu <- if (intercept) colMeans(Y) else numeric(m)
  YY <- sweep(Y,2,mu,"-")
  rows <- seq.int(r+1L,n)
  pilot <- .varma_regression(.varma_lags(YY,r)[rows,,drop=FALSE],YY[rows,,drop=FALSE])
  U <- pilot$residuals
  S <- crossprod(U)/nrow(U)
  cpp <- estimate_varma_gls_cpp(YY[rows,,drop=FALSE],U,chol2inv(chol(S)),p,q,type)
  coefs <- as.numeric(cpp$coefficients)
  aa <- if (p) array(0,c(m,m,p)) else NULL
  if (p) {
    phi <- matrix(coefs[seq_len(m*m*p)],m*p,m)
    for (j in seq_len(p)) aa[,,j] <- t(phi[(j-1L)*m+seq_len(m),,drop=FALSE])
  }
  bb <- array(0,c(m,m,q)); theta <- coefs[m*m*p+seq_len(length(coefs)-m*m*p)]
  for (j in seq_len(q)) {
    block <- diag(if (type=="final") rep(theta[j],m) else theta[(seq_len(m)-1L)*q+j],m,m)
    bb[,,j] <- block
  }
  cc <- diag(m)
  for (j in seq_len(p)) cc <- cc-matrix(aa[,,j],m,m)
  out <- varma(intercept=if (intercept) as.numeric(cc %*% mu) else NULL,ar=aa,ma=bb)
  .varma_finish(out,Y,name,length(coefs)+m*as.integer(intercept))
}
