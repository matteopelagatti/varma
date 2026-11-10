#' @useDynLib varma, .registration = TRUE
#' @import Rcpp
NULL


#' varma object class
#'
#' This constructor function has been created just to illustrate the structure
#' of the varma object class. The varma object class is designed to contain
#' all the information needed to use VARMA models at various levels.
#' You can just use it to store a VARMA process definition or to contain
#' an estimated VARMA model with data, goodness of fit, residuals, etc.
#' A VARMA process definition can be used to produce simulated sample paths,
#' population quantities (IRFs, characeteristic roots, identification, etc.).
#' An estimated VARMA model can be used to produce the sample counterparts of the
#' population quantities, bootstapped sample paths, information criteria, and forecasts.
#' In the following we use the following symbols: \eqn{m}, \eqn{p}, \eqn{q}, \eqn{n} represent
#' the dimension of the time series vector, the largest lag of the VAR part of the model,
#' the largest lag of the VMA part of the mode, and the number of time series points, respectively.
#'
#' @param intercept \eqn{m}-vector with the intercepts.
#' @param mean \eqn{m}-vector with means (typically in alternative to `intercept`).
#' @param ar array of dimensions \eqn{m \times m \times p} containing the VAR matrices.
#' @param ma array of dimensions \eqn{m \times m \times q} containing the VMA matrices.
#' @param cov \eqn{m\times m} covariance matrix of the innovations.
#' @param chol_cov  Cholesky factor of the covariance matrix of the innovations: this is rearly used,
#' but if some function working with varma objects produces it, it is good to store it to
#' make its reusable by other functions, that can be more efficient.
#' @param estimation_method string with the name of the function that produced the estimates.
#' @param loglik value of the log-likelihood function for the estimated VARMA.
#' @param n number of time-points of the time series (regardless of possible missing observations).
#' @param nobs number of vector observations contributing to the likelihood.
#' @param npar number of estimated parameters.
#' @param y \eqn{n\times m} or \eqn{n-p\times m} data matrix.
#' @param residuals \eqn{n\times m} or \eqn{n-p\times m} residuals matrix.
#' @param state_pred_mean mean of the one-step-ahead predicted state vector.
#' at time \eqn{t = n+1} (relevant if the VARMA was estimated in state-space form).
#' @param state_pred_cov covariance matrix of the one-step-ahead predicted state vector
#' at time \eqn{t = n+1} (relevant if the VARMA was estimated in state-space form).
#' @param transition_matrix transition matrix of the VARMA in state space
#' (relevant if the VARMA was estimated in state-space form).
#' @param disturbance_cov covariance matrix of the state disturbance vector
#' (relevant if the VARMA was estimated in state-space form).
#'
#' @export
varma <- function(intercept = NULL,
                  mean = NULL,
                  ar = NULL,
                  ma = NULL,
                  cov = NULL,
                  chol_cov = NULL,
                  estimation_method = NULL,
                  loglik = NULL,
                  n = NULL,
                  nobs = NULL,
                  npar = NULL,
                  y = NULL,
                  residuals = NULL,
                  state_pred_mean = NULL,
                  state_pred_cov = NULL,
                  transition_matrix = NULL,
                  disturbance_cov = NULL
) {
  structure(
    list(
      intercept = intercept,
      mean = mean,
      ar = ar,
      ma = ma,
      cov = cov,
      chol_cov = chol_cov,
      estimation_method = estimation_method,
      loglik = loglik,
      n = n,
      nobs = nobs,
      npar = npar,
      y = y,
      residuals = residuals,
      state_pred_mean = state_pred_mean,
      state_pred_cov = state_pred_cov,
      transition_matrix = transition_matrix,
      disturbance_cov = disturbance_cov
    ),
    class = "varma"
  )
}


# --- Functions that extract infos from varma objects ---


#' Dimensions of a varma object
#'
#' @param x a varma object.
#'
#' @returns A vector of length 3 with the dimensions of the varma object:
#' number of time seires (m), order of AR part (p), order of MA part (q).
#'
#' @method dim varma
#' @export
dim.varma <- function(x) {
  # returns the dimensions of the varma object
  out <- c(m = 0, p = 0, q = 0)
  if (!is.null(x$ar)) {
    out[1] <- dim(x$ar)[1]
    out[2] <- dim(x$ar)[3]
  }
  if (!is.null(x$ma)) {
    out[1] <- dim(x$ma)[1]
    out[3] <- dim(x$ma)[3]
  }
  if (out[1] == 0 && !is.null(x$cov)) out[1] <- nrow(x$cov)
  out
}

#' Print method for varma object
#'
#' @param x a varma object.
#' @param ... not used.
#'
#' @returns Nothing.
#'
#' @method print varma
#' @export
print.varma <- function(x, ...) {
  mpq <- dim(x)
  # prints the varma object
  cat("VARMA(", mpq[2], ", ", mpq[3], ") model with ",
      mpq[1]," time series\n", sep = "")
  if (!is.null(x$intercept)) {
    cat("Intercept vector:\n")
    print(x$intercept)
  }
  if (!is.null(x$mean)) {
    cat("Mean vector:\n")
    print(x$mean)
  }
  if (!is.null(x$ar) && mpq[2] > 0) {
    cat("\nAR part:\n")
    for (i in 1:mpq[2]) {
      M <- matrix(x$ar[,,i], mpq[1], mpq[1])
      dimnames(M) <- list(
        paste0("y", 1:mpq[1], "_t"),
        paste0("y", 1:mpq[1], "_[t-", i, "]")
      )
      print(M)
    }
  }
  if (!is.null(x$ma) && mpq[3] > 0) {
    cat("\nMA part:\n")
    for (i in 1:mpq[3]) {
      M <- matrix(x$ma[,,i], mpq[1], mpq[1])
      dimnames(M) <- list(
      paste0("y", 1:mpq[1], "_t"),
      paste0("e", 1:mpq[1], "_[t-", i, "]")
    )
      print(M)
    }
  }
  if (!is.null(x$cov)) {
    cat("\nCovariance matrix:\n")
    dimnames(x$cov) <- list(
      paste0("e", 1:mpq[1]),
      paste0("e", 1:mpq[1])
    )
    print(x$cov)
  }
  if (!is.null(x$estimation_method)) {
    cat("\nEstimated using", x$estimation_method, "\n")
  }
  if (!is.null(x$loglik)) {
    npar <- attr(logLik.varma(x),"df")
    n <- x$nobs
    cat("N. of vector observations      =", n, "\n")
    cat("N. of estimated parameters     =", npar, "\n")
    cat("log-lik =", as.numeric(x$loglik), "\n")
    for (criterion in c("AIC","AICc","BIC","HQC"))
      cat(criterion,"=",.varma_ic(as.numeric(x$loglik),npar,n,criterion),"\n")
  }
}

#' residuals method for varma object
#'
#' @param object a varma object.
#' @param Y data matrix if NULL the y slot in the varma object is used
#' @param recompute logical, if FALSE and the varma object contains residuals,
#' the function returns those residuals, if TRUE the function recomputes the
#' residuals.
#' @param ... not used.
#'
#' @returns Matrix of residuals. If recompute is TRUE, the residuals are always
#' recomputed by this function, if it is FALSE and the varma object slot
#' residuals is not NULL, then these residuals are returned.
#'
#' @method residuals varma
#' @export
residuals.varma <- function(object, Y = NULL, recompute = TRUE, ...) {
  if (!is.null(object$residuals) && !recompute) return(object$residuals)
  if (is.null(Y)) Y <- object$y
  if (is.null(Y)) stop("Y and object$y both missing")
  Y <- as.matrix(Y); mpq <- dim(object); m <- mpq[1L]
  ar <- if (mpq[2L]) object$ar else array(0,c(m,m,0))
  ma <- if (mpq[3L]) object$ma else array(0,c(m,m,0))
  if (!is.null(object$intercept)) {
    varma_residuals(Y,ar,ma,object$intercept)
  } else {
    if (!is.null(object$mean)) Y <- sweep(Y,2,object$mean,"-")
    varma_residuals(Y,ar,ma,numeric(0))
  }
}


#' Companion form of a VARMA
#'
#' It computes the companion form of the AR and of the MA part.
#'
#' @param varma a varma object.
#' @param part vector of strings: can be `"ar"`, `"ma"` or both `c("ar", "ma")`,
#' indicating which part should ne cast in companion form.
#'
#' @export
companion_form <- function(varma, part = c("ar", "ma")) {
  mpq <- dim.varma(varma)
  if (any(part == "ar") & mpq[2] > 0) {
    # returns the companion form of the AR part of the varma object
    ar <- matrix(0, mpq[1]*mpq[2], mpq[1]*mpq[2])
    ar[1:mpq[1], ] <- as.numeric(varma$ar)
    if (mpq[2] > 1) {
      ar[(mpq[1]+1):(mpq[1]*mpq[2]), 1:(mpq[1]*(mpq[2]-1))] <-
        diag(1, mpq[1]*(mpq[2]-1), mpq[1]*(mpq[2]-1))
    }
  } else {ar <- NULL}
  if (any(part == "ma") & mpq[3] > 0) {
    # returns the companion form of the MA part of the varma object
    ma <- matrix(0, mpq[1]*mpq[3], mpq[1]*mpq[3])
    ma[1:mpq[1], ] <- -as.numeric(varma$ma)
    if (mpq[3] > 1) {
      ma[(mpq[1]+1):(mpq[1]*mpq[3]), 1:(mpq[1]*(mpq[3]-1))] <-
        diag(1, mpq[1]*(mpq[3]-1), mpq[1]*(mpq[3]-1))
    }
  } else {ma <- NULL}
  list(
    ar = ar,
    ma = ma
  )
}


#' VARMA in state space form
#'
#' It casts the VARMA model in state space form:
#'
#' \eqn{y_t = Z_t \alpha_t}
#' \eqn{\alpha_{t+1} = T \alpha_t + R \eta_t}
#' where \eqn{\eta_t} is a white noise process with covariance matrix \eqn{Q}.
#'
#' @param varma varma object.
#' @param form string with the name of the state-space form to use: "harvey" or "pearlman".
#'
#' @returns A list with the following matrices of the state-space form:
#' T, R, Q, Z.
#'
#' @export
varma_to_ss <- function(varma, form = c("harvey", "pearlman")) {
  form <- match.arg(form)
  mpq <- dim(varma)
  if ((mpq[1] == 0) | (mpq[2] + mpq[3] == 0)) stop("VARMA object is empty")
  if (is.null(varma$cov)) stop("Covariance matrix is missing")
  m <- mpq[1]
  r <- if (form == "harvey") max(mpq[2], mpq[3]+1) else max(mpq[2], mpq[3])
  mr <- m*r
  mT <- matrix(0, mr, mr)
  if (r > 1) mT[1:(mr-m), (m+1):mr] <- diag(mr-m)
  mR <- if (form == "harvey") rbind(diag(m), matrix(0, mr-m, m)) else matrix(0, mr, m)
  mQ <- varma$cov
  mZ <- cbind(diag(1, m, m), matrix(0, m, mr-m))

  if (mpq[2] > 0) {
    mT[1:(m*mpq[2]), 1:m] <- aperm(varma$ar, c(1, 3, 2))
  }
  if (form == "pearlman" && mpq[2] > 0)
    mR[1:(m*mpq[2]), ] <- aperm(varma$ar, c(1, 3, 2))
  if (mpq[3] > 0) {
    if (form == "harvey") {
      mR[(m+1):(m+m*mpq[3]), 1:m] <- aperm(varma$ma, c(1, 3, 2))
    } else {
      if (mpq[3] > 0) mR[1:(m*mpq[3]), ] <- mR[1:(m*mpq[3]), ] + as.numeric(aperm(varma$ma, c(1, 3, 2)))
    }
  }
  list(
    T = mT,
    R = mR,
    Q = mQ,
    Z = mZ
  )
}


#' Inverse roots of the characteristic polynomials of VARMA process
#'
#' It computes the inverse roots of the AR and MA parts and produces
#' plots of the roots in the complex plane, together with the unit circle.
#'
#' @param varma a varma object.
#' @param part a character vector indicating which polynomial equation has
#' to be solved: it can be `"ar"`, `"ma"` or both `c("ar", "ma")`.
#' @param plot logical, if TRUE, the function will plot the roots in the complex plane.
#'
#' @returns A list with the AR roots (`$ar`) and the MA roots (`$ma`).
#'
#' @importFrom graphics abline lines mtext par points polypath
#' @export
inv_roots <- function(varma, part = c("ar", "ma"), plot = FALSE) {
  # returns the inverse of the roots of the characteristic polynomial
  # of the varma object
  comp <- companion_form(varma, part)
  out <- list()
  out_len <- 0
  if (!is.null(comp$ar)) {
    out$ar <- eigen(comp$ar)$values
    out_len <- out_len + 1
  }
   if (!is.null(comp$ma)) {
    out$ma <- eigen(comp$ma)$values
    out_len <- out_len + 1
  }
  if (plot) {
    if (out_len == 2) oldpar <- par(mfrow = c(1, 2))
    if (any(part == "ar") & !is.null(comp$ar)) {
      plot(Re(out$ar), Im(out$ar), xlab = "Re", ylab = "Im",
           xlim = c(min(-1, Re(out$ar)), max(1, Re(out$ar))),
           ylim = c(min(-1, Im(out$ar)), max(1, Im(out$ar))),
           main = "Inverse AR roots",
           pch = 21, bg = ifelse(Mod(out$ar) < 1, "blue", "red"))
      abline(h = 0, v = 0)
      x <- seq(0, 2*pi, length.out = 100)
      lines(cos(x), sin(x))
    }
    if (any(part == "ma") & !is.null(comp$ma)) {
      plot(Re(out$ma), Im(out$ma), xlab = "Re", ylab = "Im",
           xlim = c(min(-1, Re(out$ma)), max(1, Re(out$ma))),
           ylim = c(min(-1, Im(out$ma)), max(1, Im(out$ma))),
           main = "Inverse MA roots",
           pch = 21, bg = ifelse(Mod(out$ma) < 1, "blue", "red"))
      abline(h = 0, v = 0)
      x <- seq(0, 2*pi, length.out = 100)
      lines(cos(x), sin(x))
    }
    if (out_len == 2) par(oldpar)
  }
  out
}

#' Check the stability of the VAR and VMA operators of a varma object
#'
#' @param varma an object of varma class
#'
#' @returns a logical vector with two entries: TRUE if ar/ma stable
#'
#' @export
is_stable <- function(varma) {
  out <- inv_roots(varma, part = c("ar", "ma"))
  ar <- if (is.null(out$ar)) TRUE else all(Mod(out$ar) < 1)
  ma <- if (is.null(out$ma)) TRUE else all(Mod(out$ma) < 1)
  c(ar = ar, ma = ma)
}

#' Check if the VARMA process is identified
#'
#' Tests uniqueness of the normalized AR and MA coefficients at the supplied
#' orders using the rank of linear restrictions implied by the transfer function.
#' Common scalar polynomial roots alone do not determine multivariate identification.
#'
#' @param varma a varma object.
#' @param tol relative singular-value tolerance for the numerical rank (default: 1e-5).
#'
#' @returns TRUE if the VARMA is identified and FALSE otherwise.
#'
#' @export
is_identified <- function(varma, tol = 1e-5) {
  mpq <- dim(varma); m <- mpq[1L]; p <- mpq[2L]; q <- mpq[3L]
  if (length(tol)!=1L || !is.finite(tol) || tol<=0) stop("tol must be positive")
  if (!p || !q) return(TRUE)
  # Equal transfers satisfy Delta A(L) H(L) = Delta B(L). The coefficients
  # above degree q therefore give linear restrictions on each row of Delta A.
  # Multiplication by det A shows that degrees through q+m*p suffice.
  psi <- irf(varma,maxlag=q+m*p,orth="none")
  blocks <- lapply(seq.int(q+1L,q+m*p),function(h) {
    do.call(rbind,lapply(seq_len(p),function(i)
      if (h-i < 0L) matrix(0,m,m) else matrix(psi[,,h-i+1L],m,m)))
  })
  values <- svd(do.call(cbind,blocks),nu=0,nv=0)$d
  if (!length(values) || max(values)==0) return(FALSE)
  sum(values > tol*max(values)) == m*p
}


#' Computes the Kronecker indices of a VARMA model
#'
#' @param varma an object of varma class (only the slots ar and ma are used)
#' @param tol relative singular-value threshold for numerical rank (default \code{1e-8})
#'
#' @returns A vector with number of elements equal to the number of time series
#' with the Kronecker indices.
#'
#' @export
kronecker_indices <- function(varma, tol = 1e-8) {
  mpq <- dim(varma)
  if ((mpq[2] < 1) && (mpq[3] < 1)) return(rep(0L, mpq[1]))
  phi   <- if (mpq[2] < 1) array(dim = c(mpq[1], mpq[1], 0)) else varma$ar
  theta <- if (mpq[3] < 1) array(dim = c(mpq[1], mpq[1], 0)) else varma$ma
  kronecker_indices_cpp(phi, theta, tol)
}

#' Impulse response function of VARMA process
#'
#' @param varma a varma object.
#' @param maxlag positive integer indicating the maximum lag for which the IRF
#' is computed.
#' `lagmax`.
#' @param orth a character string indicating the orthogonalization method to be used:
#' "none", "cholesky", "spectral", "generalized" are the available options.
#' @param plot logical, if TRUE, the function plots the IRFs.
#' @param title string with a title for the plot (default to NULL, that is, no title).
#' @param varnames vector of strings with the names of the variables.
#'
#'
#' @returns A 3D array with the impulse response functions of the VARMA process.
#'
#' @importFrom graphics abline lines mtext par points polypath
#' @export
irf <- function(varma, maxlag = 10,
                orth = c("none", "cholesky", "spectral", "generalized"),
                plot = FALSE,
                title = NULL,
                varnames = NULL) {
  orth = match.arg(orth)
  mpq <- dim.varma(varma)
  if (is.null(varnames)) {
    varnames <- colnames(varma$y)
    if (is.null(varnames)) {
      varnames <- paste0("y", 1:dim(varma)[1])
    } else {
      colnames(varma$y)
    }
  }
  if (orth == "cholesky") {
    P <- t(chol(varma$cov))
  } else if (orth == "spectral") {
    eig <- eigen(varma$cov)
    P <- eig$vectors %*% diag(sqrt(eig$values), nrow = mpq[1]) %*% t(eig$vectors)
  } else if (orth == "generalized") {
    P <- varma$cov
    P <- sweep(P, 2L, sqrt(diag(varma$cov)), "/")
    diag(P) <- sqrt(diag(varma$cov))
  } else {
    P <- diag(1, mpq[1], mpq[1])
  }
  # AR part
  ar_array <- if (!is.null(varma$ar)) {
    array(varma$ar, c(mpq[1], mpq[1], mpq[2]))
  } else {
    array(0, c(mpq[1], mpq[1], 0))
  }
  ma_array <- if (!is.null(varma$ma)) {
    array(varma$ma, c(mpq[1], mpq[1], mpq[3]))
  } else {
    array(0, c(mpq[1], mpq[1], 0))
  }
  psi <- irf_varma_rcpp(ar_array, ma_array, P, maxlag)
  if (plot) {
    old <- par(mfrow = c(mpq[1], mpq[1]),
               mar = c(2, 2, 2, 0.5),
               oma = c(1, 1, 2, 1)
               )
    for (i in 1:mpq[1]) {
      for (j in 1:mpq[1]) {
        plot(0:maxlag, psi[i, j, ], type = "p", pch = 20,
             xlab = "", ylab = "",
             main = paste(varnames[j], "\u2192", varnames[i]),
             ylim = c(min(psi), max(psi)))
        lines(0:maxlag, psi[i, j, ], col = "blue")
        abline(h = 0, lty = 2)
      }
    }
    if (!is.null(title)) mtext(title, side = 3, outer = TRUE)
    par(old)
  }
  dimnames(psi) <- list(varnames, varnames, paste0(0:maxlag))
  psi
}

#' Simulate VARMA process
#'
#' @param varma a varma object.
#' @param n number of observations to be simulated.
#' @param burn_in number of initial observations to discard.
#' @param var_in_rows boolean, if TRUE, the variables are in rows,
#' otherwise they are in columns (the internal Rcpp function works with
#' variables in rows, so TRUE makes the function slightly faster).
#'
#' @returns A matrix with the simulated data.
#'
#' @export
sim_varma <- function(varma, n = 200, burn_in = 100, var_in_rows = FALSE) {
  mpq <- dim(varma)
  ar_array <- if (!is.null(varma$ar)) {
    array(varma$ar, c(mpq[1], mpq[1], mpq[2]))
  } else {
    array(0, c(mpq[1], mpq[1], 0))
  }
  ma_array <- if (!is.null(varma$ma)) {
    array(varma$ma, c(mpq[1], mpq[1], mpq[3]))
  } else {
    array(0, c(mpq[1], mpq[1], 0))
  }
  # simulate Gaussian WN with the right covariance matrix
  if (!is.null(varma$chol_cov)) {
    eps <- crossprod(varma$chol_cov, matrix(stats::rnorm((n+burn_in)*mpq[1]), mpq[1], n+burn_in))
  } else {
    P <- chol(varma$cov)
    eps <- crossprod(P, matrix(stats::rnorm((n+burn_in)*mpq[1]), mpq[1], n+burn_in))
  }
  # simulate the VARMA process (zero-mean; mean is added below)
  Y <- sim_varma_rcpp(ar_array, ma_array, eps)
  # add the process mean: mu = (I - A_1 - ... - A_p)^{-1} intercept
  if (!is.null(varma$intercept)) {
    I_minus_A <- diag(mpq[1])
    if (mpq[2] > 0) {
      for (i in seq_len(mpq[2])) I_minus_A <- I_minus_A - varma$ar[,,i]
    }
    Y <- Y + solve(I_minus_A, varma$intercept)
  } else if (!is.null(varma$mean)) {
    Y <- Y + varma$mean
  }
  # return simulated data in requested form
  if (burn_in > 0) { # if there is a burn-in perido
    if (var_in_rows) { # if vars are in rows
      return(Y[, -(1:burn_in), drop = FALSE])
    } else { # if vars are in columns
      return(t(Y[, -(1:burn_in), drop = FALSE]))
    }
  } else { # if there is no burn-in period
    if (var_in_rows) { # if vars are in rows
      return(Y)
    } else { # if vars are in columns
      return(t(Y))
    }
  }
}

#' Autocovariance function of VARMA process by vectorization and Kronecker product
#'
#' It computes the autocovariances using the state-space representation of the VARMA
#' process and the following formula for the variance of the state variable
#' \eqn{vec(\Gamma_0) = (I - T \otimes T)^{-1} vec(RQR')}
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_vk <- function(varma, maxlag = 10) { # based on vectorization and Kronecker product
  if (maxlag < 0) maxlag <- 0
  ss <- varma_to_ss(varma)
  mpq <- dim(varma)
  k <- dim(ss$T)[1]
  G <- array(0, c(k, k, maxlag + 1))
  G[,,1] <- matrix(
      solve(diag(1, k*k, k*k) - (ss$T %x% ss$T), as.numeric(ss$R %*% ss$Q %*% t(ss$R))),
    k, k)
  if (maxlag > 0) for (i in 1:maxlag) {
    G[,,i+1] <- ss$T %*% G[,,i]
  }
  G[1:mpq[1], 1:mpq[1], ]
}

#' Autocovariance function of VARMA process by vectorization, Kronecker product and
#' sparse matrices
#'
#' It computes the autocovariances using the state-space representation of the VARMA
#' process and the following formula for the variance of the state variable
#' \eqn{vec(\Gamma_0) = (I - T \otimes T)^{-1} vec(RQR')}
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_sp <- function(varma, maxlag = 10) { # based on vectorization and Kronecker product
  if (maxlag < 0) maxlag <- 0
  ss <- varma_to_ss(varma)
  mT <- Matrix::Matrix(ss$T, sparse = TRUE)
  mpq <- dim(varma)
  k <- dim(ss$T)[1]
  G <- array(0, c(k, k, maxlag + 1))
  Id <- Matrix::Diagonal(k*k)
  G[,,1] <- matrix(
    solve(Id - Matrix::kronecker(mT, mT), as.numeric(ss$R %*% tcrossprod(ss$Q, ss$R))),
    k, k
  )
  if (maxlag > 0) for (i in 1:maxlag) {
    G[,,i+1] <- as.matrix(mT %*% G[,,i])
  }
  G[1:mpq[1], 1:mpq[1], ]
}


# Autocovariance function of VARMA process
#
# It computes the autocovariances using the state-space representation of the VARMA process and
# the following formula for the variance of the state variable using the Lyapunov equation solution
# of the package netcontrol:
# \eqn{vec(\Gamma_0) = (I - T \otimes T)^{-1} vec(RQR')}
#
# @param varma a varma object.
# @param maxlag positive integer with the maximum lag for the autocovariance function.
#
# @returns Array with autocovariance matrices.
#
# @export
# autocov_ly <- function(varma, maxlag = 10) { # based on Lyapunov equation solution in package netcontrol
#   if (maxlag < 0) maxlag <- 0
#   ss <- varma_to_ss(varma)
#   mpq <- dim(varma)
#   k <- dim(ss$T)[1]
#   G <- array(0, c(k, k, maxlag + 1))
#   G[,,1] <- netcontrol::dlyap(t(ss$T), ss$R %*% ss$Q %*% t(ss$R))
#   if (maxlag > 0) for (i in 1:maxlag) {
#     G[,,i+1] <- ss$T %*% G[,,i]
#   }
#   G[1:mpq[1], 1:mpq[1], ]
# }

#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using the state-space representation of the VARMA process and
#' the formula of Giacomo Sbrana
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_gs <- function(varma, maxlag = 10) { # based on the method of Giacomo Sbrana
  if (maxlag < 0) maxlag <- 0
  autocov_gs_cpp(varma$ar, varma$ma, varma$cov, maxlag)
}

#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using the state-space representation of the VARMA process and
#' the formula of Mittnik (1990)
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_mi <- function(varma, maxlag = 10) { # based on the method of Mittnik
  if (maxlag < 0) maxlag <- 0
  mpq <- dim(varma)
  autocov_mi_cpp(varma$ar, varma$ma, varma$cov, maxlag)
}


#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using McElroy's (2017) algorithm.
#' The function is based on the original R code made available by the author.
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_mc <- function(varma, maxlag = 10) {
  if (maxlag < 0) maxlag <- 0
  autocov_mc_cpp(varma$ar, varma$ma, varma$cov, maxlag)
}

#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using Ansley's (1980) algorithm.
#' The function is based on the original R code made available by the author.
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_an <- function(varma, maxlag = 10) {
  if (maxlag < 0) maxlag <- 0
  autocov_an_cpp(varma$ar, varma$ma, varma$cov, maxlag)
}


#' VARMA estimation using the iterated Hannan-Rissanen method
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept logical, if TRUE, the model includes an intercept.
#' @param maxdiff relative tolerance for changes in fitted regression errors.
#' Convergence does not require a unique vector of AR and MA coefficients.
#' @param maxit maximum number of iterations.
#' @param r order of the pilot VAR. The default grows logarithmically with sample
#' size and is capped to leave residual degrees of freedom.
#' @param ret a character string indicating the type of output: "varma" or "regression".
#'
#' @returns If ret = "varma" a varma object, otherwise a list with the following elements:
#' coefficients, regression residuals, fitted.values, effects, weights, rank,
#' df.residual, qr, n_iter, converged, rows, and res_cov. A varma object instead
#' stores recursive innovations computed from its final coefficients, with initial
#' zeros marked by residuals_start, and convergence/iterations diagnostics.
#'
#' @export
fit_varma_ihr <- function(Y, p = 1, q = 1, intercept = TRUE,
                          maxdiff = 1e-5, maxit = 100,
                          r = .varma_long_order(as.matrix(Y), p, q, intercept),
                          ret = c("varma", "regression")) {
  fn_name <- "fit_varma_ihr"
  Y <- .varma_data(Y, p, q)
  ret <- match.arg(ret)
  maxit <- .varma_integer(maxit, "maxit")
  if (length(maxdiff) != 1L || !is.finite(maxdiff) || maxdiff <= 0) stop("maxdiff must be positive")
  m <- ncol(Y); n <- nrow(Y)
  Ylag <- .varma_lags(Y, p)
  reg_names <- c(if (intercept) "cnst",
                 if (p) paste0(rep(paste0("y", seq_len(m)), p), "_", rep(seq_len(p), each=m)),
                 if (q) paste0(rep(paste0("e", seq_len(m)), q), "_", rep(seq_len(q), each=m)))
  if (q == 0L) {
    rows <- seq.int(p+1L, n)
    X <- if (intercept) cbind(1, Ylag[rows, , drop=FALSE]) else Ylag[rows, , drop=FALSE]
    if (nrow(X) <= ncol(X)) stop("Not enough observations for the requested VAR")
    reg <- .varma_regression(X, Y[rows, , drop=FALSE])
    it <- 0L; converged <- TRUE
  } else {
    r <- .varma_integer(r, "r")
    if (r < max(p,q) || n-r <= m*r + as.integer(intercept))
      stop("r must be at least max(p,q) and leave residual degrees of freedom in the pilot VAR")
    pilot_rows <- seq.int(r+1L,n)
    XP <- .varma_lags(Y, r)[pilot_rows, , drop=FALSE]
    if (intercept) XP <- cbind(1, XP)
    pilot <- .varma_regression(XP, Y[pilot_rows, , drop=FALSE])
    E <- matrix(0,n,m); E[pilot_rows, ] <- pilot$residuals
    start <- r + max(p,q)
    if (n-start <= m*(p+q)+as.integer(intercept)) stop("Not enough observations after pilot trimming")
    rows <- seq.int(start+1L,n)
    converged <- FALSE
    response <- Y[rows,,drop=FALSE]
    X <- matrix(0,length(rows),as.integer(intercept)+m*(p+q))
    if (intercept) X[,1L] <- 1
    if (p) X[,as.integer(intercept)+seq_len(m*p)] <- Ylag[rows,,drop=FALSE]
    for (it in seq_len(maxit)) {
      for (j in seq_len(q))
        X[,as.integer(intercept)+m*p+(j-1L)*m+seq_len(m)] <- E[rows-j,,drop=FALSE]
      reg <- .varma_regression(X,response)
      # Compare fitted errors, which remain defined when coefficients are aliased.
      change <- max(abs(reg$residuals-E[rows, , drop=FALSE]))/(1+max(abs(E[rows, , drop=FALSE])))
      E[rows, ] <- reg$residuals
      if (it > 1L && change < maxdiff) { converged <- TRUE; break }
    }
  }
  dimnames(reg$coefficients) <- list(reg_names, colnames(Y))
  reg$n_iter <- it; reg$converged <- converged; reg$rows <- rows
  reg$res_cov <- crossprod(reg$residuals)/nrow(reg$residuals)
  if (ret == "regression") { reg$design <- X; return(reg) }
  fit <- .varma_from_coefficients(reg$coefficients,p,q,intercept,m)
  .varma_finish(fit,Y,fn_name,length(reg$coefficients),rank=reg$rank,
                convergence=if (converged) 0L else 1L,iterations=it)
}


#' VARMA estimation using the iterated Hannan-Rissanen method with regularization
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept boolean, if TRUE, the model will include an intercept.
#' @param lambda NULL for a cross-validated path, a non-negative scalar for a fixed
#' penalty, or a vector of non-negative candidate penalties for cross-validation.
#' @param alpha value in the interval [0, 1]; 1 for pure LASSO and 0 for pure RIDGE.
#' @param relax logical, if FALSE the classical elastic net is used, if TRUE (default)
#' the estimations and cross validation are done also using regression coefficients
#' computed as \eqn{(1-\gamma)\beta^{OLS}_i + \gamma \beta_i}, where \eqn{i} ranges
#' over the set of non-zero coefficients and \eqn{\beta^{OLS}_i} denotes the
#' regression coefficient computed by OLS (using only the active regressors).
#' @param gamma vector of weights in the range \eqn{[0, 1]} to be used when.
#' `relax = TRUE`. The default is 0, corresponding to Post-LASSO OLS.
#' @param sur logical, if TRUE a final seemingly unrelated regressions
#' step is computed, otherwise the equation by equation estimate is provided.
#' @param maxdiff convergence tolerance passed to `fit_varma_ihr`.
#' @param maxit maximum number of iterations to be used in the first step
#' based on `fit_varma_ihr`.
#' @param r integer order of the first step VAR(r) model used in `fit_varma_ihr`.
#' @param ret a character string indicating the type of output: "varma" or "regression".
#' @param parallel TRUE if parallelization is required
#' @param n_cores integer number of cores to use, if parallel = TRUE
#'
#' @returns A varma object with refreshed recursive innovations and covariance, or
#' a list of glmnet fits (fixed lambda) or cv.glmnet fits (a penalty path).
#' SUR uses the selected support. The parameter count includes innovation covariance
#' parameters and is a nominal count, not a general effective-df formula for penalized fits.
#' @importFrom stats coef
#' @export
fit_varma_net <- function(Y, p=1, q=1, intercept = TRUE,
                          lambda = NULL, alpha = 1, relax = TRUE, gamma = 0,
                          sur = TRUE, maxdiff = 1e-5, maxit = 100,
                          r = .varma_long_order(as.matrix(Y), p, q, intercept),
                          parallel = FALSE, n_cores = NULL,
                          ret = c("varma", "regression")) {
  Y <- .varma_data(Y,p,q); ret <- match.arg(ret)
  m <- ncol(Y)
  ihr <- fit_varma_ihr(Y,p,q,intercept,maxdiff,maxit,r,ret="regression")
  X <- ihr$design
  if (intercept) X <- X[, -1L, drop=FALSE]
  YY <- Y[ihr$rows, , drop=FALSE]
  fixed <- !is.null(lambda) && length(lambda)==1L
  if (!is.null(lambda) && (!length(lambda) || any(!is.finite(lambda)) || any(lambda < 0)))
    stop("lambda must contain non-negative finite values")
  if (fixed && relax && length(gamma)!=1L) stop("A fixed lambda requires a scalar gamma")
  Xnet <- if (ncol(X)==1L) cbind(X,0) else X
  fit_equation <- function(i) {
    args <- list(x=Xnet,y=YY[,i],family="gaussian",lambda=lambda,
                 alpha=alpha,intercept=intercept,relax=relax)
    # glmnet 4.1-10 drops coefficient dimensions on a one-lambda relaxed path.
    # Retain two path points, but always extract the user's exact fixed lambda.
    if (fixed && relax) args$lambda <- c(lambda+max(1e-8,lambda*1e-6),lambda)
    if (!fixed) args$gamma <- gamma
    # glmnet's relaxed refit re-evaluates its stored call in another frame.
    # Supply values, so it does not depend on symbols in this local closure.
    do.call(if (fixed) glmnet::glmnet else glmnet::cv.glmnet,args)
  }
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores()-1L)
    n_cores <- .varma_integer(n_cores,"n_cores")
    old_plan <- future::plan(); on.exit(future::plan(old_plan),add=TRUE)
    future::plan(future::multisession,workers=n_cores)
    regreg <- future.apply::future_lapply(seq_len(m),fit_equation,future.seed=TRUE)
  } else regreg <- lapply(seq_len(m),fit_equation)
  if (ret == "regression") return(regreg)
  coeff <- do.call(cbind,lapply(regreg,function(reg) {
    cc <- if (fixed && relax) stats::coef(reg,s=lambda,gamma=gamma) else
      if (fixed) stats::coef(reg,s=lambda) else stats::coef(reg)
    as.numeric(cc)[seq_len(ncol(X)+1L)]
  }))
  if (!intercept) coeff <- coeff[-1L, , drop=FALSE]
  dimnames(coeff) <- dimnames(ihr$coefficients)
  XX <- if (intercept) cbind(1,X) else X
  resid <- YY-XX %*% coeff
  if (sur) {
    active <- coeff != 0
    if (intercept) active[1L, ] <- TRUE
    xlist <- lapply(seq_len(m),function(i) XX[,active[,i],drop=FALSE])
    ylist <- lapply(seq_len(m),function(i) YY[,i,drop=FALSE])
    if (any(active)) coeff[active] <- sur_cpp(xlist,ylist,crossprod(resid)/nrow(resid))
  }
  fit <- .varma_from_coefficients(coeff,p,q,intercept,m)
  fit <- .varma_finish(fit,Y,"fit_varma_net",sum(coeff != 0),
                       convergence=if (ihr$converged) 0L else 1L,iterations=ihr$n_iter)
  fit$regression_residuals <- YY-XX %*% coeff
  fit$regression_rows <- ihr$rows
  fit
}


#' VARMA estimation using Maximum-Likelihood in State-Space form using KFAS
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept logical, if TRUE, the model will include an intercept.
#' @param maxit maximum iterations of the pilot Hannan--Rissanen estimator.
#' @param control list passed to stats::optim; control$maxit controls likelihood iterations.
#' @details The Gaussian likelihood uses a stationary initial covariance recomputed
#' at each candidate. Optimisation is restricted to stable AR and invertible MA
#' polynomials. The returned optim.out records numerical convergence; a local optimum
#' is not a guarantee of global maximisation. No echelon restrictions are imposed.
#' @param ret a character string indicating the type of output: "varma" or
#' "kfas" (a list with model, the fitted SSModel, and optim.out).
#'
#' @returns A varma object, or a list containing model and optim.out. Stored errors
#' are multivariate prediction errors from KFAS::mvInnovations.
#'
#' @importFrom KFAS SSMcustom SSModel fitSSM KFS
#' @export
fit_varma_kfas <- function(Y, p = 1, q = 1, intercept = TRUE, maxit = 100,
                           ret = c("varma", "kfas"), control = list()) {
  ret <- match.arg(ret)
  out <- .varma_ml(Y,p,q,intercept,maxit,control,"kfas")
  if (ret == "varma") out$fit else out$raw
}


#' VARMA estimation using maximum likelihood in state-space form (FKF.SP)
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept boolean, if TRUE, the model will include an intercept.
#' @param maxit maximum iterations of the pilot Hannan--Rissanen estimator.
#' @param control list passed to stats::optim; control$maxit controls likelihood iterations.
#' @details The Gaussian likelihood uses a stationary initial covariance recomputed
#' at each candidate. Optimisation is restricted to stable AR and invertible MA
#' polynomials. The returned optim.out records numerical convergence; a local optimum
#' is not a guarantee of global maximisation. No echelon restrictions are imposed.
#' @param ret a character string indicating the type of output: "varma" (default)
#' or "ss" (a list with the state-space form matrices and fkf function output).
#'
#' @returns A varma object or a list with the matrices of the state-space form,
#' according to the choice of the `ret` parameter.
#'
#' @export
fit_varma_fkf <- function(Y, p = 1, q = 1, intercept = TRUE, maxit = 100,
                           ret = c("varma", "ss"), control = list()) {
  ret <- match.arg(ret)
  out <- .varma_ml(Y,p,q,intercept,maxit,control,"fkf")
  if (ret == "varma") out$fit else out$raw
}


#' Fit a VARMA(p,q) Model using Sparse Kalman Filter
#'
#' Estimates the parameters of a VARMA(p,q) model by maximizing the
#' likelihood computed with an efficient C++/Armadillo-based Kalman filter.
#'
#' @param Y a numeric matrix or ts object with d columns (series) and N rows (observations).
#' @param p the autoregressive order.
#' @param q the moving average order.
#' @param intercept logical: if TRUE (default) a vector of intercepts is computed.
#' @param maxit maximum iterations of the pilot Hannan--Rissanen estimator.
#' @param control list passed to stats::optim; control$maxit controls likelihood iterations.
#' @details The Gaussian likelihood uses a stationary initial covariance recomputed
#' at each candidate. Optimisation is restricted to stable AR and invertible MA
#' polynomials. The returned optim.out records numerical convergence; a local optimum
#' is not a guarantee of global maximisation. No echelon restrictions are imposed.
#' @return A varma object.
#' @export
fit_varma_skf <- function(Y, p = 1, q = 1, intercept = TRUE, maxit = 100,
                           control = list()) {
  out <- .varma_ml(Y,p,q,intercept,maxit,control,"skf")
  out$fit
}


#' Fit a VARMA(p,q) Model using the iterated GLS method
#'
#' Fits a conditional Gaussian criterion using covariance-weighted Gauss--Newton
#' steps and a backtracking line search. The full innovation derivative retains
#' cross-equation terms and the order of matrix multiplication. The intercept is
#' estimated jointly. The returned convergence code is zero when tolerance is met.
#'
#' @param Y a numeric matrix or ts object with d columns (series) and N rows (observations).
#' @param p the autoregressive order.
#' @param q the moving average order.
#' @param intercept logical: if TRUE (default) a vector of intercepts is computed.
#' @param maxit integer, maximum number of GLS iterations.
#' @param tol tolerance to determine convergence.
#' @param r order of the long VAR(r) for approximating the residuals at first step
#' @param verbose logical, if TRUE some textual output is printed during execution
#' @return A varma object.
#' @export
fit_varma_rby <- function(Y, p = 1, q = 1, intercept = TRUE,
                          maxit = 500, tol = 1e-5,
                          r = .varma_long_order(as.matrix(Y),p,q,intercept),
                          verbose = FALSE) {
  Y <- .varma_data(Y,p,q); maxit <- .varma_integer(maxit,"maxit")
  if (!is.finite(tol) || length(tol)!=1L || tol<=0) stop("tol must be positive")
  n <- nrow(Y); m <- ncol(Y)
  if (!q) {
    out <- fit_varma_ihr(Y,p,0,intercept)
    out$estimation_method <- "fit_varma_rby"
    return(out)
  }
  pilot <- fit_varma_ihr(Y,p,q,intercept,maxit=1L,r=r,ret="regression")
  fit <- .varma_from_coefficients(pilot$coefficients,p,q,intercept,m)
  fit$cov <- pilot$res_cov
  for (part in c("ar","ma")) if (length(fit[[part]])) {
    for (i in seq_len(1000L)) {
      if (all(Mod(inv_roots(fit,part)[[part]]) < .98)) break
      fit[[part]] <- fit[[part]]*.9
    }
  }
  aa <- if (p) fit$ar else array(0,c(m,m,0))
  bb <- fit$ma; cc <- if (intercept) fit$intercept else numeric(0)
  width <- m*(p+q)+as.integer(intercept)
  step <- function(a,b,c) rby_optimization_step(Y,a,b,p,q,intercept,c)
  res <- step(aa,bb,cc)
  objective <- function(e) {
    e <- e[seq.int(max(p,q)+1L,n),,drop=FALSE]
    .varma_gaussian_ll(e,crossprod(e)/nrow(e))
  }
  value <- objective(res$E); history <- value
  mu <- 0; convergence <- 1L
  for (iter in seq_len(maxit)) {
    accepted <- FALSE
    # ZE is the full covariance-weighted score, with cross-equation terms.
    if (max(abs(res$ZE)) < tol) { convergence <- 0L; break }
    for (attempt in seq_len(12L)) {
      delta <- tryCatch(solve(res$ZZ+diag(mu,nrow(res$ZZ)),res$ZE),error=function(e) NULL)
      if (!is.null(delta) && all(is.finite(delta))) {
        delta <- matrix(delta,width,m)
        scale <- 1
        for (ls in seq_len(25L)) {
          at <- aa; bt <- bb
          for (i in seq_len(p)) at[,,i] <- aa[,,i]+scale*t(delta[(i-1L)*m+seq_len(m),,drop=FALSE])
          for (j in seq_len(q)) bt[,,j] <- bb[,,j]+scale*t(delta[m*p+(j-1L)*m+seq_len(m),,drop=FALSE])
          ct <- if (intercept) cc+scale*delta[width,] else numeric(0)
          candidate <- varma(ar=if (p) at else NULL,ma=bt,cov=fit$cov)
          if (all(is_stable(candidate))) {
            # Only accepted steps need a new derivative/normal-equation system.
            trial_errors <- varma_residuals(Y,at,bt,ct)
            val <- objective(trial_errors)
            if (is.finite(val) && val > value) { accepted <- TRUE; break }
          }
          scale <- scale/2
        }
      }
      if (accepted) break
      mu <- if (mu==0) max(1e-8,1e-4*mean(diag(res$ZZ))) else mu*10
    }
    if (!accepted) break
    aa <- at; bb <- bt; cc <- ct; res <- step(aa,bb,cc); value <- val
    history <- c(history,value); mu <- mu/10
    if (verbose) cat("Iteration",iter,"conditional log-likelihood",value,"\n")
    if (max(abs(scale*delta)) < tol) { convergence <- 0L; break }
  }
  fit <- varma(intercept=if (intercept) cc else NULL,ar=if (p) aa else NULL,ma=bb)
  fit <- .varma_finish(fit,Y,"fit_varma_rby",m*width,
                       convergence=convergence,iterations=iter)
  fit$loglik_history <- history
  fit
}


#' Fit a VARMA(p,q) in diagonal MA form using feasible GLS
#'
#' Computes a pilot long VAR followed by one feasible GLS regression with diagonal
#' or scalar MA restrictions. This implementation is not an iterated likelihood
#' maximiser. With an intercept, both stages use centred observations.
#'
#' @param Y a numeric matrix or ts object with d columns (series) and N rows (observations).
#' @param p the autoregressive order.
#' @param q the moving average order.
#' @param intercept logical: if TRUE (default) a vector of intercepts is computed.
#' @param r order of the long VAR(r) for approximating the residuals at first step
#' @param verbose logical, if TRUE some textual output is printed during execution
#'
#' @return A varma object.
#' @export
fit_varma_dpd <- function(Y, p = 1, q = 1, intercept = TRUE,
                          r = .varma_long_order(as.matrix(Y),p,q,intercept),
                          verbose = FALSE) {
  .varma_dp(Y,p,q,intercept,r,"diagonal")
}


#' Fit a VARMA(p,q) in final MA form using feasible GLS
#'
#' Computes a pilot long VAR followed by one feasible GLS regression with diagonal
#' or scalar MA restrictions. This implementation is not an iterated likelihood
#' maximiser. With an intercept, both stages use centred observations.
#'
#' @param Y a numeric matrix or ts object with d columns (series) and N rows (observations).
#' @param p the autoregressive order.
#' @param q the moving average order.
#' @param intercept logical: if TRUE (default) a vector of intercepts is computed.
#' @param r order of the long VAR(r) for approximating the residuals at first step
#' @param verbose logical, if TRUE some textual output is printed during execution
#'
#' @return A varma object.
#' @export
fit_varma_dpf <- function(Y, p = 1, q = 1, intercept = TRUE,
                          r = .varma_long_order(as.matrix(Y),p,q,intercept),
                          verbose = FALSE) {
  .varma_dp(Y,p,q,intercept,r,"final")
}


#' Cast VARMA estimates from MTS::VARMA into a varma object
#'
#' It takes the output of the function VARMA of the MTS package
#' and turns it into a varma object.
#'
#' @param mts_varma list returned by MTS::VARMA.
#'
#' @returns A varma object.
#'
#' @export
convert_mts_varma <- function(mts_varma) {
  m <- dim(mts_varma$Sigma)[1]
  p <- mts_varma$ARorder
  q <- mts_varma$MAorder
  structure(
    list(
      intercept = mts_varma$Ph0,
      ar = if (p>0) array(mts_varma$Phi, c(m, m, p)) else NULL,
      ma = if (q>0) -array(mts_varma$Theta, c(m, m, q)) else NULL,
      cov = mts_varma$Sigma
    ),
    class = "varma"
  )
}


#' Cast VARMA estimates from fable::VARIMA into a varma object
#'
#' It takes the output of the function VARMA of the fable package
#' and turns it into a varma object.
#'
#' @param object object returned by a code like
#' `model(Y, VARIMA(vars(V1, V2, V3) ~ pdq(1, 0, 2), identification = "kronecker_indices"))`.
#'
#' @returns A varma object.
#'
#' @export
convert_fable_varma <- function(object) {
  mod <- object[[1]][[1]]
  m <- nrow(mod$fit$Phi)
  vnames <- colnames(mod$data)
  vnames <- setdiff(vnames, attr(mod$data, "index"))
  y <- mod$fit$data
  usable <- stats::complete.cases(mod$fit$residuals)
  errors <- as.matrix(mod$fit$residuals)[usable,,drop=FALSE]
  nobs <- nrow(errors)

  structure(
    list(
      intercept = if (mod$fit$cnst) mod$fit$const else NULL,
      ar = if (is.null(mod$fit$Phi)) NULL else array(mod$fit$Phi,   c(m, m, ncol(mod$fit$Phi)/m)),
      ma = if (is.null(mod$fit$Theta)) NULL else array(-mod$fit$Theta, c(m, m, ncol(mod$fit$Theta)/m)),
      cov = mod$fit$Sigma,
      estimation_method = paste0("fable_varima_", mod$fit$identification),
      loglik = .varma_gaussian_ll(errors,mod$fit$Sigma),
      n = nrow(y),
      nobs = nobs,
      npar = length(mod$fit$coef)+m*(m+1L)/2L,
      y = mod$fit$data,
      residuals = mod$fit$residuals
    ),
    class = "varma"
  )

}


#' Compute the distance between two IRFs
#'
#' It uses the Lp norm to compute the distance between
#' two IRFs (the order-0 IRF can be omitted)
#'
#' @param irf1 an array of dimensions m x m x max_lags containing an IRF.
#' @param irf2 an array of dimensions m x m x max_lags containing an IRF.
#' @param r positive number, order of the norm to be used (default is 2 -> Euclidean norm).
#' @param omit_lag0 logical, do not consider the IRFs at lag-0 (default TRUE).
#' applied to irf1
#'
#' @returns Distance between two IRFs (a positive number).
#'
#' @export
irf_distance <- function(irf1, irf2, r = 2, omit_lag0 = TRUE) {
  if (r <= 0) stop("r must be positive")
  if (any(dim(irf1) != dim(irf2))) stop("the dimensions of the two IRF are different")
  if (length(dim(irf1)) != 3 || length(dim(irf2)) != 3) stop("irf1 and irf2 must be 3-dimensional arrays")
  di <- if (omit_lag0) irf1[,, -1] - irf2[,, -1] else irf1 - irf2
  mean(abs(di)^r)^(1/r)
}

# Generate random stable VARMA parametrizations

#' Generate Coefficients for a Stationary and Invertible VARMA(p, q) Process
#'
#' This function generates random coefficient matrices for a m-dimensional
#' Vector Autoregressive Moving Average (VARMA) process of specified orders p and q.
#' The user can control the stability of the process by setting the maximum
#' modulus of the eigenvalues of the AR and MA companion matrices, ensuring
#' stationarity and invertibility respectively.
#'
#' @param m integer with the dimension of the time series process (m > 0).
#' @param p integer with the order of the autoregressive (AR) part (p >= 0).
#' @param q integer with the order of the moving average (MA) part (q >= 0).
#' @param max_eig_ar numeric value between 0 and 1 indicating the desired maximum
#'   eigenvalue modulus for the AR companion matrix to ensure stationarity.
#' @param max_eig_ma numeric value between 0 and 1 with the desired maximum
#'   eigenvalue modulus for the MA companion matrix to ensure invertibility.
#' @param dist function to generate the random coefficients.
#' @param max_try maximum number of trials to get a stationary VARMA model
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{ar}}{An m x m x p 3-dimensional array of AR coefficient matrices.}
#'   \item{\code{ma}}{An m x m x q 3-dimensional array of MA coefficient matrices.}
#' }
#'
#' @details
#' The function ensures stationarity (for the AR part) and invertibility (for
#' the MA part) by generating random coefficient matrices, constructing the
#' corresponding companion matrix, and then rescaling the coefficients. The
#' rescaling is done such that the largest eigenvalue modulus of the companion
#' matrix matches the user-specified value (`max_eig_ar` or `max_eig_ma`).
#' The method relies on the property that scaling the i-th coefficient matrix (A_i)
#' by s^i results in scaling the eigenvalues of the companion matrix by s.
#'
#' @export
#'
#' @examples
#' mod <- rvarma(m = 3, p = 2, q = 2, max_eig_ar = 0.95, max_eig_ma = 0.95)
#'
rvarma <- function(m, p, q, max_eig_ar = 0.9, max_eig_ma = 0.9,
                   dist = function(n) stats::runif(n, -1, 1),
                   max_try = 5) {
  if (max_try <= 0) stop("VARMA generation failed")

  #-----------------------------------------------------------------------------
  # Internal helper function to generate coefficients for one part (AR or MA)
  #-----------------------------------------------------------------------------
  .generate_coeffs <- function(m_dim, order, max_eig_mod, dist) {
    # Handle order = 0 case
    if (order == 0) {
      return(NULL)
    }

    # 1. Generate random initial coefficient matrices
    coeffs <- array(dist(m_dim*m_dim*order), c(m_dim, m_dim, order))

    # 2. Construct the companion matrix for the initial coefficients
    companion_matrix <- companion_form(list(ar = coeffs), part = "ar")$ar

    # 3. Calculate the maximum modulus of the eigenvalues
    current_max_mod <- max(Mod(eigen(companion_matrix, only.values = TRUE)$values))

    # 4. Rescale the coefficient matrices to match the desired max modulus
    if (current_max_mod < max_eig_mod) {
      # If process is already stable, no rescaling is needed.
      scaled_coeffs <- coeffs
    } else {
      # This scaling factor 's' is the factor by which we want to scale the eigenvalues
      s <- max_eig_mod / current_max_mod

      # To scale the eigenvalues of the companion matrix by 's',
      # we must scale the i-th coefficient matrix by 's^i'.
      scaled_coeffs <- coeffs * rep(s^(1:order), each = m_dim*m_dim)
    }

    return(scaled_coeffs)
  }
  #-----------------------------------------------------------------------------
  # End of helper function
  #-----------------------------------------------------------------------------

  # Input validation
  if (!is.numeric(m) || m < 1 || m %% 1 != 0) {
    stop("m (dimension) must be a positive integer.")
  }
  if (!is.numeric(p) || p < 0 || p %% 1 != 0) {
    stop("p (AR order) must be a non-negative integer.")
  }
  if (!is.numeric(q) || q < 0 || q %% 1 != 0) {
    stop("q (MA order) must be a non-negative integer.")
  }
  if (!is.numeric(max_eig_ar) || max_eig_ar <= 0 || max_eig_ar >= 1) {
    stop("max_eig_ar must be between 0 (exclusive) and 1 (exclusive).")
  }
  if (!is.numeric(max_eig_ma) || max_eig_ma <= 0 || max_eig_ma >= 1) {
    stop("max_eig_ma must be between 0 (exclusive) and 1 (exclusive).")
  }

  # Generate AR coefficients (for stationarity)
  ar <- if (p > 0) .generate_coeffs(m_dim = m, order = p, max_eig_mod = max_eig_ar, dist = dist) else NULL

  # Generate MA coefficients (for invertibility)
  ma <- if (q > 0) .generate_coeffs(m_dim = m, order = q, max_eig_mod = max_eig_ma, dist = dist) else NULL

  # Return a named list with the results
  out <- structure(
    list(
      ar = ar,
      ma = ma,
      cov = diag(m)
    ),
    class = "varma"
  )
  iroots <- inv_roots(out)
  if ((is.null(iroots$ar) || all(Mod(iroots$ar) < 1)) &&
      (is.null(iroots$ma) || all(Mod(iroots$ma) < 1))) {
    out
  } else {
    rvarma(m, p, q, max_eig_ar, max_eig_ma, dist, max_try = max_try - 1)
  }
}

#' logLik method for a varma object
#'
#' @param object a varma object.
#' @param ... not used.
#'
#' @returns A numeric scalar with the log-likelihood value and the
#' attribute df with the model's degrees of freedom.
#'
#' @importFrom stats logLik
#' @method logLik varma
#' @export
logLik.varma <- function(object, ...) {
  if (is.null(object$loglik)) {
    warning("There is no log-likelihood in the varma object")
    return(structure(NA_real_,df=NA_integer_,nobs=NA_integer_,class="logLik"))
  }
  k <- object[["npar",exact=TRUE]]
  if (is.null(k)) {
    m <- dim(object)[1L]
    k <- length(object$ar)+length(object$ma)+m*(m+1L)/2L+
      if (!is.null(object$intercept)) length(object$intercept) else length(object$mean)
  }
  structure(as.numeric(object$loglik),df=k,nobs=object$nobs,class="logLik")
}


#' nobs method for a varma object
#'
#' @param object a varma object.
#' @param ... not used.
#'
#' @returns Number of vector observations used in the likelihood, not the total
#' number of scalar entries across variables.
#'
#' @importFrom stats nobs
#' @method nobs varma
#' @export
nobs.varma <- function(object, ...) {
  if (is.null(object$nobs)) {
    warning("There is no information about the number of observations in the varma object")
    return(NA)
  }
  object$nobs
}


#' predict method for a varma object
#'
#' The varma object has to be estimated using some of the methods
#' in this package.
#'
#' @param object a varma object with the estimated model.
#' @param n.ahead integer forecast horizon.
#' @param cov logical, if TRUE the covariance matrices of the forecasts are produced.
#' @param ... not used.
#'
#' @returns If `cov = FALSE`, a matrix with the forecasts, if `cov = TRUE`
#' a list with the slot `mean`, containing a matrix with the point forecasts and
#' the slot `cov` with a 3D array with the cvariance matrices of the forecasts.
#'
#' @importFrom stats predict
#' @importFrom utils tail
#' @method predict varma
#' @export
predict.varma <- function(object, n.ahead = 1, cov = TRUE, ...) {
  n.ahead <- as.integer(n.ahead)
  if (n.ahead < 1) stop("n.ahead must be an integer larger than 0")
  if (is.null(object$estimation_method)) stop("This varma object was not estimated on data")
  if (is.null(object$y)) stop("No time series to predict in the varma object")
  # Check if VARMA is in state-space form: ssf TRUE or FALSE
  ssf <- (!is.null(object$state_pred_mean) &&
            !is.null(object$state_pred_cov)  &&
            !is.null(object$transition_matrix)  &&
            !is.null(object$disturbance_cov))
  mpq <- unname(dim(object))
  m <- mpq[1]
  p <- mpq[2]
  q <- mpq[3]
  # if (object$estimation_method %in% c("fit_varma_ihr", "fit_varma_net")) {
  if (!ssf && (!is.null(object$residuals))) {
    yhat <- matrix(0, n.ahead, m)
    colnames(yhat) <- colnames(object$y)
    # AR and MA contributions must be combined in a single pass: ylag must be
    # updated with the *full* h-step forecast (AR + MA + intercept) before it
    # is used as the lagged term for the (h+1)-step forecast.  Separate loops
    # would update ylag with the incomplete AR-only value, introducing an error
    # of A_1 M_1 a_n (and higher-order terms) at every step beyond h = 1.
    ny <- nrow(object$y)
    nr <- nrow(object$residuals)
    if (p > 0) {
      ylag <- as.numeric(t(object$y[ny:(ny - p + 1L), , drop = FALSE]))
      At   <- t(matrix(object$ar, m, m * p))
    }
    if (q > 0) {
      elag <- as.numeric(t(object$residuals[nr:(nr - q + 1L), , drop = FALSE]))
      Mt   <- t(matrix(object$ma, m, m * q))
    }
    for (i in seq_len(n.ahead)) {
      yi <- if (!is.null(object$intercept)) object$intercept else {
        cc <- .varma_process_mean(object)
        for (j in seq_len(p)) cc <- cc-matrix(object$ar[,,j],m,m) %*% .varma_process_mean(object)
        as.numeric(cc)
      }
      if (p > 0) yi <- yi + as.numeric(ylag %*% At)
      if (q > 0 && i <= q) {
        yi <- yi + as.numeric(elag %*% Mt)
        if (i < q) {                                    # shift for next step
          elag <- utils::head(elag, m * (q - i))
          Mt   <- utils::tail(Mt,   m * (q - i))
        }
      }
      yhat[i, ] <- yi
      if (p > 0)
        ylag <- if (p == 1L) yi else c(yi, ylag[seq_len((p - 1L) * m)])
    }
    if (cov) {
      Psi <- irf(object, maxlag = n.ahead - 1L, orth = "none")
      # Psi[,, j] = Psi_{j-1};  Psi[,, 1] = I_m
      S <- array(0, c(m, m, n.ahead))
      dimnames(S) <- list(colnames(yhat), colnames(yhat), NULL)
      S[,, 1] <- object$cov                          # Psi_0 = I, so I Sigma I' = Sigma
      for (h in seq_len(n.ahead - 1L) + 1L) {
        S[,, h] <- S[,, h - 1L] + matrix(Psi[,, h],m,m) %*% tcrossprod(object$cov, matrix(Psi[,, h],m,m))
      }
      return(list(mean = yhat, cov = S))
    } else {
      return(yhat)
    }
  }
  # if (object$estimation_method %in% c("fit_varma_fkf",
  #                                     "fit_varma_skf",
  #                                     "fit_varma_kfas")) {
  if (ssf) {
    d <- nrow(object$transition_matrix)
    ahat <- matrix(0, n.ahead, d)
    ahat[1, ] <- object$state_pred_mean
    if (n.ahead > 1){
      for (i in 2:n.ahead) {
        ahat[i, ] <- as.numeric(object$transition_matrix %*% ahat[i-1, ])
      }
    }
    yhat <- ahat[, 1:m, drop = FALSE]
    if (!is.null(object$mean)) yhat <- yhat + rep(object$mean, each = n.ahead)
    colnames(yhat) <- colnames(object$y)
    if (cov) {
      phat <- array(0, c(d, d, n.ahead))
      phat[,, 1] <- object$state_pred_cov
      if (n.ahead > 1) {
        for (i in 2:n.ahead) {
          phat[,, i] <- object$transition_matrix %*%
            tcrossprod(phat[,, i-1], object$transition_matrix) +
            object$disturbance_cov
        }
      }
      return(
        list(mean = yhat,
             cov  = phat[1:m, 1:m, , drop = FALSE]
        )
      )
    } else {
      return(yhat)
    }
  }
  warning("Not enough information in the varma object to produce forecasts")
}

#' Simulate from an estimated varma by bootstrapping residuals
#'
#' It produces sample paths of a varma process by bootstrapping the residuals
#' of an estimated varma. A function (usually a fit_varma function) is applied
#' to the simulated sample paths.
#'
#' @param varma a varma object containing residuals.
#' @param nsim number of sample paths to simulate (default: 100).
#' @param n length of the time series (default: the `n` in the varma object).
#' @param resids_to_skip non-negative integer witht he number of initial residuals to drop because
#' they may have a very large variance (default: 5). Zero skips no valid residuals;
#' rows marked as presample by residuals_start are always excluded.
#' @param burn_in number of simulated initial observations to discard (default: 100).
#' @param fit_fn function to apply to every simulated path: it is usually a
#' `fit_varma_???` function (default: the function that generated the estimates
#' in the varma object).
#' @param ... other parameters to be passed to `fit_fn`.
#' @param parallel logical, if TRUE, parallel processing is used (default: FALSE).
#' @param n_cores number of cores for parallel processing. If NULL, uses
#'   parallel::detectCores() - 1 (default: NULL).
#'
#' @details Residual vectors are centred and resampled jointly. The fitted process
#' mean is restored after simulation. The fitted model must be stable and invertible.
#' These implementation checks do not establish bootstrap validity at singular
#' or unidentified distributions.
#' @return A list with `nsim` slots containing the applications of the
#' function `fit_fn` to the simulated sample paths.
#' @export
bootstrap_varma <- function(varma, nsim = 100, n = varma$n, resids_to_skip = 5,
                           fit_fn = get(varma$estimation_method), ...,
                           parallel = FALSE, n_cores = NULL, burn_in = 100) {
  nsim <- .varma_integer(nsim,"nsim"); n <- .varma_integer(n,"n")
  burn_in <- .varma_integer(burn_in,"burn_in",0L)
  resids_to_skip <- .varma_integer(resids_to_skip,"resids_to_skip",0L)
  if (!all(is_stable(varma))) stop("Residual bootstrap requires a stable, invertible fitted model")
  res <- as.matrix(varma$residuals)
  if (!nrow(res)) stop("No residuals available")
  valid_start <- varma$residuals_start
  if (is.null(valid_start)) valid_start <- 1L
  skip <- max(resids_to_skip,valid_start-1L)
  if (skip >= nrow(res)) stop("No residuals remain after trimming")
  res <- res[seq.int(skip+1L,nrow(res)),,drop=FALSE]
  if (any(!is.finite(res))) stop("Residuals must be finite")
  res <- t(sweep(res,2,colMeans(res),"-"))
  mpq <- dim(varma); m <- mpq[1L]
  aa <- if (mpq[2L]) varma$ar else array(0,c(m,m,0))
  bb <- if (mpq[3L]) varma$ma else array(0,c(m,m,0))
  mu <- .varma_process_mean(varma)
  intercept <- !is.null(varma$intercept) || !is.null(varma$mean)
  dots <- list(...)
  refit_args <- list(p=mpq[2L],q=mpq[3L],intercept=intercept)
  if(!is.null(varma$refit_args) &&
     identical(fit_fn,get(varma$estimation_method,mode="function")))
    refit_args <- varma$refit_args
  refit_args <- utils::modifyList(refit_args,dots,keep.null=TRUE)
  run_simulation <- function(i) {
    yy <- sim_varma_rcpp(aa,bb,res[,sample.int(ncol(res),n+burn_in,replace=TRUE),drop=FALSE])
    yy <- sweep(yy[,burn_in+seq_len(n),drop=FALSE],1,mu,"+")
    do.call(fit_fn,c(list(Y=t(yy)),refit_args))
  }
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1L,parallel::detectCores()-1L)
    n_cores <- .varma_integer(n_cores,"n_cores")
    old_plan <- future::plan(); on.exit(future::plan(old_plan),add=TRUE)
    future::plan(future::multisession,workers=n_cores)
    result <- future.apply::future_lapply(seq_len(nsim),run_simulation,future.seed=TRUE)
  } else result <- lapply(seq_len(nsim),run_simulation)
  structure(result,generator=varma)
}


#' IRF from bootstrapped varma objects
#'
#' It produces IRF plots from bootstrapped varma
#' estimates produced by `bootstrap_varma()`.
#'
#' @param varma_list the list of bootstrapped varma objects
#' generated by `bootstrap_varma()`.
#' @param maxlag the maximum lag of the IRF.
#' @param orth orthogonalization method.
#' @param title title for the plot (default: NULL, that is, no title).
#' @param probs a vector with two values between 0 and 1 representing.
#' the percentiles to be represented by the fan in the plot.
#' @param varnames names of the variables (default: y1, y2, ...).
#' @param ... other parameters to be passed to the `irf()` function.
#'
#' @returns Nothing: it produces plots.
#'
#' @importFrom graphics abline lines mtext par points polypath
#' @importFrom grDevices rgb
#' @export
booted_irf <- function(varma_list, maxlag = 10,
                       orth = c("none", "cholesky", "spectral", "generalized"),
                       title = NULL,
                       probs = c(0.1, 0.9),
                       varnames = paste0("y", 1:dim(varma_list[[1]])[1]),
                       ...) {
  orth <- match.arg(orth)
  mpq <- dim(varma_list[[1]])
  arr <- array(0, c(mpq[1], mpq[1], maxlag+1))
  a4 <- vapply(varma_list, irf, arr, maxlag = maxlag, orth = orth, ...)
  gen <- attr(varma_list, "generator")
  if (!is.null(gen)) irf_gen <- irf(gen, maxlag = maxlag, orth = orth, ...)
  old <- par(mfrow = c(mpq[1], mpq[1]),
             mar = c(2, 2, 2, 0.5),
             oma = c(1, 1, 2, 1)
  )
  mycol <- rgb(0, 0, 255, maxColorValue = 255, alpha = 70)
  for (i in 1:mpq[1]) {
    for (j in 1:mpq[1]) {
      rng <- stats::quantile(a4[i, j, , ], c(min(0.01, probs[1]), max(0.99, probs[2])))
      plot(0:maxlag, rep(0, maxlag + 1), type = "l", lty = 2,
           xlab = "", ylab = "",
           main = paste(varnames[j], "\u2192", varnames[i]),
           ylim = rng)
      polypath(x = c(0:maxlag, maxlag:0),
               y = c(apply(a4[i, j, , ], 1, stats::quantile, probs = probs[1]),
                     apply(a4[i, j, , ], 1, stats::quantile, probs = probs[2])[(maxlag+1):1]),
               col = mycol,
               border = NA)
      lines(0:maxlag, apply(a4[i, j, , ], 1, stats::median), col = "blue")
      if (!is.null(gen)) {
        points(0:maxlag, irf_gen[i, j, ], pch = 20)
      }
      if (!is.null(title)) mtext(title, side = 3, outer = TRUE)
    }
  }
  par(old)
}


#' Selct the VARMA order
#'
#' Computes information criteria over AR/MA orders and intercept choices. Conditional
#' likelihoods use a common trimmed sample; stationary likelihoods use the full sample.
#' Innovation covariance parameters are counted. Counts are nominal: conventional
#' information-criterion theory need not apply at singular representations.
#'
#' @param Y matrix of time series.
#' @param p_max maximum order for the VAR part.
#' @param q_max maximum order for the VMA part.
#' @param intercept if FALSE assume no intercept is needed, otherwise.
#' computes the information criterion both with and without intercept.
#' @param ic information criterion to use among "AIC", "AICc", "BIC", "HQC".
#' @param fit_fn fit_varma_??? function to use.
#' @param ... parameters passed to the fit_varma_??? function
#'
#' @returns A data.frame with columns p, q, intercept, ic
#'
#' @export
select_varma_order <- function(Y, p_max = 5, q_max = 5, intercept = TRUE,
                               ic = c("AIC", "AICc", "BIC", "HQC"),
                               fit_fn = fit_varma_ihr, ...) {
  ic <- match.arg(ic)
  p_max <- .varma_integer(p_max,"p_max",0L); q_max <- .varma_integer(q_max,"q_max",0L)
  Y <- as.matrix(Y); n <- nrow(Y); m <- ncol(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)) || n <= max(p_max,q_max)) stop("Invalid Y or order bounds")
  gr <- expand.grid(p=0:p_max,q=0:q_max,intercept=if (intercept) c(FALSE,TRUE) else FALSE)
  fits <- vector("list",nrow(gr))
  for (i in seq_len(nrow(gr))) {
    if (gr$p[i]+gr$q[i]==0L) {
      mu <- if (gr$intercept[i]) colMeans(Y) else numeric(m)
      E <- sweep(Y,2,mu,"-"); S <- crossprod(E)/n
      f <- varma(mean=if (gr$intercept[i]) mu else NULL,cov=S,y=Y,residuals=E,
                 loglik=.varma_gaussian_ll(E,S),n=n,nobs=n,
                 npar=m*as.integer(gr$intercept[i])+m*(m+1L)/2L)
      f$likelihood_type <- "iid"; fits[[i]] <- f
    } else {
      fits[i] <- list(tryCatch(fit_fn(Y=Y,p=gr$p[i],q=gr$q[i],intercept=gr$intercept[i],...),
        error=function(e) {
          warning(sprintf("Could not fit p=%d, q=%d, intercept=%s: %s",gr$p[i],gr$q[i],gr$intercept[i],conditionMessage(e)))
          NULL
        }))
    }
  }
  conditional <- any(vapply(fits,function(f) !is.null(f) &&
    !identical(f$likelihood_type,"iid") && !identical(f$likelihood_type,"stationary"),logical(1)))
  rows <- seq.int(if (conditional) max(p_max,q_max)+1L else 1L,n)
  gr[[ic]] <- vapply(fits,function(f) {
    if (is.null(f)) return(NA_real_)
    ll <- f$loglik; nn <- f$nobs
    if (conditional) {
      E <- residuals.varma(f,Y)[rows,,drop=FALSE]
      ll <- .varma_gaussian_ll(E,f$cov); nn <- length(rows)
    }
    .varma_ic(as.numeric(ll),attr(logLik.varma(f),"df"),nn,ic)
  },numeric(1))
  gr$rank <- rank(gr[[ic]],na.last="keep")
  gr[order(gr$rank,na.last=TRUE),]
}
