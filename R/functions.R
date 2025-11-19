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
#' @param nobs total number of observations (if no missing observations are present, `nobs`\eqn{=nm}).
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
      M <- x$ar[,,i]
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
      M <- x$ma[,,i]
      dimnames(M) <- list(
      paste0("y", 1:mpq[1], "_t"),
      paste0("e", 1:mpq[1], "_[t-", i, "]")
    )
      print(M)
    }
  }
  cat("\nCovariance matrix:\n")
  dimnames(x$cov) <- list(
    paste0("e", 1:mpq[1]),
    paste0("e", 1:mpq[1])
  )
  print(x$cov)
  if (!is.null(x$estimation_method)) {
    cat("\nEstimated using", x$estimation_method, "\n")
  }
  if (!is.null(x$loglik)) {
    npar <- if (is.null(x$npar)) {
      mpq[1]*mpq[1]*mpq[2] + mpq[1]*mpq[1]*mpq[3] + mpq[1]*(mpq[1]+1)/2
    } else x$npar
    cat("N. of non-missing observations =", x$nobs, "\n")
    cat("N. of estimated parameters     =", npar, "\n")
    cat("log-lik =", x$loglik, "\n")
    cat("AIC     =", -2*x$loglik + 2*npar, "\n")
    cat("AICc    =", -2*x$loglik + 2*npar + 2*npar*(npar + 1)/(x$nobs - npar - 1), "\n")
    cat("BIC     =", -2*x$loglik + npar*log(x$nobs), "\n")
    cat("HQC     =", -2*x$loglik + npar*log(log(x$nobs)), "\n")
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
  if (r > 1) diag(mT[1:(mr-m), (m+1):mr]) <- 1
  mR <- if (form == "harvey") rbind(diag(m), matrix(0, mr-m, m)) else matrix(0, mr, m)
  mQ <- varma$cov
  mZ <- cbind(diag(1, m, m), matrix(0, m, mr-m))

  if (mpq[2] > 0) {
    mT[1:(m*mpq[2]), 1:m] <- aperm(varma$ar, c(1, 3, 2))
  }
  if (mpq[3] > 0) {
    if (form == "harvey") {
      mR[(m+1):(m+m*mpq[3]), 1:m] <- aperm(varma$ma, c(1, 3, 2))
    } else {
      if (mpq[2] > 0) mR[1:(m*mpq[2]), ] <- aperm(varma$ar, c(1, 3, 2))
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

#' Check if the VARMA process is identified
#'
#' The function checks if the VARMA process is identified by numerically checking if the
#' roots of the characteristic polynomials are distinct.
#'
#' @param varma a varma object.
#' @param tol a small numeric value indicating the tolerance for the numerical check
#' (default: 1e-5).
#'
#' @returns TRUE if the VARMA is identified and FALSE otherwise.
#'
#' @export
is_identified <- function(varma, tol = 1e-5) {
  if (is.null(varma$ar) || is.null(varma$ma) || length(varma$ar) == 0 || length(varma$ma) == 0) {
    return(TRUE)
  }
  mpq <- dim(varma)
  if (qr(
    cbind(
      matrix(varma$ar, mpq[1], mpq[1]*mpq[2]),
      matrix(varma$ma, mpq[1], mpq[1]*mpq[3])
    )
  )$rank < mpq[1]) return(FALSE)
  ir <- inv_roots(varma, c("ar", "ma"))
  dif <- outer(ir$ar, ir$ma, function(x, y) Mod(x - y))
  if (all(dif > tol)) TRUE else FALSE
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
                varnames = paste0("y", 1:dim(varma)[1])) {
  orth = match.arg(orth)
  mpq <- dim.varma(varma)
  if (orth == "cholesky") {
    P <- t(chol(varma$cov))
  } else if (orth == "spectral") {
    eig <- eigen(varma$cov)
    P <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  } else if (orth == "generalized") {
    P <- varma$cov
    P <- P / rep(sqrt(diag(varma$cov)), mpq[1])
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
  dimnames(psi) <- list(varnames, varnames, paste0("lag", 0:maxlag))
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
  if (!is.null(varma$intercept)) {
    eps <- eps + varma$intercept
  }
  # simulate the VARMA process
  Y <- sim_varma_rcpp(ar_array, ma_array, eps)
  # return simulated data in requested form
  if (burn_in > 0) { # if there is a burn-in perido
    if (var_in_rows) { # if vars are in rows
      return(Y[, -(1:burn_in)])
    } else { # if vars are in columns
      return(t(Y[, -(1:burn_in)]))
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


#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using the state-space representation of the VARMA process and
#' the following formula for the variance of the state variable using the Lyapunov equation solution
#' of the package netcontrol:
#' \eqn{vec(\Gamma_0) = (I - T \otimes T)^{-1} vec(RQR')}
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_ly <- function(varma, maxlag = 10) { # based on Lyapunov equation solution in package netcontrol
  if (maxlag < 0) maxlag <- 0
  ss <- varma_to_ss(varma)
  mpq <- dim(varma)
  k <- dim(ss$T)[1]
  G <- array(0, c(k, k, maxlag + 1))
  G[,,1] <- netcontrol::dlyap(t(ss$T), ss$R %*% ss$Q %*% t(ss$R))
  if (maxlag > 0) for (i in 1:maxlag) {
    G[,,i+1] <- ss$T %*% G[,,i]
  }
  G[1:mpq[1], 1:mpq[1], ]
}

#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using the state-space representation of the VARMA process and
#' the algorithm og Giacomo Sbrana
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_gs <- function(varma, maxlag = 10) { # based on the method of Giacomo Sbrana
  if (maxlag < 0) maxlag <- 0
  ss <- varma_to_ss(varma)
  mpq <- dim(varma)
  G <- array(0, c(dim(ss$T)[1], dim(ss$T)[2], maxlag + 1))
  if (mpq[2] > 0) {
    ei <- eigen(ss$T)
    if (all(ei$values == 0)) {
      G[,,1] <- solve_dlyap_iter(Matrix::Matrix(ss$T, sparse = TRUE), ss$R %*% ss$Q %*% t(ss$R))
    } else {
      K <- ss$T[, 1:mpq[1]] + rbind(ss$R[-(1:mpq[1]), ], matrix(0, mpq[1], mpq[1])) # ss$T %*% ss$R
      V <- ei$vectors
      iV <- solve(V)
      iVK <- iV %*% K
      A <- 1 - outer(ei$values, ei$values)
      G[,,1] <- ss$R %*% ss$Q %*% t(ss$R) + V %*% (iVK %*% ss$Q %*% t(iVK) / A) %*% t(V)}
  } else {
    G[,,1] <- ss$R %*% ss$Q %*% t(ss$R)
  }
  if (maxlag > 0) {
    for (i in 1:maxlag) {
      G[,,i+1] <- ss$T %*% G[,,i]
    }
  }
  Re(G[1:mpq[1], 1:mpq[1], ])
}


#' Autocovariance function of VARMA process
#'
#' It computes the autocovariances using Tucker McElroy's algorithm.
#' The function is based on the original R code made available by the author.
#'
#' @param varma a varma object.
#' @param maxlag positive integer with the maximum lag for the autocovariance function.
#'
#' @returns Array with autocovariance matrices.
#'
#' @export
autocov_mc <- function(varma, maxlag = 10) { # based on Tucker McElroy's code
  # My modification of the original function
  ##################################################################################
  #
  #  VARMAauto
  #	Copyright (2015) Tucker McElroy
  #
  #	This program is free software; you can redistribute it and/or
  #	modify it under the terms of the GNU General Public License
  #	as published by the Free Software Foundation; either version 2
  #	of the License, or (at your option) any later version.
  #
  #	This program is distributed in the hope that it will be useful,
  #	but WITHOUT ANY WARRANTY; without even the implied warranty of
  #	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #	GNU General Public License for more details.
  #
  #	You should have received a copy of the GNU General Public License
  #	along with this program; if not, write to the Free Software
  #	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
  #
  ##############################################################################

  #######################################################################
  #	DOCUMENTATION:
  # Function computes autocovariances of VARMA (p,q) from lag zero
  #	to maxlag, with inputs phi and theta.
  #	(1 - phi[1]z ... - phi[p]z^p) X_t = (1 + theta[1]z ...+ theta[q]z^q) WN
  #  output: autocovariance string of length maxlag
  #  for absent AR or MA portions, pass in NULL
  #  phi and theta should be arrays of m x m matrices
  #	sigma should be an m x m matrix
  #  e.g. phi <- array(cbind(phi1,phi2,...,phip),c(m,m,p))
  ################################################################

  polymulMat <- function(amat,bmat)
  {
    p <- dim(amat)[3]
    q <- dim(bmat)[3]
    m <- dim(amat)[2]
    amatd <- amat[,,p:1]
    if(m==1) amatd <- array(amatd,c(m,m,p))
    if(q > 1) amatd <- array(c(matrix(0,m,m*(q-1)),amatd),c(m,m,p+q-1))
    bigmat <- NULL
    for(i in 1:(p+q-1))
    {
      nextmat <- matrix(amatd[,,1:(p+q-1)],m,m*(p+q-1))
      bigmat <- rbind(nextmat,bigmat)
      amatd <- amatd[,,-1]
      amatd <- array(c(amatd,matrix(0,m,m)),c(m,m,p+q-1))
    }
    bigmat <- bigmat[,1:(m*q)]
    out <- bigmat %*% t(matrix(bmat[,,q:1],m,m*q))
    out <- array(out,c(m,p+q-1,m))
    temp <- NULL
    for(i in 1:(p+q-1))
    {
      temp <- cbind(temp,out[,i,])
    }
    out <- array(temp,c(m,m,p+q-1))

    return(out)
  }

  Kcommut <- function(vect,m,n)
  {
    return(matrix(t(matrix(vect,nrow=m,ncol=n)),ncol=1))
  }

  mpq <- dim(varma)
  m <- mpq[1]
  p <- mpq[2]
  q <- mpq[3]
  Kmat <- apply(diag(m^2),1,Kcommut,m,m)

  if (q == 0) { gamMA <- array(varma$cov,c(m,m,1)) } else
  {
    temp <- polymulMat(array(cbind(diag(m),matrix(varma$ma,m,m*q)),c(m,m,q+1)),
                       array(varma$cov,c(m,m,1)))
    gamMA <- polymulMat(temp,array(cbind(diag(m),matrix(varma$ma,m,m*q)),c(m,m,q+1)))
  }
  gamMA <- gamMA[,,(q+1):(2*q+1)]
  if(m==1) gamMA <- array(gamMA,c(m,m,q+1))
  gamMAvec <- matrix(gamMA,m^2*(q+1),1)

  if (p > 0)
  {
    Amat <- matrix(0,nrow=m^2*(p+1),ncol=m^2*(2*p+1))
    Amat <- array(Amat,c(m^2,p+1,m^2,2*p+1))
    Arow <- diag(m^2)
    for(i in 1:p)
    {
      Arow <- cbind(-1*diag(m) %x% varma$ar[,,i],Arow)
    }
    for(i in 1:(p+1))
    {
      Amat[,i,,i:(i+p)] <- Arow
    }
    newA <- array(matrix(Amat[,1:(p+1),,1:p],m^2*(p+1),m^2*(p)),c(m^2,p+1,m^2,p))
    for(i in 1:(p+1))
    {
      for(j in 1:p)
      {
        newA[,i,,j] <- newA[,i,,j] %*% Kmat
      }
    }
    Amat <- cbind(matrix(Amat[,,,p+1],m^2*(p+1),m^2),
                  matrix(Amat[,,,(p+2):(2*p+1)],m^2*(p+1),m^2*(p)) +
                    matrix(newA[,,,p:1],m^2*(p+1),m^2*(p)))

    Bmat <- matrix(0,nrow=m^2*(q+1),ncol=m^2*(p+q+1))
    Bmat <- array(Bmat,c(m^2,q+1,m^2,p+q+1))
    Brow <- diag(m^2)
    for(i in 1:p)
    {
      Brow <- cbind(Brow,-1*varma$ar[,,i] %x% diag(m))
    }
    for(i in 1:(q+1))
    {
      Bmat[,i,,i:(i+p)] <- Brow
    }
    Bmat <- Bmat[,,,1:(q+1)]
    Bmat <- matrix(Bmat,m^2*(q+1),m^2*(q+1))
    Binv <- solve(Bmat)

    gamMix <- Binv %*% gamMAvec
    if (p <= q) gamMixTemp <- gamMix[1:((p+1)*m^2)] else
      gamMixTemp <- c(gamMix,rep(0,(p-q)*m^2))
    gamARMA <- solve(Amat) %*% gamMixTemp
    gamMix <- array(matrix(gamMix,m,m*(q+1)),c(m,m,q+1))
    gamARMA <- array(matrix(gamARMA,m,m*(p+1)),c(m,m,p+1))
  } else
  {
    gamARMA <- array(gamMA[,,1],c(m,m,1))
    if (q == 0) { gamMix <- array(varma$cov,c(m,m,1)) } else
    {
      gamMix <- gamMA[,,1:(q+1)]
      if(m==1) gamMix <- array(gamMix,c(1,1,q+1))
    }
  }

  if (maxlag <= p)
  {
    gamARMA <- gamARMA[,,1:(maxlag+1)]
    if(m==1) gamARMA <- array(gamARMA,c(1,1,maxlag+1))
  } else
  {
    if (maxlag > q) gamMix <- array(cbind(matrix(gamMix,m,m*(q+1)),
                                          matrix(0,m,m*(maxlag-q))),c(m,m,(maxlag+1)))
    for(k in 1:(maxlag-p))
    {
      len <- dim(gamARMA)[3]
      acf <- gamMix[,,p+1+k]
      if (p > 0)
      {
        temp <- NULL
        for(i in 1:p)
        {
          temp <- rbind(temp,gamARMA[,,len-i+1])
        }
        acf <- acf + matrix(varma$ar,m,m*p) %*% temp
      }
      gamARMA <- array(cbind(matrix(gamARMA,m,m*len),acf),c(m,m,len+1))
    }
  }

  return(gamARMA)
}

#' VARMA estimation in Diagonal MA form using the method of Dufoura and Pelletierc
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept logical, if TRUE, the model includes an intercept.
#' @param r integer order of the first step VAR(r) model.
#' @param ret a character string indicating the type of output: "varma" or "regression".
#'
#' @returns If ret = "varma" a varma object, otherwise a list with the following elements:
#' coefficients, matrix of residuals, fitted.values, effects, weights, rank,
#' df.residual, qr, n_iter, res_cov.
#'
#' @export
fit_varma_dma <- function(Y, p=1, q=1, intercept = TRUE,
                          r = max(p+q, round(nrow(Y)/(4*ncol(Y)))),
                          ret = c("varma", "regression")) {
  fn_name <- as.character(sys.call()[[1]]) # name of the function
  # input controls
  if (p < 0) stop("p must be non-negative")
  if (q < 0) stop("q must be non-negative")
  if (p == 0 & q == 0) {"p = q = 0: no VARMA model specified"}
  ret <- match.arg(ret)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  # basic quantities and variable names
  m <- ncol(Y)
  n <- nrow(Y)
  reg_names <- NULL
  if (intercept) reg_names <- "cnst"
  if (p > 0) reg_names <- c(reg_names, paste(paste0("y", 1:m), rep(1:p, each = m), sep = "_"))
  if (q > 0) reg_names <- c(reg_names, paste(paste0("e", 1:m), rep(1:q, each = m), sep = "_"))
  # pure VAR
  if (q < 1) {
    X <- if (intercept) cbind(1, stats::embed(Y, p+1)[, -(1:m)]) else stats::embed(Y, p+1)[, -(1:m)]
    reg1 <- stats::lm.fit(X, Y[-(1:p),])
    reg1$n_iter <- 0
    dimnames(reg1$coefficients) <- list(
      reg_names,
      paste0("y", 1:m)
    )
    reg1$res_cov <- crossprod(reg1$residuals, reg1$residuals)/reg1$df.residual
    if (ret == "reg") {
      return(reg1)
    } else {
      estim <- t(reg1$coefficients)
      if (intercept) {
        cnst  <- estim[,  1]
        estim <- estim[, -1]
      } else {
        cnst <- NULL
      }
      ar_array <- array(estim[, 1:(p*m)], c(m, m, p))
      nobs <- prod(dim(reg1$residuals))
      return(
        structure(
          list(
            intercept = cnst,
            ar = ar_array,
            ma = NULL,
            cov = reg1$res_cov,
            estimation_method = fn_name,
            loglik = - (nobs / 2) * (log(2 * pi) + 1) -
              (nrow(reg1$residuals) / 2) * determinant(reg1$res_cov, logarithm = TRUE)$modulus,
            n = n,
            nobs = nobs,
            npar = length(reg1$coefficients),
            y = Y,
            residuals = reg1$residuals
          ),
          class = "varma"
        )
      )
    }
  }
  # VARMA in Diagonal MA form
  # 1st step: initial VAR(r)
  X <- if (intercept) cbind(1, stats::embed(Y, r+1)[, -(1:m)]) else stats::embed(Y, r+1)[, -(1:m)]
  m1 <- stats::lm.fit(X, Y[-(1:r),])
  E <- rbind(matrix(0, r, m), m1$residuals)
  Ylag <- if (p > 0) stats::embed(Y, p+1)[, -(1:m)] else NULL
  Elag <- rbind(matrix(0, q, m*q), stats::embed(E, q+1)[, -(1:m)])[(p+1):n, ]
  # 2nd step: SUR estimation
  Y_list <- as.list(as.data.frame(Y[-(1:p), ]))
  X_list <- if (intercept) {
    lapply(1:m, function(i) cbind(1, Ylag, Elag[, seq(i, q*m, m)]))
  } else {
    lapply(1:m, function(i) cbind(Ylag, Elag[, seq(i, q*m, m)]))
  }
  S <- crossprod(m1$residuals, m1$residuals)/n
  sur_est <- sur_cpp(X_list, Y_list, S)
  B <- matrix(sur_est, ncol = m) # each column has the estimates for each time series
  # 3rd step
  # Compute new residuals
  Etilde <- rbind(matrix(0, p, m),
    sapply(1:m, function(j) Y[-(1:p), j] - X_list[[j]] %*% B[, j])
  )
  SS <- crossprod(Etilde)/n
  # Compute filtered series
  Theta <- apply(-tail(B, q), 1, diag) |> array(c(m, m, q))
  YY <- var_filter(Y, Theta)
  W  <- var_filter(Etilde, Theta)
  YYlag <- if (p > 0) stats::embed(YY, p+1)[, -(1:m)] else NULL
  EElag <- rbind(matrix(0, q, m*q), stats::embed(Etilde, q+1)[, -(1:m)])[(p+1):n, ]
  YY_list <- as.list(as.data.frame((Etilde + YY - W)[-(1:p), ]))
  XX_list <- if (intercept) {
    lapply(1:m, function(i) cbind(1, YYlag, EElag[, seq(i, q*m, m)]))
  } else {
    lapply(1:m, function(i) cbind(YYlag, EElag[, i]))
  }
  sur_est2 <- sur_cpp(XX_list, YY_list, SS)
  BB <- matrix(sur_est2, ncol = m) # each column has the estimates for each time series
  # Compute new residuals
  Efinal <- rbind(matrix(0, p, m),
                  sapply(1:m, function(j) YY[-(1:p), j] - XX_list[[j]] %*% BB[, j])
  )
  SSS <- crossprod(Efinal)/n

  if (ret == "varma") {
    nobs <- length(Efinal) - p*m
    if (intercept) {
      ar_array <- BB[2:(1+m*p), , drop = FALSE] |> t() |> array(c(m, m, p))
      ma_array <- BB[-(1:(1+m*p)), , drop = FALSE] |> t() |> apply(1, diag) |> array(c(m, m, q))
    } else {
      ar_array <- BB[1:(m*p), , drop = FALSE] |> t() |> array(c(m, m, p))
      ma_array <- BB[-(1:(m*p)), , drop = FALSE] |> apply(1, diag) |> array(c(m, m, q))
    }
    structure(
      list(
        intercept = if (intercept) BB[1, ] else NULL,
        ar = ar_array,
        ma = ma_array,
        cov = SSS,
        estimation_method = fn_name,
        loglik = - (nobs / 2) * (log(2 * pi) + 1) -
          (nrow(Efinal) / 2) * determinant(SSS, logarithm = TRUE)$modulus,
        n = nrow(Y),
        nobs = nobs,
        npar = length(BB),
        y = Y,
        residuals = Efinal
      ),
      class = "varma"
    )
  } else {
    list(
      step2_est = B,
      step3_est = BB,
      step2_cov = SS,
      step3_cov = SSS,
      step2_resid = Etilde,
      step3_resid = Efinal,
      X = YY,
      W = W,
      Y_list = YY_list,
      X_list = XX_list
    )
  }
}

#' VARMA estimation using the iterated Hannan-Rissanen method
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept logical, if TRUE, the model includes an intercept.
#' @param maxdiff maximum difference between the coefficients of two consecutive iterations:
#' it is used as a threshold to stop the iterations.
#' @param maxit maximum number of iterations.
#' @param r integer order of the first step VAR(r) model.
#' @param ret a character string indicating the type of output: "varma" or "regression".
#'
#' @returns If ret = "varma" a varma object, otherwise a list with the following elements:
#' coefficients, matrix of residuals, fitted.values, effects, weights, rank,
#' df.residual, qr, n_iter, res_cov.
#'
#' @export
fit_varma_ihr <- function(Y, p=1, q=1, intercept = TRUE,
                          maxdiff = 1.0e-5, maxit = 100,
                          r = max(p+q, round(nrow(Y)/(4*ncol(Y)))),
                          ret = c("varma", "regression")) {
  fn_name <- as.character(sys.call()[[1]]) # name of the function
  # input controls
  if (p < 0) stop("p must be non-negative")
  if (q < 0) stop("q must be non-negative")
  if (p == 0 & q == 0) {"p = q = 0: no VARMA model specified"}
  ret <- match.arg(ret)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  # basic quantities and variable names
  m <- ncol(Y)
  n <- nrow(Y)
  reg_names <- NULL
  if (intercept) reg_names <- "cnst"
  if (p > 0) reg_names <- c(reg_names, paste(paste0("y", 1:m), rep(1:p, each = m), sep = "_"))
  if (q > 0) reg_names <- c(reg_names, paste(paste0("e", 1:m), rep(1:q, each = m), sep = "_"))
  # pure VAR
  if (q < 1) {
    X <- if (intercept) cbind(1, stats::embed(Y, p+1)[, -(1:m)]) else stats::embed(Y, p+1)[, -(1:m)]
    reg1 <- stats::lm.fit(X, Y[-(1:p),])
    reg1$n_iter <- 0
    dimnames(reg1$coefficients) <- list(
      reg_names,
      paste0("y", 1:m)
    )
    reg1$res_cov <- crossprod(reg1$residuals, reg1$residuals)/reg1$df.residual
    if (ret == "reg") {
      return(reg1)
    } else {
      estim <- t(reg1$coefficients)
      if (intercept) {
        cnst  <- estim[,  1]
        estim <- estim[, -1]
      } else {
        cnst <- NULL
      }
      ar_array <- array(estim[, 1:(p*m)], c(m, m, p))
      nobs <- prod(dim(reg1$residuals))
      return(
        structure(
          list(
            intercept = cnst,
            ar = ar_array,
            ma = NULL,
            cov = reg1$res_cov,
            estimation_method = fn_name,
            loglik = - (nobs / 2) * (log(2 * pi) + 1) -
              (nrow(reg1$residuals) / 2) * determinant(reg1$res_cov, logarithm = TRUE)$modulus,
            n = n,
            nobs = nobs,
            npar = length(reg1$coefficients),
            y = Y,
            residuals = reg1$residuals
          ),
          class = "varma"
        )
      )
    }
  }
  # initial VAR with order = (p+q)
  X <- if (intercept) cbind(1, stats::embed(Y, r+1)[, -(1:m)]) else stats::embed(Y, r+1)[, -(1:m)]
  m1 <- stats::lm.fit(X, Y[-(1:r),])
  E <- rbind(matrix(0, r, m), m1$residuals)
  # main iterations
  for (it in 1:maxit) {
    Ylag <- if (p > 0) stats::embed(Y, p+1)[, -(1:m)] else NULL
    Elag <- rbind(matrix(0, q, m*q), stats::embed(E, q+1)[, -(1:m)])[(p+1):n, ]
    if (intercept) X = cbind(1, Ylag, Elag) else X = cbind(Ylag, Elag)
    if (it > 1) reg0 <- reg1
    reg1 <- stats::lm.fit(X, Y[(p+1):n, ])
    if (it > 1) {
      if (max(abs(reg1$coefficients - reg0$coefficients), na.rm = TRUE) < maxdiff) break
    }
    reg0 <- reg1
    E <- rbind(matrix(0, p, m), reg1$residuals)
  }
  reg1$n_iter <- it
  dimnames(reg1$coefficients) <- list(
    reg_names,
    paste0("y", 1:m)
  )
  reg1$res_cov <- crossprod(reg1$residuals, reg1$residuals)/reg1$df.residual

  estim <- t(reg1$coefficients)
  if (intercept) {
    cnst  <- estim[,  1]
    estim <- estim[, -1]
  } else {
    cnst <- NULL
  }
  if (p > 0) {
    ar_array <- array(estim[, 1:(p*m)], c(m, m, p))
  } else {
    ar_array <- NULL
  }
  if (q > 0) {
    ma_array <- array(estim[, (p*m+1):(p*m+q*m)], c(m, m, q))
  } else {
    ma_array <- NULL
  }
  if (ret == "varma") {
    nobs <- prod(dim(reg1$residuals))
    structure(
      list(
        intercept = cnst,
        ar = ar_array,
        ma = ma_array,
        cov = reg1$res_cov,
        estimation_method = fn_name,
        loglik = - (nobs / 2) * (log(2 * pi) + 1) -
          (nrow(reg1$residuals) / 2) * determinant(reg1$res_cov, logarithm = TRUE)$modulus,
        n = nrow(Y),
        nobs = nobs,
        npar = length(reg1$coefficients),
        y = Y,
        residuals = reg1$residuals
      ),
      class = "varma"
    )
  } else {
    reg1
  }
}


#' VARMA estimation using the iterated Hannan-Rissanen method with regularization
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept boolean, if TRUE, the model will include an intercept.
#' @param lambda positive numeric regularization parameter.
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
#' @param maxdiff maximum difference between the coefficients of two consecutive
#' iterations in the first step based on `fit_varma_ihr`:
#' it is used as a threshold to stop the iterations.
#' @param maxit maximum number of iterations to be used in the first step
#' based on `fit_varma_ihr`.
#' @param r integer order of the first step VAR(r) model used in `fit_varma_ihr`.
#' @param ret a character string indicating the type of output: "varma" or "regression".
#'
#' @returns If ret = "varma" a varma object, otherwise a list with the results of
#' cv.glmnet for each response variable.
#'
#' @export
fit_varma_net <- function(Y, p=1, q=1, intercept = TRUE,
                          lambda = NULL, alpha = 1,
                          relax = TRUE, gamma = 0,
                          sur = TRUE,
                          maxdiff = 1.0e-5, maxit = 100,
                          r = max(p+q, round(nrow(Y)/(4*ncol(Y)))),
                          parallel = FALSE,
                          n_cores = NULL,
                          ret = c("varma", "regression")) {
  fn_name <- as.character(sys.call()[[1]]) # function name
  # input controls
  if (p < 0) stop("p must be non-negative")
  if (q < 0) stop("q must be non-negative")
  if (p == 0 & q == 0) {"p = q = 0: no VARMA model specified"}
  ret <- match.arg(ret)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  # basic quantities and variable names
  m <- ncol(Y)
  n <- nrow(Y)
  # calling ihr for first step
  ihr <-fit_varma_ihr(Y, p = p, q = q, intercept = intercept,
                      maxdiff = maxdiff, maxit = maxit,
                      r = r, ret = "regression")
  regreg <- vector("list", m)
  X <- if (intercept) qr.X(ihr$qr)[,-1] else qr.X(ihr$qr)
  for (i in 1:m) {
    regreg[[i]] <- glmnet::cv.glmnet(x = X,
                                     y = if (p>0) Y[-(1:p), i] else Y[, i],
                                     family = "gaussian",
                                     alpha = alpha,
                                     intercept = intercept,
                                     relax = relax,
                                     gamma = gamma)
  }
  coeff <- sapply(regreg, function(reg) as.numeric(coef(reg)))
  if (!intercept) coeff <- coeff[-1, ] # glmnet returns zero coeff also for the missing intercept
  dimnames(coeff) <- dimnames(ihr$coefficients)
  fitted <- sapply(regreg, function(reg) predict(reg, newx = X))
  resid  <- if (p>0) Y[-(1:p), ] - fitted else Y - fitted
  res_cov <- crossprod(resid)/nrow(X)
  nobs <- sum(!is.na(Y))
  estim <- t(coeff)
  if (sur) { # SUR step if required
    not_zero <- (estim != 0)
    if (intercept) X <- cbind(1, X)
    ylist <- lapply(1:m, function(i) if (p>0) Y[-(1:p), i, drop = FALSE] else Y[, i, drop = FALSE])
    xlist <- lapply(1:m, function(i) X[, not_zero[i, ], drop = FALSE])
    beta  <- sur_cpp(xlist, ylist, res_cov)
    for (i in 1:m) {
      coeff[t(not_zero)] <- beta
      estim <- t(coeff)
    }
  }
  if (intercept) {
    cnst  <- estim[,  1]
    estim <- estim[, -1]
  } else {
    cnst <- NULL
  }
  if (p > 0) {
    ar_array <- array(estim[, 1:(p*m)], c(m, m, p))
  } else {
    ar_array <- NULL
  }
  if (q > 0) {
    ma_array <- array(estim[, (p*m+1):(p*m+q*m)], c(m, m, q))
  } else {
    ma_array <- NULL
  }
  if (ret == "varma") {
    structure(
      list(
        intercept = cnst,
        ar = ar_array,
        ma = ma_array,
        cov = res_cov,
        estimation_method = fn_name,
        loglik = - (nobs / 2) * (log(2 * pi) + 1) -
          (nrow(resid) / 2) * determinant(res_cov, logarithm = TRUE)$modulus,
        n = nrow(Y),
        nobs = nobs,
        npar = sum(estim != 0) + if (intercept) m else 0,
        y = Y,
        residuals = resid
      ),
      class = "varma"
    )
  } else {
    regreg
  }
}

#' VARMA estimation using Maximum-Likelihood in State-Space form using KFAS
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept logical, if TRUE, the model will include an intercept.
#' @param maxit integer, maximum number of iterations for the iterated Hannan-Rissanen
#' estimator, which is used for providing starting values for numerical MLE and for
#' computing the marginal covariance matrix of the VARMA in state-space form, used
#' for initializing the Kalman filter.
#' @param ret a character string indicating the type of output: "varma" or
#' "kfas" (a KFAS' SSModel object).
#'
#' @returns If ret = "varma" a varma object, otherwise a SSModel object from the KFAS package.
#'
#' @importFrom KFAS SSMcustom SSModel fitSSM KFS
#' @export
fit_varma_kfas <- function(Y, p = 1, q = 1, intercept = TRUE,
                           maxit = 100, ret = c("varma", "kfas")) {
  fn_name <- as.character(sys.call()[[1]]) # function name
  # input controls
  ret <- match.arg(ret)
  # for initial values (see below)
  ih <- fit_varma_ihr(Y, p, q, intercept, maxit = maxit)
  ss <- varma_to_ss(ih)
  n <- nrow(Y)
  m <- ncol(Y)
  r <- max(p, q+1)
  mr <- m*r
  # matrices for the ss form
  mT <- ss$T
  mR <- ss$R
  mQ <- ss$Q
  mZ <- ss$Z
  a1 <- matrix(0, nrow(mT), 1)
  P1 <- solve_dlyap_iter(Matrix::Matrix(mT, sparse = TRUE), mR %*% tcrossprod(mQ, mR))
  P1inf <- matrix(0, mr, mr)
  mod <- SSModel(
    Y ~ 0 + SSMcustom(Z = mZ, T = mT, R = mR, Q = mQ, a1 = a1,
                      P1 = P1, P1inf = P1inf, index = 1:m),
    H = matrix(0, m, m)
  )
  # update function
  # order of pars:
  # m intercepts
  # m*m*p AR coefficients
  # m*m*q MA coefficients
  # m*(m+1)/2 elements of chol(cov_eps)
  updt <- function(pars, model) {
    if (intercept) {
      model$y <- Y - rep(pars[1:m], each = n)
      pars <- pars[-(1:m)]
    }
    if (p > 0) model$T[1:(m*p), 1:m, 1] <- pars[1:(m*m*p)]
    if (q > 0) model$R[(m+1):(m+m*q), 1:m, 1] <- pars[(m*m*p+1):(m*m*p+m*m*q)]
    Pt <- matrix(0, m, m)
    Pt[upper.tri(Pt, TRUE)] <- pars[(m*m*p+m*m*q+1):(m*m*p+m*m*q+m*(m+1)/2)]
    model$Q[,, 1] <- crossprod(Pt)
    model
  }
  # initial values from fit_varma_ihr
  if (intercept) { # intercepts
    inits <- ih$intercept
  } else {
    inits <- NULL
  }
  if (p > 0) { # VAR part
    mPhi <- matrix(0, p*m, m)
    for (i in 1:p) mPhi[((i-1)*m + 1):(i*m), 1:m] <- ih$ar[,,i] # t(ih$coefficients)[1:m, ((i-1)*m + 1):(i*m)]
    inits <- c(inits, as.numeric(mPhi))
  }
  if (q > 0) { # VMA part
    mTheta <- matrix(0, q*m, m)
    for (i in 1:q) mTheta[((i-1)*m + 1):(i*m), 1:m] <- ih$ma[,,i] # t(ih$coefficients)[1:m, ((i-1)*m + 1):(i*m)]
    inits <- c(inits, as.numeric(mTheta))
  }
  # residual covariance matrix part
  inits <- c(inits, as.numeric(chol(ih$cov)[upper.tri(ih$cov, TRUE)]))
  # estimates
  ssfit <- fitSSM(mod, inits, updt, method = "BFGS")
  ssmod <- updt(ssfit$optim.out$par, ssfit$model)
  ssflt <- KFS(ssmod, filtering = "state", smoothing = "none")
  if (ret == "kfas") return(ssfit)
  if (intercept) cnst <- ssfit$optim.out$par[1:m] else cnst <- NULL
  if (p > 0) {
    array_ar <- array(0, c(m, m, p))
    for (i in 1:p) array_ar[,,i] <- ssfit$model$T[((i-1)*m+1):(i*m), 1:m, 1]
  } else array_ar <- NULL
    if (q > 0) {
    array_ma <- array(0, c(m, m, q))
    for (i in 1:q) array_ma[,,i] <- ssfit$model$R[(i*m+1):((i+1)*m), 1:m, 1]
  } else array_ma <- NULL
  structure(
    list(
      mean = cnst,
      ar = array_ar,
      ma = array_ma,
      cov = ssfit$model$Q[,,1],
      estimation_method = fn_name,
      loglik = logLik(ssfit$model),
      n = n,
      nobs = sum(!is.na(Y)),
      y = Y,
      residuals = ssflt$v,
      state_pred_mean = ssflt$a[n+1, ],
      state_pred_cov  = ssflt$P[,, n+1],
      transition_matrix = ssmod$T[,,1],
      disturbance_cov = ssmod$R[,,1] %*% tcrossprod(ssmod$Q[,,1], ssmod$R[,,1])
    ),
    class = "varma"
  )
}


#' VARMA estimation using Maximum-Likelihoo in State-Space form (FKF.SP)
#'
#' @param Y a (n x m) matrix with the time series.
#' @param p order of the AR part.
#' @param q order of the MA part.
#' @param intercept boolean, if TRUE, the model will include an intercept.
#' @param maxit integer, maximum number of iterations for the iterated Hannan-Rissanen
#' estimator, which is used for providing starting values for numerical MLE and for
#' computing the marginal covariance matrix of the VARMA in state-space form, used
#' for initializing the Kalman filter.
#' @param ret a character string indicating the type of output: "varma" (default)
#' or "ss" (a list with the state-space form matrices and fkf function output).
#'
#' @returns A varma object or a list with the matrices of the state-space form,
#' according to the choice of the `ret` parameter.
#'
#' @export
fit_varma_fkf <- function(Y, p = 1, q = 1, intercept = TRUE,
                          maxit= 100, ret = c("varma", "ss")) {
  ret <- match.arg(ret)
  fn_name <- as.character(sys.call()[[1]]) # function name
  # for initial values (see below)
  ih <- fit_varma_ihr(Y, p, q, intercept, maxit = maxit)
  n <- nrow(Y)
  m <- ncol(Y)
  r <- max(p, q+1)
  # matrices for the ss form
  mT <- matrix(0, m*r, m*r)
  mr <- m*r
  if (q > 0) diag(mT[1:(mr - m), (m+1):mr]) <- 1
  mR <- rbind(diag(m), matrix(0, mr-m, m))
  mQ <- ih$cov
  mZ <- cbind(diag(1, m, m), matrix(0, m, mr -m))
  a1 <- numeric(mr)
  mHHt <- mR %*% tcrossprod(mQ, mR)
  mGGt <- rep(0, m)
  P1 <- solve_dlyap_iter(Matrix::Matrix(mT, sparse = TRUE), mHHt)
  Pt <- matrix(0, m, m)
  ct <- matrix(0, m, 1)
  y <- t(Y)
  obj <- function(pars) {
    if (intercept) {
      ct[] <- matrix(pars[(1:m)], m, 1)
      pars <- pars[-(1:m)]
    }
    if (p > 0) mT[1:(m*p), 1:m] <- pars[1:(m*m*p)]
    if (q > 0) mR[(m+1):(m+m*q), 1:m] <- pars[(m*m*p+1):(m*m*p+m*m*q)]
    Pt[upper.tri(Pt, TRUE)] <- pars[(m*m*p+m*m*q+1):(m*m*p+m*m*q+m*(m+1)/2)]
    mQ[] <- crossprod(Pt)
    -FKF.SP::fkf.SP(a0 = a1, P0 = P1, dt = matrix(0, dim(mT)[1], 1), ct = ct,
            Tt = mT, Zt = mZ, HHt = mR %*% tcrossprod(mQ, mR), GGt = mGGt, yt = y)
  }
  final_run <- function(pars) {
    if (intercept) {
      ct[] <- matrix(pars[(1:m)], m, 1)
      pars <- pars[-(1:m)]
    }
    if (p > 0) mT[1:(m*p), 1:m] <- pars[1:(m*m*p)]
    if (q > 0) mR[(m+1):(m+m*q), 1:m] <- pars[(m*m*p+1):(m*m*p+m*m*q)]
    Pt[upper.tri(Pt, TRUE)] <- pars[(m*m*p+m*m*q+1):(m*m*p+m*m*q+m*(m+1)/2)]
    mQ[] <- crossprod(Pt)
    HH = mR %*% tcrossprod(mQ, mR)
    c(FKF.SP::fkf.SP(a0 = a1, P0 = P1, dt = matrix(0, dim(mT)[1], 1), ct = ct,
           Tt = mT, Zt = mZ, HHt = HH, GGt = mGGt, yt = y,
           verbose = TRUE),
      list(HH = HH, R = mR, Q = mQ)
    )
  }
  # initial values from fit_varma_ihr
  if (intercept) { # intercepts
    inits <- colMeans(Y)
  } else {
    inits <- NULL
  }
  if (p > 0) { # VAR part
    mPhi <- matrix(0, p*m, m)
    for (i in 1:p) mPhi[((i-1)*m + 1):(i*m), 1:m] <- ih$ar[,,i] # t(ih$coefficients)[1:m, ((i-1)*m + 1):(i*m)]
    inits <- c(inits, as.numeric(mPhi))
  }
  if (q > 0) { # VMA part
    mTheta <- matrix(0, q*m, m)
    for (i in 1:q) mTheta[((i-1)*m + 1):(i*m), 1:m] <- ih$ma[,,i] # t(ih$coefficients)[1:m, ((i-1)*m + 1):(i*m)]
    inits <- c(inits, as.numeric(mTheta))
  }
  # residual covariance matrix part
  inits <- c(inits, as.numeric(chol(ih$cov)[upper.tri(ih$cov, TRUE)]))
  # estimates
  optout <- stats::optim(inits, obj, method = "BFGS")
  ssout  <- final_run(optout$par)
  if (ret == "ss") {
    ssout
  } else {
    if (intercept) cnst <- optout$par[1:m] else cnst <- NULL
    if (p > 0) {
      array_ar <- array(0, c(m, m, p))
      for (i in 1:p) array_ar[,,i] <- ssout$Tt[((i-1)*m+1):(i*m), 1:m, 1]
    } else array_ar <- NULL
    if (q > 0) {
      array_ma <- array(0, c(m, m, q))
      for (i in 1:q) array_ma[,,i] <- ssout$R[(i*m+1):((i+1)*m), 1:m]
    } else array_ma <- NULL
    structure(
      list(
        mean = cnst,
        ar = array_ar,
        ma = array_ma,
        cov = ssout$Q,
        estimation_method = fn_name,
        loglik = -optout$value,
        nobs = sum(!is.na(Y)),
        y = Y,
        residuals = t(ssout$v),
        state_pred_mean = ssout$at[, n+1],
        state_pred_cov  = ssout$Pt[,, n+1],
        transition_matrix = ssout$Tt[,,1],
        disturbance_cov = ssout$HH
      ),
      class = "varma"
    )
  }
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
#' @param maxit integer, maximum number of iterations for the iterated Hannan-Rissanen
#' estimator, which is used for providing starting values for numerical MLE and for
#' computing the marginal covariance matrix of the VARMA in state-space form, used
#' for initializing the Kalman filter
#' @return A varma object.
#' @export
fit_varma_cpp <- function(Y, p, q, intercept = TRUE, maxit = 100) {
  fn_name <- as.character(sys.call()[[1]]) # function name
  Yt <- t(Y)
  # for initial values use intrated Hannan-Rissanen esimates
  ih <- fit_varma_ihr(Y, p, q, intercept, maxit = maxit)
  ss <- varma_to_ss(ih)
  n <- nrow(Y)
  m <- ncol(Y)
  r <- max(p, q+1)
  # matrices for the ss form
  mT <- Matrix::Matrix(ss$T, sparse = TRUE) # sparse T
  mR <- Matrix::Matrix(ss$R, sparse = TRUE) # sparse R
  mQ <- ss$Q
  a1 <- matrix(0, nrow(mT), 1)
  P1 <- solve_dlyap_iter(mT, as.matrix(mR %*% Matrix::tcrossprod(mQ, mR)))
  P1 <- (P1 + t(P1))/2
  Pt <- matrix(0, m, m)
  # Objective function
  ## order of pars:
  ## m intercepts
  ## m*m*p AR coefficients
  ## m*m*q MA coefficients
  ## m*(m+1)/2 elements of chol(cov_eps)
  ncall <- 0
  y <- Yt
  one_to_m <- 1:m
  one_to_mp <- 1:(m*p)
  one_to_mmp <- 1:(m*m*p)
  mplus1_to_mmq <- (m+1):(m+m*q)
  mmpplus1_to_mmpplusmmq <- (m*m*p+1):(m*m*p+m*m*q)
  mmpplusmmqplus1_to_last<- (m*m*p+m*m*q+1):(m*m*p+m*m*q+m*(m+1)/2)
  obj <- function(pars) {
    if (intercept) {
      y <- Yt - pars[one_to_m]
      pars <- pars[-one_to_m]
    }
    if (p > 0) mT[one_to_mp, one_to_m] <- pars[one_to_mmp]
    if (q > 0) mR[mplus1_to_mmq, one_to_m] <- pars[mmpplus1_to_mmpplusmmq]
    Pt[upper.tri(Pt, TRUE)] <- pars[mmpplusmmqplus1_to_last]
    mQ[] <- crossprod(Pt)
    -kalmanLogLik(mT, mR, mQ, a1, P1, y)
  }
  final_run <- function(pars) {
    if (intercept) {
      y <- Yt - pars[one_to_m]
      pars <- pars[-one_to_m]
    }
    if (p > 0) mT[one_to_mp, one_to_m] <- pars[one_to_mmp]
    if (q > 0) mR[mplus1_to_mmq, one_to_m] <- pars[mmpplus1_to_mmpplusmmq]
    Pt[upper.tri(Pt, TRUE)] <- pars[mmpplusmmqplus1_to_last]
    mQ[] <- crossprod(Pt)
    kalman(mT,
           mR,
           mQ,
           mat_at,
           array_Pt,
           y,
           mat_v,
           array_F)
  }

  # initial values from fit_varma_ihr
  if (intercept) { # intercepts
    inits <- colMeans(Y)
  } else {
    inits <- NULL
  }
  if (p > 0) { # VAR part
    mPhi <- matrix(0, p*m, m)
    for (i in 1:p) mPhi[((i-1)*m + 1):(i*m), 1:m] <- ih$ar[,,i] # t(ih$coefficients)[1:m, ((i-1)*m + 1):(i*m)]
    inits <- c(inits, as.numeric(mPhi))
  }
  if (q > 0) { # VMA part
    mTheta <- matrix(0, q*m, m)
    for (i in 1:q) mTheta[((i-1)*m + 1):(i*m), 1:m] <- ih$ma[,,i] # t(ih$coefficients)[1:m, ((i-1)*m + 1):(i*m)]
    inits <- c(inits, as.numeric(mTheta))
  }
  # residual covariance matrix part
  inits <- c(inits, as.numeric(chol(ih$cov)[upper.tri(ih$cov, TRUE)]))
  # estimates
  ssfit <- stats::optim(inits, obj, method = "BFGS")
  mat_at <- matrix(a1, nrow(mT), n+1)
  array_Pt <- array(P1, c(nrow(mT), nrow(mT), n+1))
  mat_v <- matrix(0, m, n)
  array_F <- array(0, c(m, m, n))
  fin_log_lik <- final_run(ssfit$par)
  pars <- ssfit$par
  if (intercept) {
    cnst <- pars[1:m]
    pars <- pars[-(1:m)]
  } else cnst <- NULL
  if (p > 0) {
    A <- matrix(pars[1:(m*m*p)], m*p, m)
    array_ar <- array(0, c(m, m, p))
    for (i in 1:p) array_ar[,,i] <- A[((i-1)*m+1):(i*m), ]
  } else array_ar <- NULL
  if (q > 0) {
    M <- matrix(pars[(m*m*p+1):(m*m*p+m*m*q)], m*q, m)
    array_ma <- array(0, c(m, m, q))
    for (i in 1:q) array_ma[,,i] <- M[((i-1)*m+1):(i*m), ]
  } else array_ma <- NULL
  structure(
    list(
      mean = cnst,
      ar = array_ar,
      ma = array_ma,
      cov = mQ,
      estimation_method = fn_name,
      loglik = fin_log_lik,
      n = n,
      nobs = sum(!is.na(Y)),
      y = Y,
      residuals = t(mat_v),
      state_pred_mean = mat_at[, n+1],
      state_pred_cov  = array_Pt[,, n+1],
      transition_matrix = as.matrix(mT),
      disturbance_cov = as.matrix(mR %*% Matrix::tcrossprod(mQ, mR))
    ),
    class = "varma"
  )
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
  nobs <- sum(!is.na(y))

  structure(
    list(
      intercept = if (mod$fit$cnst) mod$fit$const else NULL,
      ar = if (is.null(mod$fit$Phi)) NULL else array(mod$fit$Phi,   c(m, m, ncol(mod$fit$Phi)/m)),
      ma = if (is.null(mod$fit$Theta)) NULL else array(-mod$fit$Theta, c(m, m, ncol(mod$fit$Theta)/m)),
      cov = mod$fit$Sigma,
      estimation_method = paste0("fable_varima_", mod$fit$identification),
      loglik = - (nobs / 2) * (log(2 * pi) + 1) -
        (nrow(mod$fit$residuals) / 2) * determinant(mod$fit$Sigma, logarithm = TRUE)$modulus,
      n = nrow(y),
      nobs = nobs,
      npar = length(mod$fit$coef),
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
  if (all(Mod(iroots$ar) < 1) && all(Mod(iroots$ma) < 1)) {
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
    warning("There is no log-likelihood in the varma object\n")
    return(structure(NA, df = NA))
  }
  structure(
    object$loglik,
    df = length(object$ar) + length(object$ma) +
      prod(dim(object$cov) + c(0, 1))/2
  )
}

#' nobs method for a varma object
#'
#' @param object a varma object.
#' @param ... not used.
#'
#' @returns A numeric scalar with the number of observations.
#'
#' @importFrom stats nobs
#' @method nobs varma
#' @export
nobs.varma <- function(object, ...) {
  if (is.null(object$nobs)) {
    warning("There is no informatoin about the number of observations in the varma object")
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
#' @method predict varma
#' @export
predict.varma <- function(object, n.ahead = 1, cov = TRUE, ...) {
  if (is.null(object$estimation_method)) stop("This varma object was not estimated on data")
  mpq <- dim(object)
  m <- mpq[1]
  p <- mpq[2]
  q <- mpq[3]
  # if (object$estimation_method %in% c("fit_varma_ihr", "fit_varma_net")) {
  if ((!is.null(object$y)) && (!is.null(object$residuals))) {
    yhat <- matrix(0, n.ahead, m)
    colnames(yhat) <- colnames(object$y)
    if (p > 0) { # AR part
      ny   <- nrow(object$y)
      ylag <- as.numeric(t(object$y[ny:(ny-p+1), ]))
      At   <- t(matrix(object$ar, m, m*p))
      for (i in 1:n.ahead) {
        yhat[i, ] <- ylag %*% At
        if (!is.null(object$intercept)) yhat[i, ] <- yhat[i, ] + object$intercept
        ylag <- if (p == 1) yhat[i, ] else c(yhat[i, ], ylag[1:((p-1)*m)])
      }
    }
    if (q > 0) { # MA part
      nr   <- nrow(object$residuals)
      elag <- as.numeric(t(object$residuals[nr:(nr-q+1), ]))
      Mt   <- t(matrix(object$ma, m, m*q))
      for (i in 1:min(q, n.ahead)) {
        yhat[i, ] <- yhat[i, ] + elag %*% Mt
        if (q > 1) {
          elag <- utils::head(elag, m*(q-i))
          Mt   <- utils::tail(Mt, m*(q-i))
        }
      }
      if ((p < 1) && !is.null(object$intercept)) yhat <- yhat + rep(object$intercept, each = n.ahead)
    }
    if (cov) {
      S <- array(0, c(m, m, n.ahead))
      dimnames(S) <- list(colnames(yhat), colnames(yhat), NULL)
      S[,, 1] <- object$cov
      for (i in 2:n.ahead) { # AR contribution
        if (p > 0) {
          for (j in 1:min(p, i-1)) {
            S[,, i] <- S[,, i] + object$ar[,, j] %*% tcrossprod(S[,,i-j], object$ar[,, j])
          }
        }
        S[,, i] <- S[,, i] + object$cov # new innovaton contribution
        if (q > 0) {
          for (j in 1:min(q, i-1)) {
            S[,, i] <- S[,, i] + object$ma[,, j] %*% tcrossprod(object$cov, object$ma[,, j])
          }
        }
      }
      return(list(mean = yhat, cov = S))
    } else {
      return(yhat)
    }
  }
  # if (object$estimation_method %in% c("fit_varma_fkf",
  #                                     "fit_varma_cpp",
  #                                     "fit_varma_kfas")) {
  if ((!is.null(object$state_pred_mean) &&
       !is.null(object$state_pred_cov)  &&
       !is.null(object$transition_matrix)  &&
       !is.null(object$disturbance_cov))) {
    d <- nrow(object$transition_matrix)
    ahat <- matrix(0, n.ahead, d)
    ahat[1, ] <- object$state_pred_mean
    for (i in 2:n.ahead) {
      ahat[i, ] <- as.numeric(object$transition_matrix %*% ahat[i-1, ])
    }
    yhat <- ahat[, 1:m] + rep(object$mean, each = n.ahead)
    colnames(yhat) <- colnames(object$y)
    if (cov) {
      phat <- array(0, c(d, d, n.ahead))
      phat[,, 1] <- object$state_pred_cov
      for (i in 2:n.ahead) {
        phat[,, i] <- object$transition_matrix %*%
          tcrossprod(phat[,, i-1], object$transition_matrix) +
          object$disturbance_cov
      }
      return(
        list(mean = yhat,
             cov  = phat[1:m, 1:m, ]
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
#' they may have a very large variance (default: 5).
#' @param fit_fn function to apply to every simulated path: it is usually a
#' `fit_varma_???` function (default: the function that generated the estimates
#' in the varma object).
#' @param ... other parameters to be passed to `fit_fn`.
#' @param parallel logical, if TRUE, parallel processing is used (default: FALSE).
#' @param n_cores number of cores for parallel processing. If NULL, uses
#'   parallel::detectCores() - 1 (default: NULL).
#'
#' @return A list with `nsim` slots containing the applications of the
#' function `fit_fn` to the simulated sample paths.
#' @export
bootstrap_varma <- function(varma, nsim = 100, n = varma$n, resids_to_skip = 5,
                           fit_fn = get(varma$estimation_method), ...,
                           parallel = FALSE, n_cores = NULL) {

  res <- t(varma$residuals[-(1:resids_to_skip), ])
  mpq <- dim(varma)
  intercept <- !is.null(varma$intercept) || !is.null(varma$mean)

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

  # --- Start of Parallel/Sequential Logic ---

  run_simulation <- function() {
    y <- sim_varma_rcpp(
      ar_array,
      ma_array,
      res[, sample.int(ncol(res), n, replace = TRUE)] # Fixed: should be ncol
    )
    fit_fn(Y = t(y), p = mpq[2], q = mpq[3], intercept = intercept, ...)
  }

  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is needed for parallel execution. Please install it.")
    }

    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }

    # Set the parallel plan conditionally based on OS
    if (.Platform$OS.type == "windows") {
      future::plan(future::multisession, workers = n_cores)  # Fallback for Windows
    } else {
      future::plan(future::multicore, workers = n_cores)     # Preferred for Linux/macOS
    }

    # Ensure the user's original plan is restored when the function exits
    on.exit(future::plan(future::sequential), add = TRUE)

    message(paste("Running", nsim, "simulations in parallel on", n_cores, "cores..."))

    results <- future.apply::future_replicate(
      nsim,
      expr = run_simulation(),
      simplify = FALSE,
      future.seed = TRUE # Ensures reproducibility in parallel
    )

  } else {
    message(paste("Running", nsim, "simulations sequentially..."))

    results <- replicate(
      nsim,
      expr = run_simulation(),
      simplify = FALSE
    )
  }

  structure(
    results,
    generator = varma
  )
}
# boostrap_varma <- function(varma, nsim = 100, n = varma$n, resids_to_skip = 5,
#                            fit_fn = get(varma$estimation_method), ...) {
#
#   res <- t(varma$residuals[-(1:resids_to_skip), ])
#   mpq <- dim(varma)
#   if (!is.null(varma$intercept) || !is.null(varma$mean)) {
#     intercept = TRUE
#   } else {
#     intercept = FALSE
#   }
#   ar_array <- if (!is.null(varma$ar)) {
#     array(varma$ar, c(mpq[1], mpq[1], mpq[2]))
#   } else {
#     array(0, c(mpq[1], mpq[1], 0))
#   }
#   ma_array <- if (!is.null(varma$ma)) {
#     array(varma$ma, c(mpq[1], mpq[1], mpq[3]))
#   } else {
#     array(0, c(mpq[1], mpq[1], 0))
#   }
#
#   structure(
#     replicate(nsim,
#               {
#                 y <- sim_varma_rcpp(ar_array,
#                                     ma_array,
#                                     res[, sample.int(nrow(res), n, replace = TRUE)])
#                 fit_fn(Y = t(y),
#                        p = mpq[2], q = mpq[3], intercept = intercept, ...)
#               },
#               simplify = FALSE
#     ),
#     generator = varma
#   )
# }

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
