#' Solve a discrete stationary Lyapunov equation
#'
#' Solves \eqn{P = T P T' + C} without constructing a Kronecker system.
#' @param T finite square transition matrix with eigenvalues inside the unit circle.
#' @param C finite symmetric right-hand-side matrix of the same size as T.
#' @param method "doubling", "schur", or "auto". Doubling sums the covariance
#' series geometrically. Schur uses a bilinear transformation and a Sylvester
#' solver. Auto uses Schur for dimensions at least 16 and otherwise doubling,
#' falling back to doubling if the Schur residual check fails.
#' @param tol positive numerical tolerance (default 1e-10).
#' @param max_iter maximum number of doubling iterations (default 200).
#' @return A symmetric matrix P. With a positive semidefinite C this is the
#' stationary state covariance. This differs from the steady-state conditional
#' covariance of a Kalman filter, which solves a Riccati equation.
#' @export
#' @examples
#' solve_discrete_lyapunov(matrix(0.9), matrix(1))
solve_discrete_lyapunov <- function(T,C,method=c("auto","doubling","schur"),
                                   tol=1e-10,max_iter=200) {
  method <- match.arg(method)
  max_iter <- .varma_integer(max_iter,"max_iter")
  solve_dlyap_dense_cpp(as.matrix(T),as.matrix(C),max_iter,tol,method)
}

#' VARMA autocovariances using the stationary state covariance
#'
#' Computes \eqn{\Gamma_h=\operatorname{Cov}(y_{t+h},y_t)} by solving a
#' Lyapunov equation and propagating only the needed covariance columns.
#' Unlike the vectorised Kronecker methods, the solver does not build a system
#' with squared state dimension or retain all state covariance matrices by lag.
#' @param varma a stable varma model with innovation covariance.
#' @param maxlag largest non-negative autocovariance lag.
#' @param method Lyapunov method passed to solve_discrete_lyapunov.
#' @param tol Lyapunov convergence tolerance.
#' @return An m by m by (maxlag+1) array, including the lag-zero covariance.
#' @export
autocov_ss <- function(varma,maxlag=10,method=c("auto","doubling","schur"),tol=1e-10) {
  maxlag <- .varma_integer(maxlag,"maxlag",0L)
  method <- match.arg(method)
  mpq <- dim(varma); m <- as.integer(mpq[1L])
  if (is.null(varma$cov) || !is.matrix(varma$cov) ||
      !identical(dim(varma$cov),c(m,m)) || any(!is.finite(varma$cov)))
    stop("A finite m by m innovation covariance matrix is required")
  if (sum(mpq[-1L])==0L) {
    ans <- array(0,c(m,m,maxlag+1L)); ans[,,1L] <- varma$cov
    return(ans)
  }
  ss <- varma_to_ss(varma)
  P <- solve_discrete_lyapunov(ss$T,ss$R %*% ss$Q %*% t(ss$R),method,tol)
  autocov_ss_cpp(ss$T,P,m,maxlag)
}
