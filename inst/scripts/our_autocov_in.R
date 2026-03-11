# R version of our method in R

autocov_sbrana <- function(ar, ma, cov, max_lag) {
  # Infer dimension and orders
  d <- nrow(cov)
  p <- if (is.null(ar)) 0 else dim(ar)[3]
  q <- if (is.null(ma)) 0 else dim(ma)[3]

  # m = max(p, q+1) as per Lemma 2
  m <- max(p, q + 1)
  num_lags <- m + max_lag

  # 1. Recursively compute Wold coefficients (Psi) up to num_lags
  # Indexing note: psi[,,k] corresponds to Psi_{k-1}
  psi <- array(0, dim = c(d, d, num_lags + 1))
  psi[,,1] <- diag(d) # Psi_0 = I_d

  for (k in 1:num_lags) {
    psi_k <- matrix(0, d, d)
    if (k <= q) psi_k <- psi_k + ma[,,k]

    if (p > 0) {
      for (i in 1:min(k, p)) {
        psi_k <- psi_k + ar[,,i] %*% psi[,,k - i + 1]
      }
    }
    psi[,,k + 1] <- psi_k
  }

  # 2. Pre-compute elements for the persistent component (if p > 0)
  if (p > 0) {
    # Construct Companion Matrix F
    F_mat <- matrix(0, d * p, d * p)
    F_mat[1:d, ] <- matrix(ar, nrow = d)
    if (p > 1) diag(F_mat[-(1:d), 1:((p-1)*d)]) <- 1

    # Construct Matrix G
    G <- matrix(aperm(psi[,,(m+1):(m-p+2)], c(1, 3, 2)), d * p, d)

    # Eigendecomposition of F
    eig <- eigen(F_mat)
    V <- eig$vectors
    lambda <- eig$values
    V_inv <- solve(V)

    # Compute M using standard algebraic transpose t()
    V_inv_G <- V_inv %*% G
    M <- tcrossprod(V_inv_G %*% cov, V_inv_G)

    # Compute Y = M / (1 - lambda * lambda') using vectorized outer product
    Y <- M / (1 - outer(lambda, lambda))

    # Extract the Z %*% V equivalent
    Z_V <- V[1:d, , drop = FALSE]
  }

  # 3. Compute the Autocovariance sequence Gamma_h
  Gamma <- array(0, dim = c(d, d, max_lag + 1))

  for (h in 0:max_lag) {
    # Calculate finite-memory sum
    finite_part <- matrix(0, d, d)
    for (j in 0:(m - 1)) {
      finite_part <- finite_part + # psi[,,j + h + 1] %*% cov %*% t(psi[,,j + 1])
        tcrossprod(psi[,,j + h + 1] %*% cov, psi[,,j + 1])
    }

    if (p > 0) {
      # Calculate persistent-memory matrix
      Y_h <- Y * (lambda^h)
      persistent_part <- tcrossprod(Z_V %*% Y_h, Z_V) # Z_V %*% Y_h %*% t(Z_V)

      # Combine and discard any numerical-artifact imaginary parts
      Gamma[,, h + 1] <- Re(finite_part + persistent_part)
    } else {
      # Pure VMA process
      Gamma[,, h + 1] <- finite_part
    }
  }

  return(Gamma)
}
