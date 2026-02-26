#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
arma::cube compute_varma_acgf_cpp(Rcpp::Nullable<arma::cube> ar_in,
                                  Rcpp::Nullable<arma::cube> ma_in,
                                  arma::mat cov,
                                  int max_lag) {

  // 1. Handle inputs safely to allow pure AR or pure MA models
  int p = 0, q = 0;
  arma::cube ar, ma;
  if (ar_in.isNotNull()) { ar = Rcpp::as<arma::cube>(ar_in); p = ar.n_slices; }
  if (ma_in.isNotNull()) { ma = Rcpp::as<arma::cube>(ma_in); q = ma.n_slices; }

  int d = cov.n_rows;
  int m = std::max(p, q + 1);
  int num_lags = m + max_lag;

  // 2. Compute Wold coefficients (Psi)
  // C++ is 0-indexed: psi.slice(k) precisely maps to Psi_k
  arma::cube psi(d, d, num_lags + 1, arma::fill::zeros);
  psi.slice(0).eye();

  for (int k = 1; k <= num_lags; k++) {
    arma::mat psi_k(d, d, arma::fill::zeros);
    if (k <= q) {
      psi_k += ma.slice(k - 1);
    }
    if (p > 0) {
      for (int i = 1; i <= std::min(k, p); i++) {
        psi_k += ar.slice(i - 1) * psi.slice(k - i);
      }
    }
    psi.slice(k) = psi_k;
  }

  // 3. Pre-compute elements for the persistent component
  arma::cx_mat Y;
  arma::cx_mat Z_V;
  arma::cx_vec lambda;

  if (p > 0) {
    // Construct Companion Matrix F
    arma::mat F_mat(d * p, d * p, arma::fill::zeros);
    for (int i = 0; i < p; i++) {
      F_mat.submat(0, i * d, d - 1, (i + 1) * d - 1) = ar.slice(i);
    }
    if (p > 1) {
      F_mat.submat(d, 0, d * p - 1, d * (p - 1) - 1).eye();
    }

    // Construct Matrix G (0-indexed logic maps correctly to psi.slice)
    arma::mat G(d * p, d, arma::fill::zeros);
    for (int i = 1; i <= p; i++) {
      int idx = m - i + 1;
      G.submat((i - 1) * d, 0, i * d - 1, d - 1) = psi.slice(idx);
    }

    // Eigendecomposition
    arma::cx_vec eigval;
    arma::cx_mat V;
    arma::eig_gen(eigval, V, F_mat);
    lambda = eigval;

    // Cast to complex for persistent equations
    arma::cx_mat V_inv = arma::inv(V);
    arma::cx_mat cG = arma::conv_to<arma::cx_mat>::from(G);
    arma::cx_mat cCov = arma::conv_to<arma::cx_mat>::from(cov);

    // Compute M using .st() for standard (non-conjugate) transpose
    arma::cx_mat M = V_inv * cG * cCov * cG.st() * V_inv.st();

    // Compute denominator: 1 - outer(lambda, lambda)
    arma::cx_mat denom(d * p, d * p);
    for(int r = 0; r < d * p; r++) {
      for(int c = 0; c < d * p; c++) {
        denom(r, c) = 1.0 - lambda(r) * lambda(c);
      }
    }
    Y = M / denom;
    Z_V = V.submat(0, 0, d - 1, d * p - 1);
  }

  // 4. Compute the Autocovariance sequence
  arma::cube Gamma(d, d, max_lag + 1, arma::fill::zeros);

  for (int h = 0; h <= max_lag; h++) {
    // Finite memory sum (real numbers)
    arma::mat finite_part(d, d, arma::fill::zeros);
    for (int j = 0; j <= m - 1; j++) {
      finite_part += psi.slice(j + h) * cov * psi.slice(j).t();
    }

    if (p > 0) {
      // Persistent part logic
      arma::cx_vec lam_h = arma::pow(lambda, h);
      arma::cx_mat Y_h = Y;

      // Multiply each row by lambda^h (matches R's vector-matrix recycling)
      Y_h.each_col() %= lam_h;

      arma::cx_mat persistent_part = Z_V * Y_h * Z_V.st();

      // Combine and drop imaginary artifacts
      Gamma.slice(h) = finite_part + arma::real(persistent_part);
    } else {
      Gamma.slice(h) = finite_part;
    }
  }

  return Gamma;
}
