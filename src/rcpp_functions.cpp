#include <RcppArmadillo/Lighter>
using namespace Rcpp;
using namespace arma;
const double M_LOG2PI = log(2 * arma::datum::pi);

// [[Rcpp::depends(RcppArmadillo)]]

//' VARMA residuals
//'
//' Iterates the VARMA(p, q) recursion forward to produce one-step-ahead
//' residuals.  The first \code{max(p, q)} rows of the output are zero
//' (pre-sample); valid residuals begin at row \code{max(p, q) + 1} (1-indexed).
//' The AR part conditions on the observed \code{Y[1:p, ]}; pre-sample errors
//' used by the MA part are initialised to zero.
//'
//' @param Y         \eqn{n \times m} matrix of observations (rows = time points).
//' @param Phi       \eqn{m \times m \times p} cube of VAR coefficient matrices.
//'                  Pass \code{array(0, c(m, m, 0))} when \code{p = 0}.
//' @param Theta     \eqn{m \times m \times q} cube of VMA coefficient matrices.
//'                  Pass \code{array(0, c(m, m, 0))} when \code{q = 0}.
//' @param intercept \eqn{m}-vector of intercepts. Pass \code{numeric(0)} to
//'                  omit the intercept.
//' @return \eqn{n \times m} matrix of residuals.
// [[Rcpp::export]]
arma::mat varma_residuals(const arma::mat&  Y,
                           const arma::cube& Phi,
                           const arma::cube& Theta,
                           const arma::vec&  intercept) {
   const int n         = static_cast<int>(Y.n_rows);
   const int p         = static_cast<int>(Phi.n_slices);
   const int q         = static_cast<int>(Theta.n_slices);
   const int t_start   = std::max(p, q);
   const bool has_int  = intercept.n_elem > 0;

   arma::mat E(n, Y.n_cols, arma::fill::zeros);

   for (int t = t_start; t < n; ++t) {
     arma::rowvec e_t = Y.row(t);
     if (has_int)
       e_t -= intercept.t();
     for (int i = 0; i < p; ++i)
       e_t -= Y.row(t - i - 1) * Phi.slice(i).t();
     for (int j = 0; j < q; ++j)
       e_t -= E.row(t - j - 1) * Theta.slice(j).t();
     E.row(t) = e_t;
   }

   return E;
}


//' IRF of a VARMA model
//'
//' It computes an array of matrices containing the IRF
//' of a VARMA process
//'
//' @param A_cube 3d-array with AR matrices.
//' @param M_cube 3d-array with MA matrices.
//' @param P orthogonalization matrix.
//' @param horizon positive integer with the number of lags.
//'
//' @returns A 3d-array of h+1 matrices with the IRFs
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube irf_varma_rcpp(const arma::cube& A_cube,
                          const arma::cube& M_cube,
                          const arma::mat& P,
                          int horizon) {

  // Determine dimensions and orders
  int m = P.n_rows;
  int p = A_cube.n_slices;
  int q = M_cube.n_slices;

  // Create 3D cubes (arrays) for Psi and IRF
  arma::cube Psi_cube(m, m, horizon + 1, fill::zeros);
  arma::cube IRF_cube(m, m, horizon + 1, fill::zeros);

  // Psi_0 is identity matrix
  Psi_cube.slice(0) = eye(m, m);

  // Calculate Psi matrices recursively and store in the cube
  for (int j = 1; j <= horizon; j++) {
    // Add MA term if j <= q
    if (j <= q) {
      Psi_cube.slice(j) += M_cube.slice(j-1);
    }
    // Add AR terms
    for (int i = 1; i <= std::min(j, p); i++) {
      Psi_cube.slice(j) += A_cube.slice(i-1) * Psi_cube.slice(j-i);
    }
    // Calculate the IRF for this horizon
    IRF_cube.slice(j) = Psi_cube.slice(j) * P;
  }

  // Calculate IRF for horizon 0
  IRF_cube.slice(0) = Psi_cube.slice(0) * P;

  // return result
  return IRF_cube;
}

//' Simulation of a random VARMA path
//'
//' It simulates a sample path from a given VARMA process
//'
//' @param A 3d-array with AR matrices.
//' @param M 3d-array with MA matrices.
//' @param eps \eqn{m \times n} matrix with white noise variates
//' in the rows and time in the columns.
//'
//' @returns A \eqn{m \times n} matrix with the simulated sample paths.
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sim_varma_rcpp(const arma::cube& A,
                         const arma::cube& M,
                         const arma::mat& eps) {
  uword p = A.n_slices;
  uword q = M.n_slices;
  uword n = eps.n_cols;
  uword m = eps.n_rows;

  if (p > 0) {
    if (A.n_rows != m || A.n_cols != m) {
      throw std::invalid_argument("Dimension mismatch in A");
    }
  }
  if (q > 0) {
    if (M.n_rows != m || M.n_cols != m) {
      throw std::invalid_argument("Dimension mismatch in M");
    }
  }

  arma::mat y = eps;

  for (uword t = 0; t < n; t++) {

    for (uword i = 1; i <= q; i++) {
      if (t >= i) {
        y.col(t) += M.slice(i - 1) * eps.col(t - i);
      }
    }

    for (uword i = 1; i <= p; i++) {
      if (t >= i) {
        y.col(t) += A.slice(i - 1) * y.col(t - i);
      }
    }
  }

  return y;
}




//' Solve Discrete-time Lyapunov Equation via Iteration
//'
//' Solves the equation P = T * P * T' + C for P, where P is the stationary
//' covariance of the state vector in a state-space model.
//' This is used to initialize the Kalman filter for a stationary process.
//'
//' @param T state transition matrix (sparse).
//' @param C covariance matrix of the state innovations (dense).
//' @param max_iter positive integer maximum number of iterations.
//' @param tol small positive number with the tolerance for convergence.
//' @return The solution matrix P
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat solve_dlyap_iter(const arma::sp_mat& T, const arma::mat& C, int max_iter = 200, double tol = 1e-10) {
   arma::mat P = C;
   arma::mat P_old;
   for (int i = 0; i < max_iter; ++i) {
     P_old = P;
     P = T * P_old * T.t() + C;
     if (arma::approx_equal(P, P_old, "absdiff", tol)) {
       return P;
     }
   }
   // Return the last iteration if not converged, with a warning.
   Rcpp::warning("Lyapunov equation solver did not converge within %d iterations.", max_iter);
   return P;
}

//' Compute the log-likelihood of a VARMA model in state-space form
//'
//' It is efficient for large systems since it exploit the sparsity
//' of the transition matrix (\eqn(T)) and of the matrix that premultiplies the
//' disturbances (\eqn{R}). The observation matrix (\eqn(Z))and noise variance matrix
//' are not used. Call \eqn{m} the number of time series to model,
//' the state space form for this problem is:
//' \eqn{\alpha_{t+1} = T\alpha_{t} + R\eta_t}
//' \eqn{y_t = \alpha_{t}^{(1:m)}},
//' where \eqn{\alpha_{t}^{(1:m)}} denotes the subvector of \eqn{\alpha_t}
//' with the first \eqn{m} elements.
//'
//' @param T sparse matrix (use the function Matrix(T, sparse = TRUE) of the
//' Matrix package to turn the dense matrix T into a sparse one).
//' @param R sparse matrix (see above).
//' @param Q dense covariance matrix.
//' @param a1 matrix with one column with initial mean values for the state vector.
//' @param P1 initial covariance matrix of the state vector.
//' @param Yt trasposed matrix of data \eqn{n\times m}, where \eqn{n} is the
//' number of time points.
//' @param update_state: boolean (by default is false), if true the parameters
//' a1 and P1 are overwritten with the values of the state prediction for time
//' \eqn{t = n+1}: \eqn{a_{n+1|n}} and \eqn{P_{n+1|n}}.
//'
//' @returns A scalar number with the value of the log-likelihood.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double kalmanLogLik(const arma::sp_mat& T,
                    const arma::sp_mat& R,
                    const arma::mat& Q,
                    arma::mat& a1,
                    arma::mat& P1,
                    const arma::mat& Yt,
                    const bool update_state = false) { // to write a_{n+1}, P_{n+1} in a1 P1
  const arma::uword n = Yt.n_cols; // n. of observations
  const arma::uword m = Yt.n_rows; // n. of observed variables
  const arma::uword p = T.n_rows;  // n. of state variables
  double loglik = 0;

  // --- NEW LOGIC FOR CONDITIONAL UPDATE ---

  // 1. Create local matrices to serve as temporary copies if needed.
  arma::mat a1_local;
  arma::mat P1_local;

  // 2. Use references to select the target matrices.
  // If update_state is true, 'at' will refer to 'a1' (the original R matrix).
  // If update_state is false, 'at' will refer to the local copy 'a1_local'.
  arma::mat& at = update_state ? a1 : a1_local;
  arma::mat& Pt = update_state ? P1 : P1_local;

  // 3. If we are NOT updating in place, initialize the local copies.
  if (!update_state) {
    at = a1; // This copies a1's data into a1_local
    Pt = P1; // This copies P1's data into P1_local
  }

  // --- From here, the rest of the code works on 'at' and 'Pt' ---
  // It transparently modifies either the originals or the copies.

  arma::mat v(m, 1);
  arma::mat F(m, m);
  arma::mat iF(m, m);
  arma::mat PZt(p, m);
  arma::mat TPZt(p, m);
  arma::mat K(p, m);
  const arma::sp_mat Tt = T.t();
  const arma::mat RQRt = R * Q * R.t();
  arma::colvec y0;
  arma::sp_mat Z0;
  arma::uvec notna;
  bool pass = FALSE;

  for (arma::uword t=0; t < n; ++t) {
    notna = find_finite(Yt.col(t));
    if (notna.n_elem == m) { // no missing obs in y
      v = Yt.col(t) - at.rows(0, m-1);
      PZt = Pt.cols(0, m-1);
      F = Pt(span(0, m-1), span(0, m-1));
      F = (F + F.t())/2;
      pass = inv_sympd(iF, F);
      if (!pass) return -datum::inf;
      loglik += (log_det_sympd(F) + v.t()*iF*v).eval()(0,0) + m*M_LOG2PI;
      TPZt = T*PZt;
      K = TPZt*iF;
      at = T*at + K*v;
      Pt = T*Pt*Tt - K*F*K.t() + RQRt;
      continue;
    }
    if (notna.is_empty()) { // all obs in y missing
      at = T*at;
      Pt = T*Pt*Tt + RQRt;
      continue;
    }
    // y partially missing
    y0 = Yt.col(t).eval()(notna);
    v = y0 - at.rows(notna);
    PZt = Pt.cols(notna);
    F = Pt(notna, notna);
    pass = inv_sympd(iF, F);
    if (!pass) return -datum::inf;
    loglik += (real(log_det(F)) + v.t()*iF*v).eval()(0,0) + notna.n_elem*M_LOG2PI;
    TPZt = T*PZt;
    K = TPZt*iF;
    at = T*at + K*v;
    Pt = T*Pt*Tt - K*F*K.t() + RQRt;
  }

  // The 'if (update_state)' block at the end is no longer needed,
  // as this logic is now handled at the beginning.

  return -0.5*loglik;
}


//' Compute the log-likelihood of a VARMA model in state-space form
//'
//' It is efficient for large systems since it exploit the sparsity
//' of the transition matrix (\eqn(T)) and of the matrix that premultiplies the
//' disturbances (\eqn{R}). The observation matrix (\eqn(Z))and noise variance matrix
//' are not used. Call \eqn{m} the number of time series to model,
//' the state space form for this problem is:
//' \eqn{\alpha_{t+1} = T\alpha_{t} + R\eta_t}
//' \eqn{y_t = \alpha_{t}^{(1:m)}},
//' where \eqn{\alpha_{t}^{(1:m)}} denotes the subvector of \eqn{\alpha_t}
//' with the first \eqn{m} elements.
//'
//' @param T sparse matrix (use the function Matrix(T, sparse = TRUE) of the
//' Matrix package to turn the dense matrix T into a sparse one).
//' @param R sparse matrix (see above).
//' @param Q dense covariance matrix.
//' @param a1 matrix with one column with initial mean values for the state vector.
//' @param P1 initial covariance matrix of the state vector.
//' @param Yt trasposed matrix of data \eqn{n\times m}, where \eqn{n} is the
//' number of time points.
//' @param update_state: boolean (by default is false), if true the parameters
//' a1 and P1 are overwritten with the values of the state prediction for time
//' \eqn{t = n+1}: \eqn{a_{n+1|n}} and \eqn{P_{n+1|n}}.
//'
//' @returns A scalar number with the value of the log-likelihood.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double kalman(const arma::sp_mat& T,
              const arma::sp_mat& R,
              const arma::mat& Q,
                    arma::mat& at,
                    arma::Cube<double>& Pt,
              const arma::mat& Yt,
                    arma::mat& v,
                    arma::Cube<double>& F,
              const bool update_state = false) { // to write a_{n+1}, P_{n+1} in a1 P1
   const arma::uword n = Yt.n_cols; // n. of observations
   const arma::uword m = Yt.n_rows; // n. of observed variables
   const arma::uword p = T.n_rows;  // n. of state variables
   double loglik = 0;

   arma::mat iF(m, m);
   arma::mat PZt(p, m);
   arma::mat TPZt(p, m);
   arma::mat K(p, m);
   const arma::sp_mat Tt = T.t();
   const arma::mat RQRt = R * Q * R.t();
   arma::colvec y0;
   arma::sp_mat Z0;
   arma::uvec notna;
   bool pass = FALSE;

   for (arma::uword t=0; t < n; ++t) {
     notna = find_finite(Yt.col(t));
     if (notna.n_elem == m) { // no missing obs in y
       v.col(t) = Yt.col(t) - at(span(0, m-1), t); // at.rows(0, m-1);
       PZt = Pt.slice(t).cols(0, m-1); // Pt.cols(0, m-1);
       F.slice(t) = Pt.slice(t)(span(0, m-1), span(0, m-1)); // Pt(span(0, m-1), span(0, m-1));
       F.slice(t) = (F.slice(t) + F.slice(t).t())/2;
       pass = inv_sympd(iF, F.slice(t));
       if (!pass) return -datum::inf;
       loglik += (log_det_sympd(F.slice(t)) + v.col(t).t()*iF*v.col(t)).eval()(0,0) + m*M_LOG2PI;
       TPZt = T*PZt;
       K = TPZt*iF;
       at.col(t+1) = T*at.col(t) + K*v.col(t);
       Pt.slice(t+1) = T*Pt.slice(t)*Tt - K*F.slice(t)*K.t() + RQRt;
       continue;
     }
     if (notna.is_empty()) { // all obs in y missing
       at.col(t+1) = T*at.col(t);
       Pt.slice(t+1) = T*Pt.slice(t)*Tt + RQRt;
       continue;
     }
     // y partially missing
     // .col() returns a subview with no .elem(); materialise first
     const arma::vec Yt_t(Yt.col(t));
     const arma::vec at_t(at.col(t));
     arma::vec v_t = Yt_t.elem(notna) - at_t.elem(notna);            // fix 6: proper residual
     for (arma::uword k = 0; k < notna.n_elem; ++k)                  // fix 6: store residuals
       v(notna(k), t) = v_t(k);
     PZt = Pt.slice(t).cols(notna);
     arma::mat F_obs = Pt.slice(t)(notna, notna);
     F.slice(t).zeros();
     F.slice(t).submat(0, 0, notna.n_elem-1, notna.n_elem-1) = F_obs;
     arma::mat iF_obs(notna.n_elem, notna.n_elem);
     pass = inv_sympd(iF_obs, F_obs);
     if (!pass) return -datum::inf;
     loglik += (real(log_det(F_obs)) + (v_t.t()*iF_obs*v_t).eval()(0,0)) + notna.n_elem*M_LOG2PI; // fix 5
     TPZt = T*PZt;
     K = TPZt*iF_obs;
     at.col(t+1) = T*at.col(t) + K*v_t;
     Pt.slice(t+1) = T*Pt.slice(t)*Tt - K*F_obs*K.t() + RQRt;
   }

   // The 'if (update_state)' block at the end is no longer needed,
   // as this logic is now handled at the beginning.

   return -0.5*loglik;
 }


//' Seemingly unrelated regressions
//'
//' It estimates the coefficients of the SUR based on
//' generalized least squares using a given covariance matrix
//' of regression errors.
//'
//' @param X_list list of X-matrices with dependent variables
//' @param y_list list of y-vectors with regressors
//' @param sigma covariance matrix of regression errors
//'
//' @returns A vector of regression coefficients.
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sur_cpp(const Rcpp::List& X_list,
                  const Rcpp::List& y_list,
                  const arma::mat& sigma) {

  int M = X_list.size();

  // Convert R lists to a vector of Armadillo matrices/vectors
  std::vector<arma::mat> Xs(M);
  std::vector<arma::vec> ys(M);
  std::vector<int> ks(M);
  int K_total = 0;

  for(int i = 0; i < M; ++i) {
    Xs[i] = Rcpp::as<arma::mat>(X_list[i]);
    ys[i] = Rcpp::as<arma::vec>(y_list[i]);
    ks[i] = Xs[i].n_cols;
    K_total += ks[i];
  }

  // 1. Invert the error covariance matrix
  // inv_sympd is faster and more stable for symmetric positive-definite matrices
  arma::mat sigma_inv = arma::inv_sympd(sigma);

  // 2. Initialize the final "mega" matrices
  arma::mat LHS(K_total, K_total, arma::fill::zeros);
  arma::vec RHS(K_total, arma::fill::zeros);

  int row_start_idx = 0;
  for (int i = 0; i < M; ++i) {
    int col_start_idx = 0;
    for (int j = 0; j < M; ++j) {

      // 3. Get the scalar weight from the inverted covariance matrix
      double sigma_ij_inv = sigma_inv(i, j);

      // 4. Calculate the weighted cross-product X_i' * X_j
      arma::mat weighted_XX = sigma_ij_inv * (Xs[i].t() * Xs[j]);

      // 5. Place this block into the correct position in the LHS matrix
      LHS.submat(row_start_idx, col_start_idx,
                 row_start_idx + ks[i] - 1, col_start_idx + ks[j] - 1) = weighted_XX;

      col_start_idx += ks[j];
    }

    // 6. Calculate the i-th block of the RHS vector
    // This is Sum_j(sigma_ij_inv * X_i' * y_j)
    arma::vec RHS_i(ks[i], arma::fill::zeros);
    for (int j = 0; j < M; ++j) {
      RHS_i += sigma_inv(i, j) * (Xs[i].t() * ys[j]);
    }

    RHS.subvec(row_start_idx, row_start_idx + ks[i] - 1) = RHS_i;

    row_start_idx += ks[i];
  }

  // 7. Solve the system for the coefficients
  arma::vec beta_hat = arma::solve(LHS, RHS);

  return beta_hat;
}

//' VAR filter
//'
//' It computes the VAR filter on a given \eqn{n \times n} matrix containing $n$ time
//' points of a multivariate $m$-dimensional time series:
//' \eqn{y_t = x_t + A_1 y_{t-1} + \ldots A_p y_{t-p}},
//' where \eqn{x_t} is a column vector extracted from the \eqn{t}-th row of X and transposed.
//' The recursion is started using the vectors in the \eqn{p \times m} matrix.
//'
//' @param X \eqn{n \times m} matrix with the multivariate time series to filter.
//' @param A \eqn{p}-dimensional array with the \eqn{m\times m} matrices of VAR coefficients.
//' @param Y0 \eqn{p \times m} matrix with the $p$ starting observations to initialise the recursion.
//'
//' @returns A \eqn{n \times m} matrix with the filtered multivariate time series.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat var_filter(const arma::mat& X, const arma::cube& A, Rcpp::Nullable<arma::mat> Y0 = R_NilValue) {
  int n = X.n_rows;
  int m = X.n_cols;
  int p = A.n_slices;

  if (p == 0) return X;   // no AR part: output equals input

  arma::mat Y_init;
  if (Y0.isNotNull()) {
    Y_init = Rcpp::as<arma::mat>(Y0);
    if (Y_init.n_rows != static_cast<arma::uword>(p) || Y_init.n_cols != static_cast<arma::uword>(m)) {
      Rcpp::stop("Y0 must be a p x m matrix");
    }
  } else {
    Y_init = arma::zeros(p, m);
  }

  arma::mat Y(n, m);
  Y.rows(0, p - 1) = Y_init;

  for (int t = p; t < n; ++t) {
    arma::rowvec yt = X.row(t);
    for (int k = 0; k < p; ++k) {
      yt += Y.row(t - 1 - k) * A.slice(k).t();
    }
    Y.row(t) = yt;
  }

  return Y;
}

// ----------------------------------------------------------------------------
// Helper functions for Dufour & Pelletier (2022) VARMA estimation ------------
// ----------------------------------------------------------------------------

//' Helper function for building the row of Z relative to the VAR part
//'
//' This function is common to the diagonal and final forms
//'
//' @param Y data matrix
//' @param t time index
//' @param p order of VAR part
//' @param K number of time series
//'
//' @returns A row vector
//'
// [[Rcpp::depends(RcppArmadillo)]]
arma::rowvec get_var_regressors(const arma::mat& Y, int t, int p, int K) {
  arma::rowvec reg(K * p, fill::zeros);
  int count = 0;
  for (int lag = 1; lag <= p; lag++) {
    for (int k = 0; k < K; k++) {
      reg(count) = Y(t - lag, k);
      count++;
    }
  }
  return reg;
}

//' Helper function for building the row of Z relative to the VMA part
//'
//' @param U matrix of VARMA innovations
//' @param t time index
//' @param q order of MA part
//' @param k_eq MA equation to work on
//'
//' @returns A row vector
 // [[Rcpp::depends(RcppArmadillo)]]
 arma::rowvec get_ma_regressors(const arma::mat& U, int t, int q, int k_eq) {
  rowvec reg(q);
  for (int lag = 1; lag <= q; lag++) {
    reg(lag - 1) = U(t - lag, k_eq); // it takes only the equation k_eq equation
  }
  return reg;
}

//' GLS estimation of VARMA using the diagonal MA or final MA form
//'
//' @param Y data matrix
//' @param U innovation matrix
//' @param SigmaU_inv inverse of the innovation covariance matrix
//' @param p order of VAR part
//' @param q order of VMA part
//' @param type a string with one of the two choices "diagonal", "final"
// [[Rcpp::export]]
List estimate_varma_gls_cpp(arma::mat Y, arma::mat U, arma::mat SigmaU_inv,
                            int p, int q, String type) {
  int T = Y.n_rows;
  int K = Y.n_cols;
  int t_start = std::max(p, q);

  // phi and theta parameters dimensions
  int n_params_phi = K * K * p; // VAR parameters
  int n_params_theta = 0;       // VMA parameters if no VMA present

  if (type == "diagonal") {
    n_params_theta = K * q; // q parameters per equation (tot -> K*q)
  } else if (type == "final") {
    n_params_theta = q;     // q parameters common to all equations
  }

  int n_total_params = n_params_phi + n_params_theta;

  // Matrices for the linear system A * gamma = b
  arma::mat A(n_total_params, n_total_params, fill::zeros);
  arma::vec b(n_total_params, fill::zeros);

  // for loop on the time t
  for (int t = t_start; t < T; t++) {

    // Making the matrix Z_t (dimensions K x n_total_params)
    arma::mat Z_t(K, n_total_params, fill::zeros);

    // Vector of VAR regressors (equal for all equations but in different positions)
    arma::rowvec var_regs = get_var_regressors(Y, t, p, K); // Dim K*p

    for (int k = 0; k < K; k++) {
      // 1. Insert VAR regressors in the right block
      // The Phi parameters are ordere as [Phi_row1, Phi_row2, ...]
      int start_col_phi = k * (K * p);
      Z_t(k, span(start_col_phi, start_col_phi + (K * p) - 1)) = var_regs;

      // 2. Insert the VMA regressors in the right block
      if (type == "diagonal") {
        // Each equation has its own theta parameters
        // Order: [Theta_eq1, Theta_eq2, ...] after the Phi's
        int start_col_theta = n_params_phi + k * q;
        Z_t(k, span(start_col_theta, start_col_theta + q - 1)) =
          get_ma_regressors(U, t, q, k);
      } else if (type == "final") {
        // The common theta parameters are at the end
        // Order: [Theta_common] after the Phi's
        int start_col_theta = n_params_phi;
        // Note: for the Final MA, the regressor is the innovation of equation k
        // but the coefficient is common to all equations.
        Z_t(k, span(start_col_theta, start_col_theta + q - 1)) =
          get_ma_regressors(U, t, q, k);
      }
    }

    // Cumulate GLS: A += Z_t.t() * SigmaU_inv * Z_t
    // Cumulate b:   b += Z_t.t() * SigmaU_inv * y_t
    arma::mat W = Z_t.t() * SigmaU_inv;
    A += W * Z_t;
    b += W * Y.row(t).t();
  }

  // Least squares solution
  arma::vec gamma_hat = solve(A, b);

  return List::create(Named("coefficients") = gamma_hat,
                      Named("A_matrix") = A); // Approx inverse Fisher information
}


//' Autocovariance function using Sbrana's formula
//'
//' @param ar_in a 3d-array with the p AR coefficient matrices
//' @param ma_in a 3d-array with the q MA coefficient matrices
//' @param cov covariance matrix of the innovations
//' @param max_lag non-negative integer with max lag of the autocovariance function
//'
//' @returns A 3d array with the max_lag+1 autocovariance matrices
// [[Rcpp::export]]
arma::cube autocov_gs_cpp(Rcpp::Nullable<arma::cube> ar_in,
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

    // Construct Matrix G
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

// McElroy Method
// Helper function: Equivalent to Kcommut in R
arma::mat Kcommut(const arma::vec& vect, int m, int n) {
  arma::mat temp(const_cast<double*>(vect.begin()), m, n, false, true);
  return arma::vectorise(temp.t());
}


// Helper function: Polynomial multiplication of matrices
arma::cube polymulMat(const arma::cube& amat, const arma::cube& bmat) {
  int p = amat.n_slices;
  int q = bmat.n_slices;
  int m = amat.n_rows;

  arma::cube amatd(m, m, p);
  for(int i = 0; i < p; ++i) {
    amatd.slice(i) = amat.slice(p - 1 - i);
  }

  arma::cube amatd_ext(m, m, p + q - 1, arma::fill::zeros);
  for(int i = 0; i < p; ++i) {
    amatd_ext.slice(q - 1 + i) = amatd.slice(i);
  }
  amatd = amatd_ext;

  arma::mat bigmat(m * (p + q - 1), m * (p + q - 1), arma::fill::zeros);
  arma::cube current_amatd = amatd;

  for(int i = 0; i < p + q - 1; ++i) {
    arma::mat nextmat(m, m * (p + q - 1));
    for(int j = 0; j < p + q - 1; ++j) {
      nextmat.cols(j * m, (j + 1) * m - 1) = current_amatd.slice(j);
    }

    int row_idx = p + q - 2 - i;
    bigmat.rows(row_idx * m, (row_idx + 1) * m - 1) = nextmat;

    arma::cube next_amatd(m, m, p + q - 1, arma::fill::zeros);
    for(int j = 0; j < p + q - 2; ++j) {
      next_amatd.slice(j) = current_amatd.slice(j + 1);
    }
    current_amatd = next_amatd;
  }

  arma::mat bigmat_sub = bigmat.cols(0, m * q - 1);

  arma::mat bmat_rev_mat(m, m * q);
  for(int i = 0; i < q; ++i) {
    bmat_rev_mat.cols(i * m, (i + 1) * m - 1) = bmat.slice(q - 1 - i);
  }

  // Multiply matrices
  arma::mat out_mat = bigmat_sub * bmat_rev_mat.t();

  // CORRECTED RESHAPE: Map row chunks directly to cube slices
  arma::cube out(m, m, p + q - 1);
  for(int slice = 0; slice < p + q - 1; ++slice) {
    out.slice(slice) = out_mat.rows(slice * m, (slice + 1) * m - 1);
  }

  return out;
}

//' Autocovariance function using Tucker-McElroy method
//'
//' @param ar_in a 3d-array with the p AR coefficient matrices
//' @param ma_in a 3d-array with the q MA coefficient matrices
//' @param cov covariance matrix of the innovations
//' @param max_lag non-negative integer with max lag of the autocovariance function
//'
//' @returns A 3d array with the max_lag+1 autocovariance matrices
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube autocov_mc_cpp(const arma::cube& ar, const arma::cube& ma,
                          const arma::mat& cov, int maxlag = 10) {
  int m = cov.n_rows;
  int p = ar.n_slices;
  int q = ma.n_slices;
  int m2 = m * m;

  arma::mat Kmat(m2, m2);
  arma::mat I_m2 = arma::eye(m2, m2);
  for(int i = 0; i < m2; ++i) {
    Kmat.col(i) = Kcommut(I_m2.col(i), m, m);
  }

  arma::cube gamMA;
  if (q == 0) {
    gamMA = arma::cube(m, m, 1);
    gamMA.slice(0) = cov;
  } else {
    arma::cube MA_poly(m, m, q + 1);
    MA_poly.slice(0) = arma::eye(m, m);
    for(int i = 0; i < q; ++i) MA_poly.slice(i + 1) = ma.slice(i);

    arma::cube cov_cube(m, m, 1);
    cov_cube.slice(0) = cov;

    arma::cube temp = polymulMat(MA_poly, cov_cube);
    gamMA = polymulMat(temp, MA_poly);
  }

  arma::cube gamMA_sub(m, m, q + 1);
  if(q == 0) {
    gamMA_sub.slice(0) = gamMA.slice(0);
  } else {
    for(int i = 0; i <= q; ++i) gamMA_sub.slice(i) = gamMA.slice(q + i);
  }
  arma::vec gamMAvec = arma::vectorise(gamMA_sub);

  arma::cube gamARMA_cube;
  arma::cube gamMix_cube;

  if (p > 0) {
    std::vector<arma::mat> Arow_blocks(p + 1);
    for(int i = 0; i < p; ++i) {
      Arow_blocks[p - 1 - i] = -arma::kron(arma::eye(m, m), ar.slice(i));
    }
    Arow_blocks[p] = arma::eye(m2, m2);

    auto get_B = [&](int r, int c) -> arma::mat {
      if(c >= r && c <= r + p) return Arow_blocks[c - r];
      return arma::zeros(m2, m2);
    };

    arma::mat Amat_final((p + 1) * m2, (p + 1) * m2, arma::fill::zeros);
    for(int r = 0; r <= p; ++r) {
      Amat_final.submat(r * m2, 0, (r + 1) * m2 - 1, m2 - 1) = get_B(r, p);
      for(int k = 0; k < p; ++k) {
        arma::mat block = get_B(r, p + 1 + k) + get_B(r, p - 1 - k) * Kmat;
        Amat_final.submat(r * m2, (k + 1) * m2, (r + 1) * m2 - 1, (k + 2) * m2 - 1) = block;
      }
    }

    arma::mat Bmat_final((q + 1) * m2, (q + 1) * m2, arma::fill::zeros);
    for(int r = 0; r <= q; ++r) {
      for(int c = r; c <= std::min(q, r + p); ++c) {
        arma::mat block;
        if(c - r == 0) block = arma::eye(m2, m2);
        else block = -arma::kron(ar.slice(c - r - 1), arma::eye(m, m));
        Bmat_final.submat(r * m2, c * m2, (r + 1) * m2 - 1, (c + 1) * m2 - 1) = block;
      }
    }

    arma::vec gamMix = arma::solve(Bmat_final, gamMAvec);
    arma::vec gamMixTemp;
    if (p <= q) {
      gamMixTemp = gamMix.subvec(0, (p + 1) * m2 - 1);
    } else {
      gamMixTemp = arma::zeros((p + 1) * m2);
      gamMixTemp.subvec(0, gamMix.n_elem - 1) = gamMix;
    }

    arma::vec gamARMA_vec = arma::solve(Amat_final, gamMixTemp);

    gamMix_cube = arma::cube(m, m, q + 1);
    for(int i = 0; i <= q; ++i) gamMix_cube.slice(i) = arma::reshape(gamMix.subvec(i * m2, (i + 1) * m2 - 1), m, m);

    gamARMA_cube = arma::cube(m, m, p + 1);
    for(int i = 0; i <= p; ++i) gamARMA_cube.slice(i) = arma::reshape(gamARMA_vec.subvec(i * m2, (i + 1) * m2 - 1), m, m);

  } else {
    gamARMA_cube = arma::cube(m, m, 1);
    gamARMA_cube.slice(0) = gamMA_sub.slice(0);
    if (q == 0) {
      gamMix_cube = arma::cube(m, m, 1);
      gamMix_cube.slice(0) = cov;
    } else {
      gamMix_cube = arma::cube(m, m, q + 1);
      for(int i = 0; i <= q; ++i) gamMix_cube.slice(i) = gamMA_sub.slice(i);
    }
  }

  if (maxlag <= p) {
    arma::cube temp_ARMA(m, m, maxlag + 1);
    for(int i = 0; i <= maxlag; ++i) temp_ARMA.slice(i) = gamARMA_cube.slice(i);
    gamARMA_cube = temp_ARMA;
  } else {
    if (maxlag > q) {
      arma::cube temp_Mix(m, m, maxlag + 1, arma::fill::zeros);
      for(int i = 0; i <= q; ++i) temp_Mix.slice(i) = gamMix_cube.slice(i);
      gamMix_cube = temp_Mix;
    }

    // Pre-allocate the full output cube once; fill recursively.
    int n_init = gamARMA_cube.n_slices;
    arma::cube out_ARMA(m, m, maxlag + 1, arma::fill::zeros);
    for(int i = 0; i < n_init; ++i) out_ARMA.slice(i) = gamARMA_cube.slice(i);

    arma::mat ar_mat(m, m * p, arma::fill::zeros);
    if (p > 0) {
      for(int i = 0; i < p; ++i) ar_mat.cols(i * m, (i + 1) * m - 1) = ar.slice(i);
    }

    for(int k = 1; k <= maxlag - p; ++k) {
      int cur = n_init + k - 1;
      arma::mat acf = gamMix_cube.slice(p + k);
      if (p > 0) {
        arma::mat temp_mat(m * p, m);
        for(int i = 0; i < p; ++i) {
          temp_mat.rows(i * m, (i + 1) * m - 1) = out_ARMA.slice(cur - 1 - i);
        }
        acf += ar_mat * temp_mat;
      }
      out_ARMA.slice(cur) = acf;
    }
    gamARMA_cube = out_ARMA;
  }

  return gamARMA_cube;
}

//' Autocovariance function using Ansley (1980) method
//'
//' @param ar_in a 3d-array with the p AR coefficient matrices
//' @param ma_in a 3d-array with the q MA coefficient matrices
//' @param cov covariance matrix of the innovations
//' @param maxlag non-negative integer with max lag of the autocovariance function
//'
//' @returns A 3d array with the max_lag+1 autocovariance matrices
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube autocov_an_cpp(arma::cube phi,  arma::cube theta,
                          arma::mat cov, int maxlag) {
  int p = phi.n_slices;
  int q = theta.n_slices;
  int k_max = std::max(p, q);
  int m = cov.n_rows;

  // 1. Compute Cross-Covariances E[y_t epsilon_{t-j}']
  // We use a cube to store these intermediate Tau matrices
  cube tau(m, m, k_max + 1, fill::zeros);
  tau.slice(0) = cov;

  for (int j = 1; j <= k_max; ++j) {
    mat temp = zeros(m, m);
    if (j <= q) temp += theta.slice(j-1) * cov;
    for (int i = 1; i <= std::min(j, p); ++i) {
      temp += phi.slice(i-1) * tau.slice(j-i);
    }
    tau.slice(j) = temp;
  }

  // 2. Setup the Linear System for Gamma(0) ... Gamma(p)
  // System dimension: (p+1) matrices of size m x m
  int sys_dim = m * m * (p + 1);
  arma::mat A = eye(sys_dim, sys_dim);
  arma::vec b = zeros(sys_dim);

  for (int h = 0; h <= p; ++h) {
    // Build RHS: sum_{j=h}^q Theta_j Sigma Psi_{j-h}'
    mat rhs_h = zeros(m, m);
    for (int j = h; j <= q; ++j) {
      // Note: theta.slice(j-1) represents Theta_j
      mat Theta_j = (j == 0) ? eye(m, m) : theta.slice(j-1);
      rhs_h += Theta_j * tau.slice(j-h).t();
    }

    int row_start = h * m * m;
    b.subvec(row_start, row_start + m * m - 1) = vectorise(rhs_h);

    // Build LHS Matrix mapping the relationship between Gamma lags
    for (int i = 1; i <= p; ++i) {
      mat Phi_i = phi.slice(i-1);
      int target_lag = std::abs(h - i);
      int col_start = target_lag * m * m;

      mat term;
      if (h >= i) {
        // Standard case: Phi_i * Gamma(h-i)
        term = kron(eye(m, m), Phi_i);
      } else {
        // Symmetric case: Phi_i * Gamma(i-h)'
        // Uses the property vec(A X' B) = (B \otimes A) vec(X')
        // For Gamma(k)', a commutation matrix or index swap is needed
        term = kron(Phi_i, eye(m, m));
      }
      A.submat(row_start, col_start, row_start + m * m - 1, col_start + m * m - 1) -= term;
    }
  }

  // 3. Solve the system (vec_gamma contains Gamma_0 to Gamma_p)
  arma::vec vec_gamma = solve(A, b);

  // 4. Populate output cube and extend via AR recursion
  cube gamma_out(m, m, maxlag + 1, fill::zeros);
  for (int h = 0; h <= std::min(p, maxlag); ++h) {
    gamma_out.slice(h) = reshape(vec_gamma.subvec(h * m * m, (h + 1) * m * m - 1), m, m);
  }

  // Apply recursive Yule-Walker for lags > p
  for (int h = p + 1; h <= maxlag; ++h) {
    mat gh = zeros(m, m);
    for (int i = 1; i <= p; ++i) {
      gh += phi.slice(i-1) * gamma_out.slice(h-i);
    }
    gamma_out.slice(h) = gh;
  }

  return gamma_out;
}

///////////////////////////////////////////////////////////////////////////////
// Implementation of Mittnik's algorithm to compute the autocovariance function
///////////////////////////////////////////////////////////////////////////////

// mittnik1990_acov.cpp
//
// RcppArmadillo implementation of:
//   Mittnik, S. (1990). Computation of Theoretical Autocovariance Matrices of
//   Multivariate Autoregressive Moving Average Time Series.
//   Journal of the Royal Statistical Society: Series B, 52(1), 151-155.
//
// Computes the theoretical autocovariance matrices Gamma_0, ..., Gamma_{max_lag}
// for a stationary VARMA(p,q) process:
//
//   A(L) y_t = B(L) eps_t
//
// where:
//   A(L) = I - A_1*L - ... - A_p*L^p
//   B(L) = I + B_1*L + ... + B_q*L^q       (B_0 = I assumed throughout)
//   E[eps_t eps_s'] = delta_{ts} * Sigma
//
// Input conventions (identical to Sbrana's autocov_gs_cpp):
//   ar_in  : (m x m x p) cube of AR coefficients A_1, ..., A_p, or NULL for pure VMA.
//   ma_in  : (m x m x q) cube of MA coefficients B_1, ..., B_q, or NULL for pure VAR.
//   cov    : (m x m) innovation covariance matrix Sigma.
//   max_lag: number of lags to return beyond lag 0.
//
// Returns: (m x m x (max_lag+1)) cube; slice k = Gamma_k = E[y_t y_{t-k}^T].

// ============================================================
// Helper: MA-infinity (Wold) coefficients C_0, C_1, ..., C_{maxC-1}
//
//   C_0 = I  (since B_0 = I)
//   C_i = B_i + sum_{j=1}^{min(i,p)} A_j * C_{i-j}    for 1 <= i <= q
//   C_i =       sum_{j=1}^{min(i,p)} A_j * C_{i-j}    for i > q
//
// A : m x m x max(p,1)  (zero-padded when p=0)
// B : m x m x (q+1),    B.slice(0) = I, B.slice(k) = B_k for k=1..q
// ============================================================
static cube build_C_coefs(const cube& A, const cube& B, int maxC)
{
  const int m = A.n_rows;
  const int p = A.n_slices;
  const int q = (int)B.n_slices - 1;

  cube C(m, m, maxC, fill::zeros);
  for (int i = 0; i < maxC; i++) {
    mat Ci = (i <= q) ? mat(B.slice(i)) : mat(m, m, fill::zeros);
    for (int j = 1; j <= std::min(i, p); j++)
      Ci += A.slice(j - 1) * C.slice(i - j);
    C.slice(i) = Ci;
  }
  return C;
}

// ============================================================
// Helper: NC* matrices for lags tau = 0, ..., pe
//
// From eq. (2): NC*_tau = sum_{j=tau}^{q} B_j * Sigma * C_{j-tau}^T
//
// Derivation: E[eps_{t-j} y_{t-tau}^T] = Sigma * C_{j-tau}^T  (for j >= tau)
// follows from y_{t-tau} = sum_k C_k eps_{t-tau-k} and E[eps eps^T] = Sigma.
//
// Returns an m*(pe+1) x m matrix; block row tau is NC*_tau.
// ============================================================
static mat build_NC_star(const cube& B, const cube& C, int pe, int q, const mat& Sigma)
{
  const int m = B.n_rows;
  mat NC(m * (pe + 1), m, fill::zeros);
  for (int tau = 0; tau <= pe; tau++) {
    mat block(m, m, fill::zeros);
    for (int j = tau; j <= q; j++)
      block += B.slice(j) * Sigma * C.slice(j - tau).t();
    NC.rows(tau * m, (tau + 1) * m - 1) = block;
  }
  return NC;
}

//' Compute theoretical autocovariance matrices for a VARMA(p,q) process
//' using the method of Mittnik (1990). Identical call convention to
//' \code{autocov_gs_cpp}: \code{ma_in} holds B_1,...,B_q only; B_0 = I is implicit.
//'
//' @param ar_in   3D array (m x m x p) of AR coefficients A_1,...,A_p, or NULL.
//' @param ma_in   3D array (m x m x q) of MA coefficients B_1,...,B_q, or NULL.
//' @param cov     Innovation covariance matrix Sigma (m x m).
//' @param max_lag Number of lags to return beyond lag 0.
//' @return 3D array (m x m x (max_lag+1)); slice k+1 (1-based in R) is Gamma_k.
//'
//' @examples
//' \dontrun{
//' sourceCpp("mittnik1990_acov.cpp")
//'
//' # Univariate AR(1): Gamma_k = phi^k / (1 - phi^2)
//' G <- varma_acov(array(0.7, c(1,1,1)), NULL, matrix(1), max_lag = 5)
//'
//' # Bivariate VARMA(1,1)
//' m  <- 2
//' A1 <- matrix(c(0.4, 0.2, 0.1, 0.3), m, m)
//' B1 <- matrix(c(0.2, 0.1, 0.0, 0.2), m, m)
//' G  <- varma_acov(array(A1, c(m,m,1)), array(B1, c(m,m,1)), diag(m), max_lag = 8)
//' G[,,1]  # Gamma_0 (symmetric)
//' G[,,2]  # Gamma_1
//' }
// [[Rcpp::export]]
arma::cube autocov_mi_cpp(Rcpp::Nullable<arma::cube> ar_in,
                          Rcpp::Nullable<arma::cube> ma_in,
                          arma::mat cov,
                          int max_lag = 0)
{
   // ---- Unpack inputs ----
  int p = 0, q = 0;
  arma::cube AR, MA_raw;
  if (ar_in.isNotNull()) { AR     = Rcpp::as<cube>(ar_in); p = AR.n_slices; }
  if (ma_in.isNotNull()) { MA_raw = Rcpp::as<cube>(ma_in); q = MA_raw.n_slices; }

  if (p == 0 && q == 0)
    Rcpp::stop("At least one of ar_in or ma_in must be non-NULL and non-empty.");

   const int m = (p > 0) ? (int)AR.n_rows : (int)MA_raw.n_rows;

   // A: zero-padded cube used when p=0 (satisfies build_C_coefs signature)
   arma::cube A(m, m, std::max(p, 1), fill::zeros);
   for (int k = 0; k < p; k++)
     A.slice(k) = AR.slice(k);

   // B: prepend B_0 = I; B.slice(k) = B_k for k = 0..q
   arma::cube B(m, m, q + 1, fill::zeros);
   B.slice(0).eye();
   for (int k = 1; k <= q; k++)
     B.slice(k) = MA_raw.slice(k - 1);

   // ---- MA-infinity coefficients ----
   const int pe   = std::max(p, 1);   // effective AR order (>= 1 to avoid empty system)
   const int maxC = std::max({p, q, max_lag}) + 2;
   const arma::cube C_inf = build_C_coefs(A, B, maxC);

   // ---- NC* right-hand side ----
   const arma::mat NCstar = build_NC_star(B, C_inf, pe, q, cov);

   // ---- Vec-form linear system for Gamma_0 ... Gamma_pe ----
   //
   // Unknown: x_tau = vec(Gamma_tau^T),  tau = 0..pe,  stacked into a vector of
   // length m^2*(pe+1).
   //
   // Yule-Walker recurrence (eq. 2, using Gamma_{-k} = Gamma_k^T):
   //   Gamma_tau = sum_{k<=tau} A_k Gamma_{tau-k}
   //             + sum_{k>tau}  A_k Gamma_{k-tau}^T
   //             + NC*_tau
   //
   // Vectorising x_tau = vec(Gamma_tau^T):
   //   k <= tau:  coefficient on x_{tau-k} = kron(A_k, I_m)
   //   k >  tau:  coefficient on x_{k-tau} = kron(A_k, I_m) * Wm
   // where Wm is the m^2 x m^2 commutation matrix, vec(X^T) = Wm * vec(X).
   //
   // Proofs:
   //   k <= tau: vec((A_k Gamma_{tau-k})^T) = vec(Gamma_{tau-k}^T A_k^T)
   //             = kron(A_k, I_m) * vec(Gamma_{tau-k}^T)    [vec(XM) = (M^T ox I)vec(X)]
   //
   //   k >  tau: vec((A_k Gamma_{k-tau}^T)^T) = vec(Gamma_{k-tau} A_k^T)
   //             = kron(A_k, I_m) * vec(Gamma_{k-tau})
   //             = kron(A_k, I_m) * Wm * vec(Gamma_{k-tau}^T)

   const int m2     = m * m;
   const int sz_vec = m2 * (pe + 1);

   // Commutation matrix: Wm(i*m+j, j*m+i) = 1
   arma::mat Wm(m2, m2, fill::zeros);
   for (int i = 0; i < m; i++)
     for (int j = 0; j < m; j++)
       Wm(i * m + j, j * m + i) = 1.0;

   const mat I_m(m, m, fill::eye);

   arma::mat Lv(sz_vec, sz_vec, fill::eye);
   arma::vec Rv(sz_vec, fill::zeros);

   for (int tau = 0; tau <= pe; tau++)
     Rv.subvec(tau * m2, (tau + 1) * m2 - 1) =
       vectorise(NCstar.rows(tau * m, (tau + 1) * m - 1).t());

   for (int tau = 0; tau <= pe; tau++) {
     for (int k = 1; k <= p; k++) {
       const mat kAI = kron(A.slice(k - 1), I_m);
       if (tau - k >= 0)
         Lv.submat(tau * m2, (tau - k) * m2,
                   (tau + 1) * m2 - 1, (tau - k + 1) * m2 - 1) -= kAI;
       if (k - tau > 0 && k - tau <= pe)
         Lv.submat(tau * m2, (k - tau) * m2,
                   (tau + 1) * m2 - 1, (k - tau + 1) * m2 - 1) -= kAI * Wm;
     }
   }

   // ---- Solve ----
   const arma::vec gam_vec = solve(Lv, Rv);

   // ---- Unpack Gamma_0 ... Gamma_pe ----
   const int n_init = pe + 1;
   std::vector<mat> Gamma(std::max(n_init, max_lag + 1));
   for (int tau = 0; tau <= pe; tau++) {
     const arma::vec gv = gam_vec.subvec(tau * m2, (tau + 1) * m2 - 1);
     Gamma[tau] = reshape(gv, m, m).t();   // reshape gives Gamma_tau^T, .t() inverts
   }

   // ---- Higher lags by recursion (eq. 10) ----
   // Gamma_tau = sum_{k=1}^p A_k Gamma_{tau-k} + sum_{j=tau}^q B_j Sigma C_{j-tau}^T
   for (int tau = n_init; tau <= max_lag; tau++) {
     mat G(m, m, fill::zeros);
     for (int k = 1; k <= p; k++) {
       if (tau - k >= 0) G += A.slice(k - 1) * Gamma[tau - k];
     }
     for (int j = tau; j <= q; j++) {
       G += B.slice(j) * cov * C_inf.slice(j - tau).t();
     }
     Gamma[tau] = G;
   }

   // ---- Pack output ----
   cube out(m, m, max_lag + 1);
   for (int k = 0; k <= max_lag; k++)
     out.slice(k) = Gamma[k];
   return out;
}

// ----------------------------------------------------------------------------
// Helper functions for Reinsel, Basu, Yap (1992) VARMA estimation ------------
// ----------------------------------------------------------------------------

// Filter a time series with the inverse of the MA polynomial
// Solves: (I + Theta_1 L + ... + Theta_q L^q) X = Input,
// Thus: X_t = Input_t - sum_j Theta_j X_{t-j}
mat filter_series_inverse(const arma::mat& Input, const arma::cube& Theta_cube,
                          int q, int K) {
  int T = Input.n_rows;
  arma::mat X(T, K, fill::zeros);

  for (int t = 0; t < T; t++) {
    arma::vec val = Input.row(t).t();
    for (int j = 1; j <= q; j++) {
      if (t - j >= 0) {
        val -= Theta_cube.slice(j - 1) * X.row(t - j).t();
      }
    }
    X.row(t) = val.t();
  }
  return X;
}

// Gauss-Newton one-step optimization for VARMA(p, 1)
//   Y_t = [c] + Phi_1 Y_{t-1} + ... + Phi_p Y_{t-p} + E_t + Theta_1 E_{t-1} + ... + Theta_q E_{t-q}
//
// Phi_cube  : cube K x K x p (slice i-1 = Phi_i)
// Theta_cube: cube K x K x q (slice j-1 = Theta_j)
// intercept : if true, c_vec (K x 1) is subtracted when computing the resid
//             and we add a scalar regressor, 1, at the end of z_t
//
// Layout rows of ZZ/ZE (and of Delta = ZZ^{-1} ZE in R):
//   [K*p rows for Phi | K*q rows for Theta | 1 row for c (if intercept=true)]
// [[Rcpp::export]]
List rby_optimization_step(arma::mat Y, arma::cube Phi_cube, arma::cube Theta_cube,
                           int p, int q, bool intercept, arma::vec c_vec) {
  int T = Y.n_rows;
  int K = Y.n_cols;
  int t_start = std::max(p, q);
  int n_reg = K * p + K * q + (intercept ? 1 : 0);

  // 1. Compute resid E_t: E_t = Y_t - [c] - sum Phi_i Y_{t-i} - sum Theta_j E_{t-j}
  arma::mat E(T, K, fill::zeros);

  for (int t = t_start; t < T; t++) {
    arma::vec e_t = Y.row(t).t();

    if (intercept) {
      e_t -= c_vec;
    }
    for (int i = 1; i <= p; i++) {
      e_t -= Phi_cube.slice(i - 1) * Y.row(t - i).t();
    }
    for (int j = 1; j <= q; j++) {
      e_t -= Theta_cube.slice(j - 1) * E.row(t - j).t();
    }
    E.row(t) = e_t.t();
  }

  // 2. Aux series filtered by Theta(L)^{-1}
  //    U_aux: filtered Y -> regressors for the AR part
  //    V_aux: filtered E -> regressors for the MA part
  mat U_aux = filter_series_inverse(Y, Theta_cube, q, K);
  mat V_aux = filter_series_inverse(E, Theta_cube, q, K);

  // 3. Build the Gauss-Newton matrices ZZ e ZE
  //    z_t = [U_{t-1}',...,U_{t-p}', V_{t-1}',...,V_{t-q}', 1 (scalar) if intercept]
  arma::mat ZZ(n_reg, n_reg, fill::zeros);
  arma::mat ZE(n_reg, K,     fill::zeros);

  for (int t = t_start; t < T; t++) {
    rowvec z_t(n_reg);
    int idx = 0;

    for (int i = 1; i <= p; i++) {
      z_t.subvec(idx, idx + K - 1) = U_aux.row(t - i);
      idx += K;
    }
    for (int j = 1; j <= q; j++) {
      z_t.subvec(idx, idx + K - 1) = V_aux.row(t - j);
      idx += K;
    }
    if (intercept) {
      z_t[idx] = 1.0;
    }

    ZZ += z_t.t() * z_t;
    ZE += z_t.t() * E.row(t);
  }

  return List::create(Named("ZZ") = ZZ, Named("ZE") = ZE, Named("E") = E);
}

