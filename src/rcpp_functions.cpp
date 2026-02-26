#include <RcppArmadillo/Lighter>
using namespace Rcpp;
using namespace arma;
const double M_LOG2PI = log(2 * arma::datum::pi);


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

//' Solves the matrix Lypuanov equation
//'
//' Solves \eqn{P = T P T' + R Q R'} for \eqn{P} using the
//' `sylvester` solver.
//'
//' @param T matrix, see the formula above.
//' @param R matrix, see the formula above.
//' @param Q matrix, see the formula above.
//'
//' @returns A matrix with the solution.
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat solve_riccati(const arma::mat& T,
                        const arma::mat& R,
                        const arma::mat& Q) {
  arma::mat B = pinv(T.t());
  arma::mat C = R * Q * R.t();
  arma::mat X;
  bool success = sylvester(X, T, T.t(), C);
  if (!success) stop("Riccati solver failed");
  return X;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat solve_syl(const arma::mat& A, const arma::mat& B, const arma::mat& C) {
  return sylvester(A, B, C);
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
     // y0 = Yt.col(t).eval()(notna);
     v.col(t).eval()(notna) = Yt.col(t).eval()(notna) - at.col(t).eval()(notna);
     PZt = Pt.slice(t).cols(notna);
     F.slice(t) = Pt.slice(t)(notna, notna);
     pass = inv_sympd(iF, F.slice(t));
     if (!pass) return -datum::inf;
     loglik += (real(log_det(F.slice(t))) + v.row(t).t()*iF*v.row(t)).eval()(0,0) + notna.n_elem*M_LOG2PI;
     TPZt = T*PZt;
     K = TPZt*iF;
     at.col(t+1) = T*at.col(t) + K*v.col(t);
     Pt.slice(t+1) = T*Pt.slice(t)*Tt - K*F.slice(t)*K.t() + RQRt;
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
rowvec get_var_regressors(const mat& Y, int t, int p, int K) {
  rowvec reg(K * p, fill::zeros);
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
 rowvec get_ma_regressors(const mat& U, int t, int q, int k_eq) {
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
  mat A(n_total_params, n_total_params, fill::zeros);
  vec b(n_total_params, fill::zeros);

  // for loop on the time t
  for (int t = t_start; t < T; t++) {

    // Making the matrix Z_t (dimensions K x n_total_params)
    mat Z_t(K, n_total_params, fill::zeros);

    // Vector of VAR regressors (equal for all equations but in different positions)
    rowvec var_regs = get_var_regressors(Y, t, p, K); // Dim K*p

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
    mat W = Z_t.t() * SigmaU_inv;
    A += W * Z_t;
    b += W * Y.row(t).t();
  }

  // Least squares solution
  vec gamma_hat = solve(A, b);

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

// Tucker-McElroy Method
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
arma::cube autocov_mc_cpp(const arma::cube& ar, const arma::cube& ma, const arma::mat& cov, int maxlag = 10) {
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

    for(int k = 1; k <= maxlag - p; ++k) {
      int len = gamARMA_cube.n_slices;
      arma::mat acf = gamMix_cube.slice(p + k);
      if (p > 0) {
        arma::mat temp_mat(m * p, m);
        for(int i = 0; i < p; ++i) {
          temp_mat.rows(i * m, (i + 1) * m - 1) = gamARMA_cube.slice(len - 1 - i);
        }
        arma::mat ar_mat(m, m * p);
        for(int i = 0; i < p; ++i) {
          ar_mat.cols(i * m, (i + 1) * m - 1) = ar.slice(i);
        }
        acf += ar_mat * temp_mat;
      }
      arma::cube new_ARMA(m, m, len + 1);
      for(int i = 0; i < len; ++i) new_ARMA.slice(i) = gamARMA_cube.slice(i);
      new_ARMA.slice(len) = acf;
      gamARMA_cube = new_ARMA;
    }
  }

  return gamARMA_cube;
}
