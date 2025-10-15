#include <RcppArmadillo/Lighter>
using namespace Rcpp;
using namespace arma;
// const double M_LOG2PI = 1.8378770664093454835606594728112352797227949472756;
const double M_LOG2PI = log(2 * arma::datum::pi);

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
cube irf_varma_rcpp(const cube& A_cube,
                          const cube& M_cube,
                          const mat& P,
                          int horizon) {
  
  // Determine dimensions and orders
  int m = P.n_rows;
  int p = A_cube.n_slices;
  int q = M_cube.n_slices;

  // Create 3D cubes (arrays) for Psi and IRF
  cube Psi_cube(m, m, horizon + 1, fill::zeros);
  cube IRF_cube(m, m, horizon + 1, fill::zeros);
                  
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



// [[Rcpp::export]]
mat sim_varma_rcpp(const cube& A,
                   const cube& M,
                   const mat& eps) {
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
  
  mat y = eps;
  
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

// [[Rcpp::export]]
arma::mat solve_riccati(const arma::mat& T, const arma::mat& R, const arma::mat& Q) {
  arma::mat B = pinv(T.t());
  arma::mat C = R * Q * R.t();
  arma::mat X;
  bool success = sylvester(X, T, T.t(), C);
  if (!success) stop("Riccati solver failed");
  return X;
}

// [[Rcpp::export]]
arma::mat solve_syl(const arma::mat& A, const arma::mat& B, const arma::mat& C) {
  return sylvester(A, B, C);
}


/////////////////////////////////////////////////////
//' Solve Discrete-time Lyapunov Equation via Iteration
//'
//' Solves the equation P = T * P * T' + C for P, where P is the stationary
//' covariance of the state vector in a state-space model.
//' This is used to initialize the Kalman filter for a stationary process.
//'
//' @param T State transition matrix (sparse).
//' @param C Covariance matrix of the state innovations (dense).
//' @param max_iter Maximum number of iterations.
//' @param tol Tolerance for convergence.
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
//' @param T: sparse matrix (use the function Matrix(T, sparse = TRUE) of the
//' Matrix package to turn the dense matrix T into a sparse one)
//' @param R: sparse matrix (see above)
//' @param Q: dense covariance matrix
//' @param a1: matrix with one column with initial mean values for the state vector
//' @param P1: initial covariance matrix of the state vector
//' @param Yt: trasposed matrix of data \eqn{n\times m}, where \eqn{n} is the
//' number of time points.
//' @param update_state: boolean (by default is false), if true the parameters
//' a1 and P1 are overwritten with the values of the state prediction for time
//' \eqn{t = n+1}: \eqn{a_{n+1|n}} and \eqn{P_{n+1|n}}
//' 
//' @returns a scalar number with the value of the log-likelihood.
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
 //' @param T: sparse matrix (use the function Matrix(T, sparse = TRUE) of the
 //' Matrix package to turn the dense matrix T into a sparse one)
 //' @param R: sparse matrix (see above)
 //' @param Q: dense covariance matrix
 //' @param a1: matrix with one column with initial mean values for the state vector
 //' @param P1: initial covariance matrix of the state vector
 //' @param Yt: trasposed matrix of data \eqn{n\times m}, where \eqn{n} is the
 //' number of time points.
 //' @param update_state: boolean (by default is false), if true the parameters
 //' a1 and P1 are overwritten with the values of the state prediction for time
 //' \eqn{t = n+1}: \eqn{a_{n+1|n}} and \eqn{P_{n+1|n}}
 //' 
 //' @returns a scalar number with the value of the log-likelihood.
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
//' @param X_list list of y-vectors with dependent variables
//' @param y_list list of X-matrices with regressors
//' @param sigma covariance matrix of regression errors
//' 
//' @returns A vector of regression coefficients.
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

