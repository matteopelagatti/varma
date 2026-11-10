# Reproducible diagnostic audit of the CURRENT checkout, not the installed DLL.
# Run from the repository root with Rscript --vanilla audit/2026-11-09/estimation_audit.R
# This script leaves production sources unchanged. Findings are printed, not fixed.
local_lib <- "C:/Users/Matteo/AppData/Local/R/win-library/4.5"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages(library(KFAS))
suppressPackageStartupMessages(library(Matrix))
Rcpp::sourceCpp("src/rcpp_functions.cpp", rebuild = TRUE, showOutput = FALSE)
source("R/functions.R", local = .GlobalEnv)

out_dir <- "audit/2026-11-09"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "results.txt"), split = TRUE)
options(digits = 9, width = 110)
cat("VARMA source audit, 2026-11-09\n")
print(tools::md5sum(c("R/functions.R", "src/rcpp_functions.cpp")))
cat("Native routines freshly compiled from src/rcpp_functions.cpp.\n")
run <- function(name, expr) {
  cat("\n--- ", name, " ---\n", sep = "")
  tryCatch(force(expr), error = function(e) cat("CHECK ERROR: ", conditionMessage(e), "\n", sep = ""))
  invisible(NULL)
}
diffmax <- function(x, y) max(abs(as.numeric(x) - as.numeric(y)))
coefmat <- function(f) {
  m <- dim(f)[1]
  t(cbind(if (!is.null(f$intercept)) f$intercept,
          if (length(f$ar)) matrix(f$ar, m),
          if (length(f$ma)) matrix(f$ma, m)))
}
ll_gaussian <- function(E, S) {
  -0.5 * (length(E) * log(2*pi) + nrow(E) * as.numeric(determinant(S, TRUE)$modulus) +
            sum((E %*% solve(S)) * E))
}
set.seed(20261109)
m <- 2L
A <- array(matrix(c(.5, .15, -.2, .35), 2, byrow = TRUE), c(2, 2, 1))
B <- array(matrix(c(.2, -.12, .08, .25), 2, byrow = TRUE), c(2, 2, 1))
S <- matrix(c(1, .6, .6, 2), 2)
eps <- t(chol(S)) %*% matrix(rnorm(2*440), 2)
Y <- t(sim_varma_rcpp(A, B, eps))[-seq_len(200), , drop = FALSE]
colnames(Y) <- c("y1", "y2")
Yc <- sweep(Y, 2, colMeans(Y), "-")
Yshift <- sweep(Yc, 2, c(8, -5), "+")

run("DPD and DPF: translation equivariance when intercept=TRUE", {
  for (name in c("fit_varma_dpd", "fit_varma_dpf")) {
    fn <- get(name)
    f0 <- fn(Yc, 1, 1, intercept = TRUE, r = 8)
    f1 <- fn(Yshift, 1, 1, intercept = TRUE, r = 8)
    cat(name, "max AR change:", diffmax(f0$ar, f1$ar),
        "max MA change:", diffmax(f0$ma, f1$ma), "\n")
  }
})

run("SKF: returned covariance and transition match starting rather than final coefficients", {
  ihr <- fit_varma_ihr(Y, 1, 1, intercept = FALSE, maxit = 20)
  skf <- fit_varma_skf(Y, 1, 1, intercept = FALSE, maxit = 20)
  ss0 <- varma_to_ss(ihr)
  ssf <- varma_to_ss(skf)
  cat("AR coefficients final versus initial:", diffmax(skf$ar, ihr$ar), "\n")
  cat("returned covariance versus initial covariance:", diffmax(skf$cov, ihr$cov), "\n")
  cat("returned transition versus initial transition:", diffmax(skf$transition_matrix, ss0$T), "\n")
  cat("returned transition versus final coefficient transition:", diffmax(skf$transition_matrix, ssf$T), "\n")
  orig_forecast <- predict.varma(skf, 3, cov = FALSE)
  skf$transition_matrix <- ssf$T
  cat("3-step forecasts after only aligning T with returned coefficients; max change:",
      diffmax(orig_forecast, predict.varma(skf, 3, cov = FALSE)), "\n")
})

run("Sparse likelihood: false steady state for a stable invertible VMA(2)", {
  b2 <- array(0, c(2, 2, 2)); b2[, , 2] <- diag(c(.8, .6))
  f <- varma(ma = b2, cov = S)
  ss <- varma_to_ss(f)
  TT <- Matrix(ss$T, sparse = TRUE); RR <- Matrix(ss$R, sparse = TRUE)
  P0 <- solve_dlyap_iter(TT, ss$R %*% S %*% t(ss$R))
  set.seed(392)
  y <- sim_varma_rcpp(array(0, c(2, 2, 0)), b2, t(chol(S)) %*% matrix(rnorm(2*300), 2))
  aa <- matrix(0, nrow(TT), 1)
  fast <- kalmanLogLik(TT, RR, S, aa, P0, y)
  exact <- kalmanLogLik(TT, RR, S, aa, P0, y, ss_tol = 0)
  at <- matrix(0, nrow(TT), ncol(y)+1)
  Pt <- array(0, c(nrow(TT), nrow(TT), ncol(y)+1)); Pt[, , 1] <- P0
  v <- matrix(0, 2, ncol(y)); FF <- array(0, c(2, 2, ncol(y)))
  full <- kalman(TT, RR, S, at, Pt, y, v, FF)
  cat("default fast likelihood:", fast, "disabled acceleration:", exact, "full filter:", full, "\n")
  cat("F2-F1 max:", diffmax(FF[, , 2], FF[, , 1]),
      "F3-F2 max:", diffmax(FF[, , 3], FF[, , 2]), "\n")
})

run("RBY: proposed Gauss-Newton gradient versus finite differences", {
  step <- rby_optimization_step(Y, A, B, 1L, 1L, FALSE, numeric(0))
  analytic <- -2 * step$ZE
  numerical <- matrix(0, 4, 2)
  h <- 1e-6
  ssr <- function(aa, bb) sum(varma_residuals(Y, aa, bb, numeric(0))^2)
  for (eq in 1:2) for (j in 1:2) {
    ap <- am <- A; ap[eq, j, 1] <- ap[eq, j, 1] + h; am[eq, j, 1] <- am[eq, j, 1] - h
    numerical[j, eq] <- (ssr(ap, B) - ssr(am, B))/(2*h)
    bp <- bm <- B; bp[eq, j, 1] <- bp[eq, j, 1] + h; bm[eq, j, 1] <- bm[eq, j, 1] - h
    numerical[2+j, eq] <- (ssr(A, bp) - ssr(A, bm))/(2*h)
  }
  cat("-2 ZE (claimed SSR gradient):\n"); print(analytic)
  cat("Numerical SSR gradient:\n"); print(numerical)
  cat("max gradient discrepancy:", diffmax(analytic, numerical), "\n")
})

run("NET: lambda argument is ignored", {
  set.seed(700)
  a <- fit_varma_net(Y, 1, 1, intercept = FALSE, lambda = 1e-9, relax = FALSE, sur = FALSE, r = 8, maxit = 10)
  set.seed(700)
  b <- fit_varma_net(Y, 1, 1, intercept = FALSE, lambda = 1e6, relax = FALSE, sur = FALSE, r = 8, maxit = 10)
  cat("max coefficient difference for lambda=1e-9 versus 1e6:", diffmax(coefmat(a), coefmat(b)), "\n")
  cat("largest absolute coefficient:", max(abs(coefmat(a))), "\n")
})

run("NET: SUR coefficients have pre-SUR residuals and covariance", {
  # Correlated innovations, different equation supports, retained shrinkage.
  set.seed(11)
  f <- fit_varma_net(Y, 1, 1, intercept = FALSE, relax = FALSE, sur = TRUE, r = 8, maxit = 10)
  ih <- fit_varma_ihr(Y, 1, 1, intercept = FALSE, r = 8, maxit = 10, ret = "regression")
  X <- qr.X(ih$qr)
  E <- Y[-1, , drop = FALSE] - X %*% coefmat(f)
  cat("max stored versus regression residual difference:", diffmax(f$residuals, E), "\n")
  cat("max stored versus post-SUR covariance difference:", diffmax(f$cov, crossprod(E)/nrow(E)), "\n")
})

run("IHR: stored one-step regression residuals versus fitted model's innovations", {
  for (iterations in c(1, 100)) {
    f <- fit_varma_ihr(Y, 1, 1, intercept = FALSE, maxit = iterations, r = 8)
    E <- residuals.varma(f)
    cat("maxit", iterations, "max discrepancy:", diffmax(f$residuals, E[-1, ]),
        "last residual discrepancy:", diffmax(tail(f$residuals, 1), tail(E, 1)), "\n")
  }
})

run("IHR: Gaussian likelihood value and parameter count in pure VAR case", {
  f <- fit_varma_ihr(Y, 1, 0, intercept = TRUE)
  cat("reported likelihood:", as.numeric(f$loglik),
      "Gaussian likelihood at returned covariance:", ll_gaussian(f$residuals, f$cov), "\n")
  cat("stored npar:", f$npar, "logLik df:", attr(logLik.varma(f), "df"),
      "coefs plus covariance:", 2*2+2+3, "\n")
})

run("Boundary orders and single-series input", {
  tests <- list(
    IHR_univariate_VAR = function() fit_varma_ihr(Y[, 1, drop = FALSE], 1, 0),
    IHR_univariate_VARMA = function() fit_varma_ihr(Y[, 1, drop = FALSE], 1, 1),
    RBY_pure_MA = function() fit_varma_rby(Y, 0, 1, intercept = FALSE, r = 8, maxit = 2),
    DPD_pure_MA = function() fit_varma_dpd(Y, 0, 1, intercept = FALSE, r = 8),
    DPF_pure_MA = function() fit_varma_dpf(Y, 0, 1, intercept = FALSE, r = 8),
    IHR_pure_MA_control = function() fit_varma_ihr(Y, 0, 1, intercept = FALSE, r = 8)
  )
  for (name in names(tests)) {
    msg <- tryCatch({ tests[[name]](); "OK" }, error = function(e) conditionMessage(e))
    cat(name, ": ", msg, "\n", sep = "")
  }
})

run("KFAS: sequential errors versus multivariate innovations", {
  k <- fit_varma_kfas(Y, 1, 0, intercept = FALSE, ret = "kfas")
  filt <- KFAS::KFS(k$model, filtering = "state", smoothing = "none")
  mv <- KFAS::mvInnovations(filt)
  cat("covariance of stored sequential errors:\n"); print(cov(filt$v))
  cat("covariance of joint prediction errors:\n"); print(cov(mv$v))
  cat("returned innovation covariance:\n"); print(k$model$Q[, , 1])
  cat("max residual difference:", diffmax(filt$v, mv$v), "\n")
  P <- k$model$P1
  TT <- k$model$T[, , 1]; RR <- k$model$R[, , 1]; QQ <- k$model$Q[, , 1]
  cat("stationary Lyapunov residual for fitted model's P1:", diffmax(P, TT %*% P %*% t(TT) + RR %*% QQ %*% t(RR)), "\n")
})

run("KFAS and FKF: likelihood metadata and order selection", {
  k <- fit_varma_kfas(Y, 1, 0, intercept = FALSE)
  cat("exact KFAS npar slot exists:", !is.null(k[["npar", exact = TRUE]]),
      "partial-match npar:", k$npar, "\n")
  f <- fit_varma_fkf(Y, 1, 0, intercept = FALSE)
  cat("FKF n slot exists:", !is.null(f[["n", exact = TRUE]]),
      "dollar n:", f$n, "\n")
  cat("BIC from printed convention:", -2*f$loglik + f$npar*log(f$n),
      "BIC using time points:", -2*f$loglik + f$npar*log(nrow(Y)), "\n")
})

run("Generalized IRF: innovation scaling is applied to the wrong dimension", {
  f <- varma(ar = A, ma = B, cov = S)
  got <- irf(f, 0, orth = "generalized")[, , 1]
  expected <- S %*% diag(1/sqrt(diag(S)))
  cat("package impact matrix:\n"); print(got)
  cat("Sigma times inverse marginal standard deviations:\n"); print(expected)
  cat("max difference:", diffmax(got, expected), "\n")
})

run("Bootstrap: fitted mean is lost and zero skip still drops an observation", {
  f <- fit_varma_ihr(Yshift, 1, 0, intercept = TRUE)
  capture_y <- function(Y, p, q, intercept, ...) Y
  set.seed(510)
  yy <- bootstrap_varma(f, nsim = 1, n = 2000, fit_fn = capture_y)[[1]]
  cat("fitted mean:", solve(diag(2)-f$ar[, , 1], f$intercept), "\n")
  cat("bootstrap path mean:", colMeans(yy), "\n")
  nres <- nrow(f$residuals)
  retained <- f$residuals[-(1:0), , drop = FALSE]
  cat("resids_to_skip=0: input rows", nres, "retained rows", nrow(retained), "\n")
})

run("Pearlman state-space form: pure VAR versus an equivalent zero-MA representation", {
  v <- varma(ar = A, cov = S)
  w <- varma(ar = A, ma = array(0, c(2, 2, 1)), cov = S)
  cat("same pure VAR but zero MA explicitly present: max R difference:",
      diffmax(varma_to_ss(v, "pearlman")$R, varma_to_ss(w, "pearlman")$R), "\n")
})

run("Control: raw and Cholesky IRFs agree with direct recursion", {
  f <- varma(ar = A, ma = B, cov = S)
  psi <- irf(f, 4, "none")
  expected <- array(0, dim(psi)); expected[, , 1] <- diag(2)
  expected[, , 2] <- A[, , 1]+B[, , 1]
  for (j in 3:5) expected[, , j] <- A[, , 1] %*% expected[, , j-1]
  cat("raw IRF max difference:", diffmax(psi, expected), "\n")
  orth <- irf(f, 4, "cholesky")
  for (j in 1:5) expected[, , j] <- expected[, , j] %*% t(chol(S))
  cat("Cholesky IRF max difference:", diffmax(orth, expected), "\n")
})

cat("\n--- Session information ---\n")
print(sessionInfo())
sink()
