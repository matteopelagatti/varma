# varma 0.1.4

* Added `fit_var_bic()`: selects pure VAR orders by BIC on a common sample,
  then refits by OLS using all observations available at the selected order.
  The automatic upper bound is ceiling(log(n)^1.5), with a residual-degrees-
  of-freedom cap; fixed upper bounds and VAR(0) candidates are optional.
* Paper simulations now include `VAR_BIC`, without using population VARMA
  orders. Checkpoints retain selected lags and all candidate BICs. Reporting
  includes the VAR accuracy curves, lag frequencies and boundary diagnostics.

# varma 0.1.3

* Added `fit_varma_cml()` for unrestricted conditional Gaussian ML, sharing
  the echelon estimator's HR1-style pilot, compiled likelihood, analytic
  gradient and BFGS optimizer. Full and echelon estimates coincide when their
  coefficient patterns coincide.
* State-space ML (`fit_varma_skf`, `fit_varma_kfas`, `fit_varma_fkf`) now uses
  HR1 starting values by default. `pilot="IHR"` retains the iterated option.
* `initial_cov="stationary"` recomputes the initial state covariance using
  the adaptive Lyapunov solver at every candidate; `initial_cov="fixed"`
  freezes the covariance implied by the stabilized pilot. Innovation
  covariance parameters are estimated in both cases. Bootstrap refits retain
  the selected policy and recompute their own pilot covariance.
* Paper simulations and R reporting include CML. New experiments require
  fresh checkpoints; an explicit `restart` driver mode preserves the old
  experiment before starting a new one in the same output directory.

# varma 0.1.0

## Portable paper simulations (10 November 2026)

* `varma_paper_designs()` reproduces the four full-precision simulation models.
  `run_varma_paper_simulation()` retains the original estimators, defaults to
  compiled SKF MLE, and compares oracle and data-selected echelon estimates.
  A true-parameter forecast benchmark conditions on the same observed history.
* Checkpoints support interrupted runs and worker-independent RNG streams.
  `summarize_varma_paper_simulation()` reports untrimmed IRF/forecast losses,
  Monte Carlo errors, paired risk ratios, failures, convergence, and timings.
  `benchmark_varma_paper_simulation()` projects the full runtime from pilot
  replications with the actual estimator settings. A portable driver and guide
  are supplied under `inst/scripts`.

## Internal Kronecker-index selection (10 November 2026)

* `select_kronecker_indices()` estimates indices with compiled sequential
  canonical-correlation tests. It factors the past matrix once per horizon,
  incrementally updates the future basis, and calculates canonical scores
  without solving for coefficient vectors. Test details are returned silently
  unless `verbose=TRUE`.
* `fit_varma_echelon()` and bootstrap reselection use this internal selector.
  MTS has been removed from the package dependencies and runnable examples.
  The simulation script's echelon wrapper now uses the internal coefficient
  estimator as well; historical simulation results have not been regenerated.
* Selection supports a single series, discards rejected future columns at
  horizon zero as well as later horizons, and reports rank, sample-size,
  past-dimension, and serial-correction failures explicitly. These edge cases
  can differ from MTS 1.2.1. The test definition and series order are preserved
  in regular cases; the algorithm is an independent implementation.

## Faster estimation and echelon estimation (10 November 2026)

* `fit_varma_echelon()` fits a conditional Gaussian echelon model with supplied
  Kronecker indices or indices selected from the data. It uses
  a compiled analytic gradient, preserves series ordering, and supports zero
  indices and white noise. Fits work with the existing IRF, forecast, and
  residual-bootstrap interfaces. Bootstrap refits retain the selected indices
  unless `kidx = NULL` is explicitly passed to repeat selection.
* `solve_discrete_lyapunov()` adds dense doubling and Schur solvers, avoiding R
  sparse-matrix conversion. `autocov_ss()` computes stationary VARMA covariance
  matrices without a squared-dimension Kronecker system or a state-covariance
  array for every horizon. The stationary covariance is distinct from the
  conditional covariance of the Kalman filter.
* Likelihood parameter unpacking, stability checks, and stationary covariance
  construction now run in C++. The sparse-backend objective stays in C++ until
  it returns a scalar. KFAS/FKF retain their existing likelihood backends.
* IHR reuses the fixed regression design and QR residuals; NET obtains the design
  directly from its pilot. DPD/DPF assemble GLS equations from reusable cross
  products. RBY computes derivatives only for accepted line-search candidates.
* Regression checks compare the optimizations with independent Kronecker solves,
  direct GLS assembly, previous fitted results, and numerical derivatives.

## Estimation and inference corrections (10 November 2026)

* DPD/DPF now use centred data in both stages when an intercept is requested.
  They are documented as a pilot VAR followed by one feasible GLS regression.
* RBY uses the full covariance-weighted innovation derivative, retaining matrix
  multiplication order and cross-equation effects. The intercept is estimated
  jointly, and the line search increases the conditional Gaussian likelihood.
* KFAS, FKF.SP, and the sparse Kalman likelihood share parameter packing and
  recompute the stationary initial covariance for every candidate. Stable AR and
  invertible MA polynomials are required. `control` passes optimiser settings,
  and `optim.out`/`convergence` report numerical termination.
* Sparse Kalman output uses fitted coefficients, covariance, and transition
  matrices throughout. Acceleration waits for convergence of the full state
  covariance, fixing false convergence for delayed MA responses.
* KFAS and FKF.SP fits store joint vector innovations, suitable for joint residual
  resampling, instead of componentwise sequential prediction errors.
* NET honours a fixed `lambda` or a candidate path and supports its parallel
  arguments. SUR fits refresh residuals and covariance. Empty selected equations
  and a single lagged regressor are supported.
* IHR uses a logarithmically growing default pilot order and excludes unavailable
  pilot errors from regression. Its stopping criterion compares regression
  errors; aliased regression coefficients have a defined least-squares solution.
  Returned innovations are recomputed from the final coefficients, including
  when `maxit = 1`.
* Likelihoods, covariance normalisation, parameter counts, and information
  criteria are consistent. `logLik()` has class `logLik`; `nobs()` now counts
  **vector observations**, not scalar data entries. Conditional order comparisons
  use a common trimmed sample. The printed HQC penalty is corrected.
* Single-series and pure-MA estimation paths preserve matrix dimensions.
  Prediction covariance always has three dimensions, including at horizon one.
* Generalized IRFs scale the innovation covariance by columns. The Pearlman
  state-space loading is correct for pure VARs.
* Identification is checked through transfer-function coefficient restrictions,
  rather than shared scalar roots. Kronecker indices retain series order and
  use a sufficiently wide Hankel matrix, avoiding premature rank saturation.
* Residual bootstrap centres residual vectors, retains the fitted mean, supports
  burn-in, handles `resids_to_skip = 0`, and restores the caller's future plan.
  Presample residuals are explicitly marked and excluded.

These changes can change the Monte Carlo results. Existing simulation artifacts
should be regenerated before publication. Numerical convergence, the nominal
parameter counts used by information criteria, and a functioning bootstrap do
not by themselves establish inference validity at singular distributions.

`tests/regression-estimation.R` checks the corrections against numerical
derivatives, direct recursions, translation invariance, and independent Kalman
backends. The original diagnostic evidence is retained under `audit/2026-11-09`.
