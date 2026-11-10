# Estimation, covariance computation, and echelon models

Update: automatic index selection was subsequently internalized; see
`../kronid/README.md`. The MTS references below describe the initial version
at the time of these performance measurements.

Implemented in the package checkout on 10 November 2026. This extends the
estimation corrections described in `../fixes.md`. The manuscript and existing
Monte Carlo results were not changed or regenerated.

## Measured performance

Elapsed seconds per call, medians of repeated batches on the same Windows/R
installation. The baseline is the corrected package from the earlier audit,
before these performance changes, rather than the older Git HEAD.

| Calculation | Before | After | Before / after |
|---|---:|---:|---:|
| Sparse Kalman ML, 2 variables, 200 observations, VARMA(1,1) | 0.290 | 0.015 | 19.3 |
| Sparse Kalman ML, 3 variables, 400 observations, VARMA(2,2) | 1.280 | 0.300 | 4.3 |
| KFAS ML, 2 variables, 200 observations, VARMA(1,1) | 0.330 | 0.085 | 3.9 |
| KFAS ML, 3 variables, 400 observations, VARMA(2,2) | 1.750 | 0.940 | 1.9 |
| FKF ML, 2 variables, 200 observations, VARMA(1,1) | 0.240 | 0.045 | 5.3 |
| FKF ML, 3 variables, 400 observations, VARMA(2,2) | 1.500 | 0.820 | 1.8 |
| DPF, 5 variables, 800 observations, VARMA(2,2) | 0.0170 | 0.0085 | 2.0 |
| IHR, 5 variables, 800 observations, VARMA(2,2) | 0.038 | 0.032 | 1.2 |
| RBY, 3 variables, 400 observations, VARMA(2,2) | 0.0350 | 0.0275 | 1.3 |

The likelihood fits use the same optimizer budget (`control=list(maxit=10)`)
and pilot budget in both versions. These are timings for a fixed computational
budget, **not timings to a guaranteed converged optimum**. IHR uses `maxit=30`;
RBY uses `maxit=10`. Small IHR and RBY cases show no improvement beyond timing
variation, so there is no universal speedup claim.

For stationary autocovariances of a VARMA(2,2), through lag 10, the following
methods were timed in the same optimized installation:

| Variables | Dense Kronecker `autocov_vk` | Sparse Kronecker `autocov_sp` | New `autocov_ss` |
|---|---:|---:|---:|
| 2 | 0.000234 | 0.00275 | 0.000141 |
| 4 | 0.001000 | 0.00300 | 0.000172 |
| 8 | 0.042000 | 0.04000 | 0.000375 |

In the 8-variable example the new route is about **112 times faster** than
`autocov_vk` and 107 times faster than `autocov_sp`. These are moderately
persistent simulated models, not a general performance guarantee near unit
roots or for every matrix structure.

The new public Lyapunov solver also includes input and stability validation.
For state dimensions 24 and 48 it takes approximately 0.000175 and 0.000900
seconds, respectively, versus 0.000300 and 0.002100 seconds for the existing
low-level sparse-interface doubling routine in the same installation. At
dimensions 2 and 8 the public solver's validation overhead makes it slower
than that low-level routine. The old interface remains available.

### Reproducibility

`benchmark.R` fixes the simulation seed and saves every median, range, batch
size, and repetition count in `before.csv` and `after.csv`. Calls are warmed;
short calls are repeated in batches lasting at least approximately 0.1 seconds
to reduce Windows timer granularity. Medians use five batches, or three for
ML and RBY. Timing runs were sequential and separate from package compilation.
The session information is in `before-session.txt` and `after-session.txt`.

Run from the repository root, replacing the library path with an isolated
installation of the desired version:

```r
# Command line:
# Rscript audit/2026-11-10/performance/benchmark.R <library-path> after
```

The original baseline installation was
`C:/Users/Matteo/AppData/Local/Temp/varma-validation-20261110/varma.Rcheck`.
The optimized installation was
`C:/Users/Matteo/AppData/Local/Temp/varma-performance-20261110/varma.Rcheck`.
These temporary installations may eventually be removed by the operating
system; the CSVs, fitted reference objects, and session records are retained
here. The final source uses Armadillo's current `sylvester` name in place of its
deprecated forwarding alias `syl`; this does not change the calculation timed.

## What makes the code faster

- ML parameter unpacking, root checks, state-space construction, and stationary
  initial covariance now run in C++. The sparse Kalman objective remains in
  C++ until it returns the scalar likelihood. KFAS and FKF retain their own
  likelihood implementations. The stationary covariance is still recomputed
  for every candidate parameter vector.
- IHR reuses its fixed response and AR regressors, replacing only the MA
  regressors between iterations. Least-squares fits reuse QR residuals instead
  of multiplying the design matrix by the estimated coefficients again. NET
  obtains the pilot design directly instead of reconstructing it from QR.
- DPD/DPF assemble GLS equations from reusable cross products rather than
  constructing mostly zero, large regression matrices at every observation.
- RBY evaluates residuals and likelihood during its line search; derivative
  calculations are reserved for accepted candidates.
- `autocov_ss()` solves the stationary Lyapunov equation without constructing
  the squared-dimension Kronecker system. It propagates only the state-to-output
  covariance columns required at successive horizons.

## Echelon estimation

```r
library(varma)

# Rows of Y are dates; columns are series. Indices follow that column order.
known <- fit_varma_echelon(Y, kidx=c(2,1))

# Select indices using canonical-correlation tests, then estimate coefficients.
# This option requires the suggested MTS package.
selected <- fit_varma_echelon(Y, plag=5, crit=0.05)
selected$kronecker_indices
selected$convergence
selected$optim.out

irf(known, maxlag=20)
predict(known, n.ahead=12)
```

The coefficient estimator maximizes a **conditional Gaussian likelihood**,
with the innovation covariance concentrated out. It is not the stationary
state-space likelihood used by `fit_varma_kfas()`. Presample innovations are
zero through `max(kidx)` and the likelihood uses the subsequent observations.
A long-VAR pilot supplies starting values; BFGS uses a compiled analytic
gradient. Candidate models must have stable AR and invertible MA polynomials.
The optimizer reports numerical convergence, not a guarantee of a global
maximum. Increasing `control=list(maxit=500)` can help with difficult fits.

The raw echelon representation is

```text
Phi0 y[t] = c + sum_l A[l] y[t-l] + Phi0 e[t] + sum_l B[l] e[t-l].
```

`Phi0` is unit lower triangular. The function applies the restrictions implied
by the supplied indices without sorting them or reordering the series. The
returned object's usual `ar`, `ma`, and `intercept` components are premultiplied
by `solve(Phi0)` to retain the package's lag-zero identity convention. The
original restricted coefficients are in `fit$echelon`. Positive MA signs
follow the rest of this package. The parameter count includes the free
echelon coefficients and innovation covariance parameters.

The implementation supports zero indices, an all-zero index vector (white
noise), and univariate estimation with a supplied index. With `kidx=NULL`,
`MTS::Kronid()` selects indices from standardized data. Its test output and
settings are saved in `fit$selection`. Selection errors are reported explicitly;
the code does not silently substitute another specification. This selection
step introduces uncertainty beyond coefficient estimation conditional on the
indices. It is not a claim of post-selection inference validity.

The restrictions were verified against `MTS::Kronspec()`. The coefficient
optimizer is an independent implementation, not a wrapper around `Kronfit()`;
its starting values, parameter bounds, and initialization need not reproduce
that function's numerical estimates. See the
[MTS reference manual](https://cran.rstudio.com/web/packages/MTS/refman/MTS.html)
for the canonical-correlation index-selection procedure and its controls.

### Bootstrap refits

```r
# Condition on the fitted indices (default).
boot_fixed <- bootstrap_varma(known, nsim=100)

# Repeat index selection in every bootstrap sample.
boot_selected <- bootstrap_varma(selected, nsim=100, kidx=NULL)
```

Both use the existing joint residual-vector bootstrap. Supporting these refits
does not establish bootstrap validity at singular models or after selecting
indices. That remains a theoretical issue for the paper.

## Stationary covariance and autocovariances

```r
# Gamma[,,h+1] = Cov(y[t+h], y[t]).
Gamma <- autocov_ss(known, maxlag=20)

ss <- varma_to_ss(known)
P <- solve_discrete_lyapunov(ss$T, ss$R %*% ss$Q %*% t(ss$R))
```

`solve_discrete_lyapunov()` solves `P = T P t(T) + C`. Its `method` can be
`"doubling"`, `"schur"`, or `"auto"` (the default). Auto selects Schur at state
dimension 16 or larger and doubling otherwise, falling back to doubling if
the Schur calculation fails its residual check. The Schur route uses a
bilinear transformation and a Sylvester solver. It avoids the large Kronecker
linear system but still uses dense state matrices.

This is the **unconditional stationary covariance**, also used to initialize
the Kalman filter. The filter's limiting conditional prediction covariance is
a different quantity, governed by a Riccati equation. The new public routine
does not substitute the unconditional covariance for that conditional one.

## Numerical validation

- All three Lyapunov methods agree with independent Kronecker solves, including
  scalar roots close to both +1 and -1. Autocovariances agree with the existing
  vectorized method.
- Optimized GLS agrees with direct observation-by-observation normal-equation
  assembly, including pure MA specifications and diagonal/final MA forms.
- Compiled ML parameter packing and stationary covariance agree with independent
  R constructions; its likelihood agrees with an independent KFAS calculation.
- Earlier and optimized IHR/DPF fits use identical saved data. Across the three
  benchmark designs the largest coefficient change is 1.74e-13, covariance
  change 4.27e-14, and log-likelihood change 1.37e-11. Individual comparisons
  are saved in `equivalence.csv`.
- Echelon gradients agree with finite differences, including derivatives of
  nonidentity `Phi0`. Restriction patterns agree with MTS across index orders
  and zero indices. Tests cover nonzero means, white noise, univariate fits,
  automatic selection on dynamic data, forecasts, IRFs, and bootstrap refits.
- The regression suites are `tests/regression-estimation.R` and
  `tests/regression-echelon-performance.R`. The latter checks numerical
  correctness rather than imposing machine-dependent speed thresholds.

A complete runnable example is in `inst/scripts/echelon_example.R`; it also ran
successfully against the checked installation (`example-output.txt`). The final
`R CMD check` finished with **Status: OK**, with no errors, warnings, or notes.
Package-check and regression-test logs are retained alongside this report.
The check excludes PDF-manual and vignette building. Changes are in the source
checkout; the user's normal installed package library has not been replaced.
