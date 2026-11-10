# Estimation audit and corrections — 10 November 2026

The audit covered all eight `fit_varma_*` functions, the native estimation and
Kalman routines, and the residual, prediction, identification, and bootstrap
paths that consume their output. Production fixes are in `R/functions.R`,
`R/estimation_helpers.R`, and `src/rcpp_functions.cpp`. Help files were regenerated.
The pre-fix evidence, including source hashes, is in `../2026-11-09/results.txt`.

## Confirmed problems and their resolution

| Area | Observed problem | Correction |
|---|---|---|
| DPD/DPF | Adding a constant to data changed AR/MA estimates by approximately 0.4 despite requesting an intercept. | Centre data in both stages; retain the fitted mean through the model intercept. |
| RBY | The alleged gradient differed from a numerical SSR derivative by up to 17.2 in a two-series example. Filtering regressors first incorrectly commuted coefficient perturbations with the matrix MA filter. | Recurse the full innovation derivative, including cross-equation terms; use covariance-weighted normal equations and a conditional Gaussian likelihood line search. Verify the resulting profile-likelihood score by finite differences. |
| RBY mean | Sample-mean centring was described as an exact profile maximum-likelihood intercept. | Estimate the intercept jointly in the conditional criterion. |
| Sparse Kalman output | Returned covariance and transition matrices were exactly the pilot matrices, despite changed fitted AR coefficients. Aligning the transition alone changed a three-step forecast by 0.205. | Build every returned matrix from the final parameter vector. |
| Sparse likelihood | A stable VMA(2) with a delayed response gave log-likelihood -937.915 with acceleration and -897.199 without it. | Check convergence of the entire Riccati covariance; equality of the observed covariance block is insufficient. |
| All likelihood backends | Initial covariance was fixed at starting values; FKF initially used an unpopulated transition matrix. | Recompute stationary covariance for each candidate with a doubling Lyapunov solver. Use shared packing, stable/invertible candidates, positive covariance parametrisation, and explicit optimiser diagnostics. |
| KFAS/FKF residuals | Sequential component errors were stored as if they were joint vector innovations. In the KFAS example, the innovation cross covariance was about 0.471 but the stored residual cross covariance was about -0.004. | Use multivariate innovations for KFAS and observation-minus-predicted-state errors for FKF.SP. See [KFAS's multivariate-innovation documentation](https://search.r-project.org/CRAN/refmans/KFAS/html/mvInnovations.html). |
| NET penalty | Lambda values 1e-9 and 1e6 gave identical coefficients. | Honour fixed penalties and candidate paths, including relaxed fits. Accommodate a single-regressor input and glmnet 4.1-10's single-lambda relaxed-path dimension issue. |
| NET SUR | Returned residuals differed from those implied by the post-SUR regression by as much as 0.886. | Refresh regression errors, recursive innovations, covariance, and likelihood after SUR. Permit empty selected equations. |
| IHR output | With one iteration, stored regression errors were inconsistent with the returned model's recursive innovations. | Recompute innovations from final coefficients for forecasting and bootstrap; preserve regression errors only in regression output. |
| IHR pilot/termination | Default pilot order used a fixed fraction of the sample; unavailable pilot residuals entered the regression, and stopping depended on potentially unidentified coefficients. | Use a logarithmic pilot order, trim unavailable pilot lags, compare regression errors, handle aliased least-squares columns, and return convergence diagnostics. This is an algorithmic improvement, not a proof of asymptotic validity. |
| Likelihood and IC | Degrees-of-freedom covariance normalisation was inserted in a concentrated-ML formula; covariance parameters were omitted in some counts; `nobs` conventions differed; printed HQC omitted a factor of two. | Evaluate Gaussian likelihoods consistently; count covariance parameters; use vector observation counts and a proper `logLik` object; align conditional order comparisons on common observations. |
| Returned metadata | KFAS used `npars` instead of `npar`; FKF omitted `n`. | Standardise the fitted-object metadata. |
| Boundary cases | Single-series IHR and pure-MA RBY/DPD/DPF failed; some one-series state-space assignments dropped dimensions. | Preserve matrices, guard zero-order indexing, and return forecast covariance arrays at all horizons. |
| Generalized IRFs | The innovation covariance was scaled by rows instead of columns. | Use column scaling by the shock's marginal standard deviation. |
| Pearlman form | Pure VAR loading was zero unless a zero MA polynomial was explicitly supplied. | Include the AR disturbance loading independently of MA order. |
| Identification | Shared scalar roots were treated as sufficient for non-identification in a multivariate model. | Check the rank of the exact finite linear restrictions implied by equal transfer functions at the specified orders. |
| Kronecker indices | Square-Hankel rank increments and sorting could give incorrect indices for a specified series ordering or delayed responses. | Use sufficiently wide Hankel rows and retain the output-series ordering required by echelon estimation. |
| Bootstrap | Fitted means were omitted, errors were uncentred, no burn-in was provided, and zero residual skipping still dropped the first row. | Centre and jointly resample valid errors, restore the process mean, add configurable burn-in, and correct zero-skip handling. |
| Parallel execution | NET ignored parallel arguments; bootstrap reset the future plan to sequential instead of restoring it. | Implement the requested NET parallelisation and restore the previous plan on normal or exceptional exit. |

## Changes that affect existing analyses

The IHR default pilot order, trimmed pilot regression, and refreshed innovations
can change IHR, HR1, NET, and their forecasts. RBY now optimises the intended
covariance-weighted conditional criterion. The likelihood estimators now use
stationary initialisation at every parameter value and restrict optimisation to
stable, invertible representations. Results and timings should therefore be
regenerated before using the simulation figures in a submission.

`nobs()` now counts vector observations. Covariance matrices use likelihood
normalisation, and `npar` includes covariance parameters. These are intentional
changes to the fitted-object conventions. Nominal parameter counts do not solve
the separate question of information criteria at singular representations or of
effective degrees of freedom under regularisation.

DPD and DPF are accurately documented as one feasible GLS regression after a
long-VAR pilot. They are not presented as implementations of all iterations of
the cited efficient estimator. Likewise, the repaired residual bootstrap is an
implementation of the usual construction; its validity at a singular
distribution remains a theoretical question for the paper.

## Validation

`tests/regression-estimation.R` checks numerical scores, independent likelihood
backends, direct IRF recursions, shift invariance, coefficient/output agreement,
boundary orders, fixed penalties, bootstrap means, information criteria, and
future-plan restoration. The package is built and checked in an isolated
temporary installation. Full Monte Carlo experiments and the manuscript figures
are not regenerated as part of these code corrections.

Final validation on R 4.5.2, Windows, with the package installed from its source
archive: **R CMD check: OK (zero errors, warnings, or notes)**. The command used
`--no-manual --no-vignettes --no-build-vignettes`; PDF-manual and vignette building
were therefore outside this check. All eleven regression groups passed, including
parallel-plan restoration. Additional installed-package checks passed for the
default relaxed NET cross-validation and a supplied penalty path. The full check
log and regression output are saved alongside this report.
