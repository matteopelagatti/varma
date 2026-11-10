# Paper simulation runner

The default experiment reproduces the four full-precision models from
`inst/scripts/sim2.R`, with sample sizes 100, 200 and 400, 1,000 replications
per model/sample-size cell, zero means, no fitted intercepts, 100 discarded
burn-in observations, IRF lags 1–5 and forecast horizons 1–5.

It compares eight estimators: `ECH_oracle`, `ECH_selected`, `IHR`, `HR1`, `NET`,
`RBY`, `DPF`, and `MLE`. MLE defaults to `fit_varma_skf()`, which uses the
compiled stationary Gaussian likelihood. DPF retains the expanded orders in
the original experiment. The two echelon estimators use the package's
conditional Gaussian likelihood, with population indices and indices selected
from the training data, respectively. No MTS installation is required.

The additional `DGP` benchmark knows the population parameters and filters the
same observed training history. It does not see future observations or the
simulation innovations. All methods receive the same sample in a replication.

## Install on another computer

Install R first. For the supplied Windows binary, use **64-bit R 4.5.x**. In R:

```r
install.packages(c("Rcpp", "RcppArmadillo", "KFAS", "FKF.SP", "Matrix",
                   "future", "future.apply", "glmnet"))
install.packages("varma_0.1.0.zip", repos=NULL, type="win.binary")
library(varma)
```

The source archive also works on other platforms or when rebuilding for a
different R version:

```r
install.packages("varma_0.1.0.tar.gz", repos=NULL, type="source")
```

A source installation needs the appropriate compiler toolchain: Rtools on
Windows, or the normal R/C++ build tools on Linux/macOS. The source and Windows
binary in the dated replication bundle contain the same changes. Installing
dependencies needs internet access; running the simulation does not.

## Benchmark, then run

Open a terminal in the folder containing `run_paper_simulations.R`:

```text
Rscript run_paper_simulations.R benchmark timing-test 4 5
Rscript run_paper_simulations.R full paper-results 4 1000
```

The arguments after the output directory are the number of worker processes
and replications per cell. The benchmark uses the same methods and optimizer
controls as the full run, rather than a reduced iteration budget. Five pilot
replications per cell give a quick timing estimate; 20 or more are preferable
for assessing variability. Use a new output folder for each benchmark.

The full experiment performs **96,000 fits** over **12,000 independently
simulated samples**, plus the DGP forecast calculations. On a four-core laptop,
four workers usually provide a useful speed gain, but eight logical processors
do not imply an eightfold speedup. The driver limits BLAS/OpenMP threads to one
per process to avoid nested parallelism. Leave the computer connected to power
and allow for thermal throttling and other active applications.

Alternatively, use R directly:

```r
library(varma)
timing <- benchmark_varma_paper_simulation("timing-test", pilot_nsim=5, workers=4)
result <- run_varma_paper_simulation("paper-results", workers=4)
```

When using R directly, set the thread-limit environment variables before
loading numerical libraries, as the command-line driver does. The functions
default to one worker and restore the caller's RNG state and `future` plan.

## Interruptions and reproducibility

Every completed group of 10 replications is saved as a compressed checkpoint.
After an interruption, repeat the same `full` command and output directory.
Completed chunks are validated and reused; at most an unfinished chunk per
worker is repeated. You may change the number of workers on resume. Other
settings must match the saved manifest. Run only one driver per output folder.

Independent L'Ecuyer streams are assigned to design/sample-size cells and
substreams to replications. Results are reproducible across worker counts and
resumes with the same R/package versions. NET's randomized cross-validation
uses the replication's stream. The manifest retains all population matrices,
controls, the seed, package version, and session information. Each checkpoint
also retains its replication seeds. Data can therefore be regenerated without
storing every training sample.

`keep_fits=TRUE` stores fitted model objects as well, but substantially increases
disk use. The default retains raw IRF errors, forecast errors, forecast
differences from DGP, realized outcomes, and per-fit diagnostics. These are
enough to recompute the reported accuracy measures and construct plots.

## Outputs and interpretation

- `irf-summary.csv`: untrimmed MSE and RMSE of reduced-form IRFs, overall and by
  lag. `horizon=0` pools positive lags 1–5; it does not score the identity response.
- `forecast-summary.csv`: held-out forecast MSE/RMSE, by horizon and variable.
  `variable=0` pools variables within each replication.
- `forecast_gap-summary.csv`: squared distance between estimated and DGP
  conditional forecasts, separating estimation error from future innovations.
- `failures-summary.csv`: attempted fits, fit/IRF/forecast failures, convergence
  flags, instability, noninvertibility, and warnings. Failed replications stay
  in these denominators.
- `timing-summary.csv`: elapsed estimation times by method and simulation cell,
  including failed attempts, with separate evaluation times.
- `indices-summary.csv`: frequencies of returned Kronecker index vectors.
  Indices unavailable after a failed fit are labeled explicitly; these counts
  must not be interpreted as a pure index-selection error rate.
- `diagnostics-summary.csv`: replication-level diagnostics and error messages.
- `manifest.rds`, `chunk-*.rds`, and `execution-*.rds`: configuration, raw errors,
  seeds, and execution/session records.

Rows labeled `subset="available"` use all finite results, including finite
outputs from nonconverged or unstable estimates. Rows labeled `"qualified"`
add stability, invertibility, and reported-convergence requirements. HR1 is
deliberately one-step and has no iterative convergence requirement; NET's
convergence field concerns its IHR pilot. Neither subset conceals failures:
consult the failure table before interpreting rankings.

Primary RMSE is `sqrt(mean(squared_error))`, with no truncation of large errors.
The additional `trimmed_rmse` column reproduces a symmetric 1% trimming rule
for comparison with the old figures; it is not the primary estimate of risk.
Monte Carlo standard errors use independent replication-level losses, pooling
correlated variables or IRF entries within a replication first. RMSE standard
errors use the delta method and may be unreliable with heavy-tailed losses.

Relative MSEs compare methods on **common usable replications**. IRF references
are oracle echelon; forecast references are the DGP. `n_paired`, `paired_mse`,
and `reference_mse` expose the comparison's sample and denominator. These are
conditional-on-availability summaries, not unconditional risks for algorithms
that sometimes fail. The forecast-gap metric has no automatic risk ratio.

To refresh summaries after a partial run without estimating anything:

```text
Rscript run_paper_simulations.R summary paper-results
```

## Changing the experiment

Use a new directory when changing scientific settings:

```r
run_varma_paper_simulation(
  "larger-samples", sample_sizes=c(100,200,400,800), nsim=1000, workers=4,
  estimator_args=list(
    MLE=list(control=list(maxit=1000,reltol=1e-9)),
    ECH_selected=list(plag=6,crit=0.05)
  )
)
```

Default MLE and echelon limits are 500 optimizer iterations with `reltol=1e-8`.
Other estimators keep their corrected package defaults. Local optimizer
convergence does not certify a global maximum. No failed method is replaced
by a different estimator and no automatic retries alter its definition.

These experiments address point estimation of IRFs and forecasts. They do not
run bootstrap coverage or FEVD experiments, or establish a VARMA-versus-VAR
comparison; those remain separate additions to the paper. The historical
simulation outputs and manuscript are not overwritten.
