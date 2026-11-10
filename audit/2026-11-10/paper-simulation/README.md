# Paper simulation implementation and timing

Date: 10 November 2026. The main entry point is
`run_varma_paper_simulation()`, implemented in `R/paper_simulations.R`.
The portable command-line driver and installation guide are
`inst/scripts/run_paper_simulations.R` and `inst/scripts/PAPER_SIMULATIONS.md`.

## Experiment

The default reproduces the original four population models, sample sizes
100/200/400, 1,000 replications per cell, zero means, no estimated intercepts,
burn-in 100, positive IRF lags 1–5 and forecast horizons 1–5. Full-precision
population coefficients are regenerated using the original seeds, not copied
from rounded slide matrices. The generator preserves the caller's RNG state.

There are eight estimated procedures: oracle echelon, data-selected echelon,
IHR, one-stage HR, NET, RBY, DPF, and stationary Gaussian MLE. DPF uses the same
expanded orders as the old driver. A ninth, non-estimated DGP benchmark computes
the true-parameter conditional forecast using only the observed training sample.
Thus the complete experiment contains **12,000 samples and 96,000 estimated
models**, plus 12,000 true-model forecast calculations.

The revised driver uses `fit_varma_skf()` for MLE. It was the fastest backend
overall in the earlier package benchmarks and in the refreshed two-design
comparison. On the more expensive Model 4, n=200, the measured SKF/FKF/KFAS
times were 2.23/6.38/8.50 seconds. This does not imply it is fastest for every
individual fit: on Model 1 the first SKF call included startup overhead and
the numerical optimizers reached slightly different likelihoods. The scalar
likelihood itself has separate independent-backend correctness tests.
See `compare-mle.R` and `mle-backends.csv`.

MLE and both echelon optimizers use `maxit=500`, `reltol=1e-8`. Other procedures
retain their corrected package defaults, including IHR's 100 iterations and
RBY's 500; HR1 deliberately performs one iteration. Timing estimates below
use these actual settings, not a short fixed optimizer budget.

## Accuracy, diagnostics, and reproducibility

The output includes untrimmed IRF and forecast MSE/RMSE, replication-based Monte
Carlo standard errors, an explicitly labeled additional trimmed IRF/forecast
summary, and forecast distance from the true-model conditional forecast.
Comparative MSE ratios use common available replications and show their paired
sample counts and denominators. Failure, convergence, stability, and warning
tables accompany all accuracy summaries. Finite nonconverged estimates remain
in the primary available-output summary; the separate qualified subset applies
the declared convergence and stability requirements.

No fitted method sees held-out observations. No errors are truncated on their
magnitude, and no failed method is replaced by another estimator. The small
timing pilot is not evidence for a ranking of statistical accuracy. In that
pilot, IHR often reached its iteration limit; NET inherits that pilot flag.
These cases are recorded rather than silently treated as converged estimates.

Compressed checkpoints retain raw losses, actual future observations, DGP
forecasts, seeds, and per-fit diagnostics. Saving entire fitted objects is
optional. Cell-specific L'Ecuyer streams and replication substreams make
results independent of worker count and checkpoint resumption. Only completed
chunks are published. The saved configuration prevents resuming with different
scientific settings. Software/session information is retained with the run.

## Sequential timing

Host: **Intel Core i7-8550U, four physical cores, eight logical processors,
approximately 32 GB RAM**, Windows 11, R 4.5.2. BLAS/OpenMP thread limits were
set to one. Benchmarking was separate from package compilation and checks.

The sequential pilot used **10 replications in every one of the 12 cells**:
120 samples, 960 estimated models, and 120 DGP forecasts. Its computational
wall time was **207.76 seconds**; including summary generation, 220.71 seconds.

| Model | Full-run time for its three sample sizes, one worker |
|---|---:|
| 1: nilpotent, 2 variables | 35.9 minutes |
| 2: identified, 2 variables | 34.1 minutes |
| 3: common factor, 2 variables | 95.9 minutes |
| 4: identified, 3 variables | 179.6 minutes |
| Total | **5.76 hours** |

The projection sums the mean elapsed replication time in each cell, then
multiplies by 1,000. Its pilot Monte Carlo standard error is about **0.25 hour**.
That standard error does not include hardware differences, changes in power
settings, thermal throttling, competing applications, or rare slow optimizer
paths missed by a small pilot. A practical single-worker planning allowance is
**about 6–8 hours** on a similar machine. Summary writing adds little relative
to these hours of estimation.

NET and MLE account for most of the measured work: approximately 97.4 and
82.2 seconds in the pilot, respectively. Both echelon methods together took
6.6 seconds. Data-selected echelon can use smaller models than oracle echelon;
its lower runtime must not be interpreted as higher statistical accuracy.

The raw timing evidence is in `pilot-sequential/runtime-projection.csv`,
`runtime-projection.rds`, `timing-summary.csv`, and the checkpoint files.
`run-pilot.R` reproduces the measurement. The statistical point estimates from
these ten-replication cells should not be used as paper results.

## Four-worker timing

The same 120 samples were then processed with four independent R workers,
using the same settings and two-replication checkpoint chunks. Computation
including checkpoint writes and worker startup took **103.30 seconds**;
including summary generation and cleanup, the pilot took 117.41 seconds.
Extrapolating the measured throughput gives **2.87 hours** for the complete
experiment. A practical allowance is **about 3–4 hours** on this machine class.

| Execution | Pilot computational wall time | Projected full run | Planning allowance |
|---|---:|---:|---:|
| One worker | 207.76 seconds | 5.77 hours | 6–8 hours |
| Four workers | 103.30 seconds | 2.87 hours | 3–4 hours |

The measured gain is about twofold, not fourfold. Per-fit elapsed times rise
under concurrent CPU and memory use, and a mobile processor can reduce clock
speeds under sustained multicore load. Parallel projections therefore use
measured throughput rather than dividing the single-worker estimate by four.
The parallel pilot's sum of individual task durations is not a serial-runtime
estimate. Small pilots also over-amortize one-time startup when extrapolated;
conversely, they may miss rare slow optimizer paths. Rebenchmark on the target
computer with `benchmark_varma_paper_simulation()`.

`pilot-parallel` contains the raw records, and `run-parallel-pilot.R` reproduces
the measurement. All IRF, forecast, and forecast-gap summary values agree
exactly between the sequential and four-worker pilots, as recorded in
`parallel-agreement.txt`.

## Validation and delivery

The final `R CMD check` finished with **Status: OK**. All four regression suites
passed, with parallel testing enabled; PDF-manual and vignette building were
excluded. `R-CMD-check.log` and the `.Rout` files retain the results.

The new tests independently verify the DGP forecast using a Gaussian
block-Toeplitz covariance conditioning calculation. They also verify the
population designs, untrimmed losses and paired denominators, saved forecast
errors, absence of held-out observations in fitted models, failure accounting,
resume behavior, configuration mismatch rejection, and identical numerical
results with sequential and two-worker execution.

All four population models also agree with the independent full-precision
model export retained alongside the manuscript, with a maximum difference of
4e-15 from its JSON rounding (`design-verification.txt`). The Windows binary
was installed into a fresh temporary library, passed its MD5 check, and ran
both internal echelon selection and a simulation smoke test
(`binary-validation.txt`). The user's normal installed library was not replaced.

The source archive and a Windows R 4.5.x binary are prepared for transfer to
the other computer. They contain the updated package, not just the driver.
The manuscript and historical simulation outputs are unchanged; the full
1,000-replication experiment has not been run here.
