# varma 0.1.1: regression-innovation correction

This release corrects the IHR/HR1/NET forecast postprocessing introduced in
0.1.0. The estimator's final regression innovation proxies are retained;
they are no longer replaced by a zero-initialized inverse-MA recursion.
That recursion can explode for a noninvertible estimated MA polynomial.
AR/MA estimation, pilot orders, regression windows and stopping rules are
unchanged by this patch.

For 64-bit Windows R 4.5.x, restart R and install the binary:

```r
install.packages("varma_0.1.1.zip", repos=NULL, type="win.binary")
library(varma)
packageVersion("varma")
```

For other platforms or R versions, build the supplied source archive using
the appropriate compiler toolchain:

```r
install.packages("varma_0.1.1.tar.gz", repos=NULL, type="source")
```

`residuals(fit)` now returns the estimator's stored residuals by default.
`residuals(fit, recompute=TRUE)` explicitly requests the zero-initialized
inverse filter. Supplying new `Y` always computes residuals for those data.
IHR/NET's Gaussian residual criterion is labeled `regression_approximation`
when the model has an MA part. Stability and convergence flags are retained;
the patch does not impose stability constraints or remove extreme fits.

## Saved experiments

Installing this package does not alter existing checkpoints or summaries.
Use a new results directory to refresh a complete experiment, for example
from the repository root:

```text
Rscript inst/scripts/run_paper_simulations.R full paper-results-v011 7 1000
Rscript inst/scripts/report_paper_simulations.R paper-results-v011
```

The scientific settings and seed are the same as the previous run. The
version field prevents mixing 0.1.0 and 0.1.1 checkpoints in a resumed run.
The reporting script marks archived 0.1.0 forecasts as affected by the bug.

The correction was tested with the full R package check and replayed on all
12,000 saved IHR and 12,000 HR1 samples, plus 144 NET audit samples. The NET
audit contains ten ordinary draws per cell and two extreme archived draws;
it is not a replacement for a full NET simulation.
