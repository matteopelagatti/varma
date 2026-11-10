# Internal Kronecker-index selection

The package now selects Kronecker indices without MTS. The implementation is
in `R/echelon.R` and `src/rcpp_functions.cpp`; the exported entry point is
`select_kronecker_indices()`. `fit_varma_echelon()` and bootstrap reselection
use it automatically. MTS has also been removed from `DESCRIPTION`, executable
examples, and the simulation script's echelon wrapper.

```r
selection <- select_kronecker_indices(Y, plag=5, crit=0.05)
selection$index
selection$details

fit <- fit_varma_echelon(Y)             # select and estimate internally
fit <- fit_varma_echelon(Y, kidx=c(2,1)) # supply indices as before
bootstrap_varma(fit, nsim=100, kidx=NULL) # repeat selection in each replicate
```

The ordinary bootstrap still conditions on the fitted indices unless
`kidx=NULL` is explicitly passed. Supplying indices bypasses selection entirely.
`select_kronecker_indices()` accepts a matrix or a single numeric series;
its `index` component follows the original column order. It retains the
three-row `tests` matrix and adds a labeled `details` data frame containing
each horizon, component, aligned sample size, canonical correlation, serial
correction, test statistic, degrees of freedom, and p-value. It prints only
when `verbose=TRUE` (apart from R's normal automatic display of a returned
object at the console).

## Algorithm and implementation

The method tests the smallest canonical correlation between a fixed-length
past vector and an expanding set of future components. A significant candidate
is retained in the future vector. When its null hypothesis is not rejected,
the current horizon becomes that series' index and the candidate column is
removed. Earlier retained columns remain available at subsequent horizons.

The C++ implementation:

1. Centers and rescales the observations, and builds the lagged past matrix once.
2. Factors the aligned past matrix once per future horizon, instead of once
   for each component tested at that horizon.
3. Appends future columns to an orthonormal basis using two-pass modified
   Gram-Schmidt. Rejected columns do not require rebuilding the existing basis.
4. Computes canonical correlations from an SVD of the small cross-product of
   the two orthonormal bases.
5. Obtains canonical score series directly from the bases and singular vectors,
   avoiding coefficient recovery and its additional linear solves.
6. Computes the required score autocorrelations directly, avoiding repeated
   construction of R `acf` objects and console formatting.

Both bases are rebuilt when moving to a new horizon because alignment removes
one observation and changes the centering. The optimization does not reuse an
inappropriate factorization across different samples.

The serial-correlation correction and chi-square test retain the MTS
definition. At horizon `h > 0`, the squared correlation is divided by
`1 + 2 sum_{j=1}^h acf_x(j) acf_y(j)`. The multiplier is
`n - 1 - (past_dimension + future_dimension - 1)/2`, using the original sample
size `n`; the degrees of freedom are `past_dimension - future_dimension + 1`.
This preserves the established approximation rather than silently changing
its finite-sample calibration.

The implementation is newly written C++ and R code, with no MTS source code
vendored. The statistical reference is Tsay (2014), Chapter 4. Compatibility
was checked against MTS 1.2.1, whose published licence is Artistic License 2.0:
[MTS package metadata](https://raw.githubusercontent.com/cran/MTS/master/DESCRIPTION).
See also the
[official MTS reference manual](https://stat.ethz.ch/CRAN/web/packages/MTS/refman/MTS.html)
and [base R canonical-correlation documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cancor.html).

## Deliberate differences from MTS 1.2.1

- **Zero indices:** MTS leaves the final component of the horizon-zero future
  vector in place if that component's index has just been set to zero and
  other indices remain unresolved. The internal version removes every rejected
  candidate, including this one. On the saved `(1,0)` example, the subsequent
  test has 9 degrees of freedom; MTS uses 8 because it retains the extra column.
- **Single series:** matrix dimensions are maintained throughout. The internal
  version selects index 1 on the saved univariate ARMA example where MTS stops
  with `incorrect number of dimensions`.
- **Unsupported tests:** constant or numerically rank-deficient data, insufficient
  observations, exhaustion of the past dimension, and an invalid serial
  correction produce explicit errors. They do not silently become an estimated
  index. `plag` must be a positive integer; this interface does not use MTS's
  undocumented nonpositive-`plag` convention.
- **Diagnostics:** selection is quiet by default and returns labeled test
  details. Small chi-square tail probabilities use the direct upper tail.

The default rank tolerance is `1e-7`. Selection retains the sample-size guard
from the earlier echelon wrapper, `n - plag > 2*m*plag + 1`, and also checks
the aligned sample available for each test. A larger past dimension requires
more data. Series ordering remains relevant to the echelon representation.

These are approximate sequential tests at the selected significance level.
This change does not prove consistent index selection, robustness to arbitrary
conditional heteroskedasticity, or post-selection bootstrap validity.

## Agreement and correctness

`validate-reference.R` compares six saved regular designs with MTS 1.2.1:
white noise with 2 and 3 series; dynamic VARMA models with 2, 3, and 5 series;
and a `(0,1)` mixed-index model. All selected indices agree exactly. Across all
their test entries, the largest absolute difference is **5.23e-12**.
`agreement.csv` records the results.

The package tests use saved numerical outputs, so they do not need MTS.
`tests/fixtures/kronid-reference.rds` records the reference package/version and
contains the observations and numerical reference results, not MTS source.
The developer-only script here can regenerate the fixture when MTS is available.

`tests/regression-kronid.R` also compares the optimized routine with an
independent, deliberately simple implementation using base R `cancor()` and
`acf()` for every test. It covers multiple horizons, single-series estimation,
zero indices in both positions, unit/location changes, extreme finite scaling,
preserved series names and RNG state, and explicit failure cases. Existing
estimation, echelon, forecast, and bootstrap regressions also pass in source mode.

The installed-package check finished with **Status: OK**, with no errors,
warnings, or notes. All three regression suites and package examples passed;
manual and vignette building were excluded. `R-CMD-check.log` and the three
`.Rout` files retain the results. The complete echelon example also ran against
that installation, with assertions before and after confirming that MTS was
never loaded (`example-without-mts.txt`).

## Timings

Median elapsed milliseconds per selection, five warmed batches on the same
Windows/R 4.5.2 machine, with `plag=5`:

| Series | Observations | MTS | Internal | Speed ratio |
|---|---:|---:|---:|---:|
| 2 | 300 | 2.969 | 0.254 | 11.7 |
| 3 | 1,000 | 12.500 | 1.328 | 9.4 |
| 5 | 2,000 | 33.750 | 6.250 | 5.4 |
| 8 | 1,500 | 70.000 | 10.625 | 6.6 |

These timings concern **index selection**, not the subsequent coefficient
optimization. The datasets are stable diagonal VARMA(1,1) models and all select
index 1 in every component. They do not imply the same speed ratio on every
model or sample size. MTS has no quiet argument, so its output is captured as
it was in the former wrapper. Its formatting overhead remains in the timing;
the gains include removing that overhead as well as numerical optimizations.

`benchmark.R` reproduces these measurements; `benchmark.csv` includes ranges
and batch sizes. Short calls are batched to last at least approximately 0.15
seconds to reduce Windows elapsed-timer granularity. Benchmarks ran separately
from package compilation and checks.

The simulation wrapper in `inst/scripts/sim2.R` now uses the internal echelon
coefficient estimator rather than `MTS::Kronfit()`. This changes which estimator
future simulation runs use; old simulation outputs were not regenerated. The
previous forecast-error correction in that script is preserved. The user's
normal installed library and the manuscript were not changed.
