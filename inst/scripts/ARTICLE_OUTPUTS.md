# Article figures and tables, entirely in R

`article_paper_simulations.R` creates a compact set of publication figures and
LaTeX tables. It reads the raw `manifest.rds` and `chunk-*.rds` files, using the
validated collector in `report_paper_simulations.R`. Keep both scripts together.
The existing general reporting script is also written in R.

No models are estimated again. Python and a LaTeX installation are **not**
required to export the figures and table fragments. LaTeX is only needed to
compile the optional preview document.

## Run from the repository root

For an editable RStudio driver, open `inst/scripts/make_article_outputs.R`.
Run `source("inst/scripts/make_article_outputs.R")` to produce all three versions:

- `paper-results/article/full/`: all finite usable outputs.
- `paper-results/article/qualified/`: returned fits with stable AR and invertible
  MA, including nonconverged fits. Unknown stability/invertibility cannot qualify.
- `paper-results/article/trimmed/`: full usable outputs with 1% of replication
  losses removed from **each** tail, separately for every reported accuracy measure.

Each folder contains plots, LaTeX tables, CSV statistics and `article_figures.pdf`.
The R objects are `articles$full`, `articles$qualified` and `articles$trimmed`.
All three reuse the same validated numerical cache in the parent folder.
Previously exported files directly in that parent folder are not refreshed;
use these three subfolders. No reinstallation of varma or new simulations is needed.

Qualification is recomputed from saved diagnostics, ignoring the older saved
qualification flag. Nonconvergence remains a separate diagnostic. Raw checkpoints
and the simulation's original root-level summary CSVs are unchanged.

Trimming ranks whole replication losses after pooling the entries/horizons
relevant to each statistic. It removes `floor(0.01 * B)` at each end, so fewer
than 100 usable replications imply no trimming. Ties use replication order.
The retained draws may differ across methods, horizons and targets. For relative
RMSE, benchmark losses use the estimator's retained draws; the benchmark is not
trimmed independently. Thus the denominator can differ from the separately
reported DGP RMSE. CSV columns `n_eligible`, `n_trimmed`, `n_available` and
`n_paired` document the denominators. Ordinary delta-method MCSEs are suppressed
for trimmed accuracy statistics. Timing and diagnostic tables still use all attempts.

The paths, estimators, sample size, figure dimensions and font size are at the
top. Source that file to regenerate everything. The default four-panel figures
are 7 x 6 inches with 14-point titles and 13-point tick and legend labels.
Print near the exported width: shrinking a 7-inch figure to 3.5 inches halves
every font. The driver includes a separate single-design, single-column example.

The first run validates and aggregates the raw checkpoints. Subsequent runs
reuse `summary-cache.rds` only if the input content hashes, trimming, collector
and package version match. Style edits do not require repeated aggregation.
`plots.rds` retains editable ggplot objects for use in a new R session.

```r
# Once, if these plotting packages are not installed:
install.packages(c("ggplot2", "scales"))

source("inst/scripts/article_paper_simulations.R")
article <- article_varma_paper_simulation("paper-results")
```

The `varma` package must also be installed. Output goes to
`paper-results/article/`. Existing simulation checkpoints and the general
`paper-results/report/` directory are not changed.

From a terminal:

```text
Rscript inst/scripts/article_paper_simulations.R paper-results
```

In RStudio, `source()` defines the function; the following function call runs
the export. Relative paths are resolved from the current working directory.
For example, first use `setwd("C:/GitHub/varma")` on this PC, or use absolute
paths for both the script and the results directory.

## VAR benchmark (version 0.1.4)

`VAR_BIC` appears as "VAR (BIC)" in figures and tables. It selects a pure VAR
by common-sample BIC, then refits on all observations available at the selected
order. It uses no held-out outcomes or population VARMA orders.

The exporter adds `tables/var_lag_selection.tex` and
`data/var_lag_selection.csv` / `data/var_lag_frequencies.csv`.
The selected VAR lag is an approximation order; it is not judged against the
AR order of the VARMA DGP. The boundary percentage uses returned fits, while
frequency CSVs use all attempts and explicitly include failed fits.

## CML and initialization settings (version 0.1.3)

The new experiment includes `CML` in all figures and tables. CML uses the
same conditional likelihood and analytic gradient as echelon, without its
coefficient restrictions. Its curve should coincide with oracle echelon for
Models 2 and 4; overlapping curves are expected. All methods are included
by default, including CML. Legends wrap to at most four columns and support
all ten estimators plus DGP without dropping plotting symbols.

The manifest records the state-space MLE's pilot and initial covariance
policy. Keep these settings fixed within an experiment. Numerical results
from earlier package versions remain readable but do not contain CML.

## Important: the archived version 0.1.0 run

The obsolete local `paper-results` outputs were deleted at the user's request;
the same folder is now reserved for the corrected experiment. The safeguards
below still apply if an older 0.1.0 experiment is supplied to the exporter.

That run contains the confirmed IHR/HR1/NET forecast-residual bug. The article
exporter withholds those methods' forecast accuracy and forecast-gap statistics,
leaving `--` in their table cells and documenting the omission in captions.
It retains their **reduced-form IRF** statistics, estimation times and recorded
diagnostics. Omission is based on the simulation manifest's version, not on the
currently installed version of `varma`.

The exporter does not repair old forecasts, silently splice in the audit
replays, or change the estimation algorithm. For a new run made with a corrected
package, all saved methods are included automatically. Use the new run's folder
as the input. The archived diagnostic counts still describe the original run,
including its forecast failures.

## Main outputs

- `figures/irf_rmse.pdf`: ordinary untrimmed IRF RMSE versus sample size.
- `figures/irf_rmse_trimmed.pdf`: the separate trimmed IRF statistic used for
  comparison with the original paper. Default: 1% removed from **each** tail.
- `figures/forecast_rmse.pdf`: forecast RMSE versus horizon at the largest
  simulated sample size, including the DGP benchmark.
- `figures/forecast_relative_rmse.pdf`: matched-sample RMSE relative to DGP.
- `figures/forecast_gap_rmse.pdf`: distance from the DGP conditional forecast,
  pooled across forecast horizons, versus sample size.
- `figures/fit_time.pdf`: mean elapsed seconds per attempted fit.
- `article_figures.pdf`: the same figures together, with complete captions.
- `tables/*.tex`: RMSE and Monte Carlo standard errors, trimmed IRF statistics,
  matched forecast ratios, forecast gaps, timings, usable counts, diagnostics
  for each design, and data-driven Kronecker-index selection rates.
- `figures.tex`, `tables.tex`: LaTeX collections with captions and unique labels.
- `article_preview.tex`: standalone A4 preview, including all figures and tables.
- `data/*.csv`: unrounded summaries, including every sample size, forecast
  horizon and variable. Withheld forecast values are blank and explicitly flagged.
- `article-metadata.rds`, `session-info.txt`, `data/input_checksums.csv`:
  provenance, settings, execution records and hashes of the input checkpoints.

The individual figures are vector PDFs, with optional 300-dpi PNG copies.
Panels have separate vertical scales. All IRF-RMSE plots use linear axes.
Forecast and timing plots use log axes to retain very large finite values;
exact zeros cannot be shown on those axes but remain
in tables and CSVs. Lines connect observed points; no root-n rate is imposed.

## Customize

```r
article <- article_varma_paper_simulation(
  "paper-results",
  output_dir = "paper-results/article-selected",
  methods = c("CML", "MLE", "ECH_oracle", "ECH_selected"),
  forecast_n = 200,       # choose any sample size in this experiment
  subset = "available", # default; finite unstable/nonconverged results retained
  trim = 0.01,
  width = 7, height = 6, font_size = 14,
  png = TRUE
)

# The result contains ordinary editable ggplot objects:
article$plots$irf_rmse
ggplot2::ggsave("my-irf-figure.pdf", article$plots$irf_rmse,
               width = 7, height = 6)
```

Colours, shapes, line types, font and panel layout can be edited in `.pa_plot()`
near the top of the script. Summary tables are returned in `article$summaries`.
For example, filter `article$summaries$forecast` on `horizon` and `variable` to
prepare additional component-specific plots.

Command-line equivalents:

```text
Rscript inst/scripts/article_paper_simulations.R paper-results paper-results/article-selected --methods=IHR,MLE,ECH_oracle,ECH_selected --forecast-n=200
Rscript inst/scripts/article_paper_simulations.R paper-results paper-results/article-qualified --subset=qualified --no-png
```

`--allow-partial` reports an unfinished run with explicit incomplete-run labels.
Without this option, incomplete or incompatible checkpoints stop the export.
Choose a separate output directory for alternative settings. Rerunning replaces
the generated files in that output directory.

## Put the output into the article

Copy the generated `article/` directory into the Overleaf project under, for
example, `simulation-output/`. In the LaTeX preamble use:

```tex
\usepackage{graphicx,booktabs}
\newcommand{\varmaresults}{simulation-output}
```

Then insert either collection:

```tex
\input{simulation-output/figures.tex}
\input{simulation-output/tables.tex}
```

For a single figure, copy its complete environment from `figures.tex`; for a
single table use, for example:

```tex
\input{simulation-output/tables/irf_rmse.tex}
```

To check everything locally, run `pdflatex article_preview.tex` from inside the
generated directory. No automatic modification of the manuscript is performed.

## Statistical definitions

IRF loss is the mean squared error over all matrix entries and positive lags
within a replication. Forecast loss similarly pools the indicated horizons
and variables within a replication. RMSE is the square root of the mean of
these replication losses, **not** the mean of replication RMSEs. The independent
unit for Monte Carlo standard errors is a replication.

Forecast ratios use the same available replications for each method and its
DGP benchmark. They are ratios of aggregate RMSEs, not averages of ratios.
`forecast_gap` measures the difference between estimated and DGP conditional
forecasts, separating that discrepancy from future shocks.

The default `available` subset includes every usable result, even when an
estimate is unstable or does not converge. `qualified` additionally requires
stability and invertibility, irrespective of optimizer convergence. Neither is an
unconditional risk estimate when outputs are missing. Usable counts and all-
attempt diagnostic counts accompany the accuracy tables. Timing and index-
selection summaries always use all attempts, irrespective of the accuracy
subset. Delta-method MCSEs can be unreliable with very heavy-tailed losses.

These checkpoints do not supply FEVD errors or bootstrap coverage. The exporter
does not infer them from the saved IRF/forecast data. VAR benchmarks are included
when VAR_BIC is present in the simulation manifest.
