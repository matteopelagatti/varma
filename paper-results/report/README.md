# VARMA simulation report

Complete run: 12,000 saved replications; 120,000 attempted estimates; 4 designs.
Sample sizes: 100, 200, 400. Planned replications per cell: 1000. Seed: 20261110.
Recorded compute/checkpoint elapsed time: 1.41 hours; worker counts: 7.
MLE backend: skf; initial covariance: stationary; pilot: HR1. CML and the echelon estimators use conditional Gaussian likelihood.

Primary results retain every finite loss, including extreme outputs and nonconverged estimates.
Qualified results require stable AR and invertible MA, irrespective of optimizer convergence.
NET's convergence flag concerns its IHR pilot; HR1, DPF and VAR (BIC) have no iterative flag.
The paper-style trimmed comparison removes 1.0% of replication losses from each tail.
Neither availability nor qualification makes these unconditional risks of an algorithm that can fail.

IRF bands are Monte Carlo sampling quantiles, not confidence intervals.
MCSEs use independent replication-level losses. Heavy tails can make them unreliable.
Relative errors pair each method with its benchmark on common usable draws.

This run evaluates IRFs and point forecasts. It does not evaluate FEVDs, bootstrap coverage,
or prediction-interval coverage. These require additional experiments.

## Files

- `paper_report.pdf`: overview and the main figures.
- `figures/`: vector PDFs and PNG previews. `tf_consist` and `log_times` match the paper's topics.
- `irf_exp1` through `irf_exp4`: y1 response to innovation 1 at the smallest sample size, as in the slides.
- `irf_all_components_design*.pdf`: every response/innovation pair at every sample size.
- `fcst_rmse*`: every forecast variable and sample size; suffix `_qualified` denotes the diagnostic subset.
- `data/`: freshly computed statistics, component quantiles, diagnostics and replication-level losses.
- `tables/`: numeric CSVs and LaTeX tables for IRFs, forecasts at first/last horizon, and diagnostics.
- `paper_tables.tex`: standalone LaTeX document; compile from this report directory (requires booktabs).
- `report-metadata.rds`: original run configuration/execution metadata plus reporting session.

## Definitions

IRF errors are reduced-form estimated coefficients minus truth. Lag zero (the identity) is not scored.
In IRF summaries, horizon=0 means pooling all positive lags. In forecast summaries, horizon=0 means
pooling all forecast horizons; variable=0 means pooling variables. Pooling always occurs within
replication before computing Monte Carlo means/standard errors, not by treating entries as independent.

Forecast errors are realized future observations minus forecasts. Forecast gaps are fitted forecasts
minus the true-parameter conditional forecast based on the same observed training history.
Forecast relative RMSE = sqrt(paired method MSE / paired DGP MSE). IRF reference = oracle echelon.
paired_comparisons.csv additionally compares pooled forecast losses with oracle echelon.
These are finite-sample procedure comparisons, not estimates of an asymptotic efficiency penalty
attributable solely to non-identification. The model restrictions and estimation algorithms also differ.

The coefficient quantiles describe the estimates (truth is added back to saved IRF errors).
The bias column is mean estimation error. Forecast component quantiles describe forecast errors.
top_1pct_loss_share is the share of total loss from the largest ceiling(0.01*B) usable losses.
A large share flags sensitivity to a small number of replications; the observations are not discarded.

Means and variances are computed after rescaling to avoid numerical overflow with very large finite
losses. The summary CSVs originally saved by the runner are not overwritten. Raw checkpoint status
flags determine availability, consistently with the runner; a failed full forecast path is not partly salvaged.

Timing uses elapsed seconds measured inside workers, not user CPU time. Summed worker times are
not wall-clock run duration. Execution logs exclude summary writing and process cleanup; across resumed
runs the reported elapsed total is the sum of logged invocations, not calendar time including pauses.

## Overall diagnostics (counts)

```
       method fit_failures irf_failures forecast_failures nonconverged
          CML            0            0                 0            1
          DGP            0            0                 0            0
          DPF            0            0                 0            0
   ECH_oracle            0            0                 0            0
 ECH_selected            0            0                 0            0
          HR1            0            0                 0            0
          IHR            0            0                 0         6650
          MLE            0            0                 0            1
          NET            0            0                 0         6650
          RBY            0            0                 0          338
      VAR_BIC            0            0                 0            0
 unstable_ar noninvertible_ma qualified
           0                0     12000
           0                0     12000
          18               12     11977
           0                0     12000
           0                0     12000
           3               25     11973
        1117             1374     10586
           0                0     12000
         183              536     11452
           0                0     12000
           1                0     11999
```

## Exact index recovery at the largest sample size

Design 1: 90.4% exact index recovery at n = 400 (true vector: 1,0).
Design 2: 93.6% exact index recovery at n = 400 (true vector: 1,1).
Design 3: 44.1% exact index recovery at n = 400 (true vector: 1,1).
Design 4: 0.8% exact index recovery at n = 400 (true vector: 2,2,2).

## Reproduce

```r
source("inst/scripts/report_paper_simulations.R")
report_varma_paper_simulation("paper-results")
```

Only the report directory is replaced on another reporting run. No models are refitted.
The historical manuscript/slide figures are untouched. Review captions before inserting refreshed figures.
