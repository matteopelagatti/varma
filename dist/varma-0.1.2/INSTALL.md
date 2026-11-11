# varma 0.1.2: FKF replacement

The source archive `varma_0.1.2.tar.gz` replaces FKF.SP with FKF (>= 0.2.6).
The function `fit_varma_fkf()` and `mle_backend="fkf"` now use `FKF::fkf()`.
The default paper simulation backend remains `fit_varma_skf()`.

On the other PC, start a fresh R session and run:

```r
install.packages("FKF", repos = "https://cloud.r-project.org")
install.packages(file.choose(), repos = NULL, type = "source")
# Select varma_0.1.2.tar.gz in the file chooser.
packageVersion("varma")  # should be 0.1.2
```

This assumes the other varma dependencies have already been installed. To
install them as well:

```r
install.packages(c("Rcpp", "RcppArmadillo", "KFAS", "FKF", "Matrix",
                   "future", "future.apply", "glmnet"))
```

A source installation on Windows needs Rtools; consult the current CRAN
instructions for your R release:
https://cran.r-project.org/bin/windows/base/howto-R-4.6.html
The archive is source, rather than a binary built for R 4.5.

For `ret="ss"`, FKF returns joint innovations as `vt`, their covariance as
`Ft`, predicted states as `at` and `Pt`, and the likelihood as `logLik`.
These are the FKF field names, rather than the old FKF.SP raw output fields.

The installed 0.1.1 package on the development PC has not been replaced.
Existing simulations should finish with their original package version.
Checkpoint manifests deliberately reject mixing package versions.
