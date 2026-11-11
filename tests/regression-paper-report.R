# Reporting regression checks use synthetic saved errors, never fitted models.
script <- system.file("scripts", "report_paper_simulations.R", package = "varma")
if (!nzchar(script)) {
  candidates <- c("inst/scripts/report_paper_simulations.R", "../inst/scripts/report_paper_simulations.R")
  script <- candidates[file.exists(candidates)][1L]
}
source(script)
near <- function(x, y, tolerance = 1e-12) stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance)))
must_fail <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}

# Analytic reference values, including a variance calculation that would
# overflow if the unscaled 1e300 losses were squared internally.
large <- .pr_loss(c(1e300, 2e300, 3e300), c(1, 2, 3))
near(large$mse / 1e300, 2)
near(large$mse_mcse / 1e300, 1 / sqrt(3))
near(large$rmse_mcse / 1e150, 1 / (2 * sqrt(6)))
near(large$relative_rmse / 1e150, 1)
near(large$relative_rmse_mcse, 0)
small <- .pr_loss(c(1, 4, NA, 16), c(1, NA, 9, 4))
near(small$mse, 7)
near(small$relative_mse, 17 / 5)
stopifnot(small$n_available == 3, small$n_paired == 2)
near(.pr_loss(rep(0, 3), rep(1, 3))$rmse_mcse, 0)
near(.pr_component(c(-1, 1), 3)[["median"]], 3)
near(.pr_component(c(-1, 1), 3)[["rmse"]], 1)

local({
  directory <- tempfile("varma-report-regression-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  model <- varma::varma_paper_designs()[1L]
  cfg <- list(schema = 1L, designs = model, methods = c("ECH_oracle", "MLE"),
    nsim = 4L, sample_sizes = 20L, irf_lags = 2L, h = 2L, seed = 1L)
  manifest <- list(config = cfg, config_id = "regression-fixture")
  manifest_file <- file.path(directory, "manifest.rds")
  saveRDS(manifest, manifest_file)
  indices <- varma::kronecker_indices(model[[1L]])
  actual <- matrix(1:4, 2L, 2L)
  results <- lapply(1:4, function(b) {
    methods <- lapply(c("DGP", cfg$methods), function(method) {
      fc_ok <- method != "MLE" || b != 3L
      irf_ok <- method != "ECH_oracle" || b != 4L
      forecast <- if (method == "DGP") actual * 0 else actual * 0 + b
      list(status = list(method = method, fit_ok = TRUE, irf_ok = irf_ok, forecast_ok = fc_ok,
        fit_seconds = .1, evaluation_seconds = .01,
        convergence = if (method == "MLE" && b == 2L) 1L else 0L,
        convergence_type = if (method == "DGP") "not_applicable" else "estimator",
        ar_stable = TRUE, ma_invertible = TRUE, qualified = method != "MLE" || b != 2L,
        indices = if (method == "ECH_oracle") paste(indices, collapse = ",") else "",
        index_correct = if (method == "ECH_oracle") TRUE else NA,
        warnings = "", error = ""),
        irf_error = array(if (!irf_ok) NA_real_ else if (method == "DGP") 0 else b / 10, c(2L, 2L, 2L)),
        forecast_error = if (fc_ok) actual - forecast else actual * NA_real_,
        forecast_gap = if (fc_ok) forecast else actual * NA_real_)
    })
    names(methods) <- c("DGP", cfg$methods)
    list(replication = b, results = methods)
  })
  chunk <- list(schema = 1L, config_id = manifest$config_id,
    cell = list(id = 1L, design = names(model), n = 20L, indices = indices),
    replications = 1:4, results = results)
  file <- file.path(directory, "chunk-001-000001-000004.rds")
  saveRDS(chunk, file)
  before <- tools::md5sum(c(file, manifest_file))
  z <- .pr_collect(directory)
  stopifnot(identical(before, tools::md5sum(c(file, manifest_file))), z$completion$missing == 0)
  ir <- subset(z$irf, method == "MLE" & subset == "available" & horizon == 0)
  near(ir$mse, mean((1:4 / 10)^2)); stopifnot(ir$n_available == 4, ir$n_paired == 3)
  near(ir$relative_rmse, 1)
  fc <- subset(z$forecast, method == "MLE" & subset == "available" & horizon == 1 & variable == 0)
  # Horizons are rows, variables columns: h=1 corresponds to actual c(1,3).
  losses <- vapply(c(1, 2, 4), function(b) mean((c(1, 3) - b)^2), numeric(1))
  near(fc$mse, mean(losses)); near(fc$rmse, sqrt(mean(losses)))
  near(fc$mse_mcse, stats::sd(losses) / sqrt(3))
  stopifnot(fc$n_available == 3, fc$n_paired == 3)
  fcq <- subset(z$forecast, method == "MLE" & subset == "qualified" & horizon == 1 & variable == 0)
  near(fcq$mse, mean(losses[c(1, 3)])); stopifnot(fcq$n_available == 2)
  comp <- subset(z$irf_components, method == "MLE" & subset == "available" &
    response == 1 & shock == 2 & horizon == 1)
  truth <- varma::irf(model[[1L]], 2L, orth = "none")[1, 2, 2]
  near(comp$median, truth + .25); near(comp$bias, .25)
  stopifnot(subset(z$failures, method == "MLE")$forecast_failures == 1,
    subset(z$failures, method == "MLE")$nonconverged == 1)
  copy <- file.path(directory, "chunk-001-000005-000008.rds")
  file.copy(file, copy)
  must_fail(.pr_collect(directory), "Duplicate")
  unlink(copy)
  chunk$config_id <- "different-experiment"; saveRDS(chunk, file)
  must_fail(.pr_collect(directory), "Incompatible")
  chunk$config_id <- manifest$config_id
  chunk$replications <- 1:3; chunk$results <- results[1:3]; saveRDS(chunk, file)
  must_fail(.pr_collect(directory), "Incomplete")
  partial <- .pr_collect(directory, allow_partial = TRUE)
  stopifnot(partial$completion$missing == 1)
  # Optional end-to-end graphics check, including a completely absent cell.
  # Ordinary package regression checks need neither ggplot2 nor a PDF device.
  if (identical(Sys.getenv("VARMA_REPORT_RENDER_TEST"), "true")) {
    manifest$config$sample_sizes <- c(20L, 40L)
    manifest$config$mle_backend <- "skf"
    saveRDS(manifest, manifest_file)
    report <- report_varma_paper_simulation(directory, allow_partial = TRUE, details = FALSE)
    stopifnot(!report$complete, file.exists(file.path(report$report_dir, "paper_report.pdf")),
      file.exists(file.path(report$report_dir, "report-summary.rds")),
      file.exists(file.path(report$report_dir, "tables", "irf_design1_available.tex")))
  }
})
cat("Paper reporting regression checks passed.\n")
