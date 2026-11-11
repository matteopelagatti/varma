# Protect publication outputs from the known version-specific forecast bug.
paths <- c("inst/scripts/article_paper_simulations.R", "../inst/scripts/article_paper_simulations.R",
  system.file("scripts", "article_paper_simulations.R", package = "varma"))
script <- paths[nzchar(paths) & file.exists(paths)][1L]
stopifnot(!is.na(script))
source(script)
helpers <- .pa_helpers()

local({
  # The same numerical results must be withheld for the affected historical
  # release, but retained for a corrected release, regardless of installed R pkg.
  methods <- c("IHR", "HR1", "NET", "MLE", "ECH_oracle", "ECH_selected")
  base <- expand.grid(method = c("DGP", methods), subset = c("available", "qualified"),
    stringsAsFactors = FALSE)
  base$design <- "model1"; base$n <- 100L; base$horizon <- 0L
  base$variable <- 0L; base$attempted <- 10L; base$n_available <- 9L
  base$mse <- 16; base$rmse <- 4; base$rmse_mcse <- .5
  base$relative_rmse <- 2; base$n_paired <- 8L
  cfg <- list(package_version = "0.1.0", methods = methods,
    designs = list(model1 = NULL), sample_sizes = 100L)
  z <- list(manifest = list(config = cfg), irf = base, forecast = base,
    forecast_gap = base, failures = base, timing = base, selection = data.frame())
  original <- z
  old <- .pa_data(z, methods, "available")
  stopifnot(identical(z, original), all(old$irf$rmse == 4),
    all(old$forecast$subset == "available"),
    all(is.na(old$forecast$rmse[old$forecast$withheld])),
    all(is.na(old$forecast$n_available[old$forecast$withheld])),
    all(old$forecast$attempted == 10),
    all(old$forecast$rmse[!old$forecast$withheld] == 4),
    all(old$failures$rmse == 4), identical(old$withheld, c("IHR", "HR1", "NET")))
  z$manifest$config$package_version <- "0.1.1"
  corrected <- .pa_data(z, methods, "available")
  stopifnot(!any(corrected$forecast$withheld), all(corrected$forecast$rmse == 4))
  selected <- .pa_data(z, c("MLE", "ECH_selected"), "qualified")
  stopifnot(setequal(selected$irf$method, c("DGP", "MLE", "ECH_selected")),
    all(selected$irf$subset == "qualified"))

  directory <- tempfile("varma-article-")
  dir.create(directory); dir.create(file.path(directory, "tables"))
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  file <- .pa_wide(old$forecast, cfg, c("IHR", "MLE"), "rmse", "rmse_mcse",
    "Forecast RMSE", "Known faulty forecasts withheld.", "fixture", directory, helpers)
  text <- paste(readLines(file), collapse = "\n")
  stopifnot(grepl("1 & 100 & -- & \\shortstack{4.000\\\\(0.500)}", text, fixed = TRUE),
    grepl("\\label{tab:sim-fixture}", text, fixed = TRUE))
})
cat("Article export regression checks passed.\n")
