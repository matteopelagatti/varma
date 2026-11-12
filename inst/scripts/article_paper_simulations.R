# Article figures and LaTeX tables from saved VARMA simulation checkpoints.
# No estimation, simulation, Python, or LaTeX installation is needed to run this.
# Terminal: Rscript inst/scripts/article_paper_simulations.R paper-results
# R: source("inst/scripts/article_paper_simulations.R")
#    article <- article_varma_paper_simulation("paper-results")
# Keep report_paper_simulations.R beside this script: its validated collector
# supplies the statistics. See ARTICLE_OUTPUTS.md for definitions and options.

# Resolve the companion relative to this file, also when sourced from another
# working directory or from system.file("scripts", ..., package="varma").
.pa_companion <- local({
  frames <- sys.frames()
  paths <- unlist(lapply(frames, function(f) f$ofile), use.names = FALSE)
  if (!length(paths)) {
    paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  }
  if (!length(paths)) stop("Cannot locate this script; use source() or Rscript")
  file.path(dirname(normalizePath(tail(paths, 1L), winslash = "/", mustWork = TRUE)),
    "report_paper_simulations.R")
})

.pa_helpers <- function() {
  if (!file.exists(.pa_companion))
    stop("Keep report_paper_simulations.R in the same folder as this script")
  env <- new.env(parent = globalenv())
  sys.source(.pa_companion, envir = env)
  env
}

# Cache numerical aggregation, not graphics. Any change to inputs, trimming,
# collector code or installed package invalidates the cache automatically.
.pa_collect_cached <- function(input_dir, output_dir, helpers, trim, allow_partial) {
  files <- c(file.path(input_dir, "manifest.rds"),
    sort(list.files(input_dir, "^(chunk-|execution-).*[.]rds$", full.names=TRUE)))
  signature <- list(hashes=tools::md5sum(files), trim=trim,
    collector=tools::md5sum(.pa_companion), version=as.character(utils::packageVersion("varma")))
  path <- file.path(output_dir, "summary-cache.rds")
  cached <- if (file.exists(path)) tryCatch(readRDS(path), error=function(e) NULL) else NULL
  if (!is.null(cached) && identical(cached$signature, signature) &&
      (allow_partial || all(cached$summary$completion$missing == 0L))) {
    message("Using verified numerical cache; regenerating graphics and tables")
    return(cached$summary)
  }
  z <- helpers$.pr_collect(input_dir, trim, allow_partial)
  saveRDS(list(signature=signature, summary=z), path)
  z
}

# This only suppresses the known corrupted forecasts in the archived release.
# It does not remove unstable/nonconverged fits from the available subset.
.pa_data <- function(z, methods, subset) {
  cfg <- z$manifest$config
  withheld <- if (identical(cfg$package_version, "0.1.0"))
    intersect(methods, c("IHR", "HR1", "NET")) else character()
  targets <- c("irf", "forecast", "forecast_gap", "failures", "timing", "selection")
  out <- lapply(targets, function(nm) {
    d <- z[[nm]]
    if (!nrow(d)) return(d)
    d <- d[d$method %in% c("DGP", methods), , drop = FALSE]
    if ("subset" %in% colnames(d)) d <- d[d$subset == subset, , drop = FALSE]
    if (nm %in% c("forecast", "forecast_gap")) {
      d$withheld <- d$method %in% withheld
      value_cols <- setdiff(names(d)[vapply(d, is.numeric, logical(1))],
        c("n", "horizon", "variable", "attempted"))
      d[d$withheld, value_cols] <- NA_real_
    }
    d
  })
  names(out) <- targets
  out$withheld <- withheld
  out
}

.pa_prepare <- function(d, cfg) {
  d$Design <- factor(d$design, levels = names(cfg$designs),
    labels = paste("Design", seq_along(cfg$designs)))
  d
}

# Graphical settings are intentionally gathered here for easy journal styling.
.pa_plot <- function(d, cfg, helpers, x, y, xlab, ylab,
    methods, font_size, log_y = FALSE, reference = NULL) {
  d <- .pa_prepare(d, cfg)
  d$x <- d[[x]]; d$value <- d[[y]]
  d$value[!is.finite(d$value) | (log_y & d$value <= 0)] <- NA_real_
  labels <- helpers$.pr_labels[methods]
  labels[is.na(labels)] <- methods[is.na(labels)]
  p <- ggplot2::ggplot(d, ggplot2::aes(x, value, colour = method,
    linetype = method, shape = method, group = method))
  if (!is.null(reference)) p <- p + ggplot2::geom_hline(yintercept = reference,
    colour = "grey55", linewidth = .5, linetype = "dotted")
  if (length(unique(d$x)) > 1L) p <- p + ggplot2::geom_line(linewidth = .65, na.rm = TRUE)
  p <- p + ggplot2::geom_point(size = 2.2, stroke = .5, na.rm = TRUE) +
    ggplot2::facet_wrap(~Design, ncol = 2L, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$x))) +
    ggplot2::scale_colour_manual(values = helpers$.pr_colors, breaks = methods,
      labels = labels, name = NULL) +
    ggplot2::scale_shape_manual(values = stats::setNames(rep(c(16, 17, 15, 18, 0, 1, 2, 5, 6, 8),
      length.out = length(methods)), methods),
      breaks = methods, labels = labels, name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(rep(c("solid", "dashed", "dotdash"),
      length.out = length(methods)), methods), breaks = methods, labels = labels, name = NULL) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw(base_size = font_size, base_family = "serif") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = .2),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
      axis.text = ggplot2::element_text(size = font_size - 1, colour = "black"),
      axis.title = ggplot2::element_text(size = font_size),
      strip.text = ggplot2::element_text(face = "bold", size = font_size),
      legend.position = "bottom", legend.key.width = grid::unit(1.5, "lines"),
      legend.text = ggplot2::element_text(size = font_size - 1),
      panel.spacing = grid::unit(.9, "lines"),
      plot.caption = ggplot2::element_text(size = font_size - 1, hjust = 0),
      plot.margin = ggplot2::margin(7, 10, 7, 7)) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 4, byrow = TRUE),
      shape = ggplot2::guide_legend(ncol = 4, byrow = TRUE),
      linetype = ggplot2::guide_legend(ncol = 4, byrow = TRUE))
  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_number(
    scale_cut = scales::cut_short_scale()))
  p
}

.pa_cell <- function(x, se = NULL, helpers) {
  value <- helpers$.pr_num(x)
  if (!is.null(se) && is.finite(x) && is.finite(se))
    paste0("\\shortstack{", value, "\\\\(", helpers$.pr_num(se), ")}") else value
}

# Table fragments include their own float, caption, label and notes. A path
# macro lets the same fragments work in the preview and after copying to Overleaf.
.pa_table <- function(rows, headers, caption, note, name, output_dir, helpers) {
  file <- file.path(output_dir, "tables", paste0(name, ".tex"))
  lines <- c("% Generated by article_paper_simulations.R; requires booktabs and graphicx.",
    "\\begin{table}[p]\\centering", paste0("\\caption{", helpers$.pr_tex(caption), "}"),
    paste0("\\label{tab:sim-", gsub("_", "-", name), "}"),
    "\\small\\setlength{\\tabcolsep}{4pt}\\renewcommand{\\arraystretch}{1.15}",
    "\\resizebox{\\linewidth}{!}{%", paste0("\\begin{tabular}{ll",
      paste(rep("r", length(headers) - 2L), collapse = ""), "}\\toprule"),
    paste0(paste(headers, collapse = " & "), " \\\\\\midrule"),
    rows, "\\bottomrule\\end{tabular}}", "\\par\\medskip",
    paste0("\\begin{minipage}{\\linewidth}\\footnotesize ", helpers$.pr_tex(note),
      "\\end{minipage}"), "\\end{table}")
  writeLines(lines, file)
  file
}

.pa_wide <- function(d, cfg, methods, value, se, caption, note, name, output_dir, helpers) {
  headers <- c("Design", "$n$", vapply(methods, function(m) {
    if (m == "ECH_oracle") "\\shortstack{ECH\\\\oracle}" else
      if (m == "ECH_selected") "\\shortstack{ECH\\\\selected}" else helpers$.pr_tex(m)
  }, character(1)))
  rows <- character()
  for (design in names(cfg$designs)) {
    for (n in cfg$sample_sizes) {
      vals <- vapply(methods, function(method) {
        row <- d[d$design == design & d$n == n & d$method == method, , drop = FALSE]
        if (nrow(row) > 1L) stop("More than one summary for a table cell: ", name)
        if (!nrow(row)) return("--")
        if (value == "n_available") return(if (is.na(row[[value]])) "--" else as.character(row[[value]]))
        .pa_cell(row[[value]], if (is.null(se)) NULL else row[[se]], helpers)
      }, character(1))
      rows <- c(rows, paste0(paste(c(match(design, names(cfg$designs)), n, vals), collapse = " & "), " \\\\"))
    }
    rows <- c(rows, "\\addlinespace")
  }
  .pa_table(rows, headers, caption, note, name, output_dir, helpers)
}

#' Export figures and tables suitable for inclusion in an article.
#' Returns editable ggplot objects, the selected summaries, and file paths.
#' The PDF graphics contain axes/legends only; complete captions are in
#' figures.tex and the companion preview. No confidence bands are invented.
article_varma_paper_simulation <- function(input_dir = "paper-results",
    output_dir = file.path(input_dir, "article"), methods = NULL,
    subset = c("available", "qualified"), trim = .01, forecast_n = NULL,
    allow_partial = FALSE, width = 7, height = 6, font_size = 14, png = TRUE) {
  subset <- match.arg(subset)
  if (length(trim) != 1L || !is.finite(trim) || trim < 0 || trim >= .5)
    stop("trim must lie in [0, 0.5)")
  for (v in list(width, height, font_size)) if (length(v) != 1L || !is.finite(v) || v <= 0)
    stop("width, height and font_size must be positive numbers")
  for (pkg in c("varma", "ggplot2", "scales")) if (!requireNamespace(pkg, quietly = TRUE))
    stop("Please install the required R package: ", pkg)
  helpers <- .pa_helpers()
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(input_dir, "manifest.rds"))) stop("No manifest.rds in ", input_dir)
  cfg <- readRDS(file.path(input_dir, "manifest.rds"))$config
  if (is.null(methods)) methods <- cfg$methods
  if (!length(methods) || anyNA(methods) || anyDuplicated(methods) ||
      any(!methods %in% cfg$methods)) stop("methods must be distinct estimators present in the manifest")
  if (is.null(forecast_n)) forecast_n <- max(cfg$sample_sizes)
  if (length(forecast_n) != 1L || !forecast_n %in% cfg$sample_sizes)
    stop("forecast_n must be a sample size in the manifest")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  if (tolower(input_dir) == tolower(output_dir)) stop("Use a separate output directory")
  z <- .pa_collect_cached(input_dir, output_dir, helpers, trim, allow_partial)
  d <- .pa_data(z, methods, subset)
  complete <- all(z$completion$missing == 0L)
  status_note <- if (complete) "" else "INCOMPLETE RUN: only completed replications are included. "
  correction_note <- if (length(d$withheld)) paste0("Forecast results for ",
    paste(d$withheld, collapse = ", "), " are withheld: this run used varma 0.1.0 with the confirmed residual-handling bug. Rerun those estimates with a corrected release.") else ""
  if (nzchar(correction_note)) warning(correction_note, call. = FALSE)
  common <- paste0(status_note, if (subset == "available")
    "All usable estimates are included, even when nonconverged or unstable. " else
    "The qualified subset requires stable AR, invertible MA and applicable reported convergence. ",
    "Accuracy is conditional on usable outputs; consult usable counts and diagnostics. ")
  mcse <- "Parentheses contain delta-method Monte Carlo standard errors based on independent replication losses; extreme losses can make them unreliable. "
  trim_note <- sprintf("The trimmed statistic removes %.3g%% of replication losses from EACH tail before averaging and taking the square root; it is not ordinary RMSE. ", 100 * trim)
  axis_note <- "Vertical axes are logarithmic and vary by design; exact zeros appear only in the tables. Lines connect simulated points without a fitted convergence rate. "
  for (subdir in c("figures", "tables", "data")) dir.create(file.path(output_dir, subdir), showWarnings = FALSE)
  for (nm in setdiff(names(d), "withheld")) utils::write.csv(d[[nm]],
    file.path(output_dir, "data", paste0(nm, ".csv")), row.names = FALSE, na = "")
  utils::write.csv(z$completion, file.path(output_dir, "data", "completion.csv"), row.names = FALSE)
  utils::write.csv(z$indices, file.path(output_dir, "data", "index_frequencies.csv"), row.names = FALSE)
  input_files <- c(file.path(input_dir, "manifest.rds"),
    list.files(input_dir, "^chunk-.*[.]rds$", full.names = TRUE))
  utils::write.csv(data.frame(file = basename(input_files), md5 = unname(tools::md5sum(input_files))),
    file.path(output_dir, "data", "input_checksums.csv"), row.names = FALSE)
  meta <- list(generated = Sys.time(), config_id = z$manifest$config_id, config = cfg,
    execution = z$execution, input_dir = input_dir, methods = methods, subset = subset,
    trim = trim, forecast_n = forecast_n, width=width, height=height, font_size=font_size,
    complete = complete, withheld_forecasts = d$withheld,
    correction_note = correction_note, session = utils::sessionInfo())
  saveRDS(meta, file.path(output_dir, "article-metadata.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "session-info.txt"))

  ir <- d$irf[d$irf$horizon == 0L & d$irf$method != "DGP", , drop = FALSE]
  fc <- d$forecast[d$forecast$horizon == 0L & d$forecast$variable == 0L, , drop = FALSE]
  gap <- d$forecast_gap[d$forecast_gap$horizon == 0L & d$forecast_gap$variable == 0L &
    d$forecast_gap$method != "DGP", , drop = FALSE]
  fch <- d$forecast[d$forecast$n == forecast_n & d$forecast$horizon > 0L &
    d$forecast$variable == 0L & !d$forecast$withheld, , drop = FALSE]
  eligible <- setdiff(methods, d$withheld)
  plots <- list(); captions <- character()
  add <- function(name, data, x, y, xlab, ylab, caption, shown = methods,
      log_y = TRUE, reference = NULL) {
    if (!nrow(data)) return(invisible(NULL))
    plots[[name]] <<- .pa_plot(data, cfg, helpers, x, y, xlab, ylab, shown,
      font_size, log_y, reference)
    captions[[name]] <<- paste0(caption, common, if (log_y) axis_note else
      "Vertical axes vary by design. Lines connect simulated points. ")
  }
  add("irf_rmse", ir, "n", "rmse", "Sample size", "IRF RMSE (log scale)",
    sprintf("Reduced-form IRF accuracy, pooling lags 1-%d and all matrix entries within each replication before averaging losses across replications. ", cfg$irf_lags))
  add("irf_rmse_trimmed", ir, "n", "trimmed_rmse", "Sample size", "Trimmed IRF RMSE",
    paste0("Reduced-form IRF accuracy. ", trim_note), log_y = FALSE)
  add("forecast_rmse", fch, "horizon", "rmse", "Forecast horizon", "Forecast RMSE (log scale)",
    paste0("Forecast accuracy at n = ", forecast_n, ". Variables are pooled within replication. DGP uses the true parameters with the same observed training history. ", correction_note, " "), shown = c("DGP", eligible))
  add("forecast_relative_rmse", fch[fch$method != "DGP", , drop = FALSE], "horizon", "relative_rmse",
    "Forecast horizon", "RMSE / DGP RMSE (log scale)",
    paste0("Forecast accuracy relative to DGP at n = ", forecast_n,
      ". Each ratio divides two RMSEs computed on the same usable replications; it is not the mean of replication-wise ratios. The dotted line is one. ", correction_note, " "),
    shown = eligible, reference = 1)
  add("forecast_gap_rmse", gap[!gap$withheld, , drop = FALSE], "n", "rmse", "Sample size",
    "Forecast gap RMSE (log scale)", paste0("Distance from the DGP conditional forecast, pooling all horizons and variables within replication. This excludes the realized future shock contribution. ", correction_note, " "), shown = eligible)
  add("fit_time", d$timing[d$timing$method != "DGP", , drop = FALSE], "n", "mean_fit_seconds",
    "Sample size", "Seconds per attempted fit (log scale)",
    "Mean elapsed estimation time, including failed attempts, measured inside the simulation workers. All attempts are used irrespective of the accuracy subset. Worker times do not measure overall wall-clock time. ")

  device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  saveRDS(plots, file.path(output_dir, "plots.rds"))
  message("Writing ", length(plots), " article figures ...")
  for (name in names(plots)) {
    ggplot2::ggsave(file.path(output_dir, "figures", paste0(name, ".pdf")), plots[[name]],
      device = device, width = width, height = height, units = "in", bg = "white")
    if (png) ggplot2::ggsave(file.path(output_dir, "figures", paste0(name, ".png")), plots[[name]],
      width = width, height = height, units = "in", dpi = 300, bg = "white")
  }
  device(file.path(output_dir, "article_figures.pdf"), width = width + .5, height = height + 2.5)
  tryCatch({
    for (name in names(plots)) print(plots[[name]] + ggplot2::labs(
      title = gsub("_", " ", name), caption = paste(strwrap(captions[[name]], 90), collapse = "\n")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size=10,hjust=0)))
  }, finally = grDevices::dev.off())

  tables <- character()
  wide <- function(data, value, se = NULL, caption, note, name, shown = methods) {
    tables[[name]] <<- .pa_wide(data, cfg, shown, value, se, caption,
      paste0(common, note, "Missing or withheld values are shown as --."), name, output_dir, helpers)
  }
  wide(ir, "rmse", "rmse_mcse", "Reduced-form IRF accuracy.",
    paste0(mcse, "Losses pool all entries and positive response lags within replication. "), "irf_rmse")
  wide(ir, "trimmed_rmse", caption = "Trimmed reduced-form IRF accuracy.", note = trim_note, name = "irf_rmse_trimmed")
  wide(fc, "rmse", "rmse_mcse", "Forecast accuracy, pooled over horizons and variables.",
    paste0(mcse, correction_note, " "), "forecast_rmse", c("DGP", methods))
  wide(fc[fc$method != "DGP", ], "relative_rmse", "relative_rmse_mcse", "Forecast RMSE relative to DGP.",
    paste0(mcse, "Ratios use common usable replications for each estimator and DGP. ", correction_note, " "), "forecast_relative_rmse")
  wide(gap, "rmse", "rmse_mcse", "Distance from the DGP conditional forecast.",
    paste0(mcse, "All forecast horizons and variables are pooled within replication. ", correction_note, " "), "forecast_gap_rmse")
  wide(d$timing[d$timing$method != "DGP", ], "mean_fit_seconds", caption = "Mean elapsed seconds per attempted fit.",
    note = "Timings include all attempted fits irrespective of the accuracy subset, with failures and concurrent worker load. They exclude forecast evaluation and cannot be summed to obtain wall-clock duration. ", name = "fit_time")
  wide(ir, "n_available", caption = "Usable IRF replications.", note = "Counts correspond to the selected accuracy subset. ", name = "irf_counts")
  wide(fc[fc$method != "DGP", ], "n_available", caption = "Usable forecast replications.",
    note = paste0("Counts correspond to the selected accuracy subset. ", correction_note, " "), name = "forecast_counts")

  for (i in seq_along(cfg$designs)) {
    diag <- d$failures[d$failures$design == names(cfg$designs)[i] & d$failures$method != "DGP", ]
    if (!nrow(diag)) next
    diag <- diag[order(diag$n, match(diag$method, methods)), ]
    fields <- c("attempted", "fit_failures", "irf_failures", "forecast_failures",
      "nonconverged", "unstable_ar", "noninvertible_ma")
    rows <- vapply(seq_len(nrow(diag)), function(j) paste0(paste(c(diag$n[j],
      helpers$.pr_tex(helpers$.pr_labels[diag$method[j]]), unlist(diag[j, fields])), collapse = " & "), " \\\\"), character(1))
    name <- paste0("diagnostics_design", i)
    tables[[name]] <- .pa_table(rows, c("$n$", "Method", "Attempts", "Fit fail", "IRF fail",
      "Fcst fail", "Nonconv.", "AR unstable", "MA noninv."), paste0("Design ", i, ": numerical diagnostics."),
      paste0(status_note, "Counts use all attempts, including methods whose forecast accuracy is withheld. Categories overlap. Forecast failures here describe the archived calculation, not the publication exclusion. Nonconvergence is counted only where applicable; NET refers to its IHR pilot. Unknown convergence is recorded in data/failures.csv. ", correction_note), name, output_dir, helpers)
  }
  sel <- d$selection
  if (nrow(sel)) sel <- sel[sel$method == "ECH_selected", , drop = FALSE]
  if (nrow(sel)) {
    sel <- sel[order(match(sel$design, names(cfg$designs)), sel$n), ]
    rows <- vapply(seq_len(nrow(sel)), function(j) paste0(paste(c(
      match(sel$design[j], names(cfg$designs)), sel$n[j], helpers$.pr_tex(sel$true_indices[j]),
      sel$attempted[j], sel$returned[j], .pa_cell(100 * sel$correct_fraction[j],
        100 * sel$fraction_mcse[j], helpers)), collapse = " & "), " \\\\"), character(1))
    tables[["index_selection"]] <- .pa_table(rows, c("Design", "$n$", "True indices", "Attempts",
      "Returned", "Correct (\\%)"), "Data-driven Kronecker-index selection.",
      paste0(status_note, "The correct-selection percentage uses all attempts, including failed fits. Parentheses contain binomial Monte Carlo standard errors in percentage points. Returned counts refer to index vectors retained in saved diagnostics; failure may prevent recording a selected vector."),
      "index_selection", output_dir, helpers)
  }

  figure_tex <- c("% Requires graphicx. Define \\varmaresults as the path to this export folder.",
    "\\providecommand{\\varmaresults}{.}")
  for (name in names(plots)) figure_tex <- c(figure_tex, "\\begin{figure}[p]\\centering",
    paste0("\\includegraphics[width=\\linewidth]{\\varmaresults/figures/", name, ".pdf}"),
    paste0("\\caption{", helpers$.pr_tex(captions[[name]]), "}"),
    paste0("\\label{fig:sim-", gsub("_", "-", name), "}"), "\\end{figure}")
  writeLines(figure_tex, file.path(output_dir, "figures.tex"))
  writeLines(c("\\providecommand{\\varmaresults}{.}",
    paste0("\\input{\\varmaresults/tables/", names(tables), ".tex}")), file.path(output_dir, "tables.tex"))
  writeLines(c("\\documentclass[10pt]{article}", "\\ifdefined\\pdfminorversion\\pdfminorversion=7\\fi",
    "\\usepackage[a4paper,margin=20mm]{geometry}",
    "\\usepackage{graphicx,booktabs}", "\\newcommand{\\varmaresults}{.}", "\\begin{document}",
    "\\section*{VARMA simulation figures and tables}", helpers$.pr_tex(paste0(status_note,
      "Simulation package version: ", cfg$package_version, ". Accuracy subset: ", subset, ". ", correction_note)),
    "\\par\\medskip Figures use the sample size stated in each caption; tables cover all sample sizes.",
    "\\input{figures.tex}\\clearpage", "\\input{tables.tex}\\clearpage", "\\end{document}"),
    file.path(output_dir, "article_preview.tex"))
  writeLines(c("# Article exports from R", "", paste0("Source: `", input_dir, "`."),
    paste0("Simulation version: ", cfg$package_version, "; subset: ", subset, "; complete: ", complete, "."),
    "", status_note, correction_note, "", "## Files", "",
    "- figures/: vector PDFs and optional 300-dpi PNGs, with editable ggplot objects returned by the function.",
    "- article_figures.pdf: all figures with their explanatory captions.",
    "- tables/: LaTeX table floats, labels and notes, using booktabs and graphicx.",
    "- figures.tex and tables.tex: ready-to-include collections.",
    "- article_preview.tex: optional standalone preview; compile here with pdflatex article_preview.tex.",
    "- data/: numeric summaries, usable counts, diagnostics, selection frequencies and source checksums.",
    "- article-metadata.rds and session-info.txt: experiment configuration, execution records and export session.",
    "", "## Interpretation", "", common, trim_note, mcse,
    "RMSE is sqrt(mean(replication MSE)); each replication pools the relevant entries/horizons first.",
    paste0("Figure forecast curves use n=", forecast_n, "; CSV summaries retain every sample size, horizon and variable."),
    "IRFs are reduced form, without covariance orthogonalization. No FEVD or interval coverage is inferred.",
    "No data or models are refitted, no extreme finite losses are discarded from untrimmed summaries, and raw checkpoints are unchanged.",
    "", "## Include in Overleaf", "", "Copy this folder to the manuscript as simulation-output/.",
    "Add \\usepackage{graphicx,booktabs} and \\newcommand{\\varmaresults}{simulation-output} to the preamble.",
    "Then use \\input{simulation-output/figures.tex} and/or \\input{simulation-output/tables.tex}.",
    "For individual items, copy a figure environment from figures.tex or input a chosen tables/*.tex file.",
    "Captions describe the actual estimators and exclusions; existing manuscript claims are not automatically updated."),
    file.path(output_dir, "README.md"))
  message("Article exports saved in ", output_dir)
  invisible(list(output_dir = output_dir, plots = plots, captions = captions, tables = tables,
    summaries = d, metadata = meta))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat("Usage: Rscript article_paper_simulations.R [input_dir] [output_dir]\n",
      "  --methods=IHR,MLE,ECH_oracle,ECH_selected --subset=available|qualified\n",
      "  --trim=0.01 --forecast-n=400 --allow-partial --no-png\n",
      "Defaults: paper-results -> paper-results/article; all methods; available subset.\n", sep = "")
  } else {
    positional <- args[!startsWith(args, "--")]
    flags <- args[startsWith(args, "--")]
    known <- grepl("^--(methods=|subset=|trim=|forecast-n=|allow-partial$|no-png$)", flags)
    if (any(!known) || length(positional) > 2L) stop("Unknown arguments; use --help")
    option <- function(name, default = NULL) {
      x <- flags[startsWith(flags, paste0("--", name, "="))]
      if (length(x) > 1L) stop("Duplicate option: ", name)
      if (length(x)) sub(paste0("^--", name, "="), "", x) else default
    }
    input <- if (length(positional)) positional[1L] else "paper-results"
    method_option <- option("methods")
    article_varma_paper_simulation(input,
      output_dir = if (length(positional) > 1L) positional[2L] else file.path(input, "article"),
      methods = if (is.null(method_option)) NULL else strsplit(method_option, ",", fixed = TRUE)[[1L]],
      subset = option("subset", "available"), trim = as.numeric(option("trim", ".01")),
      forecast_n = if (is.null(option("forecast-n"))) NULL else as.numeric(option("forecast-n")),
      allow_partial = "--allow-partial" %in% flags, png = !"--no-png" %in% flags)
  }
}
