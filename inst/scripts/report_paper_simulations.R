# Report saved paper simulations; no model is fitted by this script.
# Terminal: Rscript inst/scripts/report_paper_simulations.R paper-results
# R: source("inst/scripts/report_paper_simulations.R")
#    report_varma_paper_simulation("paper-results")
# Requires the simulation version of varma and ggplot2. See PAPER_SIMULATIONS.md.

.pr_mean <- function(x) {
  s <- max(abs(x))
  if (s == 0) 0 else s * mean(x / s)
}

.pr_mcse <- function(x) {
  if (length(x) < 2L) return(NA_real_)
  s <- max(abs(x))
  if (s == 0) 0 else s * (stats::sd(x / s) / sqrt(length(x)))
}

# Scaling before computing variances matters here: some finite forecast losses
# are near 1e300, and sd(loss) would square them and overflow unnecessarily.
.pr_loss <- function(loss, reference = NULL, trim = .01) {
  x <- loss[is.finite(loss)]
  n <- length(x)
  mse <- if (n) .pr_mean(x) else NA_real_
  rmse <- sqrt(mse)
  se <- if (n) .pr_mcse(x) else NA_real_
  ans <- list(n_available = n, mse = mse, rmse = rmse, mse_mcse = se,
    rmse_mcse = if (isTRUE(rmse == 0)) 0 else (se / rmse) / 2,
    trimmed_rmse = if (n) sqrt(mean(x, trim = trim)) else NA_real_,
    median_replication_rmse = if (n) stats::median(sqrt(x)) else NA_real_,
    max_replication_rmse = if (n) sqrt(max(x)) else NA_real_,
    top_1pct_loss_share = if (!n) NA_real_ else if (max(x) == 0) 0 else {
      u <- sort(x / max(x), decreasing = TRUE)
      sum(head(u, max(1L, ceiling(.01 * n)))) / sum(u)
    },
    n_paired = 0L, paired_mse = NA_real_, reference_mse = NA_real_,
    relative_mse = NA_real_, relative_mse_mcse = NA_real_,
    relative_rmse = NA_real_, relative_rmse_mcse = NA_real_,
    paired_mse_difference = NA_real_, paired_mse_difference_mcse = NA_real_)
  if (!is.null(reference)) {
    paired <- is.finite(loss) & is.finite(reference)
    a <- loss[paired]; b <- reference[paired]; np <- length(a)
    ans$n_paired <- np
    if (np) {
      ma <- .pr_mean(a); mb <- .pr_mean(b)
      ans$paired_mse <- ma; ans$reference_mse <- mb
      ans$paired_mse_difference <- .pr_mean(a - b)
      ans$paired_mse_difference_mcse <- .pr_mcse(a - b)
      if (mb > 0) {
        ratio <- ma / mb
        ans$relative_mse <- ratio
        ans$relative_rmse <- sqrt(ma) / sqrt(mb)
        if (np > 1L) {
          relative_se <- if (ma == 0) 0 else
            stats::sd(a / ma - b / mb) / sqrt(np)
          ans$relative_mse_mcse <- ratio * relative_se
          ans$relative_rmse_mcse <- ans$relative_rmse * (relative_se / 2)
        }
      }
    }
  }
  ans
}

.pr_component <- function(x, truth = 0) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(n_available = 0, truth = truth,
    q05 = NA, q25 = NA, median = NA, q75 = NA, q95 = NA,
    bias = NA, rmse = NA, sd = NA))
  q <- stats::quantile(x, c(.05, .25, .5, .75, .95), names = FALSE)
  c(n_available = length(x), truth = truth,
    stats::setNames(q + truth, c("q05", "q25", "median", "q75", "q95")),
    bias = .pr_mean(x), rmse = sqrt(.pr_mean(x^2)),
    sd = .pr_mcse(x) * sqrt(length(x)))
}

.pr_row_loss <- function(x) {
  # Divide before summing to avoid overflowing the sum of finite squares.
  rowSums(x^2 / ncol(x))
}

.pr_collect <- function(input_dir, trim = .01, allow_partial = FALSE) {
  manifest <- readRDS(file.path(input_dir, "manifest.rds"))
  cfg <- manifest$config
  if (!identical(cfg$schema, 1L) || !length(cfg$designs) ||
      !length(cfg$methods) || anyDuplicated(cfg$methods))
    stop("Unsupported or malformed simulation manifest")
  methods <- c("DGP", cfg$methods)
  files <- list.files(input_dir, "^chunk-.*[.]rds$", full.names = TRUE)
  if (!length(files)) stop("No simulation checkpoints were found")
  used <- character(); cells <- list(); cell_id <- 0L
  out <- stats::setNames(vector("list", 12L), c("irf", "forecast", "forecast_gap",
    "failures", "timing", "indices", "diagnostics", "irf_components",
    "forecast_components", "replication_losses", "selection", "paired_comparisons"))
  append_row <- function(target, row) {
    out[[target]][[length(out[[target]]) + 1L]] <<- row
  }
  for (n in cfg$sample_sizes) for (design in names(cfg$designs)) {
    cell_id <- cell_id + 1L
    model <- cfg$designs[[design]]; m <- as.integer(dim(model)[1L])
    J <- cfg$irf_lags; H <- cfg$h; B <- cfg$nsim
    truth <- varma::irf(model, J, orth = "none")[, , -1L, drop = FALSE]
    true_indices <- paste(varma::kronecker_indices(model), collapse = ",")
    these <- files[grepl(sprintf("^chunk-%03d-[0-9]+-[0-9]+[.]rds$", cell_id), basename(files))]
    used <- c(used, these)
    E <- lapply(methods, function(z) matrix(NA_real_, B, m * m * J))
    F <- G <- lapply(methods, function(z) matrix(NA_real_, B, H * m))
    names(E) <- names(F) <- names(G) <- methods
    qualified <- matrix(FALSE, B, length(methods), dimnames = list(NULL, methods))
    present <- rep(FALSE, B); diagnostics <- list()
    for (file in these) {
      chunk <- readRDS(file)
      ids <- chunk$replications
      if (!identical(chunk$schema, 1L) || !identical(chunk$config_id, manifest$config_id) ||
          !identical(chunk$cell$design, design) || chunk$cell$n != n ||
          chunk$cell$id != cell_id || !length(ids) || anyNA(ids) ||
          any(ids != as.integer(ids) | ids < 1L | ids > B) || anyDuplicated(ids) ||
          length(chunk$results) != length(ids) ||
          !identical(as.integer(ids), vapply(chunk$results, function(z) as.integer(z$replication), integer(1))))
        stop("Incompatible or invalid checkpoint: ", file)
      if (paste(chunk$cell$indices, collapse = ",") != true_indices)
        stop("Population index mismatch in ", file)
      for (rep in chunk$results) {
        b <- rep$replication
        if (present[b]) stop("Duplicate replication in ", file)
        present[b] <- TRUE
        if (!setequal(names(rep$results), methods)) stop("Method mismatch in ", file)
        for (method in methods) {
          z <- rep$results[[method]]; s <- z$status
          if (!identical(s$method, method) || length(s$qualified) != 1L)
            stop("Invalid method diagnostics in ", file)
          qualified[b, method] <- isTRUE(s$qualified)
          if (isTRUE(s$irf_ok)) {
            if (!identical(dim(z$irf_error), c(m, m, as.integer(J))) ||
                any(!is.finite(z$irf_error^2))) stop("Invalid available IRF in ", file)
            E[[method]][b, ] <- as.numeric(z$irf_error)
          }
          if (isTRUE(s$forecast_ok)) {
            if (!identical(dim(z$forecast_error), c(as.integer(H), m)) ||
                !identical(dim(z$forecast_gap), c(as.integer(H), m)) ||
                any(!is.finite(c(z$forecast_error, z$forecast_gap)^2)))
              stop("Invalid available forecast in ", file)
            F[[method]][b, ] <- as.numeric(z$forecast_error)
            G[[method]][b, ] <- as.numeric(z$forecast_gap)
          }
          diagnostics[[length(diagnostics) + 1L]] <- as.data.frame(
            c(list(design = design, n = n, replication = b), s), stringsAsFactors = FALSE)
        }
      }
    }
    attempted <- sum(present)
    cells[[cell_id]] <- data.frame(design = design, n = n, planned = B,
      attempted = attempted, missing = B - attempted, chunks = length(these),
      true_indices = true_indices)
    if (attempted != B && !allow_partial)
      stop("Incomplete cell ", design, ", n=", n, ": ", attempted, "/", B,
        ". Use allow_partial=TRUE (CLI --allow-partial) to label and report an incomplete run.")
    if (!attempted) next
    message("Reading ", design, ", n=", n, ": ", attempted, "/", B, " replications")
    diag <- do.call(rbind, diagnostics); append_row("diagnostics", diag)
    # Every pooled loss uses one observation per independent replication.
    L <- list(irf = lapply(E, .pr_row_loss), forecast = lapply(F, .pr_row_loss),
              forecast_gap = lapply(G, .pr_row_loss))
    for (method in methods) {
      d <- diag[diag$method == method, , drop = FALSE]
      applicable <- d$convergence_type != "not_applicable" & d$fit_ok
      append_row("failures", data.frame(design = design, n = n, method = method,
        planned = B, attempted = attempted, fit_failures = sum(!d$fit_ok),
        irf_failures = sum(!d$irf_ok), forecast_failures = sum(!d$forecast_ok),
        nonconverged = sum(applicable & !is.na(d$convergence) & d$convergence != 0L),
        convergence_unknown = sum(applicable & is.na(d$convergence)),
        unstable_ar = sum(d$fit_ok & !is.na(d$ar_stable) & !d$ar_stable),
        noninvertible_ma = sum(d$fit_ok & !is.na(d$ma_invertible) & !d$ma_invertible),
        qualified = sum(d$qualified), with_warnings = sum(nzchar(d$warnings))))
      append_row("timing", data.frame(design = design, n = n, method = method,
        attempted = attempted, mean_fit_seconds = mean(d$fit_seconds),
        median_fit_seconds = stats::median(d$fit_seconds),
        p90_fit_seconds = as.numeric(stats::quantile(d$fit_seconds, .9)),
        total_fit_seconds = sum(d$fit_seconds),
        mean_evaluation_seconds = mean(d$evaluation_seconds)))
      if (method %in% c("ECH_oracle", "ECH_selected")) {
        labels <- ifelse(nzchar(d$indices), d$indices, "unavailable_after_failed_fit")
        counts <- table(labels)
        append_row("indices", data.frame(design = design, n = n, method = method,
          indices = names(counts), count = as.integer(counts), attempted = attempted,
          fraction = as.integer(counts) / attempted, true_indices = true_indices))
        correct <- sum(d$index_correct %in% TRUE); returned <- sum(nzchar(d$indices))
        p <- correct / attempted
        append_row("selection", data.frame(design = design, n = n, method = method,
          true_indices = true_indices, attempted = attempted, returned = returned,
          correct = correct, correct_fraction = p, fraction_mcse = sqrt(p * (1 - p) / attempted),
          correct_given_returned = if (returned) correct / returned else NA_real_))
      }
      append_row("replication_losses", data.frame(design = design, n = n, method = method,
        replication = which(present), qualified = qualified[present, method],
        irf_mse = L$irf[[method]][present], forecast_mse = L$forecast[[method]][present],
        forecast_gap_mse = L$forecast_gap[[method]][present]))
      for (subset in c("available", "qualified")) {
        mask <- if (subset == "qualified") present & qualified[, method] else present
        ref <- L$irf$ECH_oracle
        if (!is.null(ref) && subset == "qualified") ref[!qualified[, "ECH_oracle"]] <- NA_real_
        for (j in 0:J) {
          cols <- if (j == 0) seq_len(m * m * J) else (j - 1L) * m * m + seq_len(m * m)
          loss <- .pr_row_loss(E[[method]][, cols, drop = FALSE]); loss[!mask] <- NA_real_
          refj <- if (is.null(ref)) NULL else .pr_row_loss(E$ECH_oracle[, cols, drop = FALSE])
          if (!is.null(refj) && subset == "qualified") refj[!qualified[, "ECH_oracle"]] <- NA_real_
          append_row("irf", as.data.frame(c(list(design = design, n = n, method = method,
            subset = subset, horizon = j, attempted = attempted, reference = "ECH_oracle"),
            .pr_loss(loss, refj, trim))))
        }
        for (h in 0:H) for (v in 0:m) for (target in c("forecast", "forecast_gap")) {
          hs <- if (h == 0) seq_len(H) else h
          vs <- if (v == 0) seq_len(m) else v
          cols <- as.vector(outer(hs, (vs - 1L) * H, "+"))
          errs <- if (target == "forecast") F else G
          loss <- .pr_row_loss(errs[[method]][, cols, drop = FALSE]); loss[!mask] <- NA_real_
          reference <- if (target == "forecast") .pr_row_loss(F$DGP[, cols, drop = FALSE]) else NULL
          append_row(target, as.data.frame(c(list(design = design, n = n, method = method,
            subset = subset, horizon = h, variable = v, attempted = attempted,
            reference = if (target == "forecast") "DGP" else NA_character_), .pr_loss(loss, reference, trim))))
        }
        if (method != "DGP") {
          for (j in seq_len(J)) for (shock in seq_len(m)) for (response in seq_len(m)) {
            col <- response + (shock - 1L) * m + (j - 1L) * m * m
            append_row("irf_components", as.data.frame(c(list(design = design, n = n,
              method = method, subset = subset, horizon = j, response = response, shock = shock),
              as.list(.pr_component(E[[method]][mask, col], truth[response, shock, j])))))
          }
        }
        for (v in seq_len(m)) for (h in seq_len(H)) {
          append_row("forecast_components", as.data.frame(c(list(design = design, n = n,
            method = method, subset = subset, horizon = h, variable = v),
            as.list(.pr_component(F[[method]][mask, h + (v - 1L) * H])))))
        }
        for (target in names(L)) for (reference in intersect(c("DGP", "ECH_oracle"), methods)) {
          if (target == "irf" && reference == "DGP") next
          if (target == "forecast_gap") next
          a <- L[[target]][[method]]; a[!mask] <- NA_real_
          b <- L[[target]][[reference]]
          if (subset == "qualified") b[!qualified[, reference]] <- NA_real_
          append_row("paired_comparisons", as.data.frame(c(list(design = design, n = n,
            method = method, subset = subset, target = target, reference = reference),
            .pr_loss(a, b, trim))))
        }
      }
    }
  }
  if (length(setdiff(files, used))) stop("Unrecognized checkpoints: ", paste(setdiff(files, used), collapse = ", "))
  out <- lapply(out, function(x) if (length(x)) do.call(rbind, x) else data.frame())
  out$completion <- do.call(rbind, cells)
  out$manifest <- manifest
  out$execution <- lapply(list.files(input_dir, "^execution-.*[.]rds$", full.names = TRUE), readRDS)
  out
}

.pr_labels <- c(DGP = "DGP", ECH_oracle = "ECH oracle", ECH_selected = "ECH selected",
  IHR = "IHR", HR1 = "HR1", NET = "NET", RBY = "RBY", DPF = "DPF", MLE = "MLE")
.pr_colors <- c(DGP = "#222222", ECH_oracle = "#0072B2", ECH_selected = "#56B4E9",
  IHR = "#D55E00", HR1 = "#CC79A7", NET = "#E69F00", RBY = "#009E73", DPF = "#797979", MLE = "#5C3D99")

.pr_num <- function(x, digits = 3L) {
  vapply(x, function(z) if (is.na(z)) "--" else if (!is.finite(z)) as.character(z) else
    if (z != 0 && (abs(z) >= 1e4 || abs(z) < 10^(-digits))) formatC(z, format = "e", digits = 2L)
    else formatC(z, format = "f", digits = digits), character(1))
}

.pr_heat_num <- function(x) {
  label <- .pr_num(x, 2L)
  scientific <- is.finite(x) & x != 0 & (abs(x) >= 1e4 | abs(x) < .01)
  label[scientific] <- formatC(x[scientific], format = "e", digits = 1L)
  label <- gsub("e[+]0*", "e", label)
  gsub("e-0*", "e-", label)
}

.pr_theme <- function() {
  ggplot2::theme_bw(base_size = 11, base_family = "sans") + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "#E5E7EB", linewidth = .25),
    strip.background = ggplot2::element_rect(fill = "#EEF2F6", colour = "white"),
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    plot.title = ggplot2::element_text(size = 17, face = "bold", colour = "#17354D"),
    plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 12)),
    plot.caption = ggplot2::element_text(size = 9, hjust = 0, colour = "#414B55", margin = ggplot2::margin(t = 12)),
    plot.title.position = "plot", plot.caption.position = "plot",
    legend.position = "bottom", panel.spacing = grid::unit(.8, "lines"),
    plot.margin = ggplot2::margin(14, 18, 14, 14))
}

.pr_prepare <- function(d, cfg) {
  d$Method <- factor(d$method, levels = c("DGP", cfg$methods), labels =
    unname(.pr_labels[c("DGP", cfg$methods)]))
  d$Design <- factor(d$design, levels = names(cfg$designs), labels =
    paste("Design", seq_along(cfg$designs)))
  d$Sample <- factor(d$n, levels = cfg$sample_sizes, labels = paste0("n = ", cfg$sample_sizes))
  d
}

.pr_tf_plot <- function(d, cfg, subset = "available", trimmed = FALSE, trim = .01) {
  d <- .pr_prepare(d[d$subset == subset & d$horizon == 0 & d$method != "DGP", ], cfg)
  d$value <- if (trimmed) d$trimmed_rmse else d$rmse
  d$lower <- pmax(0, d$rmse - d$rmse_mcse); d$upper <- d$rmse + d$rmse_mcse
  p <- ggplot2::ggplot(d, ggplot2::aes(n, value))
  if (!trimmed) p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
    width = diff(range(cfg$sample_sizes)) / 25, colour = "#8BA3B5", na.rm = TRUE)
  p + ggplot2::geom_line(colour = "#0072B2", linewidth = .6, na.rm = TRUE) +
    ggplot2::geom_point(colour = "#17354D", size = 1.8, na.rm = TRUE) +
    ggplot2::facet_grid(Design ~ Method, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = cfg$sample_sizes) +
    ggplot2::scale_y_continuous(labels = function(x) .pr_num(x, 2)) +
    ggplot2::labs(x = "Training sample size", y = "IRF RMSE",
      title = if (trimmed) "Transfer-function accuracy: trimmed comparison" else "Transfer-function accuracy",
      subtitle = paste(if (subset == "qualified") "Qualified estimates" else "All available estimates",
        sprintf("| reduced-form IRFs, lags 1-%d | %d planned replications per cell", cfg$irf_lags, cfg$nsim)),
      caption = if (trimmed) sprintf("Square root of the mean replication loss after %.1f%% trimming in each tail. This differs from ordinary RMSE.\nLines connect observed sample sizes; no convergence-rate curve is fitted.", 100 * trim) else
        "Untrimmed RMSE; whiskers show one Monte Carlo standard error (lower ends floored at zero).\nEntries and lags are pooled within each replication. Lines do not impose a convergence rate. MCSEs can be unstable with extreme losses.") + .pr_theme()
}

.pr_ratio_plot <- function(d, cfg, target = "forecast", subset = "available") {
  d <- d[d$subset == subset & d$method != "DGP", ]
  if (target == "irf") d <- d[d$horizon == 0, ] else d <- d[d$horizon > 0 & d$variable == 0, ]
  d <- .pr_prepare(d, cfg)
  d$log_ratio <- log10(d$relative_rmse)
  d$label <- .pr_heat_num(d$relative_rmse)
  d$Horizon <- factor(d$horizon, levels = seq_len(cfg$h), labels = paste0("h = ", seq_len(cfg$h)))
  p <- ggplot2::ggplot(d, ggplot2::aes(Method, Sample, fill = log_ratio)) +
    ggplot2::geom_tile(colour = "white", linewidth = .5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7) +
    ggplot2::scale_fill_gradient(low = "#F4F7FA", high = "#E69F00", na.value = "grey88",
      name = "log10(ratio)") +
    ggplot2::scale_y_discrete(limits = rev(levels(d$Sample))) +
    ggplot2::labs(x = NULL, y = NULL,
      title = paste(if (target == "irf") "IRF" else "Forecast", "RMSE relative to", if (target == "irf") "oracle echelon" else "the DGP forecast"),
      subtitle = paste(if (subset == "qualified") "Qualified estimates" else "All available estimates",
        "| each ratio uses common usable replications | 1 = equal RMSE"),
      caption = "Numbers are RMSE ratios, not averages of replication-wise ratios. Colour represents log10(ratio).\nMissing outputs remain in failure counts. Large finite errors are retained. See CSV tables for paired counts and Monte Carlo standard errors.") +
    .pr_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 8))
  if (target == "irf") p + ggplot2::facet_wrap(~Design, ncol = 2) else
    p + ggplot2::facet_grid(Horizon ~ Design)
}

.pr_forecast_plot <- function(d, cfg, design, subset = "available") {
  d <- .pr_prepare(d[d$design == design & d$subset == subset & d$horizon > 0 & d$variable > 0, ], cfg)
  d$Variable <- factor(d$variable, labels = paste0("y", sort(unique(d$variable))))
  p <- ggplot2::ggplot(d, ggplot2::aes(horizon, rmse, colour = method, group = method)) +
    ggplot2::geom_line(linewidth = .65, na.rm = TRUE) + ggplot2::geom_point(size = 1.3, na.rm = TRUE) +
    ggplot2::facet_grid(Sample ~ Variable, scales = "free_y") +
    ggplot2::scale_colour_manual(values = .pr_colors, breaks = c("DGP", cfg$methods),
      labels = .pr_labels[c("DGP", cfg$methods)], name = NULL, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_len(cfg$h)) +
    ggplot2::labs(x = "Forecast horizon", y = "Forecast RMSE",
      title = paste0("Forecast accuracy: ", levels(d$Design)[match(design, names(cfg$designs))]),
      subtitle = paste(if (subset == "qualified") "Qualified estimates | linear vertical axis" else
        "All available estimates | logarithmic vertical axis", "| untrimmed RMSE"),
      caption = "Errors are realized future observations minus forecasts. DGP uses true parameters and the same observed training history.\nQualified means stable AR, invertible MA and reported convergence where applicable; NET uses its IHR pilot's flag. Counts are in the diagnostic tables.") + .pr_theme()
  if (subset == "available") p + ggplot2::scale_y_log10(labels = scales::label_scientific(digits = 2)) else p
}

.pr_irf_plot <- function(d, cfg, design, response = 1L, shock = 1L, n = NULL) {
  d <- d[d$design == design & d$subset == "available" & d$response == response & d$shock == shock, ]
  if (!is.null(n)) d <- d[d$n == n, ]
  d <- .pr_prepare(d, cfg)
  p <- ggplot2::ggplot(d, ggplot2::aes(horizon, median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95), fill = "#D6E6F2", na.rm = TRUE) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), fill = "#86B6D7", na.rm = TRUE) +
    ggplot2::geom_line(colour = "#006BA4", linewidth = .6, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = truth), colour = "#222222", linetype = "dashed", linewidth = .55) +
    ggplot2::scale_x_continuous(breaks = seq_len(cfg$irf_lags)) +
    ggplot2::labs(x = "Response lag", y = "Reduced-form response coefficient",
      title = sprintf("Design %d: response y%d to innovation %d", match(design, names(cfg$designs)), response, shock),
      subtitle = paste("Sampling distribution of IRF estimates", if (!is.null(n)) paste0("| n = ", n) else "| all sample sizes"),
      caption = "Blue line: median. Dark band: 25th-75th percentiles. Light band: 5th-95th percentiles. Dashed line: population response.\nBands describe variation across Monte Carlo replications; they are not confidence intervals. Extreme losses remain in the RMSE summaries.") + .pr_theme()
  if (is.null(n)) p + ggplot2::facet_grid(Sample ~ Method) else p + ggplot2::facet_wrap(~Method, ncol = 4)
}

.pr_time_plot <- function(d, cfg) {
  t <- .pr_prepare(d[d$method != "DGP", ], cfg)
  ggplot2::ggplot(t, ggplot2::aes(Method, mean_fit_seconds, colour = method)) +
    ggplot2::geom_point(size = 2.7, show.legend = FALSE) + ggplot2::facet_grid(Sample ~ Design) +
    ggplot2::scale_colour_manual(values = .pr_colors) + ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "Mean elapsed seconds per attempted fit (log scale)",
      title = "Estimation time", subtitle = "Elapsed times measured inside simulation workers, including failed attempts",
      caption = "Points show mean per-fit elapsed times under the recorded parallel workload. Summing worker times does not give wall-clock duration.\nDGP filtering, forecast evaluation and checkpoint overhead are excluded. See timing.csv for medians, 90th percentiles and total worker times.") +
    .pr_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 8))
}

.pr_text_page <- function(title, lines, subtitle = NULL) {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(title, .045, .92, just = c("left", "top"),
    gp = grid::gpar(fontsize = 25, fontface = "bold", col = "#17354D"))
  if (!is.null(subtitle)) grid::grid.text(subtitle, .045, .85, just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, col = "#414B55"))
  grid::grid.text(paste(lines, collapse = "\n"), .045, .77, just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, lineheight = 1.5, col = "#222222"))
}

.pr_tex <- function(x) {
  vapply(as.character(x), function(s) {
    chars <- strsplit(s, "", fixed = TRUE)[[1L]]
    map <- c("\\" = "\\textbackslash{}", "_" = "\\_", "%" = "\\%", "&" = "\\&",
      "#" = "\\#", "$" = "\\$", "{" = "\\{", "}" = "\\}", "^" = "\\textasciicircum{}", "~" = "\\textasciitilde{}")
    hit <- chars %in% names(map); chars[hit] <- unname(map[chars[hit]])
    paste(chars, collapse = "")
  }, character(1))
}

.pr_table <- function(data, title, note, file, csv_data = data) {
  utils::write.csv(csv_data, sub("[.]tex$", ".csv", file), row.names = FALSE, na = "")
  lines <- c("\\begin{table}[p]\\centering", paste0("\\caption{", .pr_tex(title), "}"),
    "\\small", paste0("\\begin{tabular}{", paste(rep("r", ncol(data)), collapse = ""), "}\\toprule"),
    paste0(paste(.pr_tex(names(data)), collapse = " & "), " \\\\\\midrule"),
    apply(data, 1L, function(row) paste0(paste(.pr_tex(row), collapse = " & "), " \\\\")),
    "\\bottomrule\\end{tabular}", "\\par\\medskip", paste0("\\begin{minipage}{0.98\\linewidth}\\footnotesize ",
      .pr_tex(note), "\\end{minipage}"), "\\end{table}")
  writeLines(lines, file)
}

#' Produce figures, tables and a PDF report from saved paper checkpoints.
#' No simulation or estimation routines are called. Existing raw checkpoints
#' and root-level summary CSVs are read only. Re-running replaces report files.
report_varma_paper_simulation <- function(input_dir = "paper-results",
    report_dir = file.path(input_dir, "report"), trim = .01,
    allow_partial = FALSE, details = TRUE) {
  if (length(trim) != 1L || !is.finite(trim) || trim < 0 || trim >= .5)
    stop("trim must lie in [0, 0.5)")
  for (pkg in c("varma", "ggplot2")) if (!requireNamespace(pkg, quietly = TRUE))
    stop("Please install the required package: ", pkg)
  if (!"irf" %in% getNamespaceExports("varma")) stop("Install the simulation version of varma")
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(input_dir, "manifest.rds"))) stop("No manifest.rds in ", input_dir)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  report_dir <- normalizePath(report_dir, winslash = "/", mustWork = TRUE)
  if (tolower(report_dir) == tolower(input_dir)) stop("Use a separate report subdirectory")
  z <- .pr_collect(input_dir, trim, allow_partial)
  .pr_write_report(z, input_dir, report_dir, trim, details)
}

# Kept separate from checkpoint reading so the display can also be regenerated
# from the exact saved report-summary.rds after making graphical adjustments.
.pr_write_report <- function(z, input_dir, report_dir, trim = .01, details = TRUE) {
  # Loading the namespace also registers dim.varma for deserialized models.
  if (!requireNamespace("varma", quietly = TRUE)) stop("Please install varma")
  cfg <- z$manifest$config; complete <- all(z$completion$missing == 0L)
  affected_version <- identical(cfg$package_version,"0.1.0") &&
    any(cfg$methods %in% c("IHR","HR1","NET"))
  correction_note <- if (affected_version) c(
    "CORRECTION: this run used varma 0.1.0. Its IHR/HR1/NET forecasts can be corrupted",
    "by an inverse-MA residual recursion. Version 0.1.1 fixes this; old checkpoints are not repaired.",
    "The figures below describe the archived run. Refresh affected estimates before using forecast rankings.", "") else character()
  if (affected_version) warning(paste(correction_note[1:3],collapse=" "),call.=FALSE)
  for (subdir in c("figures", "tables", "data")) dir.create(file.path(report_dir, subdir), showWarnings = FALSE)
  saveRDS(z, file.path(report_dir, "report-summary.rds"))
  for (nm in names(z)[!names(z) %in% c("manifest", "execution")])
    utils::write.csv(z[[nm]], file.path(report_dir, "data", paste0(nm, ".csv")), row.names = FALSE, na = "")
  saveRDS(list(manifest = z$manifest, execution = z$execution, trim = trim,
    generated = Sys.time(), input_dir = input_dir, report_session = utils::sessionInfo()),
    file.path(report_dir, "report-metadata.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(report_dir, "session-info.txt"))
  message("Validated checkpoints. Writing figures and tables ...")
  plots <- list(); sizes <- list()
  add <- function(name, plot, width = 16, height = 9) {
    if (!complete) plot <- plot + ggplot2::labs(tag = "INCOMPLETE RUN")
    plots[[name]] <<- plot; sizes[[name]] <<- c(width, height)
    for (ext in c("pdf", "png")) ggplot2::ggsave(file.path(report_dir, "figures", paste0(name, ".", ext)),
      plot = plot, width = width, height = height, dpi = 160, bg = "white",
      device = if (ext == "pdf") grDevices::cairo_pdf else "png", limitsize = FALSE)
  }
  add("tf_consist", .pr_tf_plot(z$irf, cfg))
  add("tf_consist_trimmed", .pr_tf_plot(z$irf, cfg, trimmed = TRUE, trim = trim))
  add("tf_consist_qualified", .pr_tf_plot(z$irf, cfg, subset = "qualified"))
  add("irf_relative_rmse", .pr_ratio_plot(z$irf, cfg, "irf"), height = 7)
  add("forecast_relative_rmse", .pr_ratio_plot(z$forecast, cfg), height = 12)
  add("forecast_relative_rmse_qualified", .pr_ratio_plot(z$forecast, cfg, subset = "qualified"), height = 12)
  add("log_times", .pr_time_plot(z$timing, cfg))
  f <- z$failures[z$failures$method != "DGP", ]
  metrics <- c(nonconverged = "Nonconverged", forecast_failures = "Forecast unavailable",
    unstable_ar = "AR unstable", noninvertible_ma = "MA noninvertible", qualified = "Qualified")
  fd <- do.call(rbind, lapply(names(metrics), function(nm) {
    x <- f; x$rate <- 100 * x[[nm]] / x$attempted; x$diagnostic <- metrics[[nm]]; x
  }))
  fd <- .pr_prepare(fd, cfg); fd$diagnostic <- factor(fd$diagnostic, levels = unname(metrics))
  add("fit_diagnostics", ggplot2::ggplot(fd, ggplot2::aes(Method, Sample, fill = rate)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", rate)), size = 2.7) +
    ggplot2::facet_grid(diagnostic ~ Design) + ggplot2::scale_y_discrete(limits = rev(levels(fd$Sample))) +
    ggplot2::scale_fill_gradient(low = "#F4F7FA", high = "#E69F00", limits = c(0, 100), name = "Percent") +
    ggplot2::labs(x = NULL, y = NULL, title = "Numerical diagnostics", subtitle = "Percent of attempted replications; categories can overlap",
      caption = "Qualified requires stable AR, invertible MA and reported convergence when applicable. NET's convergence flag refers to its IHR pilot.\nHR1 and DPF have no iterative convergence requirement. Forecast availability also requires finite squared errors. These are diagnostics, not proof of a global optimum.") +
    .pr_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 8)), height = 12)
  if (nrow(z$selection) && "ECH_selected" %in% z$selection$method) {
    s <- .pr_prepare(z$selection[z$selection$method == "ECH_selected", ], cfg)
    s$panel <- paste0(s$Design, " | k = (", s$true_indices, ")")
    add("index_selection", ggplot2::ggplot(s, ggplot2::aes(n, correct_fraction)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, correct_fraction - 1.96 * fraction_mcse),
        ymax = pmin(1, correct_fraction + 1.96 * fraction_mcse)), width = 10, colour = "#8BA3B5") +
      ggplot2::geom_line(colour = "#0072B2") + ggplot2::geom_point(size = 2.3, colour = "#0072B2") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", 100 * correct_fraction)), vjust = -1.5, size = 3.2) +
      ggplot2::facet_wrap(~panel, nrow = 1) + ggplot2::scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, .2), labels = scales::label_percent()) +
      ggplot2::scale_x_continuous(breaks = cfg$sample_sizes) +
      ggplot2::labs(x = "Training sample size", y = "Correct vector / attempted replications",
        title = "Data-driven Kronecker index selection", subtitle = "Exact recovery of the population index vector",
        caption = "Whiskers: +/-1.96 binomial Monte Carlo standard errors (clipped to [0,1]). Unreturned indices remain in the denominator.\nAn index vector may be unavailable because the subsequent fit failed; this is not a pure selection-error probability. Full vector frequencies are in indices.csv.") + .pr_theme(), height = 5.5)
  }
  gap <- .pr_prepare(z$forecast_gap[z$forecast_gap$method != "DGP" &
    z$forecast_gap$subset == "qualified" & z$forecast_gap$horizon > 0 & z$forecast_gap$variable == 0, ], cfg)
  add("forecast_gap_qualified", ggplot2::ggplot(gap, ggplot2::aes(horizon, rmse, colour = Sample, group = Sample)) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::facet_grid(Design ~ Method, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq_len(cfg$h)) +
    ggplot2::scale_colour_manual(values = rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), length.out = length(cfg$sample_sizes))) +
    ggplot2::labs(x = "Forecast horizon", y = "RMSE of forecast minus DGP forecast",
      title = "Forecast estimation error", subtitle = "Qualified estimates | variables pooled within each replication | untrimmed",
      caption = "This measures distance from the true-parameter conditional forecast using the same observed history.\nIt separates forecast estimation error from unpredictable future innovations. Available-subset values, including extreme errors, are retained in forecast_gap.csv.") + .pr_theme())
  for (i in seq_along(cfg$designs)) {
    design <- names(cfg$designs)[i]
    if (!any(z$irf_components$design == design)) next
    add(paste0("fcst_rmse", i), .pr_forecast_plot(z$forecast, cfg, design))
    add(paste0("fcst_rmse", i, "_qualified"), .pr_forecast_plot(z$forecast, cfg, design, "qualified"))
    first_n <- min(z$irf_components$n[z$irf_components$design == design])
    add(paste0("irf_exp", i), .pr_irf_plot(z$irf_components, cfg, design, n = first_n), width = 13, height = 7)
    if (details) {
      grDevices::cairo_pdf(file.path(report_dir, "figures", paste0("irf_all_components_design", i, ".pdf")), width = 16, height = 8)
      tryCatch({
        m <- dim(cfg$designs[[design]])[1L]
        for (response in seq_len(m)) for (shock in seq_len(m)) {
          p <- .pr_irf_plot(z$irf_components, cfg, design, response, shock)
          if (!complete) p <- p + ggplot2::labs(tag = "INCOMPLETE RUN")
          print(p)
        }
      }, finally = grDevices::dev.off())
    }
  }
  table_files <- character()
  for (i in seq_along(cfg$designs)) for (subset in c("available", "qualified")) {
    design <- names(cfg$designs)[i]
    d <- z$irf[z$irf$design == design & z$irf$subset == subset & z$irf$horizon == 0 & z$irf$method != "DGP", ]
    if (!nrow(d)) next
    tab <- data.frame(n = d$n, Method = .pr_labels[d$method], B = d$n_available,
      RMSE = .pr_num(d$rmse), MCSE = .pr_num(d$rmse_mcse), Trimmed = .pr_num(d$trimmed_rmse),
      `Relative RMSE` = .pr_num(d$relative_rmse), `Paired B` = d$n_paired, check.names = FALSE)
    file <- paste0("irf_design", i, "_", subset, ".tex")
    .pr_table(tab, paste0("Design ", i, ": IRF accuracy (", subset, ")."),
      sprintf("Positive lags 1-%d and all response entries are pooled within replication. RMSE is untrimmed; Trimmed removes %.1f%% of losses in each tail. MCSE is a delta-method Monte Carlo standard error, potentially unstable for heavy-tailed losses. Relative RMSE compares with oracle echelon on common usable draws. B counts usable draws. See diagnostic tables for exclusions.", cfg$irf_lags, 100 * trim),
      file.path(report_dir, "tables", file), d)
    table_files <- c(table_files, file)
    for (h in unique(c(1L, cfg$h))) {
      d <- z$forecast[z$forecast$design == design & z$forecast$subset == subset &
        z$forecast$horizon == h & z$forecast$variable == 0, ]
      tab <- data.frame(n = d$n, Method = .pr_labels[d$method], B = d$n_available,
        RMSE = .pr_num(d$rmse), MCSE = .pr_num(d$rmse_mcse),
        `Relative RMSE` = .pr_num(d$relative_rmse), `Ratio MCSE` = .pr_num(d$relative_rmse_mcse),
        `Paired B` = d$n_paired, check.names = FALSE)
      file <- paste0("forecast_design", i, "_h", h, "_", subset, ".tex")
      .pr_table(tab, paste0("Design ", i, ": forecast accuracy at horizon ", h, " (", subset, ")."),
        "Variables are pooled within replication. RMSE is untrimmed. Relative RMSE uses DGP forecasts on the same usable replications. MCSEs use independent replication losses and a delta approximation; they can be unstable with extreme losses. Qualified requires stability, invertibility and applicable convergence. B and Paired B disclose the denominators.",
        file.path(report_dir, "tables", file), d)
      table_files <- c(table_files, file)
    }
    if (subset == "available") {
      d <- z$failures[z$failures$design == design & z$failures$method != "DGP", ]
      tab <- data.frame(n = d$n, Method = .pr_labels[d$method], Attempts = d$attempted,
        `Fit fail` = d$fit_failures, `IRF fail` = d$irf_failures, `Fcst fail` = d$forecast_failures,
        `Nonconv.` = d$nonconverged, `AR unstable` = d$unstable_ar,
        `MA noninv.` = d$noninvertible_ma, Qualified = d$qualified, check.names = FALSE)
      file <- paste0("diagnostics_design", i, ".tex")
      .pr_table(tab, paste0("Design ", i, ": numerical diagnostics."),
        "Counts refer to attempted replications. Categories overlap. NET's convergence field concerns its IHR pilot. HR1 and DPF have no iterative convergence requirement. Forecast failures include nonfinite squared errors; finite extreme outputs remain available. Convergence does not certify a global maximum.",
        file.path(report_dir, "tables", file), d)
      table_files <- c(table_files, file)
    }
  }
  tex <- c("% Standalone tables. Requires booktabs. Copy this file and tables/ together.",
    "\\documentclass[11pt]{article}", "\\usepackage[a4paper,landscape,margin=16mm]{geometry}", "\\usepackage{booktabs}",
    "\\begin{document}", "\\section*{VARMA simulation tables}",
    paste0("Run status: ", if (complete) "complete" else "INCOMPLETE", ". ",
      "Results are conditional on available outputs or on the explicitly qualified subset."),
    paste0("\\input{tables/", table_files, "}"), "\\end{document}")
  writeLines(tex, file.path(report_dir, "paper_tables.tex"))
  total <- sum(z$completion$attempted)
  wall <- sum(vapply(z$execution, function(e) e$seconds, numeric(1)))
  workers <- unique(vapply(z$execution, function(e) e$workers, integer(1)))
  aggregates <- stats::aggregate(cbind(fit_failures, irf_failures, forecast_failures,
    nonconverged, unstable_ar, noninvertible_ma, qualified) ~ method, z$failures, sum)
  utils::write.csv(aggregates, file.path(report_dir, "tables", "diagnostics_overall.csv"), row.names = FALSE)
  highlights <- c(correction_note,
    sprintf("%s: %s saved replications; %s attempted estimates; %d designs.",
      if (complete) "Complete run" else "INCOMPLETE RUN", format(total, big.mark = ","),
      format(total * length(cfg$methods), big.mark = ","), length(cfg$designs)),
    sprintf("Sample sizes: %s. Planned replications per cell: %d. Seed: %s.",
      paste(cfg$sample_sizes, collapse = ", "), cfg$nsim, cfg$seed),
    sprintf("Recorded compute/checkpoint elapsed time: %.2f hours; worker counts: %s.", wall / 3600,
      paste(workers, collapse = ", ")),
    paste0("MLE backend: ", cfg$mle_backend, "; the echelon estimators use conditional Gaussian likelihood."),
    "", "Primary results retain every finite loss, including extreme outputs and nonconverged estimates.",
    "Qualified results require stable AR, invertible MA and reported convergence where applicable.",
    "NET's convergence flag concerns its IHR pilot; HR1 and DPF have no iterative flag.",
    sprintf("The paper-style trimmed comparison removes %.1f%% of replication losses from each tail.", trim * 100),
    "Neither availability nor qualification makes these unconditional risks of an algorithm that can fail.",
    "", "IRF bands are Monte Carlo sampling quantiles, not confidence intervals.",
    "MCSEs use independent replication-level losses. Heavy tails can make them unreliable.",
    "Relative errors pair each method with its benchmark on common usable draws.",
    "", "This run evaluates IRFs and point forecasts. It does not evaluate FEVDs, bootstrap coverage,",
    "prediction-interval coverage, or fitted VAR benchmarks. These require additional experiments.")
  grDevices::cairo_pdf(file.path(report_dir, "paper_report.pdf"), width = 16, height = 10)
  tryCatch({
    .pr_text_page("VARMA simulation results", highlights, paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")))
    selection <- z$selection[z$selection$method == "ECH_selected" & z$selection$n == max(cfg$sample_sizes), ]
    diagnostic_lines <- vapply(seq_len(nrow(aggregates)), function(i) {
      a <- aggregates[i, ]
      sprintf("%s: %d forecast failures; %d nonconverged; %d qualified / %d attempted.",
        .pr_labels[a$method], a$forecast_failures, a$nonconverged, a$qualified, total)
    }, character(1))
    selection_lines <- if (nrow(selection)) vapply(seq_len(nrow(selection)), function(i) {
      s <- selection[i, ]
      sprintf("Design %d: %.1f%% exact index recovery at n = %d (true vector: %s).",
        match(s$design, names(cfg$designs)), 100 * s$correct_fraction, s$n, s$true_indices)
    }, character(1)) else "No data-driven index-selection results in this run."
    .pr_text_page("Diagnostics and index selection", c(diagnostic_lines, "", selection_lines, "",
      "Counts refer to all saved replications; categories overlap. NET uses its IHR pilot's convergence flag.",
      "The figures and tables include both all available outputs and the separately qualified subset."))
    for (name in names(plots)) print(plots[[name]])
  }, finally = grDevices::dev.off())
  lines <- c("# VARMA simulation report", "", highlights, "", "## Files", "",
    "- `paper_report.pdf`: overview and the main figures.",
    "- `figures/`: vector PDFs and PNG previews. `tf_consist` and `log_times` match the paper's topics.",
    "- `irf_exp1` through `irf_exp4`: y1 response to innovation 1 at the smallest sample size, as in the slides.",
    "- `irf_all_components_design*.pdf`: every response/innovation pair at every sample size.",
    "- `fcst_rmse*`: every forecast variable and sample size; suffix `_qualified` denotes the diagnostic subset.",
    "- `data/`: freshly computed statistics, component quantiles, diagnostics and replication-level losses.",
    "- `tables/`: numeric CSVs and LaTeX tables for IRFs, forecasts at first/last horizon, and diagnostics.",
    "- `paper_tables.tex`: standalone LaTeX document; compile from this report directory (requires booktabs).",
    "- `report-metadata.rds`: original run configuration/execution metadata plus reporting session.",
    "", "## Definitions", "",
    "IRF errors are reduced-form estimated coefficients minus truth. Lag zero (the identity) is not scored.",
    "In IRF summaries, horizon=0 means pooling all positive lags. In forecast summaries, horizon=0 means",
    "pooling all forecast horizons; variable=0 means pooling variables. Pooling always occurs within",
    "replication before computing Monte Carlo means/standard errors, not by treating entries as independent.",
    "", "Forecast errors are realized future observations minus forecasts. Forecast gaps are fitted forecasts",
    "minus the true-parameter conditional forecast based on the same observed training history.",
    "Forecast relative RMSE = sqrt(paired method MSE / paired DGP MSE). IRF reference = oracle echelon.",
    "paired_comparisons.csv additionally compares pooled forecast losses with oracle echelon.",
    "These are finite-sample procedure comparisons, not estimates of an asymptotic efficiency penalty",
    "attributable solely to non-identification. The model restrictions and estimation algorithms also differ.",
    "", "The coefficient quantiles describe the estimates (truth is added back to saved IRF errors).",
    "The bias column is mean estimation error. Forecast component quantiles describe forecast errors.",
    "top_1pct_loss_share is the share of total loss from the largest ceiling(0.01*B) usable losses.",
    "A large share flags sensitivity to a small number of replications; the observations are not discarded.",
    "", "Means and variances are computed after rescaling to avoid numerical overflow with very large finite",
    "losses. The summary CSVs originally saved by the runner are not overwritten. Raw checkpoint status",
    "flags determine availability, consistently with the runner; a failed full forecast path is not partly salvaged.",
    "", "Timing uses elapsed seconds measured inside workers, not user CPU time. Summed worker times are",
    "not wall-clock run duration. Execution logs exclude summary writing and process cleanup; across resumed",
    "runs the reported elapsed total is the sum of logged invocations, not calendar time including pauses.",
    "", "## Overall diagnostics (counts)", "", "```", capture.output(print(aggregates, row.names = FALSE)), "```",
    "", "## Exact index recovery at the largest sample size", "", selection_lines,
    "", "## Reproduce", "", "```r", "source(\"inst/scripts/report_paper_simulations.R\")",
    "report_varma_paper_simulation(\"paper-results\")", "```", "",
    "Only the report directory is replaced on another reporting run. No models are refitted.",
    "The historical manuscript/slide figures are untouched. Review captions before inserting refreshed figures.")
  writeLines(lines, file.path(report_dir, "README.md"))
  message("Report written to ", report_dir)
  invisible(list(report_dir = report_dir, summary = z, figures = names(plots), complete = complete))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || "--help" %in% args) {
    cat("Usage: Rscript report_paper_simulations.R <results-dir> [report-dir] [--allow-partial] [--no-details]\n",
      "Example: Rscript inst/scripts/report_paper_simulations.R paper-results\n",
      "No simulation or model estimation is run. Requires varma and ggplot2.\n", sep = "")
  } else {
    flags <- args[startsWith(args, "--")]
    if (any(!flags %in% c("--allow-partial", "--no-details"))) stop("Unknown option: ", paste(flags, collapse = " "))
    pos <- args[!startsWith(args, "--")]
    if (!length(pos) || length(pos) > 2L) stop("Supply a results directory and, optionally, a report directory")
    report_varma_paper_simulation(pos[1L], if (length(pos) == 2L) pos[2L] else file.path(pos[1L], "report"),
      allow_partial = "--allow-partial" %in% flags, details = !"--no-details" %in% flags)
  }
}
