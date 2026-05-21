library(tidyverse)
library(varma)
library(MTS)


################################################################################
# FUNCTIONS FOR THE SIMULATION EXPERIMENTS                                     #
################################################################################

#' Estimate VARMA in Echelon form exploiting the package MTS
#'
#' @param Y matrix with time series to predict (time in rows).
#' @param intercept estimate intercept to allow non-zero mean processes.
#' @param plag number of lags used to represent the past vector. Default is 5.
#' @param crit type-I error rate used in testing for zero canonical correlations. Deafult is 0.05
#' @param kind Kronecker indices: if NULL Kronid is called with \code{plag} and \code{crit}
#'
#' @returns An object of class \code{varma}.
fit_varma_ech <- function(Y, intercept = TRUE, plag = 5, crit = 0.05, kind = NULL) {
  fn_name <- as.character(sys.call()[[1]])
  if (is.null(kind)) kind <- MTS::Kronid(Y, plag, crit)$index
  kfit <- MTS::Kronfit(Y, kind, intercept)
  m <- dim(kfit$Sigma)[1]
  nobs <- prod(dim(kfit$residuals))
  iPhi0 <- solve(kfit$Ph0)
  out <-  varma(intercept = kfit$const,
                ar = array(iPhi0 %*% kfit$Phi, c(m, m, ncol(kfit$Phi)/m)),
                ma = array(-iPhi0 %*% kfit$Theta, c(m, m, ncol(kfit$Theta)/m)),
                cov = kfit$Sigma,
                estimation_method = fn_name,
                loglik = - (nobs / 2) * (log(2 * pi) + 1) -
                  (nrow(kfit$residuals) / 2) * determinant(kfit$Sigma, logarithm = TRUE)$modulus,
                n = nrow(kfit$residuals),
                nobs = nobs,
                npar = length(kfit$coef),
                y = Y,
                residuals = kfit$residuals)
  out$residuals <- residuals(out)
  out
}



#' Simulation experiment
#'
#' It simulates VARMA time series, estimate the parameters of the
#' VARMA using different methods and evaluate the distance of the
#' estimated transfer function from the population tranfer function,
#' and the accuracy of the forecasts
#'
#' @param varma a VARMA model passed as \code{varma} object
#' @param n length of the time series to simulate
#' @param h number of steps for the forecasts and of lags for the transfer fn
#' @param nsim number of simulations (replicates)
#'
#'
experiment <- function(varma, n = 100, h = 5, nsim = 1000) {
  mpq <- dim(varma)
  id  <- is_identified(varma)
  descr <- paste0(
    if (id) "Identifiable " else "Unidentifiable ",
    mpq[1], "D VARMA(",
    mpq[2], ", ",
    mpq[3], "), n = ", n, ", h = ", h
  )

  # simulating sample paths
  list_y <- replicate(nsim, sim_varma(varma, n = n+h), simplify = FALSE)
  list_est <- vector("list", 7)
  names(list_est) <- c("ECH", "IHR", "HR1", "NET", "RBY", "DPF", "MLE")
  list_time <- vector("list", 7)
  names(list_time) <- names(list_est)

  # estimating a VARMA using seven methods on each simulated path
  # also storing the execution time
  cat("ECH estimations\n")
  kind <- kronecker_indices(varma)
  list_time[["ECH"]] <- system.time(
    list_est[["ECH"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_ech(Y = x[1:n,],
                                                                   intercept = F,
                                                                   kind = kind),
                                                     error = function(e) NA)
    )
  )
  cat("IHR estimations\n")
  list_time[["IHR"]] <- system.time(
    list_est[["IHR"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_ihr(Y = x[1:n,],
                                                                   p = mpq[2],
                                                                   q = mpq[3],
                                                                   intercept = F),
                                                     error = function(e) NA)
    )
  )
  cat("HR1 estimations\n")
  list_time[["HR1"]] <- system.time(
    list_est[["HR1"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_ihr(Y = x[1:n,],
                                                                   p = mpq[2],
                                                                   q = mpq[3],
                                                                   intercept = F,
                                                                   maxit = 1),
                                                     error = function(e) NA)
    )
  )
  cat("NET estimations\n")
  list_time[["NET"]] <- system.time(
    list_est[["NET"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_net(Y = x[1:n,],
                                                                   p = mpq[2],
                                                                   q = mpq[3],
                                                                   intercept = F),
                                                     error = function(e) NA)
    )
  )
  cat("RBY estimations\n")
  list_time[["RBY"]] <- system.time(
    list_est[["RBY"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_rby(Y = x[1:n,],
                                                                   p = mpq[2],
                                                                   q = mpq[3],
                                                                   intercept = F),
                                                     error = function(e) NA)
    )
  )
  cat("DPF estimations\n")
  list_time[["DPF"]] <- system.time(
    list_est[["DPF"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_dpf(Y = x[1:n,],
                                                                   p = mpq[2]+(mpq[1]-1)*mpq[3],
                                                                   q = mpq[1]*mpq[3],
                                                                   intercept = F),
                                                     error = function(e) NA)
    )
  )
  cat("MLE estimations\n")
  list_time[["MLE"]] <- system.time(
    list_est[["MLE"]] <- lapply(list_y,
                                function(x) tryCatch(fit_varma_fkf(Y = x[1:n,],
                                                                   p = mpq[2],
                                                                   q = mpq[3],
                                                                   intercept = F),
                                                     error = function(e) NA)
    )
  )

  ## Computing the transfer function for each estimated VARMA
  list_irf <- lapply(
    X = list_est,
    FUN = function(y) {
      lapply(X = y,
             FUN = function(x) if (!anyNA(x)) irf(x, maxlag = h, orth = "none") else NULL)
    }
  )
  ## Computing the MSE of the transfer functions
  irf_mod <- irf(varma, h, "none", plot = FALSE)
  list_irf_mse <- lapply(
    X = list_irf,
    FUN = function(y) {
      sapply(X = y,
             FUN = function(x) if (!is.null(x)) {
               tryCatch(irf_distance(x, irf_mod), error = function(e) NA)
             }
      )
    }
  )
  ## Computing the h-step-ahead forecasts
  list_pre <- lapply(
    X = list_est,
    FUN = function(y) {
      lapply(X = y,
             FUN = function(x) {
               if (!anyNA(x)) {
                 tryCatch(predict(x, n.ahead = h, cov = FALSE), error = function(e) NA)
               } else NULL
             }
      )
    }
  )
  ## Computing the h-step-ahead forecast errors
  list_pre_err <- lapply(
    X = list_pre,
    FUN = function(x) {
      if (is.null(x)) {
        NULL
      } else {
        mapply(function(y, x) {y[(n+1):(n+h), ] - x}, list_y, x, SIMPLIFY = "array")
      }
    }
  )

  # Boxplot of the MSE of the Transfer Function estimation
  tryCatch({
    list_irf_mse |>
      lapply(unlist) |>
      boxplot(outline = FALSE, ylab = "MSE",
              main = paste(descr, ": TF MSE"))
  }, error = function(e) print(e))
  # Boxplot of the one-step-ahead forecast error for y1
  old <- par(mfcol = c(mpq[1], 2),
             mar = c(2, 4, 2, 0.5),
             oma = c(1, 2, 2, 1))
  on.exit(par(old))
  tryCatch({
    for (i in 1:mpq[1]){
      list_pre_err |>
        lapply(function(x) if (is.null(x)) NA else x[1, i, ]) |>
        boxplot(outline = FALSE, ylab = "1-step", main = paste0("y", i))
      abline(h = 0, lty = 2)
      list_pre_err |>
        lapply(function(x) if (is.null(x)) NA else x[h, i, ]) |>
        boxplot(outline = FALSE, ylab = paste0(h, "-step"),
                main = paste0("y", i))
      abline(h = 0, lty = 2)
    }
    mtext(paste(descr, ": forecast errors"), outer = TRUE)
  }, error = function(e) print(e))

  list(
    description = descr,
    dgp = varma,
    identified = id,
    n = n,
    h = h,
    nsim = nsim,
    y = list_y,
    estimates = list_est,
    irfs = list_irf,
    irf_mse = list_irf_mse,
    pre = list_pre,
    pre_err = list_pre_err,
    time = list_time
  )
}


###################################################################
# SIMULATION EXPERIMENTS                                          #
###################################################################

# The results will be organized in a list of lists with the following structure
# The first-level slot is a list with the results of each experiment: experiment_#.
# Each method in each experiment is structured as follows
# results$experiment_# <- list(
#   description = string that describes the experiment
#   dgp         = varma object data generating process
#   y           = list of simulated time series: each slot is a sample path,
#   estimates   = list of est. method with list of estimated varma objects: one for each simulated path (can be NA),
#   irfs        = list of est. method with list of IRFs: one for each estimated varma object
#   irf_mse     = list of est. method with list of MSE for the IRF,
#   pre         = list of est. method with list of forecasts
#   pre_err     = list of est. method with 3D-array of forecast errors
#                 dim: forecast horizons x variables x simulations
# )


# Models definition
## mod1: 2-dim unidentifiable VARMA(1, 1) with zero roots (nilpotent)
set.seed(1)
mod1 <- rvarma(2, 1, 1)
mod1$ar[c(1, 2, 4)] <- 0
mod1$ma[c(1, 2, 4)] <- 0
mod1$cov[1, 2] <- mod1$cov[2, 1] <- 0.3

is_identified(mod1)
inv_roots(mod1, plot = TRUE)
irf_mod1 <- irf(mod1, maxlag = 5, plot = TRUE)
mod1

## mod2: 2-dim identifiable VARMA(1, 1)
set.seed(2)
mod2 <- rvarma(2, 1, 1)
mod2$cov[1, 2] <- mod2$cov[2, 1] <- 0.3

is_identified(mod2)
inv_roots(mod2, plot = TRUE)
irf_mod2 <- irf(mod2, maxlag = 5, plot = TRUE)
mod2

## mod3: 2-dim unidentifiable VARMA(2, 2)
## building a VARMA(2, 2) with a common root
set.seed(3)
tmp1 <- rvarma(2, 1, 0, dist = runif)
tmp2 <- rvarma(2, 1, 1, dist = runif)
mod3 <- rvarma(2, 2, 2)
mod3$ar[,, 1] <- (tmp1$ar[,,1] + tmp2$ar[,,1])
mod3$ar[,, 2] <- -(tmp1$ar[,,1] %*% tmp2$ar[,,1])
mod3$ma[,, 1] <- -(tmp1$ar[,,1] + tmp2$ma[,,1])
mod3$ma[,, 2] <- tmp1$ar[,,1] %*% tmp2$ma[,,1]
mod3$cov[1, 2] <- mod3$cov[2, 1] <- 0.3

is_identified(mod3)
inv_roots(mod3, plot = TRUE)
irf_mod3 <- irf(mod3, maxlag = 5, plot = TRUE)
mod3

## mod4: 3-dim identifiable VARMA(2, 2)
## building a VARMA(2, 2)
set.seed(4)
mod4 <- rvarma(3, 2, 2)
mod4$cov[1, 2] <- mod4$cov[2, 1] <- 0.3
mod4$cov[1, 3] <- mod4$cov[3, 1] <- 0.2
mod4$cov[2, 3] <- mod4$cov[3, 2] <- -0.1

is_identified(mod4)
inv_roots(mod4, plot = TRUE)
irf_mod4 <- irf(mod4, maxlag = 5, plot = TRUE)
mod4

# Simulations
results <- list()

Sys.time()
results$experiment_1 <- experiment(mod1, n = 100, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_2 <- experiment(mod2, n = 100, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_3 <- experiment(mod3, n = 100, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_4 <- experiment(mod4, n = 100, h = 5, nsim = 1000)
save(results, file = "results.Rdata")

Sys.time()
results$experiment_5 <- experiment(mod1, n = 200, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_6 <- experiment(mod2, n = 200, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_7 <- experiment(mod3, n = 200, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_8 <- experiment(mod4, n = 200, h = 5, nsim = 1000)
save(results, file = "results.Rdata")

Sys.time()
results$experiment_9  <- experiment(mod1, n = 400, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_10 <- experiment(mod2, n = 400, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_11 <- experiment(mod3, n = 400, h = 5, nsim = 1000)
save(results, file = "results.Rdata")
Sys.time()
results$experiment_12 <- experiment(mod4, n = 400, h = 5, nsim = 1000)
save(results, file = "results.Rdata")


# functions to extract contents from the simulations
tf_rmse_plot <- function(experiment, outliers = FALSE, col = "lightblue") {
  irf_mse <- lapply(experiment$irf_mse, unlist)
  boxplot(irf_mse,
          ylab = "MSE", main = experiment$description,
          outline = outliers,
          col = col)
  cat("Number of succesful estimates:\n")
  print(sapply(irf_mse, length))
  sapply(irf_mse, summary)
}



tf_rmse_ggplot <- function(res, outliers = FALSE, models = NA) {
  dt <- data.frame()
  for (i in 1:length(res)) {
    tmp <- lapply(res[[i]]$irf_mse, unlist)
    wdt <- data.frame()
    for (j in 1:length(tmp)) {
      wdt <- rbind(wdt, data.frame(Estimate = names(tmp)[j],
                                   RMSE = unlist(tmp[[j]]))
      )
    }
    wdt$n <- res[[i]]$n
    wdt$DGP <- if (is.null(models)) res[[i]]$description else models[i]
    dt <- rbind(dt, wdt)
  }

  dt |>
    group_by(n, DGP, Estimate) |>
    summarise(MedianRMSE = median(RMSE, na.rm = TRUE)) |>
    group_by(n, DGP) |>
    mutate(Rank = min_rank(MedianRMSE)) |>
    ungroup() ->
    dt_median

  dt |>
    left_join(dt_median, by = c("DGP" = "DGP", "Estimate" = "Estimate", "n" = "n")) |>
    ggplot(aes(x = Estimate, y = RMSE, fill = DGP, alpha = Rank)) +
    scale_alpha(range = c(1, 0.2)) +
    geom_boxplot(outliers = outliers, show.legend = FALSE) +
    facet_grid(rows = vars(DGP), cols = vars(n), scales = "free_y") ->
    gg1

  dt |>
    group_by(Estimate, n, DGP) |>
    summarise(RMSE = sqrt(mean(RMSE^2))) |>
    ggplot(aes(x = n, y = RMSE, color = DGP)) +
    geom_smooth(formula = y~I(sqrt(1/x)), method = "lm", se = FALSE, show.legend = FALSE) +
    geom_point(show.legend = FALSE, color = "black") +
    facet_grid(rows = vars(DGP), cols = vars(Estimate), scales = "free_y") ->
    gg2
  return(list(gg1, gg2))
}


fcst_rmse_ggplot <- function(res, outliers = FALSE, models = NA) {
  dt <- data.frame()
  for (i in 1:length(res)) {
    tmp <- res[[i]]$pre_err
    wdt <- data.frame()
    for (j in 1:length(tmp)) {
      tmp_arr <- if (is.array(tmp[[j]])) {
        tmp[[j]]
      } else { # for some reasons sometimes the fcsts over the sims are in a list of matrices
        tmp[[j]] |>
          unlist() |>
          array(c(res[[i]]$h, dim(res[[i]]$dgp)[1], length(tmp[[j]])))
      }
      tmp_mse <- apply(tmp_arr, 2, function(x) apply(x^2, 1, mean, trim = 0.01, na.rm = TRUE))
      colnames(tmp_mse) <- paste0("y", 1:ncol(tmp_mse))
      wdt <- rbind(wdt,
                   data.frame(Estimate = names(tmp)[j],
                              MSE = as.numeric(tmp_mse),
                              Horizon = 1:nrow(tmp_mse),
                              Variable = rep(colnames(tmp_mse), each = nrow(tmp_mse))
                   )
      )
    }
    wdt$n <- res[[i]]$n
    wdt$DGP <- if (is.null(models)) res[[i]]$description else models[i]
    dt <- rbind(dt, wdt)
  }
  dt$RMSE <- sqrt(dt$MSE)
  dt |>
    ggplot(aes(x = Horizon, y = RMSE, color = Variable)) +
    geom_line() +
    facet_grid(rows = vars(DGP, n), cols = vars(Estimate), scales = "free_y") +
    scale_y_log10() ->
    gg1

  dt |>
    ggplot(aes(x = Estimate, y = RMSE, fill = Variable)) +
    geom_col() +
    facet_grid(rows = vars(DGP, n), cols = vars(Horizon), scales = "free_y") ->
    gg2

  dt |>
    filter(Horizon == 1) |>
    group_by(Variable, n, Horizon, DGP) |>
    mutate(Rank = min_rank(RMSE)) |>
    ungroup() |>
    ggplot(aes(x = Estimate, y = RMSE, fill = DGP, alpha = Rank)) +
    geom_col(color = "black") +
    scale_alpha(range = c(1, 0.2)) +
    ggh4x::facet_nested(rows = vars(DGP, n), cols = vars(Variable), scales = "free_y") +
    ggtitle(paste0("Forecast Horizon: ", 1, "-step-ahead")) ->
    gg3

  dt |>
    filter(Horizon == max(Horizon)) |>
    group_by(Variable, n, Horizon, DGP) |>
    mutate(Rank = min_rank(RMSE)) |>
    ungroup() |>
    ggplot(aes(x = Estimate, y = RMSE, fill = DGP, alpha = Rank)) +
    geom_col() +
    scale_alpha(range = c(1, 0.2)) +
    ggh4x::facet_nested(rows = vars(DGP, n), cols = vars(Variable), scales = "free_y") +
    ggtitle(paste0("Forecast Horizon: ", max(dt$Horizon, na.rm = TRUE), "-step-ahead")) ->
    gg4

  dt |>
    filter(Horizon == 1) |>
    group_by(Variable, n, Horizon, DGP) |>
    mutate(`Rel. RMSE` = RMSE/RMSE[Estimate == "MLE"],
           Rank = min_rank(RMSE)) |>
    ungroup() |>
    ggplot(aes(x = Estimate, y = `Rel. RMSE`, fill = DGP, alpha = Rank)) +
    geom_col() +
    scale_alpha(range = c(1, 0.2)) +
    ggh4x::facet_nested(rows = vars(DGP, n), cols = vars(Variable), scales = "free_y") +
    ggtitle(paste0("Forecast Horizon: ", 1, "-step-ahead")) ->
    gg5

  dt |>
    filter(Horizon == max(Horizon)) |>
    group_by(Variable, n, Horizon, DGP) |>
    mutate(`Rel. RMSE` = RMSE/RMSE[Estimate == "MLE"],
           Rank = min_rank(RMSE)) |>
    ungroup() |>
    ggplot(aes(x = Estimate, y = `Rel. RMSE`, fill = DGP, alpha = Rank)) +
    geom_col() +
    scale_alpha(range = c(1, 0.2)) +
    ggh4x::facet_nested(rows = vars(DGP, n), cols = vars(Variable), scales = "free_y") +
    ggtitle(paste0("Forecast Horizon: ", max(dt$Horizon, na.rm = TRUE), "-step-ahead")) ->
    gg6

  list(gg1, gg2, gg3, gg4, gg5, gg6)
}

load("results.Rdata")

tf_plts <- tf_rmse_ggplot(results, models = rep(paste("Model", 1:4), 3))
tf_plts[[1]]
tf_plts[[2]]

fcst_plts <- fcst_rmse_ggplot(results, models = rep(paste("Model", 1:4), 3))
fcst_plts[[1]]
fcst_plts[[2]]
fcst_plts[[3]]
fcst_plts[[4]]
fcst_plts[[5]]
fcst_plts[[6]]


irf_quantiles <- function(experiment, h, perc = c(0.1, 0.5, 0.9)) {
  lapply(experiment$irfs, function(x)
         vapply(x, identity,
                FUN.VALUE = array(0, c(dim(experiment$irfs[[1]][[1]]), length(dim(experiment$irfs[[1]]))))
         )
  )
}



# NEW



#' Plot the simulated distribution of IRF estimates around the true IRF
#'
#' Returns a named list of \code{m * m} ggplot objects, one per entry \code{(i,j)}
#' of the IRF matrix.  Within each plot the facets correspond to estimation
#' methods, making cross-method comparison straightforward.
#'
#' @param sim Named list, one entry per estimation method.  Each entry
#'   is a list of \code{n_sim} arrays of dimension \code{(m, m, n_lags)}, where
#'   \code{m} is the number of variables and \code{n_lags} the number of lags
#'   (e.g. 6 for lags 0..5).
#' @param probs       Length-4 strictly-increasing numeric vector in \code{(0,1)}.
#'   \code{[probs[1], probs[4]]} is the outer (light) band;
#'   \code{[probs[2], probs[3]]} is the inner (dark) band.
#'   Default: \code{c(0.05, 0.25, 0.75, 0.95)}.
#' @param alpha_inner Opacity of the inner quantile band.  Default 0.35.
#' @param alpha_outer Opacity of the outer quantile band.  Default 0.15.
#' @param colors      Optional named character vector of colors, one per method.
#'   Each method's facet is drawn in its own color, giving a consistent palette
#'   across all returned plots.  If \code{NULL} ggplot2's default palette is used.
#' @param true_color  Color of the true IRF line.  Default \code{"black"}.
#' @param true_lwd    Line width of the true IRF line.  Default 0.9.
#' @param var_names   Character vector of length \code{m} labelling the
#'   responding variables.  Defaults to \code{"y1", "y2", ...}.
#' @param shock_names Character vector of length \code{m} labelling the shocks.
#'   Defaults to \code{var_names}.
#' @param ncol        Number of columns in \code{facet_wrap} (one facet per
#'   method).  \code{NULL} lets ggplot2 choose automatically.
#'
#' @return A named list of \code{m * m} ggplot objects.  Names are
#'   \code{"<response>_<shock>"} (e.g. \code{"y1_y2"}).  Index naturally as
#'   \code{plots[["y1_y2"]]} or iterate with \code{lapply}.
#'
#' @examples
#' \dontrun{
#' plots <- plot_irf_distribution(sim_results, pop_irf,
#'                                var_names = c("GDP", "Infl", "Rate"))
#'
#' # Print a single panel
#' print(plots[["GDP_Rate"]])
#'
#' # Arrange all panels in a grid (requires patchwork)
#' library(patchwork)
#' wrap_plots(plots, ncol = 3)
#'
#' # Save every panel
#' for (nm in names(plots))
#'   ggsave(paste0("irf_", nm, ".pdf"), plots[[nm]], width = 8, height = 4)
#' }
plot_irf_distribution <- function(sim,
                                  probs       = c(0.05, 0.25, 0.75, 0.95),
                                  alpha_inner = 0.6,
                                  alpha_outer = 0.3,
                                  colors      = NULL,
                                  true_color  = "black",
                                  true_lwd    = 0.9,
                                  var_names   = NULL,
                                  shock_names = NULL,
                                  ncol        = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  # ── Validate inputs ───────────────────────────────────────────────────────────
  if (length(probs) != 4L || !all(diff(probs) > 0) ||
      probs[1L] < 0 || probs[4L] > 1)
    stop("`probs` must be a strictly increasing length-4 vector in (0, 1).")

  sim_results <- sim$irfs
  pop_irf     <- irf(sim$dgp, sim$h)
  title       <- sim$description
  m      <- dim(pop_irf)[1L]
  n_lags <- dim(pop_irf)[3L]
  lags   <- seq_len(n_lags) - 1L          # 0, 1, ..., n_lags - 1

  methods <- names(sim_results)
  if (is.null(methods)) methods <- paste0("Method ", seq_along(sim_results))
  names(sim_results) <- methods

  if (is.null(var_names))   var_names   <- paste0("y", seq_len(m))
  if (is.null(shock_names)) shock_names <- var_names
  if (length(var_names)   != m) stop("`var_names` must have length m.")
  if (length(shock_names) != m) stop("`shock_names` must have length m.")

  # ── 1. Compute quantiles for every (method, row, col, lag) ───────────────────
  # expand.grid: row varies fastest → matches array column-major order
  base_grid <- expand.grid(
    row = seq_len(m), col = seq_len(m), lag = lags,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  q_probs <- c(probs, 0.5)                        # lo2, lo1, hi1, hi2, median
  q_names <- c("lo2", "lo1", "hi1", "hi2", "med")

  chunks <- vector("list", length(methods))

  for (k in seq_along(methods)) {
    sims  <- sim_results[[methods[k]]]
    n_sim <- length(sims)

    # Stack into (m, m, n_lags, n_sim); unlist preserves column-major order
    irf_4d <- array(unlist(sims), dim = c(m, m, n_lags, n_sim))

    # Quantiles over simulation dim → (5, m, m, n_lags)
    q_arr <- apply(irf_4d, MARGIN = 1:3,
                   FUN = quantile, probs = q_probs, na.rm = TRUE)

    chunk <- base_grid
    for (qi in seq_along(q_names))
      chunk[[q_names[qi]]] <- as.vector(q_arr[qi, , , ])
    chunk$method <- methods[k]
    chunks[[k]]  <- chunk
  }

  df        <- do.call(rbind, chunks)
  df$method <- factor(df$method, levels = methods)

  # ── 2. Population IRF in long form ────────────────────────────────────────────
  pop_df       <- base_grid
  pop_df$value <- as.vector(pop_irf)   # column-major: row fastest → matches base_grid

  # ── 3. Colors (one per method, consistent across all returned plots) ──────────
  if (!is.null(colors)) {
    colors <- setNames(as.character(colors), methods)
  } else {
    # Use ggplot2's default discrete hues so fill/color are consistent
    gg_colors <- scales::hue_pal()(length(methods))
    colors <- setNames(gg_colors, methods)
  }

  # ── 4. Shared caption ─────────────────────────────────────────────────────────
  pct <- function(p) paste0(round(p * 100L), "%")
  caption_txt <- sprintf(
    paste0("Shaded bands: [%s–%s] inner and [%s–%s] outer quantile ",
           "intervals. Solid lines: medians. Dashed: true IRF."),
    pct(probs[2L]), pct(probs[3L]), pct(probs[1L]), pct(probs[4L])
  )

  # ── 5. Shared theme ───────────────────────────────────────────────────────────
  base_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position   = "none",           # facet label already identifies method
      strip.background  = ggplot2::element_rect(fill = "grey94", color = "grey70"),
      strip.text        = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.caption      = ggplot2::element_text(hjust = 0, size = 7.5,
                                                color = "grey40"),
      plot.title        = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle     = ggplot2::element_text(size = 9, color = "grey30")
    )

  # ── 6. One plot per (i, j) IRF entry ─────────────────────────────────────────
  plots      <- vector("list", m * m)
  plot_names <- character(m * m)
  k <- 0L

  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      k <- k + 1L

      sub_df  <- df[df$row == i & df$col == j, , drop = FALSE]
      sub_pop <- pop_df[pop_df$row == i & pop_df$col == j, , drop = FALSE]

      p <- ggplot2::ggplot(
        sub_df,
        ggplot2::aes(x = lag, group = method, color = method, fill = method)
      ) +
        # Outer quantile band
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lo2, ymax = hi2),
          alpha = alpha_outer, color = NA
        ) +
        # Inner quantile band
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lo1, ymax = hi1),
          alpha = alpha_inner, color = NA
        ) +
        # Median line
        ggplot2::geom_line(ggplot2::aes(y = med), linewidth = 0.55) +
        # Zero reference
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
        # True IRF (same in every method facet)
        ggplot2::geom_line(
          data        = sub_pop,
          mapping     = ggplot2::aes(x = lag, y = value),
          color       = true_color,
          linewidth   = true_lwd,
          linetype    = "dashed",
          inherit.aes = FALSE
        ) +
        # Facets = methods
        ggplot2::facet_wrap(
          ggplot2::vars(method),
          ncol   = ncol,
          scales = "fixed"
        ) +
        # Method-consistent colors (facet fill matches band fill even without legend)
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_fill_manual(values  = colors) +
        ggplot2::scale_x_continuous(breaks = lags) +
        ggplot2::labs(
          title = title,
          subtitle    = sprintf("Response: %s  —  Shock: %s",
                             var_names[i], shock_names[j]),
          # subtitle = sprintf("IRF[%d, %d]", i, j),
          x        = "Lag",
          y        = "IRF"#,
          # caption  = caption_txt
        ) +
        base_theme

      plots[[k]]     <- p
      plot_names[[k]] <- paste0(var_names[i], "_", shock_names[j])
    }
  }

  names(plots) <- plot_names
  plots
}


plot_fcst_distribution <- function(sim,
                                   probs       = c(0.05, 0.25, 0.75, 0.95),
                                   alpha_inner = 0.6,
                                   alpha_outer = 0.3,
                                   colors      = NULL,
                                   true_color  = "black",
                                   true_lwd    = 0.9,
                                   var_names   = NULL,
                                   ncol        = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  # ── Validate inputs ───────────────────────────────────────────────────────────
  if (length(probs) != 4L || !all(diff(probs) > 0) ||
      probs[1L] < 0 || probs[4L] > 1)
    stop("`probs` must be a strictly increasing length-4 vector in (0, 1).")

  sim_results <- sim$pre_err
  title       <- sim$description
  m      <- dim(sim$y[[1]])[2]
  n_lags <- sim$h
  lags   <- seq_len(n_lags)          # 1, ..., n_lags

  methods <- names(sim_results)
  if (is.null(methods)) methods <- paste0("Method ", seq_along(sim_results))
  names(sim_results) <- methods

  if (is.null(var_names))   var_names   <- paste0("y", seq_len(m))
  if (length(var_names)   != m) stop("`var_names` must have length m.")

  # ── 1. Compute quantiles for every (method, row, col, lag) ───────────────────
  # expand.grid: row varies fastest → matches array column-major order
  base_grid <- expand.grid(
    col = seq_len(m), lag = lags,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  q_probs <- c(probs, 0.5)                        # lo2, lo1, hi1, hi2, median
  q_names <- c("lo2", "lo1", "hi1", "hi2", "med")

  chunks <- vector("list", length(methods))
  rmse <- matrix(0, n_lags, length(methods), dimnames = list(1:n_lags, methods))

  for (k in seq_along(methods)) {
    cat("method", methods[k], "\n")
    sims  <- sim_results[[methods[k]]]
    if (is.list(sims)) sims <- array(unlist(sims), dim = c(n_lags, m, length(sims)))
    n_sim <- dim(sims)[3]

    # Quantiles over simulation dim → (5, m, n_lags)
    q_arr <- apply(sims, MARGIN = 1:2,
                   FUN = quantile, probs = q_probs, na.rm = TRUE)

    # RMSE
    rmse[, k] <- apply(sims^2, 1, median, na.rm = TRUE)

    chunk <- base_grid
    for (qi in seq_along(q_names))
      chunk[[q_names[qi]]] <- as.vector(q_arr[qi, , ])
    chunk$method <- methods[k]
    chunks[[k]]  <- chunk
  }
  print(rmse)

  df        <- do.call(rbind, chunks)
  df$method <- factor(df$method, levels = methods)

  # ── 3. Colors (one per method, consistent across all returned plots) ──────────
  if (!is.null(colors)) {
    colors <- setNames(as.character(colors), methods)
  } else {
    # Use ggplot2's default discrete hues so fill/color are consistent
    gg_colors <- scales::hue_pal()(length(methods))
    colors <- setNames(gg_colors, methods)
  }

  # ── 4. Shared caption ─────────────────────────────────────────────────────────
  pct <- function(p) paste0(round(p * 100L), "%")
  caption_txt <- sprintf(
    paste0("Shaded bands: [%s–%s] inner and [%s–%s] outer quantile ",
           "intervals. Solid lines: medians."),
    pct(probs[2L]), pct(probs[3L]), pct(probs[1L]), pct(probs[4L])
  )

  # ── 5. Shared theme ───────────────────────────────────────────────────────────
  base_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position   = "none",           # facet label already identifies method
      strip.background  = ggplot2::element_rect(fill = "grey94", color = "grey70"),
      strip.text        = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor  = ggplot2::element_blank(),
      plot.caption      = ggplot2::element_text(hjust = 0, size = 7.5,
                                                color = "grey40"),
      plot.title        = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle     = ggplot2::element_text(size = 9, color = "grey30")
    )

  # ── 6. One plot per (i, j) IRF entry ─────────────────────────────────────────
  plots      <- vector("list", m)
  plot_names <- character(m)
  k <- 0L

  for (i in seq_len(m)) {
      k <- k + 1L

      sub_df  <- df[df$col == i, , drop = FALSE]

      p <- ggplot2::ggplot(
        sub_df,
        ggplot2::aes(x = lag, group = method, color = method, fill = method)
      ) +
        # Outer quantile band
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lo2, ymax = hi2),
          alpha = alpha_outer, color = NA
        ) +
        # Inner quantile band
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lo1, ymax = hi1),
          alpha = alpha_inner, color = NA
        ) +
        # Median line
        ggplot2::geom_line(ggplot2::aes(y = med), linewidth = 0.55) +
        # Zero reference
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
        # Facets = methods
        ggplot2::facet_wrap(
          ggplot2::vars(method),
          ncol   = ncol,
          scales = "fixed"
        ) +
        # Method-consistent colors (facet fill matches band fill even without legend)
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_fill_manual(values  = colors) +
        ggplot2::scale_x_continuous(breaks = lags) +
        ggplot2::labs(
          title = title,
          subtitle    = sprintf("Time Series: %s",
                                var_names[i]),
          # subtitle = sprintf("IRF[%d, %d]", i, j),
          x        = "Horizon",
          y        = "Forecast error"#,
          # caption  = caption_txt
        ) +
        base_theme

      plots[[k]]     <- p
      plot_names[[k]] <- var_names[i]
  }

  names(plots) <- plot_names
  plots
}



# USE IT

g <- plot_irf_distribution(results$experiment_12)
print(g$y1_y2)

gf <- plot_fcst_distribution(results$experiment_8)

for (i in seq_along(results)) {
  g <- plot_irf_distribution(results[[i]])
  ggsave(filename = paste0("experiment", i, ".pdf"), plot = g$y1_y2)
}

for (i in seq_along(results)) {
  gf <- plot_fcst_distribution(results[[i]])
  ggsave(filename = paste0("experiment", i, "_fcst_err.pdf"), plot = gf$y1)
}
