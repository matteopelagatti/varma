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
  iPhi0 <- solve(iPhi0)
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
  out$residuals <- residuals.varma(out)
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
    pre_err = list_pre_err
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
