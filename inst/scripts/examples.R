set.seed(2026)
mod1 <- rvarma(3, 2, 1)
P <- matrix(rnorm(9), 3, 3)
mod1$cov <- crossprod(P)
mod1
mod1$cov |> cov2cor()
inv_roots(mod1, plot = TRUE)

Y1 <- sim_varma(mod1, 500)
plot(ts(Y1))

library(microbenchmark)
microbenchmark(
  est_rby <- fit_varma_rby(Y1, 3, 2, FALSE),
  est_ihr <- fit_varma_ihr(Y1, 3, 2, FALSE),
  est_hr1 <- fit_varma_ihr(Y1, 3, 2, FALSE, maxit = 1),
  est_net <- fit_varma_net(Y1, 3, 2, FALSE),
  est_fkf <- fit_varma_fkf(Y1, 3, 2, FALSE),
  est_dpd <- fit_varma_dpd(Y1, 4, 4, FALSE),
  est_dpf <- fit_varma_dpf(Y1, 4, 4, FALSE),
  times = 5
)

mod1
est_rby
est_ihr
est_hr1
est_net
est_fkf
est_dpd
est_dpf

irf(mod1, plot = TRUE, title = "Population")
irf(est_rby, plot = TRUE, title = "RBY")
irf(est_ihr, plot = TRUE, title = "IHR")
irf(est_hr1, plot = TRUE, title = "IHR1")
irf(est_fkf, plot = TRUE, title = "FKF")
irf(est_net, plot = TRUE, title = "NET")
irf(est_dpd, plot = TRUE, title = "DPD")
irf(est_dpf, plot = TRUE, title = "DPF")

c(
  RBY = irf_distance(irf(mod1), irf(est_rby)),
  IHR = irf_distance(irf(mod1), irf(est_ihr)),
  HR1 = irf_distance(irf(mod1), irf(est_hr1)),
  NET = irf_distance(irf(mod1), irf(est_net)),
  FKF = irf_distance(irf(mod1), irf(est_fkf)),
  DPD = irf_distance(irf(mod1), irf(est_dpd)),
  DPF = irf_distance(irf(mod1), irf(est_dpf))
)

kronecker_indices(mod1$ar, mod1$ma, 3)
MTS::Kronid(Y1)
