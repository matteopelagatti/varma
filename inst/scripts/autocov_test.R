# extra function, since the package lyapunov is not in the CRAN we define this
# functio here

autocov_ly <- function(varma, maxlag = 10) { # based on Lyapunov equation solution in package netcontrol
  if (maxlag < 0) maxlag <- 0
  ss <- varma_to_ss(varma)
  mpq <- dim(varma)
  k <- dim(ss$T)[1]
  G <- array(0, c(k, k, maxlag + 1))
  G[,,1] <- lyapunov::dlyap(ss$T, tcrossprod(ss$R %*% ss$Q, ss$R))
  if (maxlag > 0) for (i in 1:maxlag) {
    G[,,i+1] <- ss$T %*% G[,,i]
  }
  G[1:mpq[1], 1:mpq[1], ]
}


# function to compare times of autocov_xx functions
library(bench)

autocov_comp_time <- function(fn_list, vdim = c(5, 10, 15, 20, 25, 30),
                              vp = c(1, 5), vq = c(1, 5),
                              n_times = 10, maxlag = 10) {

  if (is.null(names(fn_list))) {
    names(fn_list) <- paste0("fn_", seq_along(fn_list))
  }

  param_grid <- expand.grid(d = vdim, p = vp, q = vq)

  results_list <- lapply(seq_len(nrow(param_grid)), function(i) {
    d <- param_grid$d[i]
    p <- param_grid$p[i]
    q <- param_grid$q[i]

    cat("Working on case d =", d, "p =", p, "q =", q, "\n")
    print(Sys.time())

    # Generate the model exactly ONCE for this specific combination
    model <- rvarma(d, p, q)

    fn_results <- lapply(names(fn_list), function(fn_name) {
      f <- fn_list[[fn_name]]

      # Run the benchmark
      bm <- bench::mark(
        expr = f(model, maxlag),
        iterations = n_times,
        min_time = 0,       # Ensures it runs exactly n_times (overrides bench's default minimum time)
        filter_gc = FALSE,  # Keeps times even if a garbage collection event occurred during the run
        check = FALSE
      )

      # Extract the vector of ALL times (in seconds)
      times_vec <- as.numeric(bm$time[[1]])

      # Extract the single memory allocation value (in bytes)
      mem_val <- as.numeric(bm$mem_alloc)

      # Create a data frame with one row per iteration
      data.frame(
        function_name = fn_name,
        d = d,
        p = p,
        q = q,
        iteration = seq_along(times_vec),
        time_seconds = times_vec,
        memory_bytes = mem_val,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, fn_results)
  })

  final_df <- do.call(rbind, results_list)
  return(final_df)
}

# timing the autocov_xx functions

set.seed(2026)
dfbench <- autocov_comp_time(fn_list = list(Ansley  = autocov_an,
                                            dLyap   = autocov_ly,
                                            McElroy = autocov_mc,
                                            Mittnik = autocov_mi,
                                            Sbrana  = autocov_gs),
                             vdim = c(5, 10, 15, 20, 25, 30),
                             vp = c(1, 5),
                             vq = c(1, 5),
                             n_times = 10, maxlag = 10)

# plotting the times
library(tidyverse)

# lines based on the median

dfbench |>
  mutate(function_name = ifelse(function_name == "Sbrana", "Theorem 1", function_name)) |>
  select(-iteration) |>
  group_by(function_name, d, p, q) |>
  summarise(seconds = median(time_seconds), bytes = first(memory_bytes)) |>
  rename(method = function_name) |>
  ggplot(aes(x = d, y = seconds, color = method, shape = method)) +
  geom_line() +
  geom_point() +
  facet_grid(rows = vars(p), cols = vars(q), labeller = \(x) label_both(x, sep = " = ")) +
  scale_y_log10()

ggsave(filename = "speed.pdf", width = 7, height = 5)



dfbench |>
  select(-iteration, -memory_bytes) |>
  group_by(function_name, d, p, q) |>
  summarise(seconds = median(time_seconds)) |>
  rename(method = function_name) |>
  pivot_wider(names_from = method, values_from = seconds) |>
  mutate(ratio = dLyap/Sbrana) |>
  ggplot(aes(x=d, y=ratio)) +
  geom_line() +
  facet_grid(rows = vars(p), cols = vars(q), labeller = \(x) label_both(x, sep = " = "),
             scales = "free")



# dfbench |>
#   rename(method = function_name, seconds = time_seconds) |>
#   ggplot(aes(x = method, y = seconds)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(p, q), cols = vars(d)) +
#   scale_y_log10()
#
#
# dfbench |>
#   rename(method = function_name, seconds = time_seconds) |>
#   mutate(d = as.factor(d)) |>
#   ggplot(aes(x = method, y = seconds, fill = d)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(p), cols = vars(q)) +
#   scale_y_log10()
#
# # d as x and method as columns
#
# dfbench |>
#   rename(method = function_name, seconds = time_seconds) |>
#   mutate(d = as.factor(d)) |>
#   ggplot(aes(x = d, y = seconds)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(p, q), cols = vars(method)) +
#   scale_y_log10()
#
#


# set.seed(2026)
# mod1 <- rvarma(5, 5, 5)
# mod1$cov[1,2] <- mod1$cov[2, 1] <- 0.7
# inv_roots(mod1, plot = TRUE)
# is_identified(mod1)
# maxlags <- 10
#
# library(microbenchmark)
#
# mPhi <- mod1$ar |> array(c(dim(mod1)[1], dim(mod1)[1]*dim(mod1)[2]))
# mTheta <- mod1$ma |> array(c(dim(mod1)[1], dim(mod1)[1]*dim(mod1)[3]))
#
# foo <- microbenchmark(
#   acf_gs <- autocov_gs(mod1, maxlags),
#   acf_mc <- autocov_mc(mod1, maxlags),
#   # acf_ts <- MTS::VARMAcov(mPhi, -mTheta, mod1$cov, 3),
#   acf_an <- autocov_mc(mod1, maxlags),
#   # acf_sp <- autocov_sp(mod1, 3),
#   acf_ly <- autocov_ly(mod1, maxlags),
#   acf_mi <- autocov_mi(mod1, maxlags),
#   times = 10)
#
# cbind(
#   GS = acf_gs[1, 2, ],
#   MC = acf_mc[1, 2, ],
#   AN = acf_an[1, 2, ],
#   LY = acf_ly[1, 2, ],
#   MI = acf_mi[1, 2, ]
# )
#
# levels(foo$expr) <- c("Sbrana", "McElroy", #"MTS",
#                       "Ansley", #"Sp-Kron",
#                       "dLyap", "Mittnik")
#
# plot(foo, log = "y")
#
# times <- 10
# sizes <- seq(3, 18, 3)
# df <- data.frame(expr = character(2*times*length(sizes)),
#                  time = integer(2*times*length(sizes)),
#                  dim = integer(2*times*length(sizes)))
# for (i in 1:length(sizes)) {
#   set.seed(2026)
#   mod1 <- rvarma(sizes[i], 5, 5)
#   foo <- microbenchmark(
#     acf5 <- autocov_gs(mod1, 3),
#     acf5 <- autocov_mc(mod1, 3),
#     times = times)
#   df[((i-1)*times*2 + 1):(i*times*2), 1] <- as.character(foo[, 1])
#   df[((i-1)*times*2 + 1):(i*times*2), 2] <- foo[, 2]
#   df[((i-1)*times*2 + 1):(i*times*2), 3] <- sizes[i]
# }
#
# library(tidyverse)
# df |>
#   mutate(algorithm = ifelse(expr == "acf5 <- autocov_mc(mod1, 3)", "McElroy", "Our")) ->
#   df
#
# df |>
#   group_by(dim, algorithm) |>
#   summarise(millisec = mean(time)/1000000) |>
#   ggplot(aes(x = dim, y = millisec, col = algorithm)) +
#   geom_line() +
#   geom_point() +
#   scale_y_log10()
#
#
# df |>
#   group_by(dim, algorithm) |>
#   summarise(millisec = mean(time)/1000000) |>
#   ggplot(aes(x = dim, y = millisec, col = algorithm)) +
#   geom_point() +
#   geom_smooth(formula = y ~ s(x, bs = "cs"), method = "gam") +
#   scale_y_log10()
