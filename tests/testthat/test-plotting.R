# Plotting tests: check that functions return ggplot objects without error

params <- list(beta = 0.3, gamma = 0.1, sigma = 0.2, delta = 0.005)
init   <- c(S = 0.9989, E = 0.0006, I = 0.0005, R = 0, D = 0)
traj   <- solve_model(params, init, 0:60, model = "SEIRD")

obs_Y  <- traj$I[1:15] + rnorm(15, 0, 1e-5)
obs_R  <- traj$R[1:15] + rnorm(15, 0, 1e-5)

test_that("plot_trajectory returns a ggplot", {
  p <- plot_trajectory(traj, obs_Y = obs_Y, obs_R = obs_R, T_obs_end = 14)
  expect_s3_class(p, "ggplot")
})

test_that("plot_forecast returns a ggplot", {
  fc <- data.frame(time = 0:60,
                   I_median = traj$I,
                   I_lower  = traj$I * 0.8,
                   I_upper  = traj$I * 1.2)
  p <- plot_forecast(fc, obs_Y = obs_Y, T_obs_end = 14)
  expect_s3_class(p, "ggplot")
})

test_that("plot_scenarios returns a ggplot", {
  pi_a <- build_pi_step(c(10), c(1, 0.5))
  pi_b <- function(t) 1

  sc_df <- do.call(rbind, lapply(
    list("Lockdown" = pi_a, "No NPI" = pi_b),
    function(pi_fn) {
      tr <- solve_model(params, init, 0:60,
                        model = "SEIRD", pi_fn = pi_fn)
      data.frame(time=tr$time, I_median=tr$I,
                 I_lower=tr$I*0.8, I_upper=tr$I*1.2,
                 scenario=deparse(substitute(pi_fn)))
    }
  ))
  sc_df$scenario <- rep(c("Lockdown","No NPI"), each = 61)
  p <- plot_scenarios(sc_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_Rt returns a ggplot", {
  Rt_df <- data.frame(t_end    = 5:20,
                       Rt_mean  = runif(16, 0.8, 2.0),
                       Rt_lower = runif(16, 0.5, 0.9),
                       Rt_upper = runif(16, 1.8, 2.5))
  p <- plot_Rt(Rt_df, change_times = c(8, 15))
  expect_s3_class(p, "ggplot")
})

test_that("plot_multipatch_snapshot returns a ggplot", {
  M <- matrix(c(0, 0.01, 0.01, 0), 2, 2)
  ode_fn <- build_multipatch_SEIR(
    n_patches = 2, M = M,
    beta_vec  = c(0.3, 0.25),
    gamma_vec = c(0.1, 0.1),
    sigma_vec = c(0.2, 0.2)
  )
  init_mat <- matrix(c(0.9999, 0.99995, 5e-5, 2e-5, 1e-4, 1e-5, 0, 0), 2)
  mp_df <- solve_multipatch(ode_fn, init_mat, times = 0:30, n_patches = 2)
  p <- plot_multipatch_snapshot(mp_df, t_snapshot = 15,
                                 patch_names = c("Region A", "Region B"))
  expect_s3_class(p, "ggplot")
})
