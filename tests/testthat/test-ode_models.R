test_that("solve_model returns correct structure for SIR", {
  params <- list(beta = 0.3, gamma = 0.1)
  init   <- c(S = 0.999, I = 0.001, R = 0)
  times  <- 0:50
  out    <- solve_model(params, init, times, model = "SIR")

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "S", "I", "R"))
  expect_equal(nrow(out), length(times))
  # Populations must stay in [0,1]
  expect_true(all(out[, -1] >= -1e-9))
  expect_true(all(out[, -1] <= 1 + 1e-9))
})

test_that("solve_model conserves population for SEIR", {
  params <- list(beta = 0.3, gamma = 0.1, sigma = 0.2)
  init   <- c(S = 0.998, E = 0.001, I = 0.001, R = 0)
  times  <- 0:100
  out    <- solve_model(params, init, times, model = "SEIR")

  row_sums <- rowSums(out[, c("S","E","I","R")])
  expect_true(all(abs(row_sums - 1) < 1e-6),
              label = "SEIR compartments sum to 1")
})

test_that("solve_model conserves population for SEIRD (S+E+I+R+D=1)", {
  params <- list(beta = 0.35, gamma = 0.07, sigma = 0.2, delta = 0.005)
  init   <- c(S = 0.9989, E = 0.0006, I = 0.0005, R = 0, D = 0)
  times  <- 0:100
  out    <- solve_model(params, init, times, model = "SEIRD")

  row_sums <- rowSums(out[, c("S","E","I","R","D")])
  expect_true(all(abs(row_sums - 1) < 1e-5))
})

test_that("pi_fn = 0 halts transmission", {
  params  <- list(beta = 0.5, gamma = 0.1)
  init    <- c(S = 0.99, I = 0.01, R = 0)
  times   <- 0:30
  zero_pi <- function(t) 0
  out     <- solve_model(params, init, times, model = "SIR",
                         pi_fn = zero_pi)
  # With no transmission, I should not increase
  expect_true(all(diff(out$I) <= 1e-9),
              label = "Infected should not grow when pi(t)=0")
})

test_that("SVEIRD model runs without error", {
  params <- list(beta = 0.3, gamma = 0.08, sigma = 0.2,
                 delta = 0.003, vax_rate = 0.005, wane = 0.01, ve = 0.85)
  init   <- c(S = 0.994, V = 0, E = 0.003, I = 0.002, R = 0.001, D = 0)
  times  <- 0:60
  expect_no_error(solve_model(params, init, times, model = "SVEIRD"))
})

test_that("solve_model dispatches correctly to all model types", {
  base_p <- list(beta = 0.3, gamma = 0.1, sigma = 0.2,
                 delta = 0.005, vax_rate = 0.002, wane = 0.01, ve = 0.8)
  models <- c("SIR","SEIR","SEIRD","SVEIR","SVEIRD")
  inits  <- list(
    SIR    = c(S=0.999, I=0.001, R=0),
    SEIR   = c(S=0.998, E=0.001, I=0.001, R=0),
    SEIRD  = c(S=0.998, E=0.001, I=0.001, R=0, D=0),
    SVEIR  = c(S=0.998, V=0, E=0.001, I=0.001, R=0),
    SVEIRD = c(S=0.998, V=0, E=0.001, I=0.001, R=0, D=0)
  )
  for (m in models) {
    out <- solve_model(base_p, inits[[m]], 0:10, model = m)
    expect_true(inherits(out, "data.frame"), info = paste("model:", m))
    expect_gt(nrow(out), 1L)
  }
})
