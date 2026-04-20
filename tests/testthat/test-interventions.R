test_that("build_pi_step returns correct values", {
  pi_fn <- build_pi_step(change_times = c(10, 20),
                          pi_values    = c(1.0, 0.5, 0.2))
  expect_equal(pi_fn(5),  1.0)
  expect_equal(pi_fn(10), 0.5)   # at the knot, steps to new value
  expect_equal(pi_fn(15), 0.5)
  expect_equal(pi_fn(20), 0.2)
  expect_equal(pi_fn(99), 0.2)
})

test_that("build_pi_spline returns values in [0,1]", {
  pi_fn <- build_pi_spline(
    knot_times  = c(0, 10, 20, 60),
    knot_values = c(1, 0.8, 0.4, 0.3)
  )
  vals <- sapply(c(0, 5, 10, 15, 20, 40, 60), pi_fn)
  expect_true(all(vals >= 0 & vals <= 1),
              label = "spline pi(t) in [0,1]")
})

test_that("build_pi_spline requires at least 4 knots", {
  expect_error(
    build_pi_spline(c(0, 10, 20), c(1, 0.5, 0.3)),
    regexp = "ord"   # interpSpline error message
  )
})

test_that("build_pi_exp decays correctly", {
  pi_fn <- build_pi_exp(lambda = 0.1, t0 = 0)
  expect_equal(pi_fn(0),  1,   tolerance = 1e-9)
  expect_lt(pi_fn(10), pi_fn(5))   # strictly decreasing
  expect_gte(pi_fn(100), 0)        # non-negative
})

test_that("compose_pi multiplies correctly", {
  f1 <- function(t) 0.5
  f2 <- function(t) 0.4
  f3 <- function(t) 0.5
  combined <- compose_pi(f1, f2, f3)
  expect_equal(combined(0), 0.5 * 0.4 * 0.5, tolerance = 1e-12)
})

test_that("compose_pi with step and spline stays in [0,1]", {
  step   <- build_pi_step(c(10), c(1, 0.5))
  spline <- build_pi_spline(c(0, 5, 15, 30), c(1, 0.9, 0.7, 0.6))
  combo  <- compose_pi(step, spline)
  vals   <- sapply(0:30, combo)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("build_phi_pulse integrates approximately to phi_values", {
  phi_fn <- build_phi_pulse(change_times = c(10),
                             phi_values   = c(0.3),
                             bandwidth    = 0.5)
  # Peak should be close to 0.3
  expect_equal(phi_fn(10), 0.3, tolerance = 1e-6)
  # Away from the pulse it should be near 0
  expect_lt(phi_fn(20), 0.01)
})

test_that("gp_cov_sqexp returns symmetric PD matrix", {
  times <- seq(0, 50, by = 5)
  K     <- gp_cov_sqexp(times, l = 14, sigma = 0.3)
  expect_equal(dim(K), c(length(times), length(times)))
  expect_true(isSymmetric(K, tol = 1e-10))
  # All eigenvalues positive (positive definite)
  eigs <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})
