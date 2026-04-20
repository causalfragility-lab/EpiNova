test_that("gravity_mobility returns valid matrix", {
  N_vec    <- c(1e6, 5e6, 2e6)
  dist_mat <- matrix(c(0,100,200, 100,0,150, 200,150,0), 3, 3)
  M <- gravity_mobility(N_vec, dist_mat, max_travel = 0.05)

  expect_equal(dim(M), c(3, 3))
  # Diagonal should be 0
  expect_equal(diag(M), c(0, 0, 0))
  # Row sums should not exceed max_travel
  expect_true(all(rowSums(M) <= 0.05 + 1e-9))
  # All values non-negative
  expect_true(all(M >= 0))
})

test_that("build_multipatch_SEIR and solve_multipatch run without error", {
  n <- 2
  M <- matrix(c(0, 0.01, 0.01, 0), 2, 2)

  ode_fn <- build_multipatch_SEIR(
    n_patches  = n,
    M          = M,
    beta_vec   = c(0.3, 0.25),
    gamma_vec  = c(0.1, 0.1),
    sigma_vec  = c(0.2, 0.2)
  )

  init_mat <- matrix(
    c(0.9999, 0.99995, 5e-5, 2e-5, 1e-4, 1e-5, 0, 0),
    nrow = 2
  )

  expect_no_error({
    mp_df <- solve_multipatch(ode_fn, init_mat,
                               times = 0:30, n_patches = n)
  })

  mp_df <- solve_multipatch(ode_fn, init_mat,
                             times = 0:30, n_patches = n)
  expect_s3_class(mp_df, "data.frame")
  expect_true("I_1" %in% names(mp_df))
  expect_true("I_2" %in% names(mp_df))
  # Infected proportions must remain non-negative
  expect_true(all(mp_df$I_1 >= -1e-9))
  expect_true(all(mp_df$I_2 >= -1e-9))
})

test_that("isolation (M=0) gives independent single-patch dynamics", {
  # With zero mobility, each patch should evolve exactly like solve_model
  M_zero <- matrix(0, 2, 2)
  ode_fn <- build_multipatch_SEIR(
    n_patches = 2, M = M_zero,
    beta_vec  = c(0.3, 0.3),
    gamma_vec = c(0.1, 0.1),
    sigma_vec = c(0.2, 0.2)
  )
  init_mat <- matrix(
    c(0.999, 0.999, 5e-4, 5e-4, 5e-4, 5e-4, 0, 0),
    nrow = 2
  )
  mp_df <- solve_multipatch(ode_fn, init_mat,
                             times = 0:50, n_patches = 2)

  single <- solve_model(
    list(beta = 0.3, gamma = 0.1, sigma = 0.2),
    c(S = 0.999, E = 5e-4, I = 5e-4, R = 0),
    times = 0:50, model = "SEIR"
  )
  # Patch 1 I should match single-patch I closely
  expect_equal(mp_df$I_1, single$I, tolerance = 1e-4)
})
