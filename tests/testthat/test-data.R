test_that("hubei_covid is well-formed", {
  expect_true(is.list(hubei_covid))
  expect_named(hubei_covid,
               c("NI","RI","N","begin_date","end_date","description"),
               ignore.order = TRUE)
  expect_length(hubei_covid$NI, 30)
  expect_length(hubei_covid$RI, 30)
  expect_true(all(hubei_covid$NI >= hubei_covid$RI),
              label = "Cumulative confirmed >= cumulative removed")
  expect_gt(hubei_covid$N, 1e6)
})

test_that("prep_proportions returns valid proportions", {
  prop <- prep_proportions(hubei_covid)
  expect_named(prop, c("Y", "R"))
  expect_length(prop$Y, 30)
  expect_true(all(prop$Y >= 0))
  expect_true(all(prop$R >= 0))
  expect_true(all(prop$Y + prop$R <= 1 + 1e-9))
})

test_that("hubei_covid works end-to-end with solve_model", {
  p  <- prep_proportions(hubei_covid)
  params <- list(beta = 0.35, gamma = 0.07, sigma = 0.2, delta = 0.005)
  init   <- c(S = 1 - p$Y[1] * 2 - p$Y[1],
              E = p$Y[1] * 2, I = p$Y[1], R = p$R[1], D = 0)
  expect_no_error(
    solve_model(params, init, 0:30, model = "SEIRD")
  )
})
