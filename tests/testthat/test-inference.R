test_that("estimate_Rt_simple returns correct structure", {
  incidence <- c(1,2,3,5,8,13,21,18,14,10,7,5,3,2,1)
  out <- estimate_Rt_simple(incidence, mean_si = 5, sd_si = 2,
                             window = 5L, n_boot = 50L)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("t_end","Rt_mean","Rt_lower","Rt_upper"))
  expect_true(all(out$Rt_lower <= out$Rt_mean + 1e-9))
  expect_true(all(out$Rt_upper >= out$Rt_mean - 1e-9))
  expect_true(all(out$Rt_mean >= 0, na.rm = TRUE))
})

test_that("estimate_Rt_simple Rt > 1 during growth phase", {
  # Exponentially growing epidemic
  incidence <- round(2^(0:14))
  out <- estimate_Rt_simple(incidence, mean_si = 3, sd_si = 1,
                             window = 5L, n_boot = 50L)
  # Most Rt estimates in growth phase should exceed 1
  expect_gt(mean(out$Rt_mean, na.rm = TRUE), 1)
})

test_that("estimate_Rt falls back gracefully without EpiEstim", {
  skip_if(requireNamespace("EpiEstim", quietly = TRUE),
          "EpiEstim is installed — fallback not triggered")
  incidence <- c(1,2,4,8,10,9,7,5,3,2)
  # Should not error; should return a data frame via fallback
  out <- suppressMessages(
    estimate_Rt(incidence, mean_si = 5, sd_si = 2)
  )
  expect_s3_class(out, "data.frame")
})

test_that("score_forecast returns correct columns", {
  fc_df <- data.frame(
    time     = 1:10,
    I_median = seq(0.01, 0.1, length.out = 10),
    I_lower  = seq(0.005, 0.08, length.out = 10),
    I_upper  = seq(0.015, 0.12, length.out = 10)
  )
  actual <- seq(0.012, 0.095, length.out = 10)
  out <- score_forecast(fc_df, actual)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("CRPS","coverage_95","MAE"))
  expect_gte(out$coverage_95, 0)
  expect_lte(out$coverage_95, 1)
  expect_gte(out$MAE, 0)
})
