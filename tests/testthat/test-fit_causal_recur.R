# tests/testthat/test-fit_causal_recur.R

library(testthat)
library(BayCauRETM)

test_that("fit_causal_recur fails early when required cols are missing", {
  df <- data.frame(sid = 1:3, period = 1:3, event_count = 0, death_flag = 0)
  expect_error(
    fit_causal_recur(
      data = df, K = 2,
      id_col = "sid", time_col = "period", treat_col = "trt_arm",
      formula_T = death_flag ~ 1, formula_Y = event_count ~ 1,
      stan_model_file = tempfile(fileext = ".rds")
    ),
    "Columns not found"
  )
})

test_that("fit_causal_recur errors if Stan model file is missing", {
  df <- data.frame(
    sid         = rep(1:2, each = 2),
    period      = rep(1:2, 2),
    event_count = 0,
    death_flag  = 0,
    trt_arm     = 0
  )

  expect_error(
    fit_causal_recur(
      data = df, K = 2,
      id_col = "sid", time_col = "period", treat_col = "trt_arm",
      formula_T = death_flag  ~ A,
      formula_Y = event_count ~ A,
      num_chains = 1, iter = 50,
      stan_model_file = tempfile(fileext = ".rds")   # file does not exist
    ),
    "Stan model file doesn't exist"
  )
})
