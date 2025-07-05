# tests/testthat/test-preprocess.R

library(testthat)
library(BayCauRETM)

test_that("preprocess_data returns expected structure", {
  set.seed(1)
  df <- data.frame(
    pat_id = rep(1:4, each = 3),
    k_idx  = rep(1:3, 4),
    Y_obs  = rpois(12, 1),
    T_obs  = rbinom(12, 1, 0.4),
    A      = rbinom(12, 1, 0.4)
  )

  out <- preprocess_data(df, K = 3)

  expect_type(out, "list")
  expect_true(all(c("processed_df", "n_pat") %in% names(out)))
  expect_equal(out$n_pat, 4)
  expect_true(all(c("Y_prev", "T_prev") %in% names(out$processed_df)))
})
