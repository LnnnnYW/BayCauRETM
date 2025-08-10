# tests/testthat/test-preprocess.R

library(testthat)
library(BayCauRETM)

test_that("preprocess_data sorts, remaps IDs, and can create lag column", {
  set.seed(1)
  df <- data.frame(
    pat_id = c("S2","S2","S1"),
    k_idx  = c(2,1,1),
    Y_obs  = c(1,0,2),
    T_obs  = 0,
    A      = c(1,1,0)
  )

  out <- suppressWarnings(preprocess_data(df, K = 3, lag_col = "lagYK"))
  expect_s3_class(out, "preprocess_data")
  pd <- out$processed_df

  expect_equal(out$n_pat, 2L)
  expect_true(all(pd$pat_id %in% 1:2))

  ord <- order(pd$pat_id, pd$k_idx)
  expect_identical(ord, seq_along(ord))

  expect_true("lagYK" %in% names(pd))
  by_first_k <- ave(pd$k_idx, pd$pat_id, FUN = function(x) seq_along(x) == 1)
  expect_true(all(pd$lagYK[by_first_k] == 0))
})
