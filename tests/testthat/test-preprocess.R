# tests/testthat/test-preprocess.R

library(testthat)
library(BayCauRETM)

test_that("preprocess_data returns expected structure", {
  set.seed(1)
  df <- data.frame(
    id    = rep(1:4, each = 3),
    k     = rep(1:3, 4),
    Y     = rpois(12, 1),
    T     = rbinom(12, 1, 0.4),
    C     = rbinom(12, 1, 0.02),
    A     = rbinom(12, 1, 0.4),
    regexp = "No variation in T_obs"
  )
  out <- preprocess_data(df, id_col = "id", k_col = "k",
                         y_col = "Y", t_col = "T", c_col = "C",
                         a_col = "A", x_cols = NULL, K = 3)
  expect_type(out, "list")
  expect_true(all(c("processed_df", "n_pat") %in% names(out)))
  expect_equal(out$n_pat, 4)
  expect_true(all(c("Y_prev", "T_prev", "C_prev") %in%
                    names(out$processed_df)))
})
