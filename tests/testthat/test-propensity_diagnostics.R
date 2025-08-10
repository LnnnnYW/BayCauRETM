# tests/testthat/test-propensity_diagnostics.R

library(testthat)
library(BayCauRETM)

test_that("propensity_score_diagnostics fits model and builds plots", {
  skip_if_not_installed("ggplot2")

  set.seed(1)
  df <- data.frame(
    A      = rbinom(200, 1, 0.5),
    Y_prev = rpois(200, 1),
    X1     = rnorm(200)
  )

  ps <- propensity_score_diagnostics(
    data       = df,
    treat_col  = "A",
    covariates = c("Y_prev", "X1"),
    trim       = c(0.01, 0.99),
    plot_types = c("histogram", "density"),
    bins       = 20
  )

  expect_s3_class(ps, "ps_diag")
  expect_true(is.data.frame(ps$summary))
  expect_true(all(c("model","ps","trim_flag","plots") %in% names(ps)))
  expect_true(all(c("histogram","density") %in% names(ps$plots)))
})
