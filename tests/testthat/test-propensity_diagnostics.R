# tests/testthat/test-propensity_diagnostics.R

library(testthat)
library(BayCauRETM)

test_that("propensity_score_diagnostics returns ps_diag object", {
  dat <- data.frame(
    A      = rbinom(120, 1, 0.4),
    Y_prev = rpois(120, 1),
    X1     = rnorm(120),
    X2     = rnorm(120)
  )

  psd <- propensity_score_diagnostics(
    data       = dat,
    treat_col  = "A",
    covariates = c("Y_prev", "X1", "X2"),
    plot_types = "histogram"
  )
  expect_s3_class(psd, "ps_diag")
  expect_equal(nrow(psd$summary), 6)
  expect_true("histogram" %in% names(psd$plots))
})
