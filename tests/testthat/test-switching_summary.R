# tests/testthat/test-switching_summary.R

library(testthat)
library(BayCauRETM)

test_that("switching_probability_summary computes hazards and summaries", {
  set.seed(11)
  df <- data.frame(
    pat_id = rep(1:15, each = 4),
    k_idx  = rep(1:4, times = 15),
    A      = rbinom(60, 1, 0.2),
    Y_obs  = rpois(60, 1),
    age    = rnorm(60, 60, 8)
  )
  sw <- suppressWarnings(switching_probability_summary(df, covariates = "age"))
  expect_s3_class(sw, "switching_summary")
  expect_true(all(c("df2","summary_df","model") %in% names(sw)))
  expect_true(all(c("hazard","surv_prob","switch_prob") %in% names(sw$df2)))
})

test_that("plot.switching_summary returns ggplot objects", {
  skip_if_not_installed("ggplot2")
  set.seed(7)
  df <- data.frame(
    pat_id = rep(1:8, each = 3),
    k_idx  = rep(1:3, times = 8),
    A      = rbinom(24, 1, 0.5),
    Y_obs  = rpois(24, 1)
  )
  sw <- suppressWarnings(switching_probability_summary(df))

  p1 <- plot(sw, type = "boxplot")
  expect_s3_class(p1, "ggplot")
  p2 <- plot(sw, type = "ribbon")
  expect_s3_class(p2, "ggplot")
})
