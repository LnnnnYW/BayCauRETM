# tests/testthat/test-switching_summary.R

library(testthat)
library(BayCauRETM)

test_that("switching_probability_summary returns class and plot", {
  df <- data.frame(
    pat_id = rep(1:5, each = 4),
    k_idx  = rep(1:4, 5),
    A      = rbinom(20, 1, 0.3),
    Y_prev = rpois(20, 0.6)
  )
  res <- switching_probability_summary(df)
  expect_s3_class(res, "switching_summary")
  expect_true(all(c("df2", "summary_df", "model") %in% names(res)))
  expect_equal(nrow(res$summary_df), 4)

  gp <- plot(res)
  expect_s3_class(gp, "ggplot")
})
