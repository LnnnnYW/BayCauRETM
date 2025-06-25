# tests/testthat/test_result_summary_table.R

library(testthat)
library(BayCauRETM)

test_that("result_summary_table merges two summaries", {
  fit_stub <- list(stan_fit = NULL); class(fit_stub) <- "causal_recur_fit"
  gfake <- list(delta = list("s=1" = list(mean = 0, CI_lower = -1, CI_upper = 1)))
  class(gfake) <- "gcomp_out"

  res <- result_summary_table(fit_stub, gfake, format = "data.frame")
  expect_s3_class(res, "result_summary_table")
  expect_true(all(c("param_summary", "delta_summary") %in% names(res)))
})
