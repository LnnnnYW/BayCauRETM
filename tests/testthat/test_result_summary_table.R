# tests/testthat/test_result_summary_table.R

library(testthat)
library(BayCauRETM)

make_fake_gcomp <- function(K = 3) {
  delta <- lapply(1:K, function(s) {
    list(mean = 0.1 * s, CI_lower = 0.1 * s - 0.02, CI_upper = 0.1 * s + 0.02, draws = rnorm(20))
  })
  names(delta) <- paste0("s=", 1:K)
  structure(
    list(R_mat = matrix(runif(20 * (K + 1)), ncol = K + 1), delta = delta),
    class = "gcomp_out"
  )
}

test_that("result_summary_table merges parameter and delta summaries", {
  g <- make_fake_gcomp(3)
  res <- result_summary_table(
    fit_out = list(), gcomp_out = g,
    pars_to_report = "beta1", s_vec = 1:3, format = "data.frame"
  )
  expect_s3_class(res, "result_summary_table")
  expect_true(is.data.frame(res$param_summary))
  expect_true(is.data.frame(res$delta_summary))
  expect_true(all(c("param_table","delta_table") %in% names(res)))
})
