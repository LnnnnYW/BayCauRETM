# tests/testthat/test-summary_visualization.R
library(testthat)
library(BayCauRETM)

mock_stanfit <- methods::new("stanfit")


test_that("summarize_results works on fake output", {
  fake_out <- list(
    stan_fit = list(),   # minimal stanfit mock
    delta    = list("s=1" = list(mean=1, CI_lower=0, CI_upper=2),
                    "s=2" = list(mean=2, CI_lower=1, CI_upper=3))
  )
  df <- summarize_results(fit_out = list(stan_fit = mock_stanfit),
                          gcomp_out = list(delta = fake_out$delta))
  expect_s3_class(df$param_summary, "data.frame")
  expect_s3_class(df$delta_summary, "data.frame")
  expect_equal(nrow(df$delta_summary), 2)
})

test_that("plot_delta_vs_s returns ggplot", {
  fake_delta <- list("s=1"=list(mean=1,CI_lower=0,CI_upper=2),
                     "s=2"=list(mean=2,CI_lower=1,CI_upper=3))
  p <- plot_delta_vs_s(gcomp_out = list(delta = fake_delta), s_vec = 1:2)
  expect_s3_class(p, "ggplot")
})

