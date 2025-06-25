# tests/testthat/test_plot_wrappers.R

library(testthat)
library(BayCauRETM)

test_that("static & interactive wrapper functions return ggplot/plotly", {
  fake <- list(
    delta = list(
      "s=1" = list(mean = 0, CI_lower = -1, CI_upper = 1),
      "s=2" = list(mean = -.1, CI_lower = -1.2, CI_upper = 0.9)
    ),
    R_mat = matrix(0, 1, 2)
  )
  class(fake) <- "gcomp_out"

  p_static <- plot_posterior_causal_contrast_static(fake)
  expect_s3_class(p_static, "ggplot")

  p_wrap <- plot_posterior_causal_contrast_interactive(fake)
  expect_s3_class(p_wrap, "ggplot")

  if (requireNamespace("plotly", quietly = TRUE)) {
    p_int <- plot_posterior_causal_contrast_interactive(fake, interactive = TRUE)
    expect_true(inherits(p_int, "plotly"))
  }
})
