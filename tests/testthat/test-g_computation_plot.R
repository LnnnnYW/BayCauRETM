# tests/testthat/test_g_computation_plot.R

library(testthat)
library(BayCauRETM)

test_that("gcomp_out plotting via S3 works", {
  fake <- list(
    R_mat = matrix(0, 1, 3),
    delta = list(
      "s=1" = list(mean = 0,   CI_lower = -1, CI_upper = 1),
      "s=2" = list(mean = 0.2, CI_lower = -0.8, CI_upper = 1.4)
    )
  )
  class(fake) <- "gcomp_out"

  p <- plot(fake)
  expect_s3_class(p, "ggplot")

  if (requireNamespace("plotly", quietly = TRUE)) {
    pint <- plot(fake, interactive = TRUE)
    expect_true(inherits(pint, "plotly"))
  }
})

