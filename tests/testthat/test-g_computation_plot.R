# tests/testthat/test_g_computation_plot.R

library(testthat)
library(BayCauRETM)

make_fake_gcomp <- function(K = 4, draws = 50) {
  delta <- lapply(1:K, function(s) {
    list(
      draws    = rnorm(draws, 0.05 * s, 0.01),
      mean     = 0.05 * s,
      CI_lower = 0.05 * s - 0.02,
      CI_upper = 0.05 * s + 0.02
    )
  })
  names(delta) <- paste0("s=", 1:K)
  R_mat <- matrix(runif(draws * (K + 1)), ncol = K + 1)
  structure(list(R_mat = R_mat, delta = delta), class = "gcomp_out")
}

test_that("plot.gcomp_out returns ggplot", {
  skip_if_not_installed("ggplot2")
  g <- make_fake_gcomp(4)
  p <- plot(g, ref_line = 0)
  expect_s3_class(p, "ggplot")
})
