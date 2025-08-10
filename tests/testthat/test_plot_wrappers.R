# tests/testthat/test_plot_wrappers.R

library(testthat)
library(BayCauRETM)

make_fake_gcomp <- function(K = 3, draws = 40) {
  delta <- lapply(1:K, function(s) {
    list(
      draws    = rnorm(draws, 0.03 * s, 0.01),
      mean     = 0.03 * s,
      CI_lower = 0.03 * s - 0.02,
      CI_upper = 0.03 * s + 0.02
    )
  })
  names(delta) <- paste0("s=", 1:K)
  R_mat <- matrix(runif(draws * (K + 1)), ncol = K + 1)
  structure(list(R_mat = R_mat, delta = delta), class = "gcomp_out")
}

test_that("static wrapper returns ggplot", {
  skip_if_not_installed("ggplot2")
  g <- make_fake_gcomp(3)
  p <- plot_posterior_causal_contrast_static(g, s_vec = 1:3, ref_line = 0)
  expect_s3_class(p, "ggplot")
})

test_that("interactive wrapper returns plotly when requested", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("plotly")
  g <- make_fake_gcomp(3)
  p <- plot_posterior_causal_contrast_interactive(g, interactive = TRUE)
  expect_true("plotly" %in% class(p))
})
