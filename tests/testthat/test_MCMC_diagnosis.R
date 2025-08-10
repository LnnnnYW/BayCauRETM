# tests/testthat/test_MCMC_diagnosis.R

library(testthat)
library(BayCauRETM)

test_that("mcmc_diagnosis returns stats and trace plots", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")

  set.seed(123)
  sc  <- "parameters{real y;} model{y ~ normal(0,1);} "
  fit <- suppressWarnings(rstan::stan(model_code = sc, iter = 200, chains = 1, refresh = 0))

  fake_fit <- list(stan_fit = fit, data_preprocessed = data.frame(A = 0))
  diag <- mcmc_diagnosis(fake_fit, pars_to_check = "y",
                         save_plots = FALSE, positivity = FALSE)

  expect_s3_class(diag, "mcmc_diag")
  expect_true(is.data.frame(diag$stats))
  expect_equal(nrow(diag$stats), 1)
  expect_true("y" %in% names(diag$plots))
})

test_that("mcmc_diagnosis positivity diagnostics run when data present", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")

  set.seed(42)
  sc  <- "parameters{real y;} model{y ~ normal(0,1);} "
  fit <- suppressWarnings(rstan::stan(model_code = sc, iter = 150, chains = 1, refresh = 0))

  n <- 200
  df_pos <- data.frame(
    A      = rbinom(n, 1, 0.5),
    Y_obs  = rpois(n, 1),
    k_idx  = sample(1:3, n, TRUE)
  )

  fit_out <- list(stan_fit = fit, data_preprocessed = df_pos)

  diag_pos <- mcmc_diagnosis(fit_out, pars_to_check = "y",
                             save_plots = FALSE, positivity = TRUE,
                             ps_covariates = c("Y_obs","k_idx"))

  expect_s3_class(diag_pos, "mcmc_diag")
  expect_true(is.data.frame(diag_pos$stats))
  expect_equal(nrow(diag_pos$stats), 1)
})
