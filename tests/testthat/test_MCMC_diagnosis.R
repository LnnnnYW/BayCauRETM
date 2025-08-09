# tests/testthat/test_MCMC_diagnosis.R

library(testthat)
library(BayCauRETM)

test_that("mcmc_diagnosis returns stats list", {
  skip_on_cran()
  skip_if_not_installed("rstan")

  sc <- "parameters{real y;} model{y~normal(0,1);}"
  fit <- rstan::stan(model_code = sc, iter = 100, chains = 1, refresh = 0)

  fake_fit <- list(stan_fit = fit, data_preprocessed = data.frame(A = 0))
  diag <- mcmc_diagnosis(fake_fit, pars_to_check = "y",
                         save_plots = FALSE, positivity = FALSE)
  expect_s3_class(diag, "mcmc_diag")
  expect_equal(nrow(diag$stats), 1)
})
