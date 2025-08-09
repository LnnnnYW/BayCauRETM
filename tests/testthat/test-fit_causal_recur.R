# tests/testthat/test-fit_causal_recur.R

library(testthat)
library(BayCauRETM)

test_that("fit_causal_recur runs on tiny toy data (no-variation warnings suppressed)", {

  skip_on_cran()
  skip_if_not_installed("rstan")

  stan_code <- "parameters{real y;} model{y ~ normal(0,1);}"
  tmp_stan  <- tempfile(fileext = ".stan")
  writeLines(stan_code, tmp_stan)
  withr::defer(unlink(tmp_stan), teardown_env())

  set.seed(1)
  df <- expand.grid(pat_id = 1:2, k_idx = 1:2)
  n  <- nrow(df)
  df$Y_obs <- rpois(n, 1)
  df$T_obs <- rbinom(n, 1, 0.30)
  df$A     <- rbinom(n, 1, 0.40)

  prior <- list(
    eta_beta  = 0, sigma_beta  = 1, rho_beta  = 0.5,
    eta_gamma = 0, sigma_gamma = 1, rho_gamma = 0.5
  )

  fit <- suppressWarnings(
    fit_causal_recur(
      data            = df,
      K               = 2,
      formula_T       = T_obs ~ Y_prev + A + k_idx,
      formula_Y       = Y_obs ~ Y_prev + A + k_idx,
      prior           = prior,
      stan_model_file = tmp_stan,
      num_chains      = 1,
      iter            = 20,
      verbose         = FALSE
    )
  )

  expect_s3_class(fit, "causal_recur_fit")
  expect_true(inherits(fit$stan_fit, "stanfit"))
})

