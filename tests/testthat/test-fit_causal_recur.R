# tests/testthat/test-fit_causal_recur.R

library(testthat)
library(BayCauRETM)

test_that("fit_causal_recur runs on small dataset", {
  df <- data.frame(
    patient_id = rep(1:2, each = 2),
    k_idx      = rep(1:2, times = 2),
    Y          = rpois(4, 1),
    T          = rbinom(4, 1, 0.2),
    C          = rbinom(4, 1, 0.05),
    A          = rbinom(4, 1, 0.5)
  )

  prior <- list(
    eta_beta    = 0, sigma_beta  = 1, rho_beta   = 0.5,
    eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = 0.5
  )

  expect_silent({
    suppressMessages({
      suppressWarnings({
        invisible(capture.output({
          fit <- fit_causal_recur(
            data       = df,
            K          = 2,
            id_col     = "patient_id",
            k_col      = "k_idx",
            y_col      = "Y",
            t_col      = "T",
            c_col      = "C",
            a_col      = "A",
            x_cols     = NULL,
            formula_T  = T_obs ~ Y_prev + A + k_idx,
            formula_Y  = Y_obs ~ Y_prev + A + k_idx,
            prior      = prior,
            num_chains = 1,
            iter       = 50
          )
        }))
      })
    })
  })

  expect_true(inherits(fit$stan_fit, "stanfit"))
  expect_type(fit$data_preprocessed, "list")
})
