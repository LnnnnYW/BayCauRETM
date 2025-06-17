# tests/testthat/test-g_computation.R

test_that("g_computation errors with missing input", {
  expect_error(g_computation(list()), "stan_fit")
})


# tests/testthat/test-mcmc_diagnosis.R

test_that("mcmc_diagnosis runs without positivity", {
  skip_on_cran()
  # Create minimal fake stan_fit and data_preprocessed
  dummy_fit <- list(
    stan_fit           = NULL,
    data_preprocessed  = data.frame(
      pat_id = 1,
      k_idx  = 1,
      Y_prev = 0,
      T_prev = 0,
      C_prev = 0,
      A      = 0
    )
  )
  # Expect error because stan_fit NULL
  expect_error(mcmc_diagnosis(dummy_fit, positivity = FALSE), "stan_fit")
})
