#' Run full BayCauRETM pipeline: prepare, fit, and g-compute
#'
#' @param data Long-format dataset
#' @param i Subject ID column name
#' @param K Time index column name
#' @param Y_k Recurrent event column
#' @param T_k Terminal event column
#' @param C_k Censor indicator column
#' @param A Treatment column name
#' @param covariates Vector of baseline covariate column names
#' @param formula_T Formula for terminal event model
#' @param formula_Y Formula for recurrent event model
#' @param treatment_vals Treatment values to compare (e.g., c(0, 1))
#' @param iter Iterations for Stan
#' @param chains Number of chains
#' @param stan_model_path Path to Stan model
#' @param n_draws Number of posterior samples for g-computation
#' @param seed Random seed
#'
#' @return A list with model fit, posterior draws, phi results
#' @export
run_baycause <- function(data, i, K, Y_k, T_k, C_k, A, covariates,
                         formula_T, formula_Y,
                         treatment_vals = c(0, 1),
                         iter = 2000, chains = 4,
                         stan_model_path = system.file("stan/baycause_model.stan", package = "BayCauRETM"),
                         n_draws = 1000,
                         seed = 2025) {

  message("✅ Preparing Stan data...")
  stan_data <- make_stan_data(
    data = data,
    i = i, K = K,
    Y_k = Y_k, T_k = T_k, C_k = C_k, A = A,
    formula_T = formula_T, formula_Y = formula_Y
  )

  message("✅ Fitting Bayesian model in Stan...")
  fit <- rstan::stan(
    file = stan_model_path,
    data = stan_data,
    iter = iter,
    chains = chains,
    seed = seed
  )

  message("✅ Running g-computation...")
  g_result <- g_computation(
    fit = fit,
    data = data,
    i = i,
    K = K,
    A = A,
    covariates = covariates,
    treatment_vals = treatment_vals,
    n_draws = n_draws
  )

  message("✅ Done.")
  list(
    fit = fit,
    phi = g_result$phi,
    mean_phi = g_result$mean,
    ci = g_result$ci
  )
}
