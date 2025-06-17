#' Fit Bayesian Causal Recurrent- and Terminal-Event Model
#'
#' Fit a discrete-time Bayesian model for recurrent event counts and a terminal event (death),
#' with gAR(1) smoothing priors on time-varying intercepts.
#'
#' @param data A data.frame in long format, one row per subject and interval.
#' @param K Integer. Total number of discrete intervals in the study.
#' @param id_col Character. Name of the subject identifier column.
#' @param k_col Character. Name of the interval index column (1..K).
#' @param y_col Character. Name of the recurrent event count outcome column.
#' @param t_col Character. Name of the terminal event indicator column (0/1).
#' @param c_col Character. Name of the censoring indicator column (0/1).
#' @param a_col Character. Name of the treatment indicator column (0/1 or factor).
#' @param x_cols Character vector of additional covariate names, or NULL.
#' @param formula_T Formula for the death submodel, e.g. \code{T_obs ~ Y_prev + A + k_idx}.
#' @param formula_Y Formula for the count submodel, e.g. \code{Y_obs ~ Y_prev + A + k_idx}.
#' @param prior Named list of gAR(1) hyperparameters:
#'   \describe{
#'     \item{eta_beta}{Mean for death intercepts.}
#'     \item{sigma_beta}{SD for death intercepts.}
#'     \item{rho_beta}{AR(1) correlation for death intercepts.}
#'     \item{eta_gamma}{Mean for count intercepts.}
#'     \item{sigma_gamma}{SD for count intercepts.}
#'     \item{rho_gamma}{AR(1) correlation for count intercepts.}
#'   }
#' @param num_chains Integer. Number of MCMC chains (default: 4).
#' @param iter Integer. Total iterations per chain, including warmup (default: 2000).
#' @param stan_model_file Optional character. Path to a pre-compiled Stan model file.
#' @param control Optional list. Passed to Stan sampling (see \code{cmdstanr} docs).
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return A list with components:
#'   \describe{
#'     \item{stan_fit}{An rstan::stanfit object with posterior samples.}
#'     \item{data_preprocessed}{Processed data.frame with lagged columns.}
#'     \item{n_pat}{Number of unique subjects.}
#'     \item{design_info}{List with formula and covariate names.}
#'     \item{prior}{Prior hyperparameters as supplied.}
#'     \item{stan_data_list}{List passed to Stan.}
#'   }
#'
#' @details
#' Internally, this function calls \code{preprocess_data()}, builds design matrices,
#' compiles or loads a Stan model, and runs MCMC sampling using \code{rstan::sampling()} or \code{cmdstanr}.
#'
#' @examples
#' # Minimal example with fake data
#' df <- data.frame(
#'   patient_id = rep(1:2, each = 2),
#'   k_idx      = rep(1:2, 2),
#'   Y          = rpois(4, 1),
#'   T          = rbinom(4, 1, 0.2),
#'   C          = rbinom(4, 1, 0.05),
#'   A          = rbinom(4, 1, 0.5),
#'   X1         = rnorm(4)
#' )
#'
#' prior <- list(
#'   eta_beta    = 0, sigma_beta  = 1, rho_beta   = 0.5,
#'   eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = 0.5
#' )
#'
#' \dontrun{
#' fit <- fit_causal_recur(
#'   data       = df,
#'   K          = 2,
#'   id_col     = "patient_id",
#'   k_col      = "k_idx",
#'   y_col      = "Y",
#'   t_col      = "T",
#'   c_col      = "C",
#'   a_col      = "A",
#'   x_cols     = "X1",
#'   formula_T  = T_obs ~ Y_prev + A + k_idx,
#'   formula_Y  = Y_obs ~ Y_prev + A + k_idx,
#'   prior      = prior,
#'   num_chains = 1,
#'   iter       = 100
#' )
#' print(fit)
#' }
#'
#' @importFrom rstan stan_model sampling
#' @importFrom stats model.matrix
#' @importFrom dplyr mutate
#' @export







fit_causal_recur <- function(data,
                             K,
                             id_col,
                             k_col,
                             y_col,
                             t_col,
                             c_col,
                             a_col,
                             x_cols = NULL,
                             formula_T,
                             formula_Y,
                             prior,
                             num_chains = 4,
                             iter = 2000,
                             stan_model_file = NULL,
                             control = list(adapt_delta = 0.95, max_treedepth = 15),
                             verbose = TRUE) {
  # Preprocess data
  prep <- preprocess_data(data, id_col, k_col, y_col, t_col, c_col, a_col, x_cols, K)
  df <- prep$processed_df
  n_pat <- prep$n_pat
  N     <- nrow(df)

  # Ensure numeric history and treatment
  df <- df %>% mutate(Y_prev = as.numeric(Y_prev), A = as.numeric(A))

  # Build design matrices for each submodel
  X_T_mat <- model.matrix(formula_T, data = df)
  X_Y_mat <- model.matrix(formula_Y, data = df)
  # Drop intercept columns
  X_T <- as.matrix(X_T_mat[, colnames(X_T_mat) != "(Intercept)"])
  X_Y <- as.matrix(X_Y_mat[, colnames(X_Y_mat) != "(Intercept)"])
  p_T <- ncol(X_T); p_Y <- ncol(X_Y)

  # Record design info
  design_info <- list(
    formula_T = formula_T,
    formula_Y = formula_Y,
    covariate_names_T = colnames(X_T),
    covariate_names_Y = colnames(X_Y)
  )

  # Locate built in Stan model file
   if (is.null(stan_model_file)) {
     stan_model_file <- system.file("stan", "causal_recur_model.stan", package = "BayCauRETM")
     if (!nzchar(stan_model_file) || !file.exists(stan_model_file)) {
       stop("Stan model file not found in inst/stan/causal_recur_model.stan")
     }
     }
   message("Compiling Stan model (this may take a moment)...")
   stan_mod <- rstan::stan_model(file = stan_model_file, verbose = FALSE)

  # Assemble data for Stan
  stan_data <- list(
    N = N, n_pat = n_pat, K = K,
    p_T = p_T, p_Y = p_Y,
    pat_id = df$pat_id, k_idx = df$k_idx,
    T_prev = df$T_prev, C_prev = df$C_prev, A = df$A,
    X_T = X_T, Y_prev = df$Y_prev,
    X_Y = X_Y, T_obs = df$T_obs, Y_obs = df$Y_obs,
    eta_beta = prior$eta_beta,
    sigma_beta = prior$sigma_beta,
    rho_beta = prior$rho_beta,
    eta_gamma = prior$eta_gamma,
    sigma_gamma = prior$sigma_gamma,
    rho_gamma = prior$rho_gamma
  )

  if (verbose) message("Starting MCMC sampling (chains=", num_chains, ", iter=", iter, ")...")
  stan_fit <- rstan::sampling(
    object = stan_mod,
    data = stan_data,
    chains = num_chains,
    iter = iter,
    seed = 1234,
    control = control
  )

  list(
    stan_fit = stan_fit,
    data_preprocessed = df,
    n_pat = n_pat,
    K = K,
    design_info = design_info,
    prior = prior,
    stan_data_list = stan_data
  )
}
