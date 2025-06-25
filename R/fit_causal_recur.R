# fit_causal_recur

#' Fit Bayesian Causal Recurrent- and Terminal-Event Model
#'
#' @description
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
#' @return An object of class \code{causal_recur_fit}, a list with components:
#'   \describe{
#'     \item{stan_fit}{An \code{rstan::stanfit} object with posterior samples.}
#'     \item{data_preprocessed}{Processed data.frame with lagged columns.}
#'     \item{n_pat}{Number of unique subjects.}
#'     \item{K}{Number of intervals.}
#'     \item{design_info}{List with formula and covariate names for submodels.}
#'     \item{prior}{Prior hyperparameters as supplied.}
#'     \item{stan_data_list}{List passed to Stan.}
#'   }
#'   This object has class \code{causal_recur_fit}, so methods \code{print()},
#'   \code{summary()}, and \code{plot()} are available.
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
#' prior <- list(
#'   eta_beta    = 0, sigma_beta  = 1, rho_beta   = 0.5,
#'   eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = 0.5
#' )
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
#' summary(fit)
#' plot(fit)
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
  # Drop intercept columns if present
  if ("(Intercept)" %in% colnames(X_T_mat)) {
    X_T <- X_T_mat[, setdiff(colnames(X_T_mat), "(Intercept)"), drop = FALSE]
  } else {
    X_T <- X_T_mat
  }
  if ("(Intercept)" %in% colnames(X_Y_mat)) {
    X_Y <- X_Y_mat[, setdiff(colnames(X_Y_mat), "(Intercept)"), drop = FALSE]
  } else {
    X_Y <- X_Y_mat
  }
  p_T <- ncol(X_T); p_Y <- ncol(X_Y)

  # Record design info
  design_info <- list(
    formula_T = formula_T,
    formula_Y = formula_Y,
    covariate_names_T = colnames(X_T),
    covariate_names_Y = colnames(X_Y)
  )

  # Locate built-in Stan model file
  if (is.null(stan_model_file)) {
    pkg_dir <- system.file(package = "BayCauRETM")
    stan_model_file <- file.path(pkg_dir, "stan", "causal_recur_model.stan")
    if (!file.exists(stan_model_file)) {
      stop("Cannot find inst/stan/causal_recur_model.stan in package BayCauRETM")
    }
  }
  if (verbose) message("Compiling Stan model (this may take a moment)...")
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

  if (verbose) message(sprintf("Starting MCMC sampling (chains=%d, iter=%d)...", num_chains, iter))
  stan_fit <- rstan::sampling(
    object = stan_mod,
    data = stan_data,
    chains = num_chains,
    iter = iter,
    seed = 1234,
    control = control
  )

  # Assemble result and assign class
  result <- list(
    stan_fit = stan_fit,
    data_preprocessed = df,
    n_pat = n_pat,
    K = K,
    design_info = design_info,
    prior = prior,
    stan_data_list = stan_data
  )
  class(result) <- "causal_recur_fit"
  return(result)
}

#' Print method for causal_recur_fit
#'
#' @description
#' Print a brief summary of a causal recurrent-event model fit.
#'
#' @param x An object of class \code{causal_recur_fit}, output of \code{fit_causal_recur()}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the object itself.
#' @rdname causal_recur_fit-print
#' @method print causal_recur_fit
#' @export
print.causal_recur_fit <- function(x, ...) {
  cat("Causal recurrent-event model fit\n")
  if (!is.null(x$n_pat)) {
    cat("  Number of subjects:", x$n_pat, "\n")
  }
  if (!is.null(x$K)) {
    cat("  Number of intervals K:", x$K, "\n")
  }
  cat("Use summary() for posterior parameter estimates.\n")
  cat("Use mcmc_diagnosis() for MCMC convergence diagnostics.\n")
  invisible(x)
}

#' Summary method for causal_recur_fit
#'
#' @description
#' Summarize posterior estimates for a causal recurrent-event model fit.
#'
#' @param object An object of class \code{causal_recur_fit}.
#' @param pars_to_report Character vector of parameter names or patterns to report
#'   (default: c("beta_Y","beta_A","gamma_Y","gamma_A")).
#' @param ... Additional arguments passed to \code{rstan::summary()}.
#' @return Invisibly returns a data.frame of posterior means and 95% intervals.
#' @rdname causal_recur_fit-summary
#' @method summary causal_recur_fit
#' @export
summary.causal_recur_fit <- function(object,
                                     pars_to_report = c("beta_Y","beta_A","gamma_Y","gamma_A"),
                                     ...) {
  stan_fit <- object$stan_fit
  sum_obj <- tryCatch(rstan::summary(stan_fit, pars = pars_to_report, ...),
                      error = function(e) NULL)
  if (is.null(sum_obj) || is.null(sum_obj$summary)) {
    cat("No summary available for specified parameters.\n")
    return(invisible(NULL))
  }
  sum_stan <- sum_obj$summary
  df <- data.frame(
    Parameter = rownames(sum_stan),
    Mean      = sum_stan[, "mean"],
    `2.5%`    = sum_stan[, "2.5%"],
    `97.5%`   = sum_stan[, "97.5%"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  print(df, row.names = FALSE)
  invisible(df)
}

#' Plot method for causal_recur_fit
#'
#' @description
#' For a \code{causal_recur_fit} object, prompt the user to run MCMC diagnostics.
#'
#' @param x An object of class \code{causal_recur_fit}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns NULL.
#' @rdname causal_recur_fit-plot
#' @method plot causal_recur_fit
#' @export
plot.causal_recur_fit <- function(x, ...) {
  cat("To check MCMC convergence, please run:\n")
  cat("  mcmc_diagnosis(fit_out, pars_to_check = ..., save_plots = ..., positivity = ...)\n")
  invisible(NULL)
}

