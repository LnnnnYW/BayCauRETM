# fit_causal_recur

#' Fit Bayesian Causal Recurrent - and Terminal - Event Model
#'
#' @description
#' Fit a discrete - time Bayesian model for recurrent event counts and a terminal
#' event (death), using gAR(1) smoothing priors on the time - varying intercepts.
#'
#' @param data A **long-format** `data.frame` that already contains the standard
#'   analysis columns: `pat_id`, `k_idx`, `Y_obs`, `T_obs`, `A`, and *optionally*
#'   extra covariates.  The function will run a light version of
#'   `preprocess_data()` to add `Y_prev` and `T_prev` (and to fill in any missing
#'   intervals) but **will not** rename columns.
#' @param K Integer. Total number of discrete intervals in the study.
#' @param x_cols Character vector of additional (static or time - varying)
#'   covariate names to keep; `NULL` if none.
#' @param formula_T A formula for the terminal-event (death) sub‑model, e.g.
#'   `T_obs ~ Y_prev + A + k_idx`.
#' @param formula_Y A formula for the recurrent-count sub‑model, e.g.
#'   `Y_obs ~ Y_prev + A + k_idx`.
#' @param prior Named list of gAR(1) hyperparameters with elements
#'   `eta_beta`, `sigma_beta`, `rho_beta`, `eta_gamma`, `sigma_gamma`,
#'   `rho_gamma`.
#' @param num_chains Integer. Number of MCMC chains (default `4`).
#' @param iter Integer. Total iterations *per* chain including warm-up (default
#'   `2000`).
#' @param stan_model_file Optional path to a pre-compiled Stan model
#'   (`.stan` -> `*.rds`).  If `NULL`, the package-internal model is used.
#' @param control List passed to Stan sampling (see **cmdstanr** docs).
#' @param verbose Logical. Print progress messages (default `TRUE`).
#' @param pars_to_report Character vector of parameter names to include in the summary/plot S3 methods.
#' @param object Object returned by this function (used by summary method).
#' @inheritParams base::print
#'
#' @return An object of class `causal_recur_fit` (list) with elements
#'   `stan_fit`, `data_preprocessed`, `n_pat`, `K`, `design_info`, `prior`, and
#'   `stan_data_list`.
#'
#' @details
#' Internally, the function first calls [preprocess_data()] (which assumes the
#' input already uses the `_obs` / `_prev` naming convention), then builds the
#' required design matrices, compiles or loads the Stan model, and finally runs
#' MCMC sampling via **rstan**.
#'
#' @examples
#' df <- data.frame(
#'   pat_id = rep(1:2, each = 2),
#'   k_idx  = rep(1:2, 2),
#'   Y_obs  = rpois(4, 1),
#'   T_obs  = rbinom(4, 1, 0.2),
#'   A      = rbinom(4, 1, 0.5)
#' )
#' prior <- list(
#'   eta_beta    = 0, sigma_beta  = 1, rho_beta   = 0.5,
#'   eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = 0.5
#' )
#' \dontrun{
#' fit <- fit_causal_recur(
#'   data       = df,
#'   K          = 2,
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
#' @name causal_recur_fit
#' @docType class
#' @export
fit_causal_recur <- function(data,
                             K,
                             x_cols = NULL,
                             formula_T,
                             formula_Y,
                             prior,
                             num_chains = 4,
                             iter = 2000,
                             stan_model_file = NULL,
                             control = list(adapt_delta = 0.95,
                                            max_treedepth = 15),
                             verbose = TRUE) {

  df <- data

  # 1) Identify the user’s event and count columns from the left-hand sides of the formulas
  event_col <- all.vars(formula_T)[1]   # e.g. "death"
  count_col <- all.vars(formula_Y)[1]   # e.g. "num_evt"

  # Copy into the internal standard names T_obs and Y_obs
  df$T_obs <- df[[event_col]]
  df$Y_obs <- df[[count_col]]

  # 2) Find the set of covariates that appear on the right-hand sides of both formulas
  rhs_T      <- setdiff(all.vars(update(formula_T, 0 ~ .)), event_col)
  rhs_Y      <- setdiff(all.vars(update(formula_Y, 0 ~ .)), count_col)
  common_rhs <- intersect(rhs_T, rhs_Y)

  # 3) Automatically detect the treatment variable A:
  if ("A" %in% names(df)) {
    # use df$A as is
  } else if (length(common_rhs) == 1) {
    df$A <- df[[common_rhs]]
    message("Auto‐detected treatment column: ", common_rhs, "A")
    common_rhs <- setdiff(common_rhs, common_rhs)
  } else {
    stop("Could not uniquely identify treatment column; please supply a column named 'A'.")
  }

  # 4) Automatically detect the time index k_idx:
  if (!"k_idx" %in% names(df) && length(common_rhs) == 1) {
    df$k_idx <- df[[common_rhs]]
    message("Auto‐detected time index: ", common_rhs, "k_idx")
    common_rhs <- setdiff(common_rhs, common_rhs)
  } else if (!"k_idx" %in% names(df)) {
    stop("Could not uniquely identify time index; please supply a column named 'k_idx'.")
  }

  # 5) Automatically detect the subject ID pat_id:
  used_cols <- c(event_col, count_col, "A", "k_idx", x_cols)
  extras    <- setdiff(names(df), c(used_cols, "T_obs", "Y_obs"))
  if (!"pat_id" %in% names(df) && length(extras) == 1) {
    df$pat_id <- df[[extras]]
    message("Auto‐detected subject ID: ", extras, "pat_id")
  } else if (!"pat_id" %in% names(df)) {
    stop("Could not uniquely identify subject ID; please supply a column named 'pat_id'.")
  }

  # 6) Drop any other columns that are not T_obs, Y_obs, A, k_idx, pat_id, or x_cols
  keep_cols <- c("pat_id", "k_idx", "T_obs", "Y_obs", "A", x_cols)
  df <- df[, keep_cols, drop = FALSE]

  #  Preprocess
  prep <- preprocess_data(df = df, K = K, x_cols = x_cols)
  df   <- prep$processed_df
  n_pat <- prep$n_pat
  N     <- nrow(df)

  # Copy back to user column names for model.matrix()
  df[[event_col]] <- df$T_obs
  df[[count_col]] <- df$Y_obs

  # Design matrices
  df <- df %>%
    dplyr::mutate(Y_prev = as.numeric(Y_prev),
                  A      = as.numeric(A))

  drop_int <- function(mat) {
    if ("(Intercept)" %in% colnames(mat))
      mat[, setdiff(colnames(mat), "(Intercept)"), drop = FALSE]
    else mat
  }
  X_T <- drop_int(model.matrix(formula_T, data = df))
  X_Y <- drop_int(model.matrix(formula_Y, data = df))
  p_T <- ncol(X_T); p_Y <- ncol(X_Y)

  design_info <- list(
    formula_T = formula_T,
    formula_Y = formula_Y,
    covariate_names_T = colnames(X_T),
    covariate_names_Y = colnames(X_Y)
  )

  # Stan model & sampling
  if (is.null(stan_model_file)) {
    pkg_dir <- system.file(package = "BayCauRETM")
    stan_model_file <- file.path(pkg_dir, "stan",
                                 "causal_recur_model.stan")
  }
  if (verbose) message("Compiling Stan model...")
  stan_mod <- rstan::stan_model(file = stan_model_file, verbose = FALSE)

  stan_data <- list(
    N       = N,
    n_pat   = n_pat,
    K       = K,
    p_T     = p_T,
    p_Y     = p_Y,
    pat_id  = df$pat_id,
    k_idx   = df$k_idx,
    T_prev  = df$T_prev,
    A       = df$A,
    X_T     = X_T,
    Y_prev  = df$Y_prev,
    X_Y     = X_Y,
    T_obs   = df$T_obs,
    Y_obs   = df$Y_obs,
    eta_beta    = prior$eta_beta,
    sigma_beta  = prior$sigma_beta,
    rho_beta    = prior$rho_beta,
    eta_gamma   = prior$eta_gamma,
    sigma_gamma = prior$sigma_gamma,
    rho_gamma   = prior$rho_gamma
  )

  if (verbose)
    message(sprintf("Sampling (%d chains × %d iter)...",
                    num_chains, iter))
  stan_fit <- rstan::sampling(
    object   = stan_mod,
    data     = stan_data,
    chains   = num_chains,
    iter     = iter,
    seed     = 1234,
    control  = control
  )

  structure(
    list(
      stan_fit          = stan_fit,
      data_preprocessed = df,
      n_pat             = n_pat,
      K                 = K,
      design_info       = design_info,
      prior             = prior,
      stan_data_list    = stan_data
    ),
    class = "causal_recur_fit"
  )
}



# print / summary / plot methods

#' @describeIn causal_recur_fit Print a brief summary
#' @export
print.causal_recur_fit <- function(x, ...) {
  cat("Causal recurrent event model fit\n")
  if (!is.null(x$n_pat)) cat("  Number of subjects:", x$n_pat, "\n")
  if (!is.null(x$K))     cat("  Number of intervals K:", x$K, "\n")
  cat("Use summary() for posterior parameter estimates.\n")
  cat("Use mcmc_diagnosis() for MCMC convergence diagnostics.\n")
  invisible(x)
}

#' @describeIn causal_recur_fit Summarise posterior estimates
#' @method summary causal_recur_fit
#' @export
summary.causal_recur_fit <- function(object,
                                     pars_to_report = c("beta_Y", "beta_A", "gamma_Y", "gamma_A"),
                                     ...) {
  stan_fit <- object$stan_fit
  sum_obj  <- tryCatch(rstan::summary(stan_fit, pars = pars_to_report, ...),
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

#' @describeIn causal_recur_fit Prompt user to run diagnostics
#' @export
plot.causal_recur_fit <- function(x, ...) {
  cat("To check MCMC convergence, please run:\n")
  cat("  mcmc_diagnosis(fit_out, pars_to_check = ..., save_plots = ..., positivity = ...)\n")
  invisible(NULL)
}
