# fit_causal_recur

#' Fit Bayesian Causal Recurrent and Terminal‑Event Model
#'
#' @description
#' Fits a discrete‑time Bayesian model for recurrent event counts and a terminal
#' event (death) using gAR(1) smoothing priors on the time‑varying intercepts.
#'
#' @param data A **long‑format** `data.frame` that contains the user‑named
#'   identifier, time index, treatment, outcome, and any covariate columns.
#'   These names are supplied via `id_col`, `time_col`, `treat_col`, plus the
#'   left‑hand sides of `formula_T` and `formula_Y`.
#' @param K Integer. Total number of discrete intervals in the study.
#' @param id_col Character scalar. Column name holding the **subject ID**.
#' @param time_col Character scalar. Column name holding the **discrete time
#'   index** (`1,...,K`).
#' @param treat_col Character scalar. Column name holding the **treatment
#'   indicator** (`0/1`).
#' @param x_cols Character vector of additional (static or time‑varying)
#'   covariate names to keep; `NULL` if none.
#' @param formula_T A formula for the terminal‑event (death) sub‑model, e.g.
#'   `death_flag ~ Y_prev + A + k_idx`.
#' @param formula_Y A formula for the recurrent‑count sub‑model, e.g.
#'   `event_count ~ Y_prev + A + k_idx`.
#' @param prior Named list of gAR(1) hyperparameters with elements
#'   `eta_beta`, `sigma_beta`, `rho_beta`, `eta_gamma`, `sigma_gamma`,
#'   `rho_gamma`.
#' @param num_chains Integer. Number of MCMC chains (default `4`).
#' @param iter Integer. Total iterations *per* chain including warm‑up
#'   (default `2000`).
#' @param stan_model_file Optional path to a pre‑compiled Stan model.
#'   If `NULL`, the package‑internal teacher‑style model is used.
#' @param control List passed to Stan sampling (see **rstan** docs).
#' @param verbose Logical. Print progress messages (default `TRUE`).
#' @inheritParams base::print
#'
#' @return An object of class `causal_recur_fit` (list) with elements
#'   `stan_fit`, `data_preprocessed`, `n_pat`, `K`, `design_info`, `prior`, and
#'   `stan_data_list`.
#'
#' @details
#' Internally the function
#' 1. Copies `id_col`, `time_col`, and `treat_col` into the canonical names
#'    `pat_id`, `k_idx`, and `A`, and copies the outcomes named on the left‑hand
#'    sides of `formula_T` / `formula_Y` into `T_obs` / `Y_obs`;
#' 2. Calls [preprocess_data()] to fill missing intervals, create lagged
#'    variables, and run basic checks;
#' 3. Builds teacher‑style standata, compiles the Stan model, and runs MCMC
#'    sampling via **rstan**.
#'
#' @examples
#' df <- data.frame(
#'   sid         = rep(1:2, each = 2),
#'   period      = rep(1:2, 2),
#'   event_count = rpois(4, 1),
#'   death_flag  = rbinom(4, 1, 0.2),
#'   trt_arm     = rbinom(4, 1, 0.5)
#' )
#' prior <- list(
#'   eta_beta  = 0, sigma_beta  = 1, rho_beta   = 0.5,
#'   eta_gamma = 0, sigma_gamma = 1, rho_gamma  = 0.5
#' )
#' \dontrun{
#' fit <- fit_causal_recur(
#'   data       = df, K = 2,
#'   id_col     = "sid",
#'   time_col   = "period",
#'   treat_col  = "trt_arm",
#'   formula_T  = death_flag  ~ Y_prev + A + k_idx,
#'   formula_Y  = event_count ~ Y_prev + A + k_idx,
#'   prior      = prior,
#'   num_chains = 1, iter = 100
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
                             id_col, time_col, treat_col,
                             x_cols = NULL,
                             formula_T, formula_Y,
                             prior,
                             num_chains = 4, iter = 2000,
                             stan_model_file = NULL,
                             compiled_model_cache = NULL,
                             control = list(adapt_delta = 0.95,
                                            max_treedepth = 15),
                             verbose = TRUE) {

  needed <- c(id_col, time_col, treat_col)
  if (any(miss <- !needed %in% names(data)))
    stop("Columns not found: ", paste(needed[miss], collapse = ", "))

  ## detect user lag variable
  rhs_vars <- union(all.vars(formula_T)[-1], all.vars(formula_Y)[-1])
  lag_var  <- rhs_vars[grepl("^lag", rhs_vars)]
  lag_var  <- if (length(lag_var) == 1) lag_var else "lagYk"

  ## copies
  df <- data
  event_col <- all.vars(formula_T)[1]
  count_col <- all.vars(formula_Y)[1]

  df$T_obs  <- df[[event_col]]
  df$Y_obs  <- df[[count_col]]
  df$pat_id <- df[[id_col]]
  df$k_idx  <- df[[time_col]]
  df$A      <- df[[treat_col]]

  keep_cols <- c("pat_id", "k_idx",
                 "T_obs", "Y_obs",
                 "A",
                 treat_col, time_col,
                 x_cols)

  df <- df[, keep_cols, drop = FALSE]

  ## preprocess
  prep <- preprocess_data(df, K = K, x_cols = x_cols, lag_col = lag_var)
  df   <- prep$processed_df
  n_pat <- prep$n_pat
  df[[event_col]] <- df$T_obs
  df[[count_col]] <- df$Y_obs
  df[[treat_col]] <- df$A
  df[[time_col]]  <- df$k_idx

  ## baseline design matrix
  baseline_df <- data |>
    dplyr::group_by(!!rlang::sym(id_col)) |>
    dplyr::slice_min(order_by = !!rlang::sym(time_col), n = 1) |>
    dplyr::ungroup()

  get_X <- function(d) {
    mm <- model.matrix(update(~ ., formula_T, 0 ~ .), data = d)
    mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  X_base <- get_X(baseline_df)
  P <- ncol(X_base)
  if (P == 0) X_base <- matrix(0, nrow(baseline_df), 0)
  X_all <- X_base[df$pat_id, , drop = FALSE]

  ## Stan data slices
  df_Y1 <- df[df$k_idx == 1, ]
  df_Yk <- df[df$k_idx > 1, ]
  NY1   <- nrow(df_Y1); NYk <- nrow(df_Yk); NTk <- nrow(df)

  stan_data <- list(
    NY1 = NY1, NYk = NYk, NTk = NTk, K = K, P = P,

    kvecT = df$k_idx,
    L_Tk  = if (P) X_all else matrix(0, NTk, 0),
    A_Tk  = df$A,
    Tk    = df$T_obs,

    kvecY = df_Yk$k_idx,
    L_Yk  = if (P) X_all[df$k_idx > 1, , drop = FALSE] else matrix(0, NYk, 0),
    lagYk = df_Yk[[lag_var]],
    A_Yk  = df_Yk$A,
    Yk    = df_Yk$Y_obs,

    L_Y1  = if (NY1 == 0) matrix(0, 0, P)
    else if (P) X_base else matrix(0, NY1, 0),
    A_Y1  = if (NY1 == 0) numeric(0) else df_Y1$A,
    Y1    = if (NY1 == 0) integer(0) else df_Y1$Y_obs
  )

  ## compile Stan model
  if (is.null(stan_model_file)) {
    pkg_dir <- system.file(package = "BayCauRETM")
    stan_model_file <- file.path(pkg_dir, "stan", "causal_recur_model.stan")
  }
  if (!is.null(compiled_model_cache) && file.exists(compiled_model_cache)) {
    if (verbose) message("Loading cached compiled Stan model …")
    stan_mod <- readRDS(compiled_model_cache)
  } else {
    if (verbose) message("Compiling Stan model …")
    stan_mod <- rstan::stan_model(stan_model_file, verbose = FALSE)
    if (!is.null(compiled_model_cache)) saveRDS(stan_mod, compiled_model_cache)
  }

  if (verbose)
    message(sprintf("Sampling (%d chains × %d iter) …", num_chains, iter))
  stan_fit <- rstan::sampling(
    stan_mod, data = stan_data,
    chains = num_chains, iter = iter,
    seed = 1234, control = control
  )

  structure(
    list(stan_fit          = stan_fit,
         data_preprocessed = df,
         n_pat             = n_pat,
         K                 = K,
         design_info       = list(formula_T = formula_T,
                                  formula_Y = formula_Y,
                                  covariate_names = colnames(X_base)),
         prior             = prior,
         stan_data_list    = stan_data),
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
                                     pars_to_report = c("beta1", "theta1", "theta_lag",
                                                        "betaL[1]", "thetaL[1]"),
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
