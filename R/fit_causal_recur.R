#' Fit Bayesian Causal Recurrent and Terminal-Event Model
#'
#' @description
#' Fits a discrete-time Bayesian model for recurrent counts and a terminal
<<<<<<< HEAD
#' event using gAR(1) smoothing priors on time-varying intercepts. This
#' implementation expects a **pre-compiled** Stan model (`.rds`) and uses
#' `rstan::sampling()` for MCMC.
#'
#' @param data A **long-format** `data.frame` that contains the user-named
#'   identifier, time index, treatment, outcome, and any covariate columns.
#'   These names are supplied via `id_col`, `time_col`, `treat_col`, plus the
#'   left-hand sides of `formula_T` and `formula_Y`.
#' @param K Integer. Total number of discrete intervals in the study.
#' @param id_col Character scalar. Column name holding the **subject ID**.
#' @param time_col Factor scalar. Column name holding the **discrete time
#'   index** (`1,...,K`).
#' @param treat_col Character scalar. Column name holding the **treatment
#'   indicator** (`0/1`).
#' @param formula_T A formula for the terminal-event (death) sub-model, e.g.
#'   `death_flag ~ Y_prev + A + k_idx`. The **right-hand side terms excluding
#'   the treatment column** are used to build a (possibly empty) design matrix
#'   for time-varying/lagged predictors in the terminal sub-model.
#' @param formula_Y A formula for the recurrent-count sub-model, e.g.
#'   `event_count ~ Y_prev + A + k_idx`. The **right-hand side terms excluding
#'   the treatment column** are used to build a (possibly empty) design matrix
#'   for the recurrent sub-model (for `k_idx > 1` rows).
=======
#' event using gAR(1) smoothing priors on time-varying intercepts. A pre-compiled
#' Stan model (.rds) is loaded from the package and MCMC is run via `rstan::sampling()`.
#'
#' @param data A long-format data.frame containing user-named identifier, time index,
#'   treatment, outcomes, and covariates. The column names for ID, time, and treatment
#'   are provided by `id_col`, `time_col`, and `treat_col`. The left-hand sides of
#'   `formula_T` and `formula_Y` give the terminal-event and recurrent-count variables.
#' @param K Integer. Total number of discrete intervals in the study.
#' @param id_col Character scalar.Column name for the subject identifier.
#' @param time_col Character scalar. Name of the column that stores the discrete
#'   time index per subject. Values in that column should be integers 1..K.
#' @param treat_col Character scalar. Column name for the binary treatment indicator (0/1).
#' @param formula_T A formula for the terminal-event sub-model, for example
#'   `death_flag ~ lagYk + A + k_idx`. Terms on the right-hand side, excluding the
#'   treatment column, form the design matrix for the terminal model.
#' @param formula_Y A formula for the recurrent-count sub-model, for example
#'   `event_count ~ lagYk + A + k_idx`. Terms on the right-hand side, excluding the
#'   treatment column, form the design matrices for the recurrent model (separately for
#'   baseline rows with `k_idx == 1` and follow-up rows with `k_idx > 1`).
>>>>>>> 759e44d2adaa77e0e2258aa398e3f620bf7ff1ce
#' @param prior Named list of gAR(1) hyperparameters. Supported elements:
#'   `eta_beta`, `sigma_beta`, `rho_beta`, `eta_gamma`, `sigma_gamma`,
#'   `rho_gamma`, `sigma_beta1`, `sigma_theta1`, `sigma_theta_lag`.
#'   Missing entries fall back to internal defaults.
#' @param num_chains Integer. Number of MCMC chains (default `4`).
<<<<<<< HEAD
#' @param iter Integer. Total iterations *per* chain including warm-up
=======
#' @param iter Integer. Total iterations per chain including warm-up
>>>>>>> 759e44d2adaa77e0e2258aa398e3f620bf7ff1ce
#'   (default `2000`).
#' @param control List passed to `rstan::sampling()` (e.g.,
#'   `list(adapt_delta = 0.95, max_treedepth = 15)`).
#' @param cores Integer. Number of CPU cores to use for sampling (passed to
#'   `rstan::sampling()`).
<<<<<<< HEAD
#' @param verbose Logical. Print progress messages (default `TRUE`).
#' @param lag_col Character scalar or `NULL`. If `NULL`, a column named
#'   `"lagYk"` is created and initialized to 0 before preprocessing. If a name
#'   is provided and that column is **absent** in `data`, `preprocess_data()`
#'   will generate it as an integer indicator based on the subject-specific lag
#'   of `I(Y_obs > 0)`. If the column already exists, it is left unchanged.
=======
#' @param verbose Logical. Index of whether to print progress messages (default `TRUE`).
#' @param lag_col Character scalar. Name of the lag indicator if lagged terms are used in either
#'   formula. If the formulas include this name but the column is absent in `data`, the function
#'   calls `preprocess_data()` to generate it within subject as an indicator based on
#'   the previous interval outcome; otherwise the input is left unchanged. Default "lagYk".
>>>>>>> 759e44d2adaa77e0e2258aa398e3f620bf7ff1ce
#' @inheritParams base::print
#'
#' @return An object of class `causal_recur_fit` (list) with elements
#'   `stan_fit`, `data_preprocessed`, `n_pat`, `K`, `design_info`, `prior`, and
#'   `stan_data_list`.
#'
#' @details
#' Internally the function:
<<<<<<< HEAD
#' 1. Copies `id_col`, `time_col`, and `treat_col` into the canonical names
#'    `pat_id`, `k_idx`, and `A`, and copies the outcomes named on the left-hand
#'    sides of `formula_T` / `formula_Y` into `T_obs` / `Y_obs`;
#' 2. Calls `preprocess_data()` for **row ordering by (`pat_id`,`k_idx`)**,
#'    **subject-ID remapping**, and **optional lag-column creation** when the
#'    requested `lag_col` is absent. It **does not** pad to a full grid,
#'    **does not** carry treatment forward, and **does not** truncate after
#'    terminal events in the current implementation;
#' 3. Constructs design matrices from the right-hand sides of `formula_T` and
#'    `formula_Y` **excluding** the treatment column; missing values are set to
#'    zero. The terminal model uses all rows; the recurrent model uses rows with
#'    `k_idx > 1`. Baseline (`k_idx == 1`) and time-varying covariate matrices
#'    unrelated to the lag terms are set to zero-column matrices (`P = 0`);
#' 4. Loads the pre-compiled Stan model from `stan_model_file` and runs MCMC via
#'    `rstan::sampling()`, returning the fitted object along with the data,
#'    design information, prior settings, and the full list of inputs to Stan.
=======
#' 1) Copies user-named columns into canonical names `pat_id`, `k_idx`, `A`, `T_obs`, `Y_obs`.
#' 2) Calls `preprocess_data()` to order rows by (`pat_id`, `k_idx`), remap subject IDs to
#'    consecutive integers, and, when lag terms are requested but missing, create the `lag_col`
#'    within subject. Rows after the first `T_obs == 1` per subject are removed (the death row
#'    itself is kept). The data are not padded to a full grid.
#' 3) Builds design matrices from the right-hand sides of `formula_T` and `formula_Y`, excluding
#'    the treatment column. Missing values in the model matrices are set to zero. The terminal
#'    model uses all rows; the recurrent model uses `k_idx == 1` for baseline and `k_idx > 1`
#'    for follow-up.
#' 4) Loads the pre-compiled Stan model from the package and runs MCMC via `rstan::sampling()`.
>>>>>>> 759e44d2adaa77e0e2258aa398e3f620bf7ff1ce
#'
#' @examples
#' df <- data.frame(
#'   sid         = rep(1:2, each = 2),
#'   period      = rep(1:2, 2),
<<<<<<< HEAD
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
#'   lag_col    = "lagYk"
#'   formula_T  = death_flag  ~ lagYk + A ,
#'   formula_Y  = event_count ~ lagYk + A ,
#'   prior      = prior,
#'   num_chains = 1, iter = 100,
#'   stan_model_file = "causal_recur_model.rds",
#' )
#' print(fit)
#' summary(fit)
#' plot(fit)
=======
#'   event_count = c(0,1, 1,0),
#'   death_flag  = c(0,0, 0,1),
#'   trt_arm     = c(0,1, 0,1)
#' )
#' prior <- list(
#'   eta_beta = 0, sigma_beta = 1, rho_beta = 0.5,
#'   eta_gamma = 0, sigma_gamma = 1, rho_gamma = 0.5
#' )
#' \dontrun{
#' fit <- fit_causal_recur(
#'   data      = df, K = 2,
#'   id_col    = "sid",
#'   time_col  = "period",
#'   treat_col = "trt_arm",
#'   formula_T = death_flag  ~ lagYk + A,
#'   formula_Y = event_count ~ lagYk + A,
#'   prior     = prior,
#'   num_chains = 1, iter = 200,
#'   cores = 1, verbose = FALSE,
#'   lag_col = "lagYk"
#' )
#' print(fit)
>>>>>>> 759e44d2adaa77e0e2258aa398e3f620bf7ff1ce
#' }
#'
#' @importFrom rstan stan_model sampling
#' @importFrom stats model.matrix
#' @importFrom dplyr mutate
#' @importFrom stats terms reformulate na.pass plogis as.formula median
#' @name causal_recur_fit
#' @docType class
#' @export

fit_causal_recur <- function(
    data, K,
    id_col, time_col, treat_col,
    formula_T, formula_Y,
    prior = NULL,
    num_chains = 4, iter = 2000,
    control = list(adapt_delta = 0.95, max_treedepth = 15),
    cores = 1,
    verbose = TRUE,
    lag_col = "lagYk"
) {
  stan_model_file <- system.file("stan", "causal_recur_model.rds", package = "BayCauRETM")
  if (!nzchar(stan_model_file) || !file.exists(stan_model_file))
    stop("Stan model file not found at: ", stan_model_file)

  if (is.null(prior)) {
    prior <- list(
      eta_beta = 0,  sigma_beta = 0.7,  rho_beta = 0.6,
      eta_gamma= 0,  sigma_gamma= 0.7,  rho_gamma= 0.6,
      sigma_beta1 = 0.5,
      sigma_theta1 = 0.5,
      sigma_theta_lag = 0.5
    )
  }

  need_cols <- c(id_col, time_col, treat_col)
  if (any(miss <- !need_cols %in% names(data)))
    stop("Columns not found: ", paste(need_cols[miss], collapse = ", "))

  event_col <- all.vars(formula_T)[1]   # Tk
  count_col <- all.vars(formula_Y)[1]   # Yk
  if (!(event_col %in% names(data))) stop("Event column '", event_col, "' not found")
  if (!(count_col %in% names(data))) stop("Count column '", count_col, "' not found")

  df <- as.data.frame(data)
  df$T_obs  <- df[[event_col]]
  df$Y_obs  <- df[[count_col]]
  df$pat_id <- df[[id_col]]
  df$k_idx  <- df[[time_col]]
  df$A      <- df[[treat_col]]

  terms_no_resp <- function(f) stats::delete.response(stats::terms(f))
  vars_in_rhs   <- function(f) unique(all.vars(terms_no_resp(f)))
  uses_lag      <- !is.null(lag_col) && lag_col %in% c(vars_in_rhs(formula_T), vars_in_rhs(formula_Y))
  need_lag      <- uses_lag && !(lag_col %in% names(df))

  prep  <- preprocess_data(df, K = K, lag_col = if (uses_lag) lag_col else NULL, need_lag = need_lag)
  df    <- prep$processed_df
  n_pat <- prep$n_pat

  k1_mask <- df$k_idx == 1
  df_Y1   <- df[k1_mask, ]
  df_Yk   <- df[!k1_mask, ]

  NTk <- nrow(df)
  NY1 <- nrow(df_Y1)
  NYk <- nrow(df_Yk)

  `%||%` <- function(x, y) if (is.null(x)) y else x

  get_labels <- function(f) attr(stats::terms(f), "term.labels") %||% character(0)
  rhs_T <- setdiff(get_labels(formula_T), treat_col)
  rhs_Y <- setdiff(get_labels(formula_Y), treat_col)

  pat <- if (!is.null(lag_col)) paste0("\\b", lag_col, "\\b") else "^$"
  T_lag_terms <- rhs_T[grepl(pat, rhs_T)]
  Y_lag_terms <- rhs_Y[grepl(pat, rhs_Y)]
  T_cov_terms <- setdiff(rhs_T, T_lag_terms)
  Y_cov_terms <- setdiff(rhs_Y, Y_lag_terms)

  build_mm <- function(terms, d, env) {
    if (!length(terms)) return(matrix(0, nrow(d), 0))
    f <- stats::reformulate(terms, intercept = FALSE)
    environment(f) <- env
    mm <- stats::model.matrix(f, d, na.action = stats::na.pass)
    mm[is.na(mm)] <- 0
    storage.mode(mm) <- "double"
    mm
  }
  env_T <- attr(stats::terms(formula_T), ".Environment")
  env_Y <- attr(stats::terms(formula_Y), ".Environment")


  MM_T_cov_all  <- build_mm(T_cov_terms, df,    env_T)
  MM_Y_cov_y1   <- build_mm(Y_cov_terms, df_Y1, env_Y)
  MM_Y_cov_yk   <- build_mm(Y_cov_terms, df_Yk, env_Y)

  cov_cols <- union(colnames(MM_T_cov_all) %||% character(0),
                    union(colnames(MM_Y_cov_y1) %||% character(0),
                          colnames(MM_Y_cov_yk) %||% character(0)))

  align_mm <- function(mm, cols) {
    if (!length(cols)) return(matrix(0, nrow(mm), 0))
    if (is.null(colnames(mm))) {
      out <- matrix(0, nrow(mm), length(cols)); colnames(out) <- cols; return(out)
    }
    miss <- setdiff(cols, colnames(mm))
    mm2  <- if (length(miss)) cbind(mm, matrix(0, nrow(mm), length(miss), dimnames = list(NULL, miss))) else mm
    mm2[, cols, drop = FALSE]
  }

  L_Tk  <- align_mm(MM_T_cov_all, cov_cols)
  L_Y1  <- if (NY1) align_mm(MM_Y_cov_y1, cov_cols) else matrix(0, 0, length(cov_cols), dimnames = list(NULL, cov_cols))
  L_Yk  <- if (NYk) align_mm(MM_Y_cov_yk, cov_cols) else matrix(0, 0, length(cov_cols), dimnames = list(NULL, cov_cols))
  P     <- length(cov_cols)

  Lag_Tk <- build_mm(T_lag_terms, df,    env_T)
  Lag_Yk <- build_mm(Y_lag_terms, df_Yk, env_Y)
  QlagT  <- ncol(Lag_Tk)
  QlagY  <- ncol(Lag_Yk)

  param_labels <- list(
    cov_terms = union(T_cov_terms, Y_cov_terms),
    cov_cols  = cov_cols,
    T_lag     = colnames(Lag_Tk) %||% character(0),
    Y_lag     = colnames(Lag_Yk) %||% character(0)
  )

  prior_def <- list(
    eta_beta = 0,  sigma_beta = 1,  rho_beta = 0,
    eta_gamma = 0, sigma_gamma = 1, rho_gamma = 0,
    sigma_beta1 = 1,
    sigma_theta1 = 1,
    sigma_theta_lag = 1
  )
  prior_use <- modifyList(prior_def, prior)

  stan_data <- list(
    NY1 = NY1, NYk = NYk, NTk = NTk,
    K = K, P = P,
    kvecT = df$k_idx,
    L_Tk  = L_Tk,
    A_Tk  = df$A,
    Tk    = df$T_obs,
    kvecY = if (NYk) df_Yk$k_idx else integer(0),
    L_Yk  = L_Yk,
    A_Yk  = if (NYk) df_Yk$A else numeric(0),
    Yk    = if (NYk) df_Yk$Y_obs else integer(0),
    L_Y1  = L_Y1,
    A_Y1  = if (NY1) df_Y1$A else numeric(0),
    Y1    = if (NY1) df_Y1$Y_obs else integer(0),
    QlagT  = QlagT,
    Lag_Tk = if (QlagT) Lag_Tk else matrix(0, NTk, 0),
    QlagY  = QlagY,
    Lag_Yk = if (QlagY) Lag_Yk else matrix(0, NYk, 0),
    eta_beta = prior_use$eta_beta,
    sigma_beta = prior_use$sigma_beta,
    rho_beta = prior_use$rho_beta,
    eta_gamma = prior_use$eta_gamma,
    sigma_gamma = prior_use$sigma_gamma,
    rho_gamma = prior_use$rho_gamma,
    sigma_beta1 = prior_use$sigma_beta1,
    sigma_theta1 = prior_use$sigma_theta1,
    sigma_theta_lag = prior_use$sigma_theta_lag
  )

  if (verbose) message("Loading pre-compiled Stan model from: ", stan_model_file)
  stan_mod <- readRDS(stan_model_file)
  if (!inherits(stan_mod, "stanmodel"))
    stop("The RDS at ", stan_model_file, " is not a 'stanmodel' object.")

  if (verbose) message(sprintf("Sampling (%d chains * %d iter, cores=%d)...", num_chains, iter, cores))
  stan_fit <- rstan::sampling(
    object = stan_mod,
    data   = stan_data,
    chains = num_chains, iter = iter,
    cores  = cores, seed = 1234,
    control = control
  )

  structure(
    list(
      stan_fit          = stan_fit,
      data_preprocessed = df,
      n_pat             = n_pat,
      K                 = K,
      param_labels      = param_labels,
      design_info       = list(
        id_col = id_col, time_col = time_col, treat_col = treat_col,
        formula_T = formula_T, formula_Y = formula_Y,
        lag_col = if (uses_lag) lag_col else NULL
      ),
      prior             = prior_use,
      stan_data_list    = stan_data
    ),
    class = c("causal_recur_fit", "list")
  )
}




# print / summary / plot methods

#' @describeIn causal_recur_fit Print a brief object summary.
#' @param x A `causal_recur_fit` object.
#' @export

print.causal_recur_fit <- function(x, ...) {
  cat("Causal recurrent event model fit\n")
  if (!is.null(x$n_pat)) cat("  Number of subjects:", x$n_pat, "\n")
  if (!is.null(x$K))     cat("  Number of intervals K:", x$K, "\n")
  cat("Use summary() for posterior parameter estimates.\n")
  cat("Use mcmc_diagnosis() for MCMC convergence diagnostics.\n")
  invisible(x)
}


#' @describeIn causal_recur_fit Summarize posterior parameter estimates.
#' @param object A `causal_recur_fit` object.
#' @param pars_to_report Character vector of parameter names (regex allowed).
#' @method summary causal_recur_fit
#' @export

summary.causal_recur_fit <- function(object,
                                     pars_to_report = c("beta1", "theta1", "thetaLag",
                                                        "betaL[1]", "thetaL[1]"),
                                     ...) {
  stan_fit <- object$stan_fit
  sum_obj  <- tryCatch(rstan::summary(stan_fit, pars = pars_to_report, ...),
                       error = function(e) {
                         warning("Error extracting summary for parameters: ",
                                 paste(pars_to_report, collapse = ", "))
                         NULL
                       })

  if (is.null(sum_obj) || is.null(sum_obj$summary)) {
    cat("No summary available for specified parameters.\n")
    return(invisible(NULL))
  }

  sum_stan <- sum_obj$summary
  df <- data.frame(
    Parameter = rownames(sum_stan),
    Mean      = round(sum_stan[, "mean"], 4),
    `2.5%`    = round(sum_stan[, "2.5%"], 4),
    `97.5%`   = round(sum_stan[, "97.5%"], 4),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  print(df, row.names = FALSE)
  invisible(df)
}

#' @describeIn causal_recur_fit Display MCMC diagnostic guidance (no plot produced).
#' @param x A `causal_recur_fit` object.
#' @export

plot.causal_recur_fit <- function(x, ...) {
  cat("To check MCMC convergence, please run:\n")
  cat("  mcmc_diagnosis(fit_out, pars_to_check = ..., save_plots = ..., positivity = ...)\n")
  invisible(NULL)
}
