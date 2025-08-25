#' Fit Bayesian Causal Recurrent and Terminal-Event Model
#'
#' @description
#' Fits a discrete-time Bayesian model for recurrent counts and a terminal
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
#' @param prior Named list of gAR(1) hyperparameters. Supported elements:
#'   `eta_beta`, `sigma_beta`, `rho_beta`, `eta_gamma`, `sigma_gamma`,
#'   `rho_gamma`, `sigma_beta1`, `sigma_theta1`, `sigma_theta_lag`.
#'   Missing entries fall back to internal defaults.
#' @param num_chains Integer. Number of MCMC chains (default `4`).
#' @param iter Integer. Total iterations *per* chain including warm-up
#'   (default `2000`).
#' @param stan_model_file Path to a **pre-compiled** Stan model object saved
#'   via `saveRDS()` (default `"inst/stan/causal_recur_model.rds"`). The file
#'   must exist; otherwise the function stops.
#' @param control List passed to `rstan::sampling()` (e.g.,
#'   `list(adapt_delta = 0.95, max_treedepth = 15)`).
#' @param cores Integer. Number of CPU cores to use for sampling (passed to
#'   `rstan::sampling()`).
#' @param verbose Logical. Print progress messages (default `TRUE`).
#' @param lag_col Character scalar or `NULL`. If `NULL`, a column named
#'   `"lagYk"` is created and initialized to 0 before preprocessing. If a name
#'   is provided and that column is **absent** in `data`, `preprocess_data()`
#'   will generate it as an integer indicator based on the subject-specific lag
#'   of `I(Y_obs > 0)`. If the column already exists, it is left unchanged.
#' @inheritParams base::print
#'
#' @return An object of class `causal_recur_fit` (list) with elements
#'   `stan_fit`, `data_preprocessed`, `n_pat`, `K`, `design_info`, `prior`, and
#'   `stan_data_list`.
#'
#' @details
#' Internally the function:
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
    lag_col = NULL
) {
  stan_model_file <- system.file("stan", "causal_recur_model.rds", package = "BayCauRETM")
  if (!nzchar(stan_model_file) || !file.exists(stan_model_file)) {
    stop("Stan model file not found at: ", stan_model_file)
  }

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

  if (is.null(lag_col)) {
    lag_col <- "lagYk"
    df[[lag_col]] <- 0
  }

  prep  <- preprocess_data(df, K = K, lag_col = lag_col)
  df    <- prep$processed_df
  n_pat <- prep$n_pat

  k1_mask <- df$k_idx == 1
  df_Y1   <- df[k1_mask, ]
  df_Yk   <- df[!k1_mask, ]

  NTk <- nrow(df)
  NY1 <- nrow(df_Y1)
  NYk <- nrow(df_Yk)

  P      <- 0
  X_all  <- matrix(0, NTk, 0)
  X_base <- matrix(0, NY1, 0)
  L_Yk   <- matrix(0, NYk, 0)

  `%||%` <- function(x, y) if (is.null(x)) y else x
  get_rhs <- function(f) attr(terms(f), "term.labels") %||% character(0)
  tv_terms_T <- setdiff(get_rhs(formula_T), treat_col)
  tv_terms_Y <- setdiff(get_rhs(formula_Y), treat_col)

  build_mm <- function(rhs, d) {
    if (length(rhs) == 0) return(matrix(0, nrow(d), 0))
    mm <- model.matrix(reformulate(rhs, intercept = FALSE), d, na.action = na.pass)
    mm[is.na(mm)] <- 0
    storage.mode(mm) <- "double"
    mm
  }

  Lag_Tk <- build_mm(tv_terms_T, df)
  Lag_Yk <- build_mm(tv_terms_Y, df_Yk)

  QlagT <- ncol(Lag_Tk)
  QlagY <- ncol(Lag_Yk)

  lag_pat <- paste0("\\b", lag_col, "\\b")

  is_lag_T <- grepl(lag_pat, colnames(Lag_Tk) %||% character(0))
  is_lag_Y <- grepl(lag_pat, colnames(Lag_Yk) %||% character(0))

  names_T_cov <- (colnames(Lag_Tk) %||% character(0))[!is_lag_T]
  names_T_lag <- (colnames(Lag_Tk) %||% character(0))[ is_lag_T]

  names_Y_cov <- (colnames(Lag_Yk) %||% character(0))[!is_lag_Y]
  names_Y_lag <- (colnames(Lag_Yk) %||% character(0))[ is_lag_Y]

  param_labels <- list(
    T_cov = names_T_cov,
    T_lag = names_T_lag,
    Y_cov = names_Y_cov,
    Y_lag = names_Y_lag
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
    L_Tk  = X_all,
    A_Tk  = df$A,
    Tk    = df$T_obs,

    kvecY = df_Yk$k_idx,
    L_Yk  = L_Yk,
    A_Yk  = df_Yk$A,
    Yk    = df_Yk$Y_obs,

    L_Y1 = X_base,
    A_Y1 = if (NY1) df_Y1$A else integer(0),
    Y1   = if (NY1) df_Y1$Y_obs else integer(0),

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
  if (!inherits(stan_mod, "stanmodel")) {
    stop("The RDS at ", stan_model_file, " is not a 'stanmodel' object. Got class: ",
         paste(class(stan_mod), collapse = ", "))
  }

  if (verbose) {
    message(sprintf("Sampling (%d chains * %d iter, cores=%d)...",
                    num_chains, iter, cores))
  }

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
        formula_T   = formula_T,
        formula_Y   = formula_Y,
        covariate_names = character(0),
        treat_col   = treat_col
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
