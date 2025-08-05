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
                             control = list(adapt_delta = 0.95, max_treedepth = 15),
                             cores = 1,
                             verbose = TRUE,
                             lag_col = NULL) {

  needed <- c(id_col, time_col, treat_col)
  if (any(miss <- !needed %in% names(data)))
    stop("Columns not found: ", paste(needed[miss], collapse = ", "))

  event_col <- all.vars(formula_T)[1]
  count_col <- all.vars(formula_Y)[1]

  if (!(event_col %in% names(data))) stop("Event column '", event_col, "' not found in data")
  if (!(count_col %in% names(data))) stop("Count column '", count_col, "' not found in data")

  df <- data.frame(
    T_obs  = data[[event_col]],
    Y_obs  = data[[count_col]],
    pat_id = data[[id_col]],
    k_idx  = data[[time_col]],
    A      = data[[treat_col]],
    stringsAsFactors = FALSE
  )

  other_cols <- setdiff(names(data), c(event_col, count_col, id_col, time_col, treat_col))
  df[other_cols] <- data[other_cols]

  prep <- preprocess_data(df, K = K, x_cols = x_cols, lag_col = lag_col)
  df   <- prep$processed_df
  n_pat <- prep$n_pat

  first_obs_idx <- match(unique(df$pat_id), df$pat_id)
  baseline_df <- df[first_obs_idx, ]

  get_X <- function(d) {
    if (is.null(x_cols) || length(x_cols) == 0) return(matrix(0, nrow(d), 0))
    mm <- model.matrix(stats::reformulate(x_cols, intercept = FALSE), data = d)
    storage.mode(mm) <- "double"
    mm
  }
  X_base <- get_X(baseline_df)
  P <- ncol(X_base)

  if (P > 0) {
    X_all <- X_base[df$pat_id, , drop = FALSE]
  } else {
    X_all <- matrix(0, nrow(df), 0)
  }

  k1_mask <- df$k_idx == 1
  kgt1_mask <- df$k_idx > 1

  df_Y1 <- df[k1_mask, ]
  df_Yk <- df[kgt1_mask, ]
  NTk <- nrow(df); NY1 <- sum(k1_mask); NYk <- sum(kgt1_mask)

  get_rhs_terms <- function(f) {
    tl <- attr(stats::terms(f), "term.labels")
    if (is.null(tl)) character(0) else tl
  }
  terms_T <- get_rhs_terms(formula_T)
  terms_Y <- get_rhs_terms(formula_Y)
  drop_set <- unique(c(x_cols, treat_col))
  tv_terms_T <- setdiff(terms_T, drop_set)
  tv_terms_Y <- setdiff(terms_Y, drop_set)

  build_mm_safe <- function(rhs_terms, d) {
    if (length(rhs_terms) == 0) {
      return(matrix(0, nrow(d), 0, dimnames = list(NULL, character(0))))
    }

    tryCatch({
      formula_obj <- stats::reformulate(rhs_terms, intercept = FALSE)
      mm <- model.matrix(formula_obj, data = d, na.action = na.pass)

      if (anyNA(mm)) {
        mm[is.na(mm)] <- 0
      }

      storage.mode(mm) <- "double"
      mm
    }, error = function(e) {
      warning("Error building design matrix for terms: ", paste(rhs_terms, collapse = ", "),
              ". Using zero matrix. Error: ", e$message)
      matrix(0, nrow(d), 0)
    })
  }

  Lag_Tk <- build_mm_safe(tv_terms_T, df)
  Lag_Yk <- build_mm_safe(tv_terms_Y, df_Yk)

  if (nrow(Lag_Tk) != NTk) {
    stop(sprintf("Dimension mismatch: Lag_Tk has %d rows but expected %d",
                 nrow(Lag_Tk), NTk))
  }
  if (nrow(Lag_Yk) != NYk) {
    stop(sprintf("Dimension mismatch: Lag_Yk has %d rows but expected %d",
                 nrow(Lag_Yk), NYk))
  }

  QlagT <- ncol(Lag_Tk)
  QlagY <- ncol(Lag_Yk)


  L_Yk <- if (P > 0) X_base[df_Yk$pat_id, , drop = FALSE] else matrix(0, NYk, 0)

  stan_data <- list(
    NY1 = as.integer(NY1),
    NYk = as.integer(NYk),
    NTk = as.integer(NTk),
    K = as.integer(K),
    P = as.integer(P),

    kvecT = as.integer(df$k_idx),
    L_Tk  = X_all,
    A_Tk  = as.integer(df$A),
    Tk    = as.integer(df$T_obs),

    kvecY = as.integer(df_Yk$k_idx),
    L_Yk  = L_Yk,
    A_Yk  = as.integer(df_Yk$A),
    Yk    = as.integer(df_Yk$Y_obs),

    L_Y1  = if (P > 0) X_base else matrix(0, NY1, 0),
    A_Y1  = if (NY1 > 0) as.integer(df_Y1$A) else integer(0),
    Y1    = if (NY1 > 0) as.integer(df_Y1$Y_obs) else integer(0),

    QlagT = as.integer(QlagT),
    Lag_Tk = if (QlagT > 0) Lag_Tk else matrix(0, NTk, 0),
    QlagY = as.integer(QlagY),
    Lag_Yk = if (QlagY > 0) Lag_Yk else matrix(0, NYk, 0)
  )

  if (is.null(stan_model_file)) {
    pkg_dir <- system.file(package = "BayCauRETM")
    stan_model_file <- file.path(pkg_dir, "stan", "causal_recur_model.stan")
  }

  if (!is.null(compiled_model_cache) && file.exists(compiled_model_cache)) {
    if (verbose) message("Loading cached compiled Stan model...")
    stan_mod <- readRDS(compiled_model_cache)
  } else {
    if (verbose) message("Compiling Stan model...")
    stan_mod <- rstan::stan_model(stan_model_file, verbose = FALSE)
    if (!is.null(compiled_model_cache)) {
      dir.create(dirname(compiled_model_cache), recursive = TRUE, showWarnings = FALSE)
      saveRDS(stan_mod, compiled_model_cache)
    }
  }

  if (verbose) {
    message(sprintf("Sampling (%d chains × %d iter, cores=%d)...",
                    num_chains, iter, cores))
  }

  stan_fit <- rstan::sampling(
    stan_mod, data = stan_data,
    chains = num_chains, iter = iter,
    cores = cores,
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

#' @describeIn causal_recur_fit Prompt user to run diagnostics
#' @export
plot.causal_recur_fit <- function(x, ...) {
  cat("To check MCMC convergence, please run:\n")
  cat("  mcmc_diagnosis(fit_out, pars_to_check = ..., save_plots = ..., positivity = ...)\n")
  invisible(NULL)
}
