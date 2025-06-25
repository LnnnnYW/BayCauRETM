utils::globalVariables(c(
  "Y_prev", "A", "s", "Mean", "CI_low", "CI_high",
  "pat_orig", "k_idx", "pat_id", "T_obs", "C_obs",
  "Y_obs", "T_prev", "C_prev", "n", "hazard", "surv_prob",
  "switch_prob", ".data", "scenario", "lower", "upper"
))

#' Preprocess Long-Format Data for Causal Recurrent-Event Modeling
#'
#' Prepares a long-format dataset for discrete-time Bayesian causal modeling
#' of recurrent and terminal events.
#'
#' @param df A data.frame in long format, one row per subject and interval.
#' @param id_col Character. Name of the column identifying subjects.
#' @param k_col Character. Name of the interval index column (should be integers from 1 to K).
#' @param y_col Character. Name of the column containing recurrent-event counts.
#' @param t_col Character. Name of the terminal-event (death) indicator column (0/1).
#' @param c_col Character. Name of the censoring indicator column (0/1).
#' @param a_col Character. Name of the treatment indicator column (0/1 or factor).
#' @param x_cols Character vector of additional covariate column names, or NULL.
#' @param K Integer. Total number of discrete intervals per subject.
#'
#' @return A list with components:
#' \describe{
#'   \item{processed_df}{A data.frame with columns \code{pat_id}, \code{k_idx}, \code{Y_obs}, \code{T_obs}, \code{C_obs}, \code{A}, lagged values \code{Y_prev}, \code{T_prev}, \code{C_prev}, and any \code{x_cols}.}
#'   \item{n_pat}{Integer. Number of unique subjects.}
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Renames columns to standardized names.
#'   \item Converts subject IDs to integer \code{pat_id}.
#'   \item Computes lagged variables (\code{Y_prev}, \code{T_prev}, \code{C_prev}).
#'   \item Replaces missing values in \code{Y_obs} with zero.
#'   \item Checks covariate existence, value ranges, duplicates, missing intervals, and binary coding.
#' }
#'
#' @examples
#' df_long <- data.frame(
#'   patient_id = rep(1:5, each = 3),
#'   k_idx      = rep(1:3, times = 5),
#'   Y          = sample(0:2, 15, TRUE),
#'   T          = rbinom(15, 1, 0.1),
#'   C          = rbinom(15, 1, 0.05),
#'   A          = rbinom(15, 1, 0.5),
#'   age        = rnorm(15, 60, 10),
#'   sex        = sample(c("M","F"), 15, TRUE)
#' )
#' \dontrun{
#' out <- preprocess_data(
#'   df     = df_long,
#'   id_col = "patient_id",
#'   k_col  = "k_idx",
#'   y_col  = "Y",
#'   t_col  = "T",
#'   c_col  = "C",
#'   a_col  = "A",
#'   x_cols = c("age", "sex"),
#'   K      = 3
#' )
#' str(out)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom dplyr rename arrange mutate group_by ungroup count filter anti_join
#' @importFrom tidyr complete
#' @importFrom forcats fct_lump
#' @export


preprocess_data <- function(df,
                            id_col,
                            k_col,
                            y_col,
                            t_col,
                            c_col,
                            a_col,
                            x_cols = NULL,
                            K) {
  # Basic parameter checks
  stopifnot(is.data.frame(df))
  stopifnot(is.character(id_col), is.character(k_col), is.character(y_col),
            is.character(t_col), is.character(c_col), is.character(a_col))
  stopifnot(is.numeric(K), length(K) == 1, K >= 1)
  if (!is.null(x_cols)) stopifnot(is.character(x_cols))

  # Rename and standardize columns
  df2 <- df %>%
    dplyr::rename(
      pat_orig = !!rlang::sym(id_col),
      k_idx     = !!rlang::sym(k_col),
      Y_obs     = !!rlang::sym(y_col),
      T_obs     = !!rlang::sym(t_col),
      C_obs     = !!rlang::sym(c_col),
      A         = !!rlang::sym(a_col)
    ) %>%
    dplyr::arrange(pat_orig, k_idx) %>%
    dplyr::mutate(
      pat_orig = as.factor(pat_orig),
      pat_id   = as.integer(pat_orig)
    ) %>%
    dplyr::ungroup()

  # Compute lagged values within each subject
  df2 <- df2 %>%
    dplyr::group_by(pat_id) %>%
    dplyr::arrange(k_idx) %>%
    dplyr::mutate(
      T_prev = dplyr::lag(T_obs, default = 0),
      C_prev = dplyr::lag(C_obs, default = 0),
      Y_prev = dplyr::lag(Y_obs, default = 0)
    ) %>%
    dplyr::ungroup()
  df2 <- df2 %>%
    dplyr::mutate(
      Y_prev = ifelse(is.na(Y_prev), 0, Y_prev),
      T_prev = ifelse(is.na(T_prev), 0, T_prev),
      C_prev = ifelse(is.na(C_prev), 0, C_prev)
    )
  # Replace NA in Y_obs with zero
  df2 <- df2 %>%
        group_by(pat_id) %>%
        tidyr::complete(
            k_idx = 1:K,
            fill = list(Y_obs = 0, T_obs = 0, C_obs = 0, A = 0)
          ) %>%
        ungroup()

  # Check existence of additional covariates
  if (!is.null(x_cols)) {
    missing_x <- setdiff(x_cols, names(df2))
    if (length(missing_x) > 0) {
      stop("Columns not found: ", paste(missing_x, collapse = ", "))
    }
  }

  # Validate k_idx range
  if (any(df2$k_idx < 1 | df2$k_idx > K)) {
    stop("k_idx values must lie between 1 and K.")
  }

  # Check duplicate (pat_id, k_idx)
  dup <- df2 %>% dplyr::count(pat_id, k_idx) %>% dplyr::filter(n > 1)
  if (nrow(dup) > 0) {
    stop("Duplicate rows for some (pat_id, k_idx).")
  }


  # Check coding of A, T_obs, C_obs
  if (!is.numeric(df2$A)) stop("A must be numeric 0/1.")
  if (!all(df2$A %in% c(0,1))) warning("A contains values outside 0/1.")
  if (any(!df2$T_obs %in% c(0,1))) stop("T_obs must be coded 0/1.")
  if (any(!df2$C_obs %in% c(0,1))) stop("C_obs must be coded 0/1.")

  # Check Y_obs integrity
  if (any(df2$Y_obs < 0 | df2$Y_obs != as.integer(df2$Y_obs))) {
    stop("Y_obs must be non-negative integers.")
  }

  # Check for variation
  if (length(unique(df2$A)) < 2) warning("No variation in A; effect estimation impossible.")
  if (length(unique(df2$T_obs)) < 2) warning("No variation in T_obs.")

  # Final NA check
  key_vars <- c("pat_id","k_idx","Y_obs","T_obs","C_obs","A", x_cols)
  if (anyNA(df2[, na.omit(key_vars)])) {
    stop("NA values found in key columns.")
  }

  # Number of subjects
  n_pat <- length(unique(df2$pat_id))

  result <- list(
    processed_df = df2,
    n_pat        = n_pat
  )
  class(result) <- c("preprocess_data", "list")
  return(result)
}
