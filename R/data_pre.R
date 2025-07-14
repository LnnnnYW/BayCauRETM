utils::globalVariables(c(
  "Y_prev", "A", "s", "Mean", "CI_low", "CI_high",
  "k_idx", "pat_id", "T_obs",
  "Y_obs", "T_prev", "n", "hazard", "surv_prob",
  "switch_prob", ".data", "scenario", "lower", "upper"
))

#' Preprocess Long-Format Data for Causal Recurrent-Event Modeling
#'
#' Prepares a long-format dataset (already using standardized column names) for
#' discrete-time Bayesian causal modeling of recurrent and terminal events.
#'
#' @param df A *long-format* `data.frame` that **already contains** the analysis
#'        columns `pat_id`, `k_idx`, `Y_obs`, `T_obs`, `A`, plus any additional
#'        covariates.
#' @param K  Integer. Total number of discrete intervals per subject.
#' @param x_cols Character vector of additional covariate column names, or `NULL`.
#'
#' @return A list with components:
#'   * **processed_df** - The input data with added lagged variables `Y_prev`,
#'     `T_prev`, and filled-in intervals.
#'   * **n_pat** - Integer, number of unique subjects.
#'
#' @details The data **must already follow** the naming convention
#' `*_obs` for observed quantities and `*_prev` for their lagged counterparts.
#' The function therefore *does not* rename columns.  Its only tasks are to
#' fill in missing intervals, create lagged variables, run basic validity
#' checks, and return the results.
#'
#' @examples
#' df_long <- data.frame(
#'   pat_id = rep(1:5, each = 3),
#'   k_idx  = rep(1:3, times = 5),
#'   Y_obs  = sample(0:2, 15, TRUE),
#'   T_obs  = rbinom(15, 1, 0.1),
#'   A      = rbinom(15, 1, 0.5),
#'   age    = rnorm(15, 60, 10)
#' )
#' \dontrun{
#' out <- preprocess_data(df_long, K = 3, x_cols = "age")
#' str(out)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange mutate group_by ungroup count filter select all_of
#' @importFrom tidyr complete
#' @keywords internal
#' @noRd


preprocess_data <- function(df,
                            K,
                            x_cols = NULL) {
  #  basic checks
  stopifnot(is.data.frame(df))
  stopifnot(is.numeric(K), length(K) == 1, K >= 1)
  if (!is.null(x_cols)) stopifnot(is.character(x_cols))

  required_cols <- c("pat_id", "k_idx", "Y_obs", "T_obs", "A")
  missing_req   <- setdiff(required_cols, names(df))
  if (length(missing_req) > 0) {
    stop("df is missing required columns: ", paste(missing_req, collapse = ", "))
  }

  # prepare data
  df2 <- df %>%
    dplyr::arrange(pat_id, k_idx) %>%
    # ensure pat_id is numeric (required by some modelling back-ends)
    dplyr::mutate(pat_id = as.integer(as.factor(pat_id))) %>%
    dplyr::select(pat_id, k_idx, Y_obs, T_obs, A, dplyr::all_of(x_cols))

  # fill missing intervals & NAs
  df2 <- df2 %>%
    dplyr::group_by(pat_id) %>%
    tidyr::complete(
      k_idx = 1:K,
      fill = list(Y_obs = 0, T_obs = 0, A = 0)
    ) %>%
    dplyr::ungroup()

  # lagged variables
  df2 <- df2 %>%
    dplyr::group_by(pat_id) %>%
    dplyr::arrange(k_idx) %>%
    dplyr::mutate(
      T_prev = dplyr::lag(T_obs, default = 0),
      Y_prev = dplyr::lag(Y_obs, default = 0)
    ) %>%
    dplyr::ungroup()

  # validations
  # k_idx range
  if (any(df2$k_idx < 1 | df2$k_idx > K)) {
    stop("k_idx values must lie between 1 and K.")
  }

  # duplicates
  dup <- df2 %>% dplyr::count(pat_id, k_idx) %>% dplyr::filter(n > 1)
  if (nrow(dup) > 0) {
    stop("Duplicate rows for some (pat_id, k_idx).")
  }

  # coding checks
  if (!is.numeric(df2$A)) stop("A must be numeric 0/1.")
  if (!all(df2$A %in% c(0, 1))) warning("A contains values outside 0/1.")
  if (any(!df2$T_obs %in% c(0, 1))) stop("T_obs must be coded 0/1.")

  # outcome integrity
  if (any(df2$Y_obs < 0 | df2$Y_obs != as.integer(df2$Y_obs))) {
    stop("Y_obs must be non-negative integers.")
  }

  # variation warnings
  if (length(unique(df2$A))   < 2) warning("No variation in A; effect estimation impossible.")
  if (length(unique(df2$T_obs)) < 2) warning("No variation in T_obs.")

  # final NA check
  key_vars <- c("pat_id", "k_idx", "Y_obs", "T_obs", "A", x_cols)
  if (anyNA(df2[, na.omit(key_vars)])) {
    stop("NA values found in key columns.")
  }

  n_pat <- length(unique(df2$pat_id))

  structure(
    list(processed_df = df2, n_pat = n_pat),
    class = c("preprocess_data", "list")
  )
}
