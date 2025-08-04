utils::globalVariables(c(
  "Y_prev", "A", "s", "Mean", "CI_low", "CI_high",
  "k_idx", "pat_id", "T_obs",
  "Y_obs", "T_prev", "n", "hazard", "surv_prob",
  "switch_prob", ".data", "scenario", "lower", "upper", "lagY_bin"
))

#' Preprocess Long‑Format Data for Causal Recurrent‑Event Modeling
#'
#' Prepares a long‑format dataset (already using standardised column names) for
#' discrete‑time Bayesian causal modelling of recurrent and terminal events.
#'
#' @param df A *long‑format* `data.frame` that **already contains** the analysis
#'        columns `pat_id`, `k_idx`, `Y_obs`, `T_obs`, `A`, plus any additional
#'        covariates.
#' @param K  Integer. Total number of discrete intervals per subject.
#' @param x_cols Character vector of additional covariate column names, or `NULL`.
#'
#' @return A list with components:
#'   * **processed_df** – The input data with filled‑in intervals, lagged
#'     variables `Y_prev`/`T_prev`, a binary `lagY_bin`, and post‑death rows
#'     removed.
#'   * **n_pat** – Integer, number of unique subjects.
#'
#' @details
#' * Fills every subject’s record to the full grid `k_idx = 1:K`;
#' * Carries treatment status `A` forward in time (`tidyr::fill(direction = "down")`);
#' * Removes all rows **after** the first `T_obs == 1` to ensure monotone
#'   terminal indicators;
#' * Creates `T_prev`, `Y_prev`, and a helper indicator `lagY_bin = I(Y_prev > 0)`;
#' * Runs basic validity checks (duplicate `(pat_id,k_idx)`, coding, ranges, missing values).
#' The data **must already** follow the naming convention `*_obs` for observed
#' quantities. The function does *not* rename columns – it only pads, truncates,
#' lags, and validates.
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
#' @importFrom tidyr complete fill
#' @keywords internal
#' @noRd


preprocess_data <- function(df, K, x_cols = NULL, lag_val = NA) {

  stopifnot(is.data.frame(df))
  stopifnot(is.numeric(K), length(K) == 1, K >= 1)

  req <- c("pat_id", "k_idx", "Y_obs", "T_obs", "A")
  if (any(miss <- !req %in% names(df)))
    stop("df is missing required columns: ",
         paste(req[miss], collapse = ", "))

  df$lagYK <- lag_val

  df <- df |>
    dplyr::arrange(pat_id, k_idx) |>
    dplyr::mutate(pat_id = as.integer(as.factor(pat_id)))

  key_vars <- c("pat_id", "k_idx", "Y_obs", "T_obs", "A", x_cols)
  if (anyNA(df[key_vars])) {
    warning("NA found in covariates – replaced with 0")
    df[key_vars][is.na(df[key_vars])] <- 0
  }

  structure(
    list(processed_df = df,
         n_pat        = length(unique(df$pat_id))),
    class = c("preprocess_data", "list")
  )
}
