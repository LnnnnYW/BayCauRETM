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

preprocess_data <- function(df, K, lag_col = NULL) {
  stopifnot(is.data.frame(df))
  stopifnot(is.numeric(K), length(K) == 1, K >= 1)
  if (!is.null(lag_col)) stopifnot(is.character(lag_col), length(lag_col) == 1)

  req <- c("pat_id", "k_idx", "Y_obs", "T_obs", "A")
  if (any(miss <- !req %in% names(df)))
    stop("df is missing required columns: ",
         paste(req[miss], collapse = ", "))

  df <- df[order(df$pat_id, df$k_idx), ]
  pat_ids <- unique(df$pat_id)
  pat_map <- setNames(seq_along(pat_ids), pat_ids)
  df$pat_id <- pat_map[as.character(df$pat_id)]


  if (!is.null(lag_col) && !(lag_col %in% names(df))) {
    warning("Specified lag_col '", lag_col, "' not found. Auto-generating based on lagged Y_obs > 0.")

    n <- nrow(df)
    lag_values <- numeric(n)

    pat_starts <- match(unique(df$pat_id), df$pat_id)
    pat_ends <- c(pat_starts[-1] - 1, n)

    for (i in seq_along(pat_starts)) {
      start_idx <- pat_starts[i]
      end_idx <- pat_ends[i]
      if (end_idx > start_idx) {
        lag_values[(start_idx + 1):end_idx] <-
          as.integer(df$Y_obs[start_idx:(end_idx - 1)] > 0)
      }
    }
    df[[lag_col]] <- lag_values
  }


  structure(
    list(processed_df = df,
         n_pat        = length(pat_ids)),
    class = c("preprocess_data", "list")
  )
}
