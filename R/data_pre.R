utils::globalVariables(c(
  "Y_prev", "A", "s", "Mean", "CI_low", "CI_high",
  "k_idx", "pat_id", "T_obs",
  "Y_obs", "T_prev", "n", "hazard", "surv_prob",
  "switch_prob", ".data", "scenario", "lower", "upper", "lagY_bin",
  "p25", "p75"
))

#' Preprocess Long-Format Data (ordering, ID remapping, optional lag creation)
#'
#' Prepares a long-format dataset that already follows the standard column
#' names for downstream discrete-time causal modeling. This lightweight
#' preprocessor currently performs ordering, subject-ID normalization, basic
#' input checks, and (optionally) the creation of a user-named lag indicator
#' derived from `Y_obs`.
#'
#' @param df A *long-format* `data.frame` that **already contains** the analysis
#'        columns `pat_id`, `k_idx`, `Y_obs`, `T_obs`, and `A`, plus any
#'        additional covariates.
#' @param K  Integer. Total number of discrete intervals per subject. Validated
#'        for correctness but **not** used to pad/expand the data in the current
#'        implementation.
#' @param lag_col Character scalar or `NULL`. If provided and the column is
#'        **absent** in `df`, the function auto-generates it as an integer
#'        indicator based on the subject-specific lag of `I(Y_obs > 0)`. If the
#'        named column already exists, it is left unchanged.
#'
#' @return A list with components:
#'   * **processed_df** - The input data sorted by `pat_id` and `k_idx`, with
#'     `pat_id` remapped to consecutive integers `1:n_pat` and, if requested and
#'     missing, a new `lag_col` created from `lag(I(Y_obs > 0))` within subject.
#'   * **n_pat** - Integer, number of unique subjects (after remapping).
#'
#' @details
#' * Validates inputs: `df` must be a `data.frame`; `K` must be a single numeric
#'   value `>= 1`; `lag_col` (if not `NULL`) must be a single character string.
#' * Checks that required columns `pat_id`, `k_idx`, `Y_obs`, `T_obs`, `A`
#'   exist; otherwise the function errors.
#' * Sorts rows by `pat_id`, then `k_idx`.
#' * Remaps (normalizes) `pat_id` to consecutive integers `1, 2, ..., n_pat`
#'   while preserving subject grouping.
#' * If `lag_col` is provided **and not found** in `df`, creates it as:
#'   `as.integer(lag(I(Y_obs > 0)))` computed **within each subject**; the first
#'   row of each subject receives `0`. Existing `lag_col` is respected.
#' * **Does not** pad to a full grid `k_idx = 1:K`, **does not** carry forward
#'   treatment, **does not** truncate rows after terminal events, and **does not**
#'   create `T_prev` or `Y_prev`. Column names are not changed.
#'
#' @examples
#' df_long <- data.frame(
#'   pat_id = rep(c("S1","S2","S3","S4","S5"), each = 3),
#'   k_idx  = rep(1:3, times = 5),
#'   Y_obs  = sample(0:2, 15, TRUE),
#'   T_obs  = rbinom(15, 1, 0.1),
#'   A      = rbinom(15, 1, 0.5),
#'   age    = rnorm(15, 60, 10)
#' )
#' \dontrun{
#' # Request creation of a missing lag column based on lagged I(Y_obs > 0)
#' out <- preprocess_data(df_long, K = 3, lag_col = "lagY_bin")
#' str(out)
#' head(out$processed_df)
#' table(out$processed_df$lagY_bin)
#' }
#'
#' @importFrom stats setNames
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
