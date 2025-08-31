utils::globalVariables(c(
  "Y_prev", "A", "s", "Mean", "CI_low", "CI_high",
  "k_idx", "pat_id", "T_obs",
  "Y_obs", "T_prev", "n", "hazard", "surv_prob",
  "switch_prob", ".data", "scenario", "lower", "upper", "lagY_bin",
  "p25", "p75"
))

`%||%` <- function(x, y) if (is.null(x)) y else x

.vars_in_formula <- function(f) {
  if (is.null(f)) return(character(0))
  unique(all.vars(stats::delete.response(stats::terms(f))))
}

.formula_uses_var <- function(f, varname) {
  varname %in% .vars_in_formula(f)
}


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


preprocess_data <- function(df, K, lag_col = NULL, need_lag = FALSE) {
  stopifnot(is.data.frame(df))
  stopifnot(is.numeric(K), length(K) == 1, K >= 1)

  req <- c("pat_id", "k_idx", "Y_obs", "T_obs", "A")
  if (any(miss <- !req %in% names(df)))
    stop("df is missing required columns: ", paste(req[miss], collapse = ", "))

  df <- df[order(df$pat_id, df$k_idx), ]
  pat_ids <- unique(df$pat_id)
  pat_map <- setNames(seq_along(pat_ids), pat_ids)
  df$pat_id <- pat_map[as.character(df$pat_id)]

  if (need_lag && !is.null(lag_col) && !(lag_col %in% names(df))) {
    message("Lag column '", lag_col, "' not found; auto-generating from I(Y_{k-1} > 0).")
    n <- nrow(df)
    lag_values <- numeric(n)

    pat_starts <- match(unique(df$pat_id), df$pat_id)
    pat_ends   <- c(pat_starts[-1] - 1, n)

    for (i in seq_along(pat_starts)) {
      s <- pat_starts[i]; e <- pat_ends[i]
      if (e > s) lag_values[(s + 1):e] <- as.integer(df$Y_obs[s:(e - 1)] > 0)
      lag_values[s] <- 0   # every patient has lag=0 at first term
    }
    df[[lag_col]] <- lag_values
  }

  if (!is.null(lag_col) && lag_col %in% names(df)) {
    idx_first <- !duplicated(df$pat_id)
    df[[lag_col]][idx_first & is.na(df[[lag_col]])] <- 0
  }

  trim_after_death <- function(kk, tt) {
    pos <- which(tt == 1)
    if (length(pos)) {
      last <- min(pos)
      keep <- seq_len(length(kk)) <= last
    } else keep <- rep(TRUE, length(kk))
    keep
  }
  keep_lst <- mapply(
    trim_after_death,
    split(df$k_idx, df$pat_id),
    split(df$T_obs, df$pat_id),
    SIMPLIFY = FALSE
  )
  df <- df[unsplit(keep_lst, df$pat_id), , drop = FALSE]

  structure(
    list(processed_df = df,
         n_pat        = length(unique(df$pat_id))),
    class = c("preprocess_data", "list")
  )
}

