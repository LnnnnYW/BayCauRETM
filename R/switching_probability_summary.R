#' Compute switching-probability summary
#'
#' @description
#' Fit a discrete-time logistic hazard model for treatment initiation
#' ("switching"), compute subject-level switching probabilities per interval,
#' and summarize them by interval.
#'
#' @param df A `data.frame` from `preprocess_data()` (or similarly structured),
#'   containing at least `pat_id`, `k_idx`, `A`, `Y_obs` (if `lagYK` is missing),
#'   and any additional covariates. If `k_idx` is missing, it is created as the
#'   within-subject row index. If `lagYK` is missing, it is auto-generated as
#'   `I(lag(Y_obs) > 0)` within subject.
#' @param covariates Character vector of covariate names to include in the
#'   hazard model, or `NULL` to use only `lagYK` and a factorized time index.
#' @inheritParams base::print
#'
#' @return An object of class `switching_summary`, a list with:
#'   * **df2** - the input data with added columns:
#'     `hazard` (predicted switching probability at each interval),
#'     `surv_prob` (probability of not yet switched before the interval),
#'     and `switch_prob` (joint probability of surviving to the interval and
#'     switching at that interval).
#'   * **summary_df** - data frame summarizing `switch_prob` by `k_idx`:
#'     mean, median, 5th/25th/75th/95th percentiles.
#'   * **model** - the fitted `glm` object (`family = binomial`).
#'
#' @details
#' The discrete-time hazard is modeled via `glm(A ~ lagYK + factor(k_idx) + covariates,
#' family = binomial)`. After prediction, for each subject the survival
#' probability is computed as the cumulative product of `(1 - hazard)` lagged by
#' one interval; the switching probability at each interval is
#' `surv_prob * hazard`. If a provided `lagYK` is not binary, it is coerced to
#' `as.numeric(lagYK > 0)`.
#'
#' @examples
#' \dontrun{
#' sw <- switching_probability_summary(processed_df, covariates = c("age", "sex"))
#' print(sw)       # prints summary_df
#' summary(sw)     # same as print
#' plot(sw)        # ribbon or boxplot summaries
#' }
#'
#' @importFrom stats glm predict quantile
#' @importFrom dplyr arrange group_by mutate summarise ungroup lag
#' @importFrom utils modifyList
#' @export

switching_probability_summary <- function(df, covariates = NULL) {
  # helpers
  .match_any <- function(cands, pool) {
    pool_l <- tolower(pool)
    for (nm in cands) {
      hit <- which(pool_l == tolower(nm))
      if (length(hit)) return(pool[hit[1]])
    }
    NA_character_
  }
  .bt <- function(nm) { # backtick if needed
    ok <- grepl("^[A-Za-z.][A-Za-z0-9_.]*$", nm)
    ifelse(ok, nm, sprintf("`%s`", nm))
  }

  id_col   <- .match_any(c("pat_id","id"), names(df))
  time_col <- .match_any(c("k_idx","k","k_idx_int","time_idx","time"), names(df))
  tr_col   <- .match_any(c("A","Ak","treat","treatment"), names(df))
  lag_col  <- .match_any(c("lagYK","lagYk","lagyk"), names(df))
  y_col    <- .match_any(c("Y_obs","Yk"), names(df))

  if (is.na(id_col) || is.na(tr_col)) {
    stop("Data must contain an ID column (pat_id/id) and a treatment column (A/Ak).")
  }

  if (is.na(time_col)) {
    message("Time index column not found; creating k_idx from row order within each patient.")
    df <- df |>
      dplyr::group_by(!!rlang::sym(id_col)) |>
      dplyr::mutate(k_idx = dplyr::row_number()) |>
      dplyr::ungroup()
    time_col <- "k_idx"
  }

  if (is.na(lag_col)) {
    if (is.na(y_col))
      stop("lagYK not found and no Y column (Y_obs/Yk) to derive it from.")
    df <- df |>
      dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col)) |>
      dplyr::group_by(!!rlang::sym(id_col)) |>
      dplyr::mutate(lagYK = as.numeric(dplyr::lag(.data[[y_col]], default = 0) > 0)) |>
      dplyr::ungroup()
    lag_col <- "lagYK"
    message("lagYK auto-generated from Y (indicator of previous > 0).")
  } else {
    if (any(df[[lag_col]] < 0, na.rm = TRUE) || any(df[[lag_col]] > 1, na.rm = TRUE)) {
      df[[lag_col]] <- as.numeric(df[[lag_col]] > 0)
    }
  }

  if (!is.null(covariates)) {
    cov_map <- vapply(covariates, function(v) .match_any(v, names(df)), character(1))
    if (anyNA(cov_map)) {
      warning("Dropped missing covariate(s): ",
              paste(covariates[is.na(cov_map)], collapse = ", "))
      cov_map <- cov_map[!is.na(cov_map)]
    }
  } else cov_map <- character(0)

  fac_col <- "__factor_time_col__"
  while (fac_col %in% names(df)) fac_col <- paste0(fac_col, "_")
  df[[fac_col]] <- factor(df[[time_col]])

  lhs  <- .bt(tr_col)
  rhsv <- c(.bt(lag_col), .bt(fac_col), .bt(cov_map))
  form <- stats::as.formula(paste(lhs, "~", paste(rhsv, collapse = " + ")))

  mf <- stats::model.frame(form, data = df, na.action = stats::na.omit)
  if (!nrow(mf)) stop("No rows available to fit the switching model (all NA after subsetting).")
  rows_fit <- as.integer(rownames(mf))

  mod_hazard <- stats::glm(form, data = mf,
                           family = stats::binomial,
                           control = stats::glm.control(maxit = 100))

  hazard <- rep(NA_real_, nrow(df))
  hazard[rows_fit] <- stats::predict(mod_hazard, type = "response")

  df2 <- df |>
    dplyr::mutate(hazard = hazard) |>
    dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col)) |>
    dplyr::group_by(!!rlang::sym(id_col)) |>
    dplyr::mutate(
      surv_prob   = cumprod(dplyr::lag(1 - hazard, default = 1)),
      switch_prob = surv_prob * hazard
    ) |>
    dplyr::ungroup()

  summary_df <- df2 |>
    dplyr::group_by(!!rlang::sym(time_col)) |>
    dplyr::summarise(
      mean   = mean(switch_prob, na.rm = TRUE),
      median = median(switch_prob, na.rm = TRUE),
      p05    = stats::quantile(switch_prob, 0.05, na.rm = TRUE),
      p25    = stats::quantile(switch_prob, 0.25, na.rm = TRUE),
      p75    = stats::quantile(switch_prob, 0.75, na.rm = TRUE),
      p95    = stats::quantile(switch_prob, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(k_idx = !!rlang::sym(time_col))

  out <- list(df2 = df2, summary_df = summary_df, model = mod_hazard)
  class(out) <- "switching_summary"
  out
}

# print / summary / plot methods

#' @describeIn switching_probability_summary Print the summary table of switching probabilities by interval.
#'
#' @export

print.switching_summary <- function(x, ...) {
  cat("Switching probability summary by interval:\n")
  print(x$summary_df)
  invisible(x$summary_df)
}

#' @describeIn switching_probability_summary Alias for `print()`.
#' @param object A `switching_summary` object.
#' @method summary switching_summary
#' @export

summary.switching_summary <- function(object, ...) {
  print(object)
}

#' @describeIn switching_probability_summary Plot switching probability over time.
#'
#' @param type One of `"ribbon"` or `"boxplot"`. `"ribbon"` shows mean with an
#'   interquartile ribbon (p25-p75); `"boxplot"` shows per-interval distributions
#'   with an optional mean overlay.
#' @param show_mean Logical. Whether to overlay the mean curve (default `TRUE`).
#' @param theme_fn A **ggplot2** theme function (default `ggplot2::theme_minimal`).
#' @param ... Passed to `ggplot2::geom_boxplot()` when `type = "boxplot"`.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_line geom_point geom_ribbon labs theme scale_x_continuous
#' @export

plot.switching_summary <- function(x,
                                   type = c("ribbon","boxplot"),
                                   show_mean = TRUE,
                                   theme_fn = ggplot2::theme_minimal,
                                   ...) {
  type <- match.arg(type)
  df2 <- x$df2
  df2 <- df2[is.finite(df2$switch_prob) & is.finite(df2$k_idx), , drop = FALSE]

  if (type == "boxplot") {
    p <- ggplot2::ggplot(df2, ggplot2::aes(x = factor(k_idx), y = switch_prob)) +
      ggplot2::geom_boxplot(na.rm = TRUE, ...) +
      ggplot2::labs(
        x = "Interval (k_idx)",
        y = "Switching probability",
        title = "Switching probability over intervals"
      ) +
      theme_fn()
    if (show_mean) {
      mean_df <- x$summary_df
      p <- p +
        ggplot2::geom_line(
          data = mean_df,
          mapping = ggplot2::aes(x = as.numeric(k_idx), y = mean, group = 1),
          linewidth = 0.8
        )
    }
  } else {
    mean_df <- x$summary_df
    p <- ggplot2::ggplot(mean_df, ggplot2::aes(x = as.numeric(k_idx), y = mean)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = p25, ymax = p75), alpha = 0.2) +
      ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(size = 1, na.rm = TRUE) +
      ggplot2::scale_x_continuous(breaks = sort(unique(as.numeric(mean_df$k_idx))),
                                  labels = sort(unique(mean_df$k_idx))) +
      ggplot2::labs(title = "Switching probability by interval",
                    x = "Interval (k_idx)", y = "Switching probability") +
      theme_fn()
  }
  print(p)
  invisible(p)
}
