utils::globalVariables(c("ps", "trim_flag"))

#' Propensity score diagnostics
#'
#' @description
#' Fit a logistic regression model for treatment assignment to estimate propensity
#' scores, provide numerical summaries, flag extreme values, and generate
#' diagnostic plots.
#'
#' @details
#' The function:
#' 1. Fits a logistic regression `glm(treat ~ covariates, family = binomial)`.
#' 2. Computes predicted propensity scores (`ps`) for all observations.
#' 3. Summarizes key statistics: minimum, lower trim quantile, 5%/95%, upper trim
#'    quantile, maximum.
#' 4. Flags observations with PS outside the trimming bounds.
#' 5. Builds diagnostic plots:
#'    - Histogram of `ps` with vertical lines at trimming and 5%/95% cutoffs.
#'    - Density curves of `ps` by treatment group with the same vertical lines.
#'
#' If `auto_lag = TRUE` and `"lagyk"` appears among `covariates` (case-insensitive)
#' but the corresponding column is missing while `pat_id`, `k_idx`, and `Y_obs`
#' exist, a `lagYK` column is auto-generated as the within-subject lag of `Y_obs`
#' ordered by `(pat_id, k_idx)` with `0` at each subject's first row.
#'
#' @param data A `data.frame` containing the treatment column and covariates.
#' @param treat_col Character. Name of the binary treatment column (`0/1` or factor).
#'   Must exist in `data`.
#' @param covariates Character vector. Names of covariate columns to include in
#'   the PS model. Must exist in `data` (see `auto_lag` note above).
#' @param trim Numeric vector of length 2. Lower and upper quantile bounds used
#'   to flag PS tails. For example, `c(0.01, 0.99)` flags observations below the
#'   1% or above the 99% quantile. Default `c(0.01, 0.99)`.
#' @param plot_types Character vector. Which plots to produce; subset of
#'   `c("histogram", "density")`. Default `c("histogram", "density")`.
#' @param bins Integer. Number of bins for the histogram (default `30`).
#' @param auto_lag Logical. If `TRUE`, attempt to auto-create `lagYK` as described
#'   in **Details** when requested by `covariates`. Default `TRUE`.
#'
#' @return An object of class `ps_diag`, a list with:
#'   * **model** - the fitted `glm` object (`family = binomial`).
#'   * **ps** - numeric vector of predicted propensity scores aligned to `data`.
#'   * **summary** - data frame with columns `stat` (e.g., `"min"`, `"p_lo"`,
#'     `"p05"`, `"p95"`, `"p_hi"`, `"max"`) and `value`.
#'   * **trim_flag** - logical vector indicating observations outside trimming bounds.
#'   * **plots** - named list of **ggplot2** objects for requested `plot_types`
#'     (e.g., `"histogram"`, `"density"`).
#'
#' @examples
#' \dontrun{
#' # Assume preprocess_data() has generated processed_df
#' df2 <- prep_out$processed_df
#' ps_diag <- propensity_score_diagnostics(
#'   data       = df2,
#'   treat_col  = "A",
#'   covariates = c("Y_prev", "X1", "X2"),
#'   trim       = c(0.01, 0.99),
#'   plot_types = c("histogram", "density"),
#'   bins       = 30
#' )
#' # Print numerical summary
#' print(ps_diag)
#' # Summary (same as print)
#' summary(ps_diag)
#' # Plot histogram
#' plot(ps_diag, type = "histogram")
#' # Plot density by treatment
#' plot(ps_diag, type = "density")
#' }
#'
#' @export
#' @importFrom stats glm predict quantile
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline labs theme element_text scale_color_manual scale_fill_manual
#' @importFrom dplyr mutate



propensity_score_diagnostics <- function(data,
                                         treat_col,
                                         covariates,
                                         trim       = c(0.01, 0.99),
                                         plot_types = c("histogram", "density"),
                                         bins       = 30,
                                         auto_lag   = TRUE) {
  # helpers
  .match_col <- function(name, pool) {
    idx <- which(tolower(pool) == tolower(name))
    if (length(idx)) pool[idx[1]] else NA_character_
  }
  .safe_vline <- function(p, xs, ...) {
    xs <- xs[is.finite(xs)]
    if (length(xs) == 0L) return(p)
    p + ggplot2::geom_vline(xintercept = xs, ...)
  }

  tcol <- .match_col(treat_col, names(data))
  if (is.na(tcol))
    stop("treat_col '", treat_col, "' not found. Available: ",
         paste(names(data), collapse = ", "))

  real_cov <- vapply(covariates, .match_col, character(1), pool = names(data))

  if (anyNA(real_cov) && auto_lag &&
      any(tolower(covariates) == "lagyk") &&
      all(c("pat_id","k_idx","Y_obs") %in% names(data))) {

    if (!"lagYK" %in% names(data)) {
      data <- data |>
        dplyr::arrange(pat_id, k_idx) |>
        dplyr::group_by(pat_id) |>
        dplyr::mutate(lagYK = dplyr::lag(Y_obs, default = 0)) |>
        dplyr::ungroup()
      message("lagYK column auto-generated from Y_obs")
    }
    real_cov <- vapply(covariates, .match_col, character(1), pool = names(data))
  }

  if (anyNA(real_cov)) {
    miss <- covariates[is.na(real_cov)]
    stop("Covariate(s) not found: ", paste(miss, collapse = ", "),
         "\nAvailable columns: ", paste(names(data), collapse = ", "))
  }

  formula_ps <- stats::as.formula(paste(tcol, "~", paste(real_cov, collapse = " + ")))
  mf         <- stats::model.frame(formula_ps, data = data, na.action = stats::na.omit)
  rows_fit   <- as.integer(rownames(mf))
  mod_ps     <- stats::glm(formula_ps, data = mf, family = stats::binomial)

  ps_vec <- rep(NA_real_, nrow(data))
  ps_vec[rows_fit] <- stats::predict(mod_ps, type = "response")

  ps_nonmiss <- ps_vec[is.finite(ps_vec)]
  stats_vals <- c(
    min  = min(ps_nonmiss),
    p_lo = stats::quantile(ps_nonmiss, trim[1]),
    p05  = stats::quantile(ps_nonmiss, 0.05),
    p95  = stats::quantile(ps_nonmiss, 0.95),
    p_hi = stats::quantile(ps_nonmiss, trim[2]),
    max  = max(ps_nonmiss)
  )
  ps_rng <- range(ps_nonmiss)
  eps    <- 1e-6
  x_lo   <- max(stats_vals["p_lo"] + eps, ps_rng[1])
  x_hi   <- min(stats_vals["p_hi"] - eps, ps_rng[2])

  summary_df <- data.frame(
    stat  = names(stats_vals),
    value = as.numeric(stats_vals),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  trim_flag <- !is.na(ps_vec) &
    (ps_vec < stats_vals["p_lo"] | ps_vec > stats_vals["p_hi"])

  df2 <- data
  df2$ps        <- ps_vec
  df2$trim_flag <- trim_flag

  plots <- list()
  if ("histogram" %in% plot_types) {
    p <- ggplot2::ggplot(df2, ggplot2::aes(x = ps)) +
      ggplot2::geom_histogram(bins = bins, fill = "grey80", color = "white") +
      ggplot2::labs(title = "Propensity Score Histogram",
                    x = "Propensity Score", y = "Count") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p <- .safe_vline(p, c(x_lo, x_hi), linetype = "dashed", color = "red")
    p <- .safe_vline(p, c(stats_vals["p05"], stats_vals["p95"]),
                     linetype = "dotted", color = "blue")
    plots$histogram <- p
  }

  if ("density" %in% plot_types) {
    treat_factor <- as.factor(data[[tcol]])
    p <- ggplot2::ggplot(df2, ggplot2::aes(x = ps, color = treat_factor, fill = treat_factor)) +
      ggplot2::geom_density(alpha = 0.2, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = c("0" = "grey60", "1" = "dodgerblue"),
                                  name   = tcol) +
      ggplot2::scale_fill_manual(values  = c("0" = "grey80", "1" = "lightblue"),
                                 name    = tcol) +
      ggplot2::labs(title = "Propensity Score Density by Treatment",
                    x = "Propensity Score", y = "Density") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p <- .safe_vline(p, c(stats_vals["p_lo"], stats_vals["p_hi"]),
                     linetype = "dashed", color = "red")
    p <- .safe_vline(p, c(stats_vals["p05"], stats_vals["p95"]),
                     linetype = "dotted", color = "blue")
    plots$density <- p
  }

  out <- list(
    model     = mod_ps,
    ps        = ps_vec,
    summary   = summary_df,
    trim_flag = trim_flag,
    plots     = plots
  )
  class(out) <- "ps_diag"
  out
}


#' Print method for ps_diag
#'
#' @description
#' Print a summary of propensity scores and the number of observations flagged
#' outside the trimming bounds.
#'
#' @param x An object of class `ps_diag`.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the `ps_diag` object.
#' @rdname ps_diag-print
#' @method print ps_diag
#' @export

print.ps_diag <- function(x, ...) {
  cat("Propensity Score Summary:\n")
  print(x$summary)
  cat("\nFlagged (outside trim bounds): ",
      sum(x$trim_flag), " of ", length(x$trim_flag), "\n", sep = "")
  invisible(x)
}

#' Summary method for ps_diag
#'
#' @description
#' Summary for a `ps_diag` object; same as `print.ps_diag()`.
#'
#' @param object An object of class `ps_diag`.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the `ps_diag` object.
#' @rdname ps_diag-summary
#' @method summary ps_diag
#' @export

summary.ps_diag <- function(object, ...) {
  print(object)
}

#' Plot method for ps_diag
#'
#' @description
#' Plot propensity score diagnostics. The `type` argument selects one of the
#' available plots in `x$plots` (typically `"histogram"` or `"density"`).
#'
#' @param x An object of class `ps_diag`.
#' @param type Character. One of `names(x$plots)`; typically `"histogram"` or
#'   `"density"`.
#' @param ... Additional arguments passed to the underlying ggplot layers.
#'
#' @return Invisibly returns the ggplot object.
#' @rdname ps_diag-plot
#' @method plot ps_diag
#' @export

plot.ps_diag <- function(x, type = c("histogram", "density"), ...) {
  type <- match.arg(type, choices = names(x$plots))
  p <- x$plots[[type]]
  print(p)
  invisible(p)
}
