# propensity_score_diagnostics

utils::globalVariables(c("ps", "trim_flag"))

#' Propensity Score Diagnostics
#'
#' @description
#' Fit a logistic regression model for treatment assignment to estimate propensity scores,
#' provide numerical summaries, flag extreme values, and generate diagnostic plots.
#'
#' @details
#' This function performs:
#' \enumerate{
#'   \item Fit a logistic regression via \code{glm(treat ~ covariates, family=binomial)}.
#'   \item Compute predicted propensity scores (\code{ps}) for all observations.
#'   \item Summarize key statistics: minimum, trimmed quantile, 5%/95%, trimmed upper, maximum.
#'   \item Flag observations with PS outside trimming bounds (e.g., below lower or above upper quantile).
#'   \item Generate diagnostic plots:
#'     \itemize{
#'       \item Histogram of \code{ps} with vertical lines at trimming and 5%/95% cutoffs.
#'       \item Density curves of \code{ps} by treatment group with the same vertical lines.
#'     }
#' }
#'
#' @param data A \code{data.frame} containing the treatment column and covariates.
#' @param treat_col Character. Name of the binary treatment column (0/1 or factor). Must exist in \code{data}.
#' @param covariates Character vector. Names of covariate columns to include in the PS model. Must exist in \code{data}.
#' @param trim Numeric vector of length 2. Lower and upper quantile bounds for flagging PS tails.
#'   For example, \code{c(0.01, 0.99)} flags observations with PS below 1% or above 99%. Default is \code{c(0.01, 0.99)}.
#' @param plot_types Character vector. Which plots to produce; subset of \code{c("histogram", "density")}.
#'   Default is \code{c("histogram", "density")}.
#' @param bins Integer. Number of bins for the histogram (default: 30).
#'
#' @return An object of class \code{ps_diag}, a list with components:
#'   \describe{
#'     \item{\code{model}}{The fitted \code{glm} object (family = binomial).}
#'     \item{\code{ps}}{Numeric vector of predicted propensity scores.}
#'     \item{\code{summary}}{Data frame with \code{stat} (e.g., "min","p_lo","p05","p95","p_hi","max") and \code{value}.}
#'     \item{\code{trim_flag}}{Logical vector indicating observations with PS outside trimming bounds.}
#'     \item{\code{plots}}{Named list of \code{ggplot2} objects per requested \code{plot_types}.}
#'   }
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
                                         bins       = 30) {
  # Input checks
  if (!treat_col %in% names(data)) {
    stop("treat_col '", treat_col, "' not found in data")
  }
  if (!all(covariates %in% names(data))) {
    missing <- setdiff(covariates, names(data))
    stop("Covariate(s) not found: ", paste(missing, collapse = ", "))
  }


  formula_ps <- as.formula(paste(treat_col, "~", paste(covariates, collapse = " + ")))
  mf         <- stats::model.frame(formula_ps, data = data, na.action = stats::na.omit)
  rows_fit   <- as.integer(rownames(mf))

  mod_ps <- stats::glm(formula_ps, data = mf, family = stats::binomial)

  ps_vec <- rep(NA_real_, nrow(data))
  ps_vec[rows_fit] <- stats::predict(mod_ps, type = "response")

  # Summary statistics
  ps_nonmiss <- ps_vec[!is.na(ps_vec)]
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

  # Flag observations outside trimming bounds
  trim_flag <- !is.na(ps_vec) &
    (ps_vec < stats_vals["p_lo"] | ps_vec > stats_vals["p_hi"])

  # Prepare data for plotting
  df2 <- data
  df2$ps        <- ps_vec
  df2$trim_flag <- trim_flag

  plots <- list()

  if ("histogram" %in% plot_types) {
    p_hist <- ggplot2::ggplot(df2, ggplot2::aes(x = ps)) +
      ggplot2::geom_histogram(bins = bins, fill = "grey80", color = "white") +
      ggplot2::geom_vline(xintercept = c(x_lo, x_hi),
                          linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = c(stats_vals["p05"], stats_vals["p95"]),
                          linetype = "dotted", color = "blue") +
      ggplot2::labs(title = "Propensity Score Histogram",
                    x = "Propensity Score", y = "Count") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plots$histogram <- p_hist
  }

  if ("density" %in% plot_types) {
    treat_factor <- as.factor(data[[treat_col]])
    p_den <- ggplot2::ggplot(df2, ggplot2::aes(x = ps, color = treat_factor, fill = treat_factor)) +
      ggplot2::geom_density(alpha = 0.2, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = c("0" = "grey60", "1" = "dodgerblue"),
                                  name   = treat_col) +
      ggplot2::scale_fill_manual(values  = c("0" = "grey80", "1" = "lightblue"),
                                 name    = treat_col) +
      ggplot2::geom_vline(xintercept = c(stats_vals["p_lo"], stats_vals["p_hi"]),
                          linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = c(stats_vals["p05"], stats_vals["p95"]),
                          linetype = "dotted", color = "blue") +
      ggplot2::labs(title = "Propensity Score Density by Treatment",
                    x = "Propensity Score", y = "Density") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plots$density <- p_den
  }

  # Output -------------------------------------------------------------------
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
#' Print summary of propensity scores and number of flagged observations.
#'
#' @param x An object of class \code{ps_diag}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the \code{ps_diag} object.
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
#' Summary for a \code{ps_diag} object; same as \code{print.ps_diag()}.
#'
#' @param object An object of class \code{ps_diag}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the \code{ps_diag} object.
#' @rdname ps_diag-summary
#' @method summary ps_diag
#' @export


summary.ps_diag <- function(object, ...) {
  print(object)
}

#' Plot method for ps_diag
#'
#' @description
#' Plot propensity score diagnostics. User specifies \code{type} as "histogram" or "density".
#'
#' @param x An object of class \code{ps_diag}.
#' @param type Character, either "histogram" or "density" (default: both available).
#' @param ... Additional arguments passed to the underlying ggplot layers (e.g., adjust bins).
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
