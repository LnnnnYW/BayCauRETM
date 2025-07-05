# switching_probability_summary

#' Compute Switching Probability Summary
#'
#' @description
#' Fit a discrete-time logistic hazard model for treatment initiation ("switching"),
#' compute subject-level switching probabilities per interval, and summarize by interval.
#'
#' @param df A data.frame from \code{preprocess_data()}, containing at least:
#'   \code{pat_id}, \code{k_idx}, \code{A}, \code{Y_prev}, and any additional covariates.
#' @param covariates Character vector of covariate names to include in the hazard model,
#'   or NULL for only \code{Y_prev} and \code{k_idx}.
#' @param object Object returned by this function (used by summary method).
#' @inheritParams base::print
#'
#' @return An object of class \code{switching_summary}, which is a list with components:
#'   \describe{
#'     \item{\code{df2}}{Data.frame including original columns plus:
#'       \code{hazard} (predicted probability of switching at each interval),
#'       \code{surv_prob} (probability of not yet switched before each interval),
#'       \code{switch_prob} (joint probability of surviving to interval and switching at that interval).}
#'     \item{\code{summary_df}}{Data.frame summarizing \code{switch_prob} by \code{k_idx}:
#'       mean, median, 5th percentile, 25th, 75th, and 95th percentiles.}
#'     \item{\code{model}}{The fitted \code{glm} object (family = binomial).}
#'   }
#'
#' @details
#' This function fits a logistic regression for the hazard of initiating treatment
#' at each discrete interval, then computes for each subject the survival probability
#' (not yet switched) and the switching probability at each interval. It returns
#' the augmented data.frame and a summary by interval.
#'
#' @examples
#' \dontrun{
#' sw <- switching_probability_summary(processed_df, covariates = c("age", "sex"))
#' print(sw)       # prints summary_df
#' summary(sw)     # same as print
#' plot(sw)        # draws boxplot + mean curve
#' }
#'
#' @importFrom stats glm predict quantile
#' @importFrom dplyr arrange group_by mutate summarise ungroup lag
#' @export


switching_probability_summary <- function(df, covariates = NULL) {
  preds <- c("Y_prev", "factor(k_idx)", covariates)
  preds <- preds[!vapply(preds, is.null, logical(1))]
  formula_str <- paste("A ~", paste(preds, collapse = " + "))
  mod_hazard <- stats::glm(as.formula(formula_str), data = df, family = stats::binomial)

  df2 <- df %>%
    dplyr::arrange(pat_id, k_idx) %>%
    dplyr::mutate(
      hazard      = stats::predict(mod_hazard, type = "response"),
      surv_prob   = cumprod(dplyr::lag(1 - hazard, default = 1)),
      switch_prob = surv_prob * hazard
    )

  summary_df <- df2 %>%
    dplyr::group_by(k_idx) %>%
    dplyr::summarise(
      mean   = mean(switch_prob, na.rm = TRUE),
      median = median(switch_prob, na.rm = TRUE),
      p05    = stats::quantile(switch_prob, 0.05, na.rm = TRUE),
      p25    = stats::quantile(switch_prob, 0.25, na.rm = TRUE),
      p75    = stats::quantile(switch_prob, 0.75, na.rm = TRUE),
      p95    = stats::quantile(switch_prob, 0.95, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  out <- list(df2 = df2, summary_df = summary_df, model = mod_hazard)
  class(out) <- "switching_summary"
  out
}

# print

#' @describeIn switching_probability_summary Print the summary table of switching probabilities by interval.
#' @export
print.switching_summary <- function(x, ...) {
  cat("Switching probability summary by interval:\n")
  print(x$summary_df)
  invisible(x$summary_df)
}

# summary

#' @describeIn switching_probability_summary Alias for print.
#' @method summary switching_summary
#' @export
summary.switching_summary <- function(object, ...) {
  print(object)
}

# plot

#' @describeIn switching_probability_summary Plot switching probability over time.
#' @param show_mean Logical. Whether to overlay the mean curve (default: TRUE).
#' @param theme_fn  A ggplot2 theme function.
#' @param ...       Passed to \code{ggplot2::geom_boxplot()}.
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_line labs theme
#' @export
plot.switching_summary <- function(x, show_mean = TRUE, theme_fn = ggplot2::theme_minimal, ...) {
  df2 <- x$df2
  p <- ggplot2::ggplot(df2, ggplot2::aes(x = factor(k_idx), y = switch_prob)) +
    ggplot2::geom_boxplot(...) +
    ggplot2::labs(
      x = "Interval (k_idx)",
      y = "Switching probability",
      title = "Switching probability over intervals"
    ) +
    theme_fn()

  if (show_mean) {
    mean_df <- x$summary_df %>% dplyr::mutate(k_idx = as.numeric(k_idx))
    p <- p +
      ggplot2::geom_line(
        data = mean_df,
        mapping = ggplot2::aes(x = factor(k_idx), y = mean, group = 1),
        linewidth = 0.8, color = "blue"
      )
  }
  print(p)
  invisible(p)
}

