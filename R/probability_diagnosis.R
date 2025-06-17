#' Diagnostics for Discrete-Time Switching Probability
#'
#' @description
#' Fits a discrete-time logistic-hazard model for treatment switching,
#' computes interval-specific switching probabilities, and visualises diagnostics.
#'
#' @param df A pre-processed data frame from \code{\link{preprocess_data}}
#'   containing at least \code{pat_id}, \code{k_idx}, \code{A} and any covariates.
#' @param covariates Character vector of covariate names to include, or \code{NULL}.
#' @param save_plots Logical; if \code{TRUE}, saves the plot as PNG.
#' @param plot_prefix File-name prefix when saving the plot.
#'
#' @return A list with components:
#' \describe{
#'   \item{summary}{Data frame of mean / median / 5-95 % quantiles of switching probability by interval.}
#'   \item{plot}{A \code{ggplot} box-plot object.}
#'   \item{model}{The fitted \code{glm} object.}
#' }
#'
#' @details
#' The function fits a logistic regression for the hazard of switching in each
#' discrete interval, predicts subject-level probabilities, summarises them,
#' and (optionally) writes a PNG of the diagnostic box-plot.
#'
#' @examples
#' \dontrun{
#' prep <- preprocess_data(
#'   data.frame(
#'     id = rep(1:3, each = 2), k = rep(1:2, 3),
#'     Y = sample(0:2, 6, TRUE),
#'     T = rbinom(6, 1, 0.1),
#'     C = rbinom(6, 1, 0.05),
#'     A = rbinom(6, 1, 0.5),
#'     X = rnorm(6)
#'   ),
#'   id_col = "id", k_col = "k", y_col = "Y",
#'   t_col = "T", c_col = "C", a_col = "A",
#'   x_cols = "X", K = 2
#' )
#' res <- switching_probability_diagnostics(prep$processed_df)
#' res$plot
#'}
#'
#' @importFrom stats glm as.formula quantile median
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal
#' @importFrom dplyr summarise ungroup
#' @export


switching_probability_diagnostics <- function(df,
                                              covariates = NULL,
                                              save_plots = FALSE,
                                              plot_prefix = "switch_") {
  # 1. Fit discrete-time hazard model for switching
  pred_vars <- c("Y_prev", "factor(k_idx)", covariates)
  formula_str <- paste("A ~", paste(pred_vars, collapse = " + "))
  mod_hazard  <- stats::glm(
    stats::as.formula(formula_str),
    data   = df,
    family = stats::binomial
  )

  # 2. Predict hazard (P[A=1 | ...])
  df$hazard <- stats::predict(mod_hazard, type = "response")

  # 3. Compute survival and switching probabilities per subject
  df <- df %>%
    arrange(pat_id, k_idx) %>%
    group_by(pat_id) %>%
    mutate(
      surv_prob   = cumprod(lag(1 - hazard, default = 1)),
      switch_prob = surv_prob * hazard
    ) %>%
    ungroup()

  # 4. Summarize distribution by k_idx
  summary_df <- df %>%
    group_by(k_idx) %>%
    summarise(
      mean   = mean(switch_prob),
      median = median(switch_prob),
      p05    = stats::quantile(switch_prob, 0.05),
      p25    = stats::quantile(switch_prob, 0.25),
      p75    = stats::quantile(switch_prob, 0.75),
      p95    = stats::quantile(switch_prob, 0.95)
    ) %>%
    ungroup()

  # 5. Create boxplot of switching probabilities
  p <- ggplot2::ggplot(df, ggplot2::aes(x = k_idx, y = switch_prob)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      title = "Switching Probability over Time",
      x     = "Weeks since Diagnosis, s",
      y     = "Probability of Initiating at s"
    )

  # 6. Save plot if requested
  if (save_plots) {
    ggplot2::ggsave(
      filename = paste0(plot_prefix, "switch_prob_boxplot.png"),
      plot     = p
    )
  }

  list(
    summary = summary_df,
    plot    = p,
    model   = mod_hazard
  )
}
