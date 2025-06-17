#' Summarize Stan and g-Computation Results
#'
#' Create summary tables of posterior parameter estimates and causal contrasts.
#'
#' @param fit_out Output list from \code{fit_causal_recur()}.
#' @param gcomp_out Output list from \code{g_computation()}.
#' @param pars_to_report Character vector of Stan parameter names to report (default:
#'   \code{c("beta_Y","beta_A","gamma_Y","gamma_A")}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{param_summary}{A data frame with columns: \code{Parameter}, \code{Mean},
#'     \code{2.5\%}, and \code{97.5\%}.}
#'   \item{delta_summary}{A data frame with columns: \code{s}, \code{Mean},
#'     \code{2.5\%}, and \code{97.5\%} for each treatment start time.}
#' }
#'
#' @details
#' This function produces side-by-side summaries for both the posterior model parameters
#' and the estimated treatment contrasts from g-computation.
#'
#' @examples
#' \dontrun{
#' res <- summarize_results(fit_out, gcomp_out)
#' print(res$param_summary)
#' print(res$delta_summary)
#' }
#'
#' @importFrom rstan summary
#' @export




summarize_results <- function(fit_out,
                              gcomp_out,
                              pars_to_report = c("beta_Y", "beta_A", "gamma_Y", "gamma_A")) {
  cat("===== Stan Posterior Parameter Summary =====\n")
  stan_fit <- fit_out$stan_fit

  stan_mode <- tryCatch(methods::slot(stan_fit, "mode"), error = function(e) NULL)
  if (is.null(stan_mode) || length(stan_mode) == 0) {
    df_par <- data.frame(
      Parameter = character(0),
      Mean      = numeric(0),
      `2.5%`    = numeric(0),
      `97.5%`   = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {

    sum_obj <- tryCatch(
      rstan::summary(stan_fit, pars = pars_to_report),
      error = function(e) NULL
    )
    sum_stan <- if (!is.null(sum_obj) && !is.null(sum_obj$summary)) sum_obj$summary else NULL

    if (is.null(sum_stan)) {
      df_par <- data.frame(
        Parameter = character(0),
        Mean      = numeric(0),
        `2.5%`    = numeric(0),
        `97.5%`   = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      df_par <- data.frame(
        Parameter = rownames(sum_stan),
        Mean      = sum_stan[, "mean"],
        `2.5%`    = sum_stan[, "2.5%"],
        `97.5%`   = sum_stan[, "97.5%"],
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
  }
  print(df_par)

  cat("\n===== g-computation: delta(s, K+1) Summary =====\n")
  delta_list <- gcomp_out$delta
  df_delta <- do.call(rbind, lapply(names(delta_list), function(nm) {
    x <- delta_list[[nm]]
    data.frame(
      s       = as.integer(sub("s=", "", nm)),
      Mean    = x$mean,
      `2.5%`  = x$CI_lower,
      `97.5%` = x$CI_upper,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_delta) <- NULL
  print(df_delta)

  return(list(
    param_summary = df_par,
    delta_summary = df_delta
  ))
}









#' Plot Causal Effect \eqn{\Delta(s, K+1)} vs Treatment Start
#'
#' Draw a line plot with a 95\% credible ribbon for the estimated causal
#' contrast \eqn{\Delta(s, K+1)} as a function of \eqn{s}.
#'
#' @param gcomp_out Output list from \code{g_computation()}.
#' @param s_vec Integer vector of treatment start intervals to plot.
#'
#' @return A \code{ggplot} object representing the causal contrast and its uncertainty.
#'
#' @examples
#' \dontrun{
#' p <- plot_delta_vs_s(gcomp_out, s_vec = 1:10)
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_minimal
#' @importFrom dplyr arrange
#'
#' @export





plot_delta_vs_s <- function(gcomp_out, s_vec) {
  delta_list <- gcomp_out$delta
  df_plot <- do.call(rbind, lapply(names(delta_list), function(nm) {
    x <- delta_list[[nm]]
    data.frame(
      s       = as.integer(sub("s=", "", nm)),
      Mean    = x$mean,
      CI_low  = x$CI_lower,
      CI_high = x$CI_upper,
      row.names = NULL
    )
  }))
  df_plot <- dplyr::arrange(df_plot, s)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = s, y = Mean)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_low, ymax = CI_high), alpha = 0.3) +
    ggplot2::labs(
      x = "Treatment Start Interval (s)",
      y = expression(Delta(s, K+1)),
      title = expression("Estimated Causal Effect " ~ Delta(s, K+1)),
      subtitle = "Posterior mean (line) and 95% credible interval (shaded)"
    ) +
    ggplot2::theme_minimal()

  return(p)
}

