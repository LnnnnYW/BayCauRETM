# mcmc_diagnosis

#' MCMC Convergence and Positivity Diagnostics
#'
#' @description
#' Summarize MCMC convergence for a fitted model (R‑hat and effective sample
#' size), generate trace‑plots for specified parameters, and optionally perform
#' positivity (overlap) diagnostics on treatment assignment.
#'
#' @param fit_out Output list from [fit_causal_recur()], containing at least
#'   `stan_fit` (an **rstan** `stanfit` object) and — if `positivity = TRUE` —
#'   `data_preprocessed`.
#' @param pars_to_check Character vector of parameter‑name patterns to diagnose.
#'   **Defaults have been updated for the teacher‑style model:**
#'   `c("beta0", "theta0", "beta1", "theta1", "theta_lag")`.
#'   Patterns are passed to **bayesplot**’s `mcmc_trace()`.
#' @param save_plots Logical; if `TRUE`, saves each trace‑plot as a PNG file
#'   using `plot_prefix`.
#' @param plot_prefix Filename prefix for saved plots (default: `"traceplot_"`).
#' @param positivity Logical; if `TRUE`, perform additional positivity
#'   diagnostics on treatment assignment (default `FALSE`).
#' @inheritParams base::print
#'
#' @return An object of class `mcmc_diag`, a list with components:
#'   * **stats** —data frame with `Parameter`, `n_eff`, and `Rhat`.
#'   * **plots** —named list of ggplot objects for the trace‑plots.
#'
#' @details
#' If `positivity = TRUE`, a logistic regression of treatment `A` on history
#' (`Y_prev`, `k_idx`) plus any extra covariates is fitted using the
#' pre‑processed data. Propensity‑score summaries and a histogram are produced.
#'
#' @examples
#' \dontrun{
#' fit_out <- fit_causal_recur(...)
#' diag <- mcmc_diagnosis(fit_out, positivity = TRUE)
#' print(diag)
#' plot(diag, pars = "beta1")
#' }
#'
#' @importFrom rstan summary
#' @importFrom bayesplot mcmc_trace
#' @importFrom ggplot2 ggtitle ggsave ggplot aes geom_histogram labs theme
#' @importFrom stats glm predict quantile
#'
#' @name mcmc_diag
#' @docType class
#' @export
mcmc_diagnosis <- function(fit_out,
                           pars_to_check = c("beta0", "theta0",
                                             "beta1", "theta1", "theta_lag"),
                           save_plots    = FALSE,
                           plot_prefix   = "traceplot_",
                           positivity    = FALSE) {

  # sanity
  if (is.null(fit_out$stan_fit))
    stop("stan_fit is missing in fit_out")
  stan_fit <- fit_out$stan_fit

  # convergence stats
  cat("----- MCMC Rhat & Effective Sample Size -----\n")
  sum_stats <- rstan::summary(stan_fit, pars = pars_to_check)$summary
  df_stats  <- data.frame(
    Parameter = rownames(sum_stats),
    n_eff     = sum_stats[, "n_eff"],
    Rhat      = sum_stats[, "Rhat"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  print(df_stats)
  cat("(Values close to Rhat = 1 and large n_eff indicate good convergence.)\n\n")

  # trace‑plots
  plots <- lapply(pars_to_check, function(par) {
    p <- bayesplot::mcmc_trace(as.array(stan_fit), regex_pars = par) +
      ggplot2::ggtitle(paste0("Traceplot: ", par))
    if (save_plots)
      ggplot2::ggsave(filename = paste0(plot_prefix, par, ".png"), plot = p)
    p
  })
  names(plots) <- pars_to_check

  # positivity diagnostics
  if (positivity) {
    if (is.null(fit_out$data_preprocessed)) {
      warning("data_preprocessed not found in fit_out; cannot perform positivity diagnostics.")
    } else {
      cat("Positivity Diagnostics for Treatment A \n")
      df <- fit_out$data_preprocessed
      exclude <- c("pat_id", "k_idx", "Y_obs", "T_obs",
                   "Y_prev", "T_prev", "A")
      covariates <- setdiff(names(df), exclude)

      formula_A <- stats::as.formula(paste("A ~",
                                           paste(c("Y_prev", "k_idx", covariates),
                                                 collapse = " + ")))
      mod_ps <- stats::glm(formula_A, data = df, family = stats::binomial)
      ps     <- stats::predict(mod_ps, type = "response")

      cat(sprintf("Propensity score range: [%.3f, %.3f]\n", min(ps), max(ps)))
      cat("Percentiles (1%, 5%, 95%, 99%):\n")
      print(stats::quantile(ps, c(.01, .05, .95, .99)))
      cat("\n")

      p_ps <- ggplot2::ggplot(data.frame(ps = ps),
                              ggplot2::aes(x = ps)) +
        ggplot2::geom_histogram(bins = 30, boundary = 0) +
        ggplot2::labs(
          title = "Estimated Propensity Score Distribution for A",
          x     = "P(A = 1 | history, covariates)",
          y     = "Count"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5))
      print(p_ps)
      if (save_plots)
        ggplot2::ggsave(filename = paste0(plot_prefix, "ps_hist.png"), plot = p_ps)
    }
  }

  # return
  out <- list(stats = df_stats, plots = plots)
  class(out) <- "mcmc_diag"
  invisible(out)
}


# print / summary / plot methods

#' @describeIn mcmc_diag Print the table of MCMC convergence statistics (R‑hat & n_eff).
#' @export
print.mcmc_diag <- function(x, ...) {
  cat("MCMC convergence diagnostics (R-hat & n_eff):\n")
  print(x$stats)
  invisible(x$stats)
}

#' @describeIn mcmc_diag Same as \code{print()} for a \code{mcmc_diag} object.
#' @method summary mcmc_diag
#' @export
summary.mcmc_diag <- function(object, ...) {
  print(object)
}

#' @describeIn mcmc_diag Display stored trace‑plots; optionally filter by \code{pars}.
#' @export
plot.mcmc_diag <- function(x, pars = NULL, ...) {
  plots <- x$plots
  if (!is.null(pars)) {
    plots <- plots[intersect(names(plots), pars)]
  }
  if (length(plots) == 0) {
    warning("No traceplots to display for the specified parameters.")
    return(invisible(NULL))
  }
  for (p in plots) print(p)
  invisible(plots)
}
