# mcmc_diagnosis

#' MCMC Convergence and Positivity Diagnostics
#'
#' @description
#' Summarize MCMC convergence for a fitted model (R-hat and effective sample size),
#' generate traceplots for specified parameters, and optionally perform positivity
#' (overlap) diagnostics on treatment assignment.
#'
#' @param fit_out Output list from \code{fit_causal_recur()}, containing at least
#'   \code{stan_fit} (an RStan stanfit object) and \code{data_preprocessed} if
#'   \code{positivity = TRUE}.
#' @param pars_to_check Character vector of parameter name patterns to diagnose
#'   (default: \code{c("beta0","gamma0","beta_Y","beta_A","gamma_Y","gamma_A")}).
#'   These are passed as regex to \code{bayesplot::mcmc_trace}.
#' @param save_plots Logical; if \code{TRUE}, saves each traceplot as a PNG file
#'   with prefix \code{plot_prefix}. Default FALSE.
#' @param plot_prefix Filename prefix for saved plots (default: \code{"traceplot_"}).
#' @param positivity Logical; if \code{TRUE}, perform additional positivity
#'   diagnostics on treatment assignment (default: \code{FALSE}).
#'
#' @return An object of class \code{mcmc_diag}, a list with components:
#'   \describe{
#'     \item{\code{stats}}{Data.frame with columns \code{Parameter}, \code{n_eff}, and \code{Rhat}.}
#'     \item{\code{plots}}{Named list of ggplot objects for traceplots, names are the parameter patterns.}
#'   }
#'   Use \code{print()}, \code{summary()}, or \code{plot()} on the returned object to inspect results.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Computes R-hat and effective sample size for specified parameters via \code{rstan::summary()}.
#'   \item Generates traceplots for each parameter pattern using \code{bayesplot::mcmc_trace()}.
#'   \item If \code{positivity = TRUE}, fits a logistic regression of treatment on history/covariates
#'         from \code{fit_out$data_preprocessed}, prints propensity score range and percentiles,
#'         and plots a histogram of propensity scores.
#' }
#' The returned object has class \code{mcmc_diag}, so users can call
#' \code{print()}, \code{summary()}, or \code{plot()} directly.
#'
#' @examples
#' \dontrun{
#' fit_out <- fit_causal_recur(...)
#' diag <- mcmc_diagnosis(
#'   fit_out       = fit_out,
#'   pars_to_check = c("beta_Y", "beta_A"),
#'   save_plots    = FALSE,
#'   positivity    = FALSE
#' )
#' print(diag)        # prints R-hat and n_eff table
#' summary(diag)      # same as print
#' plot(diag)         # shows all traceplots
#' plot(diag, pars = "beta_Y")  # show traceplot(s) matching "beta_Y"
#' }
#'
#' @importFrom rstan summary
#' @importFrom bayesplot mcmc_trace
#' @importFrom ggplot2 ggtitle ggsave ggplot aes geom_histogram labs theme
#' @importFrom stats glm predict quantile
#' @export


mcmc_diagnosis <- function(fit_out,
                           pars_to_check = c("beta0", "gamma0", "beta_Y", "beta_A", "gamma_Y", "gamma_A"),
                           save_plots    = FALSE,
                           plot_prefix   = "traceplot_",
                           positivity    = FALSE) {
  # Check that stan_fit is present
  if (is.null(fit_out$stan_fit)) {
    stop("stan_fit is missing in fit_out")
  }
  stan_fit <- fit_out$stan_fit

  # 1) Compute R-hat and effective sample size
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

  # 2) Generate traceplots
  plots <- lapply(pars_to_check, function(par) {
    p <- bayesplot::mcmc_trace(as.array(stan_fit), regex_pars = par) +
      ggplot2::ggtitle(paste0("Traceplot: ", par))
    if (save_plots) {
      ggplot2::ggsave(filename = paste0(plot_prefix, par, ".png"), plot = p)
    }
    p
  })
  names(plots) <- pars_to_check

  # 3) Positivity diagnostics if requested
  if (positivity) {
    if (is.null(fit_out$data_preprocessed)) {
      warning("data_preprocessed not found in fit_out; cannot perform positivity diagnostics.")
    } else {
      cat("----- Positivity Diagnostics for Treatment A -----\n")
      df <- fit_out$data_preprocessed
      # Exclude outcome/history vars; include remaining covariates
      exclude <- c("pat_id","k_idx","Y_obs","T_obs","C_obs","Y_prev","T_prev","C_prev","A")
      covariates <- setdiff(names(df), exclude)
      formula_A <- as.formula(paste("A ~", paste(c("Y_prev", "k_idx", covariates), collapse = " + ")))
      mod_ps <- stats::glm(formula_A, data = df, family = stats::binomial)
      ps     <- stats::predict(mod_ps, type = "response")
      cat(sprintf("Propensity score range: [%.3f, %.3f]\n", min(ps), max(ps)))
      cat("Propensity score percentiles (1%, 5%, 95%, 99%):\n")
      print(stats::quantile(ps, probs = c(0.01, 0.05, 0.95, 0.99)))
      cat("\n")
      df_ps <- data.frame(ps = ps)
      p_ps <- ggplot2::ggplot(df_ps, ggplot2::aes(x = ps)) +
        ggplot2::geom_histogram(bins = 30, boundary = 0) +
        ggplot2::labs(
          title = "Estimated Propensity Score Distribution for A",
          x     = "P(A=1 | history, covariates)",
          y     = "Count"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      print(p_ps)
      if (save_plots) {
        ggplot2::ggsave(filename = paste0(plot_prefix, "ps_hist.png"), plot = p_ps)
      }
    }
  }

  # 4) Return structured results as class 'mcmc_diag'
  out <- list(
    stats = df_stats,
    plots = plots
  )
  class(out) <- "mcmc_diag"
  invisible(out)
}

#' Print method for mcmc_diag
#'
#' @description
#' Print the MCMC convergence statistics (R-hat and effective sample size).
#'
#' @param x An object of class \code{mcmc_diag}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the \code{stats} data.frame.
#' @rdname mcmc_diag-print
#' @method print mcmc_diag
#' @export


print.mcmc_diag <- function(x, ...) {
  cat("MCMC convergence diagnostics (R-hat & n_eff):\n")
  print(x$stats)
  invisible(x$stats)
}

#' Summary method for mcmc_diag
#'
#' @description
#' Summary for a \code{mcmc_diag} object; same as \code{print()}.
#'
#' @param object An object of class \code{mcmc_diag}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the \code{stats} data.frame.
#' @rdname mcmc_diag-summary
#' @method summary mcmc_diag
#' @export


summary.mcmc_diag <- function(object, ...) {
  print(object)
}

#' Plot method for mcmc_diag
#'
#' @description
#' Plot traceplots stored in a \code{mcmc_diag} object. Optionally filter
#' by parameter patterns.
#'
#' @param x An object of class \code{mcmc_diag}.
#' @param pars Character vector of parameter patterns to plot (default: all).
#'   Should match names used in \code{pars_to_check}. If provided, only those
#'   plots are shown.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the list of ggplot objects.
#' @rdname mcmc_diag-plot
#' @method plot mcmc_diag
#' @export
plot.mcmc_diag <- function(x, pars = NULL, ...) {
  plots <- x$plots
  if (!is.null(pars)) {
    # keep only those matching exactly in names; user may supply subset of names
    plots <- plots[intersect(names(plots), pars)]
  }
  if (length(plots) == 0) {
    warning("No traceplots to display for the specified parameters.")
    return(invisible(NULL))
  }
  for (p in plots) {
    print(p)
  }
  invisible(plots)
}
