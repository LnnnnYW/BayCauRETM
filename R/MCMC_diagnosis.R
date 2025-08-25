#' MCMC convergence and positivity diagnostics
#'
#' @description
#' Summarize MCMC convergence for a fitted model (R-hat and effective sample
#' size), generate trace plots for selected parameters via
#' `bayesplot::mcmc_trace()`, and optionally perform positivity diagnostics
#' for treatment assignment.
#'
#' @param fit_out Output list from `fit_causal_recur()`, containing `stan_fit`
#'   (an rstan `stanfit` object) and, if `positivity = TRUE`, `data_preprocessed`.
#' @param pars_to_check Character vector of parameter-name patterns to diagnose
#'   (patterns are passed to `bayesplot::mcmc_trace()`).
#' @param save_plots Logical; if `TRUE`, save each trace plot as a PNG using `plot_prefix`.
#' @param plot_prefix Filename prefix for saved plots (default `"traceplot_"`).
#' @param positivity Logical; if `TRUE`, run an additional logistic-regression
#'   propensity-score check (default `FALSE`).
#' @param ps_covariates Character vector of column names used in the positivity
#'   model (required when `positivity = TRUE`).
#'
#' @return An object of class `mcmc_diag`, a list with:
#'   * `stats` — data frame with `Parameter`, `n_eff`, and `Rhat`;
#'   * `plots` — named list of ggplot objects for the trace plots.
#'
#' @examples
#' \dontrun{
#' # diag <- mcmc_diagnosis(fit_out, pars_to_check = c("beta0","theta0"))
#' # print(diag)
#' # plot(diag, pars = "beta0")
#' }
#'
#' @importFrom rstan summary
#' @importFrom bayesplot mcmc_trace
#' @importFrom ggplot2 ggtitle ggsave ggplot aes geom_histogram labs theme
#' @importFrom stats glm predict quantile
#' @name mcmc_diag
#' @docType class
#' @export



mcmc_diagnosis <- function(fit_out,
                           pars_to_check = c("beta0", "theta0", "beta1", "theta1", "thetaLag"),
                           save_plots     = FALSE,
                           plot_prefix    = "traceplot_",
                           positivity     = FALSE,
                           ps_covariates  = NULL) {

  if (inherits(fit_out, "causal_recur_fit"))
    fit_out <- unclass(fit_out)

  stan_fit <- fit_out$stan_fit

  if (!("stan_fit" %in% names(fit_out)) || !inherits(fit_out$stan_fit, "stanfit"))
    stop("input must contain a valid 'stan_fit' (rstan::stanfit object)")

  stan_pars <- names(rstan::extract(stan_fit))
  use_pars  <- intersect(pars_to_check, stan_pars)
  if (length(use_pars) == 0)
    stop("None of the specified 'pars_to_check' exist in the fitted Stan object.")

  cat("----- MCMC Rhat & Effective Sample Size -----\n")
  sum_stats <- rstan::summary(stan_fit, pars = use_pars)$summary
  df_stats  <- data.frame(
    Parameter = rownames(sum_stats),
    n_eff     = sum_stats[, "n_eff"],
    Rhat      = sum_stats[, "Rhat"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  df_stats$Parameter <- .map_param_names(fit_out, df_stats$Parameter)

  print(df_stats)
  cat("(Values close to Rhat = 1 and large n_eff indicate good convergence.)\n\n")

  # trace-plots
  .title_pretty <- function(tag) {
    switch(tag,
           "beta0"    = "T-model: covariates + lags + time-baseline",
           "theta0"   = "Y-model: covariates + lags + time-baseline",
           "beta1"    = "T-model: treatment_effect_T",
           "theta1"   = "Y-model: treatment_effect_Y",
           "thetaLag" = "Lag-kernel (theta_lag_extra)",
           tag
    )
  }

  bayesplot::color_scheme_set("blue")

  arr_all  <- as.array(stan_fit)
  all_vars <- dimnames(arr_all)[[3]]

  plots <- lapply(use_pars, function(par) {
    idx <- grep(paste0("^", par, "(\\[|$|_star$)"), all_vars)
    if (length(idx) == 0L) return(NULL)

    vars_pretty <- .map_param_names(fit_out, all_vars[idx])
    arr_sub     <- arr_all[, , idx, drop = FALSE]
    dimnames(arr_sub)[[3]] <- vars_pretty

    n_pan <- length(vars_pretty)
    ncol  <- 2L
    nrow  <- ceiling(n_pan / ncol)

    fig_w <- 11
    fig_h <- max(6.5, 3.0 + 2.4 * nrow)

    p <- bayesplot::mcmc_trace(
      arr_sub,
      pars = vars_pretty,
      facet_args = list(scales = "free_y", ncol = ncol)
    ) +
      ggplot2::labs(title = paste0("Traceplot: ", .title_pretty(par))) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
        strip.text    = ggplot2::element_text(size = 12, face = "bold"),
        axis.title    = ggplot2::element_text(size = 12),
        axis.text     = ggplot2::element_text(size = 10),
        panel.spacing = grid::unit(1.4, "lines"),
        plot.margin   = grid::unit(rep(0.6, 4), "lines")
      )

    if (save_plots) {
      ggplot2::ggsave(paste0(plot_prefix, par, ".png"),
                      plot = p, width = fig_w, height = fig_h, dpi = 150)
    }
    attr(p, "BayCauRETM_fig_width")  <- fig_w
    attr(p, "BayCauRETM_fig_height") <- fig_h
    p
  })
  plots <- Filter(Negate(is.null), plots)
  names(plots) <- use_pars[seq_along(plots)]

  # positivity diagnostics
  if (positivity) {
    if (is.null(fit_out$data_preprocessed)) {
      warning("data_preprocessed not found in fit_out; cannot perform positivity diagnostics.")
    } else {
      cat("Positivity Diagnostics for Treatment A \n")

      df <- fit_out$data_preprocessed

      if (is.null(ps_covariates))
        stop("positivity = TRUE but no 'ps_covariates' provided.")

      missing <- setdiff(ps_covariates, names(df))
      if (length(missing) > 0) {
        stop("ps_covariates not found in data: ", paste(missing, collapse = ", "))
      }

      formula_A <- as.formula(paste("A ~", paste(ps_covariates, collapse = " + ")))
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

  out <- list(stats = df_stats, plots = plots)
  class(out) <- "mcmc_diag"
  invisible(out)
}



# print / summary / plot methods

#' @describeIn mcmc_diag Print the table of MCMC convergence statistics (R-hat & n_eff).
#' @param x An `mcmc_diag` object.
#' @param ... Additional arguments (ignored).
#' @export

print.mcmc_diag <- function(x, ...) {
  cat("MCMC convergence diagnostics (R-hat & n_eff):\n")
  print(x$stats)
  invisible(x$stats)
}

#' @describeIn mcmc_diag Same as `print()` for an `mcmc_diag` object.
#' @param object An `mcmc_diag` object.
#' @param ... Additional arguments (ignored).
#' @method summary mcmc_diag
#' @export

summary.mcmc_diag <- function(object, ...) {
  print(object)
}

#' @describeIn mcmc_diag Display stored trace plots; optionally filter by `pars`.
#' @param x An `mcmc_diag` object.
#' @param pars Optional character vector of parameter names to display.
#' @param ... Additional arguments (ignored).
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
