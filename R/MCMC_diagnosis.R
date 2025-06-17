#' MCMC Convergence and Positivity Diagnostics
#'
#' Print R-hat and effective sample size, generate traceplots, and
#' optionally check positivity of the treatment mechanism.
#'
#' @param fit_out Output list from \code{fit_causal_recur()}, containing
#'   \code{stan_fit}.
#' @param pars_to_check Character vector of parameter names to diagnose (default:
#'   \code{c("beta0","gamma0","beta_Y","beta_A","gamma_Y","gamma_A")}).
#' @param save_plots Logical; if \code{TRUE}, saves each traceplot as a PNG file.
#' @param plot_prefix Filename prefix for saved plots (default: \code{"traceplot_"}).
#' @param positivity Logical; if \code{TRUE}, perform additional positivity (overlap) diagnostics on the treatment assignment (default: \code{FALSE}).
#'
#' @return Invisibly returns \code{NULL}. Side effects: prints summary and creates plots.
#'
#' @details
#' Summarizes MCMC convergence using R-hat and effective sample size.
#' Generates traceplots for specified parameters. Optionally, traceplots
#' are saved as PNG files in the working directory, with filenames
#' beginning with \code{plot_prefix}. If \code{positivity = TRUE}, the function will also check for lack of overlap in the propensity model.
#'
#' @examples
#' \dontrun{
#' # Example assumes fit_causal_recur() has been run:
#' fit_out <- fit_causal_recur(...)
#' mcmc_diagnosis(fit_out, pars_to_check = c("beta_Y", "beta_A"))
#' }
#'
#' @importFrom rstan summary
#' @importFrom bayesplot mcmc_trace
#' @importFrom ggplot2 ggtitle ggsave
#' @export





mcmc_diagnosis <- function(fit_out,
                           pars_to_check = c("beta0", "gamma0", "beta_Y", "beta_A", "gamma_Y", "gamma_A"),
                           save_plots    = FALSE,
                           plot_prefix   = "traceplot_",
                           positivity    = TRUE) {
  # — Require a stanfit object —
    if (is.null(fit_out$stan_fit)) {
      stop("stan_fit is missing in fit_out")
      }
  if (!inherits(fit_out$stan_fit, "stanfit")) {
    stop("stan_fit must be a stanfit object")
    }
  # 1) Extract stanfit and show Rhat / n_eff
  stan_fit <- fit_out$stan_fit
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
  cat("(Values close to Rhat=1 and large n_eff indicate good convergence.)\n\n")

  # 2) Traceplots
  for (par in pars_to_check) {
    plt <- bayesplot::mcmc_trace(as.array(stan_fit), regex_pars = par) +
      ggplot2::ggtitle(paste0("Traceplot: ", par))
    print(plt)
    if (save_plots) {
      ggplot2::ggsave(
        filename = paste0(plot_prefix, par, ".png"),
        plot     = plt
      )
    }
  }

  # 3) Causal‐positivity diagnostics
  if (positivity) {
    cat("----- Positivity Diagnostics for Treatment A -----\n")
    df <- fit_out$data_preprocessed

    # Identify covariates: everything except key outcome/history vars
    exclude <- c("pat_id","k_idx",
                 "Y_obs","T_obs","C_obs",
                 "Y_prev","T_prev","C_prev","A")
    covariate_cols <- setdiff(names(df), exclude)

    # Build logistic model formula: A ~ Y_prev + k_idx + other covariates
    preds     <- c("Y_prev", "k_idx", covariate_cols)
    formula_A <- as.formula(paste("A ~", paste(preds, collapse = " + ")))

    # Fit and predict propensity scores
    mod_ps <- stats::glm(formula_A, data = df, family = stats::binomial)
    ps     <- stats::predict(mod_ps, type = "response")

    # Summarize extremes & key quantiles
    cat(sprintf("Propensity score range: [%.3f, %.3f]\n", min(ps), max(ps)))
    cat("Propensity score percentiles (1%, 5%, 95%, 99%):\n")
    print(stats::quantile(ps, probs = c(0.01, 0.05, 0.95, 0.99)))
    cat("\n")

    # Plot distribution
    df_ps <- data.frame(ps = ps)
    p_ps  <- ggplot2::ggplot(df_ps, ggplot2::aes(x = ps)) +
      ggplot2::geom_histogram(bins = 30, boundary = 0) +
      ggplot2::labs(
        title = "Estimated Propensity Score Distribution for A",
        x     = "P(A=1 | history, covariates)",
        y     = "Count"
      )
    print(p_ps)
    if (save_plots) {
      ggplot2::ggsave(
        filename = paste0(plot_prefix, "ps_hist.png"),
        plot     = p_ps
      )
    }
  }

  invisible(NULL)
}
