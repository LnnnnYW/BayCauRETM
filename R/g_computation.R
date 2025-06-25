# g_computation
#' Bayesian g-Computation for Recurrent-Event Rate Contrasts
#'
#' Perform Bayesian g-computation to estimate average recurrent-event
#' rates under different treatment initiation times versus never treating.
#'
#' @param fit_out Output list from \code{\link{fit_causal_recur}}.
#' @param s_vec Integer vector of treatment-start intervals.
#' @param B     Number of Monte-Carlo replicates per draw (default 50).
#'
#' @return An object of class \code{gcomp_out}, a list with components:
#'   \describe{
#'     \item{R_mat}{Numeric matrix of dimension M * (length(s_vec)+1), where M is the
#'       number of posterior draws. Columns are the weighted mean event rates under each
#'       treatment initiation strategy (including "never treat" as the last column).}
#'     \item{delta}{Named list of length \code{length(s_vec)}, each element is a list with:
#'       \code{draws} (numeric vector of contrasts), \code{mean}, \code{CI_lower}, \code{CI_upper}.}
#'   }
#'
#' @details
#' For each posterior draw m:
#' \enumerate{
#'   \item Draw Dirichlet weights across subjects.
#'   \item For each treatment start time s in \code{s_vec}, simulate B replicate paths
#'         of death and recurrent-event counts over intervals 1..K, conditional on starting
#'         treatment at s (or never treating for s = K+1).
#'   \item Compute subject-level average event rate (total events / time at risk) for each replicate,
#'         then average over replicates and weight by Dirichlet weights to obtain the mean rate R(s).
#'   \item Store R(s) for each s, and for "never treat" (s = K+1). Compute contrast delta(s) = R(s) − R(K+1).
#' }
#' The returned object has class \code{gcomp_out}. S3 methods \code{print}, \code{summary}, and \code{plot}
#' are provided so that users can call \code{print(gcomp_out)}, \code{summary(gcomp_out)}, or \code{plot(gcomp_out, ...)}.
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is the output of fit_causal_recur(...)
#' gcomp_out <- g_computation(fit_out = fit, s_vec = 1:K, B = 20)
#' print(gcomp_out)                # print summary table of delta(s)
#' summary(gcomp_out)              # same as print
#' plot(gcomp_out)                 # static ggplot of delta(s) vs s, with reference line at 0
#' plot(gcomp_out, interactive = TRUE)  # interactive Plotly plot (requires plotly)
#' plot(gcomp_out, save_file = "delta_plot.png")  # save static plot
#' plot(gcomp_out, s_vec = 1:5, line_size = 1.2, ribbon_alpha = 0.2)
#' }
#'
#' @importFrom rstan extract
#' @importFrom stats rgamma rbinom rpois weighted.mean quantile
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar geom_line geom_point geom_text labs scale_color_manual scale_shape_manual theme position_dodge ggsave
#' @importFrom dplyr bind_rows
#' @importFrom plotly ggplotly
#' @importFrom stats as.formula median na.omit setNames
#' @importFrom utils write.csv
#' @export


g_computation <- function(fit_out, s_vec, B = 50) {
  # Validate input
  if (is.null(fit_out$stan_fit)) stop("stan_fit is missing in fit_out")
  if (!inherits(fit_out$stan_fit, "stanfit")) stop("stan_fit must be an RStan stanfit object")
  df    <- fit_out$data_preprocessed
  n_pat <- fit_out$n_pat
  K     <- fit_out$K

  # 1) Extract posterior draws
  post <- rstan::extract(
    fit_out$stan_fit,
    pars = c("beta0", "gamma0", "beta_X", "beta_Y", "beta_A",
             "gamma_X", "gamma_Y", "gamma_A"),
    permuted = TRUE
  )
  M <- length(post$beta_Y)

  # 2) Dirichlet weights for each draw
  W      <- matrix(stats::rgamma(M * n_pat, shape = 10, rate = 10),
                   nrow = M, ncol = n_pat)
  pi_mat <- W / rowSums(W)

  # 3) Baseline covariates at k = 1
  mat_cov        <- stats::model.matrix(fit_out$design_info$formula_T, data = df)
  X_cov          <- mat_cov[, colnames(mat_cov) != "(Intercept)", drop = FALSE]
  idx_baseline   <- which(df$k_idx == 1)
  X_cov_baseline <- X_cov[idx_baseline, , drop = FALSE]

  # 4) Simulation function for one posterior draw m and start time s_start
  simulate_once <- function(m, s_start) {
    # Extract parameter vectors for draw m
    beta0_m   <- post$beta0[m, ]
    beta_X_m  <- post$beta_X[m, ]
    beta_Y_m  <- post$beta_Y[m]
    beta_A_m  <- post$beta_A[m]
    gamma0_m  <- post$gamma0[m, ]
    gamma_X_m <- post$gamma_X[m, ]
    gamma_Y_m <- post$gamma_Y[m]
    gamma_A_m <- post$gamma_A[m]

    rates <- numeric(n_pat)
    for (i in seq_len(n_pat)) {
      Xi_cov   <- X_cov_baseline[i, ]
      rate_vec <- numeric(B)
      for (b in seq_len(B)) {
        T_prev_sim <- 0L
        Y_prev_sim <- 0.0
        num_Y       <- 0.0
        num_risk    <- 0L

        for (k in seq_len(K)) {
          if (T_prev_sim == 1L) break

          A_k <- as.integer(k >= s_start)

          # (a) death hazard
          logit_h <- beta0_m[k] +
            sum(Xi_cov * beta_X_m) +
            beta_Y_m * Y_prev_sim +
            beta_A_m * A_k
          p_death <- if (is.finite(logit_h)) {
            pd <- 1 / (1 + exp(-logit_h))
            ifelse(is.na(pd), 0, pmin(pmax(pd, 0), 1))
          } else {
            0
          }

          # draw death indicator safely
          T_k_sim <- suppressWarnings(stats::rbinom(1, size = 1, prob = p_death))
          if (identical(T_k_sim, NA_integer_)) T_k_sim <- 0L
          if (T_k_sim == 1L) {
            T_prev_sim <- 1L
            next
          }

          # (b) recurrent event count
          num_risk <- num_risk + 1L
          log_mu   <- gamma0_m[k] +
            sum(Xi_cov * gamma_X_m) +
            gamma_Y_m * Y_prev_sim +
            gamma_A_m * A_k
          mu_k <- if (is.finite(log_mu)) {
            m0 <- exp(log_mu)
            m0 <- max(min(m0, 1e4), 1e-6)
            m0
          } else {
            0
          }

          # draw Poisson count safely
          Y_k_sim <- suppressWarnings(stats::rpois(1, lambda = mu_k))
          if (!is.finite(Y_k_sim)) Y_k_sim <- 0
          Y_prev_sim <- Y_k_sim
          num_Y      <- num_Y + Y_k_sim
        } # end for k

        rate_vec[b] <- if (num_risk > 0L) num_Y / num_risk else 0
      } # end for b
      rates[i] <- mean(rate_vec, na.rm = TRUE)
    } # end for i
    rates
  }

  # 5) Build R_mat: rows = draws, cols = strategies (s_vec and never treat = K+1)
  S     <- length(s_vec)
  R_mat <- matrix(NA_real_, nrow = M, ncol = S + 1)
  colnames(R_mat) <- c(paste0("s=", s_vec), paste0("s=", K + 1))

  for (m in seq_len(M)) {
    # never treat case
    R_never       <- simulate_once(m, K + 1)
    R_mat[m, S+1] <- stats::weighted.mean(R_never, w = pi_mat[m, ], na.rm = TRUE)
    # each treatment start
    for (j in seq_len(S)) {
      R_s         <- simulate_once(m, s_vec[j])
      R_mat[m, j] <- stats::weighted.mean(R_s, w = pi_mat[m, ], na.rm = TRUE)
    }
  }

  # 6) Compute contrasts delta(s) = R(s) - R(never)
  delta <- setNames(vector("list", S), paste0("s=", s_vec))
  for (j in seq_len(S)) {
    raw_vals <- R_mat[, j] - R_mat[, S + 1]
    good     <- raw_vals[is.finite(raw_vals)]
    if (length(good) == 0L) {
      warning("All delta draws non-finite for s=", s_vec[j], "; setting to NA")
      delta[[j]] <- list(
        draws    = numeric(0),
        mean     = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_
      )
    } else {
      qs <- stats::quantile(good, probs = c(0.025, 0.975), na.rm = TRUE)
      delta[[j]] <- list(
        draws    = good,
        mean     = mean(good, na.rm = TRUE),
        CI_lower = as.numeric(qs[1]),
        CI_upper = as.numeric(qs[2])
      )
    }
  }

  out <- list(R_mat = R_mat, delta = delta)
  class(out) <- "gcomp_out"
  out
}

#' Print method for g_computation output
#'
#' @description
#' Print a summary table of causal contrasts delta(s, K+1) for a \code{gcomp_out} object.
#'
#' @param x An object of class \code{gcomp_out}, output of \code{g_computation()}.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns a data.frame of the summary.
#' @rdname gcomp_out-print
#' @method print gcomp_out
#' @export
print.gcomp_out <- function(x, ...) {
  df <- do.call(rbind, lapply(names(x$delta), function(nm) {
    xi <- x$delta[[nm]]
    data.frame(
      s      = as.integer(sub("^s=", "", nm)),
      Mean   = xi$mean,
      `2.5%` = xi$CI_lower,
      `97.5%` = xi$CI_upper,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
  cat("Causal contrast delta(s, K+1) summary:\n")
  print(df, row.names = FALSE)
  invisible(df)
}

#' Summary method for g_computation output
#'
#' @description
#' Summary for a \code{gcomp_out} object; equivalent to \code{print()}.
#'
#' @param object An object of class \code{gcomp_out}.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns the same data.frame as printed.
#' @rdname gcomp_out-summary
#' @method summary gcomp_out
#' @export
summary.gcomp_out <- function(object, ...) {
  print(object)
}

#' Plot method for g_computation output
#'
#' @description
#' Plot the causal contrast delta(s, K+1) versus treatment start interval s for a
#' \code{gcomp_out} object. By default returns a static ggplot. To obtain an
#' interactive Plotly plot or save to file, set \code{interactive = TRUE} or
#' provide \code{save_file}.
#'
#' @param x An object of class \code{gcomp_out}.
#' @param s_vec Integer vector of intervals to plot. If NULL, all available s are used.
#' @param ref_line Numeric or NULL. y-coordinate for horizontal reference line (e.g., 0). Default 0.
#' @param theme_fn A ggplot2 theme function (default: \code{ggplot2::theme_minimal}).
#' @param interactive Logical. If TRUE, returns a Plotly object via \code{plotly::ggplotly()}. Default FALSE.
#' @param save_file Optional character. File path to save the static ggplot (e.g., "plot.png"). Default NULL.
#' @param width Numeric. Width in inches for saving (default: 8).
#' @param height Numeric. Height in inches for saving (default: 5).
#' @param dpi Numeric. Resolution in dpi for saving (default: 300).
#' @param ... Additional arguments passed to underlying static/interactive plot functions,
#'   such as \code{line_size}, \code{ribbon_alpha}, \code{show_points}, \code{label_points}.
#'
#' @return Invisibly returns the ggplot or plotly object.
#' @rdname gcomp_out-plot
#' @method plot gcomp_out
#' @export
plot.gcomp_out <- function(x,
                           s_vec       = NULL,
                           ref_line    = 0,
                           theme_fn    = ggplot2::theme_minimal,
                           interactive = FALSE,
                           save_file   = NULL,
                           width       = 8,
                           height      = 5,
                           dpi         = 300,
                           ...) {
  # Decide static vs interactive or save
  if (interactive || !is.null(save_file)) {
    # Use interactive wrapper; requires plot_posterior_causal_contrast_interactive to exist
    p <- plot_posterior_causal_contrast_interactive(
      contrast_list = x,
      s_vec         = s_vec,
      ref_line      = ref_line,
      theme_fn      = theme_fn,
      interactive   = interactive,
      save_file     = save_file,
      width         = width,
      height        = height,
      dpi           = dpi,
      ...
    )
  } else {
    # Static plot; requires plot_posterior_causal_contrast_static to exist
    p <- plot_posterior_causal_contrast_static(
      contrast_list = x,
      s_vec         = s_vec,
      ref_line      = ref_line,
      theme_fn      = theme_fn,
      ...
    )
  }
  print(p)
  invisible(p)
}
