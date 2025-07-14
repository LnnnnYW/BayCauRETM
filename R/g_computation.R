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
  if (!inherits(fit_out$stan_fit, "stanfit"))
    stop("fit_out$stan_fit must be a 'stanfit' object")

  df    <- fit_out$data_preprocessed
  n_pat <- fit_out$n_pat
  K     <- fit_out$K

  # Posterior draws
  post <- rstan::extract(
    fit_out$stan_fit,
    pars = c("beta0", "gamma0", "beta_X", "beta_Y", "beta_A",
             "gamma_X", "gamma_Y", "gamma_A"),
    permuted = TRUE
  )
  M <- length(post$beta_Y)                 # number of posterior iterations

  # Dirichlet weights
  W      <- matrix(stats::rgamma(M * n_pat, shape = 10, rate = 10),
                   nrow = M, ncol = n_pat)
  pi_mat <- W / rowSums(W)

  # Baseline covariates
  X_cov_full <- model.matrix(fit_out$design_info$formula_T, data = df)
  X_cov      <- X_cov_full[, colnames(X_cov_full) != "(Intercept)", drop = FALSE]
  X_cov_baseline <- X_cov[df$k_idx == 1, , drop = FALSE]   # n_pat × p
  p_baseline <- ncol(X_cov_baseline)

  # helper: simulate one posterior draw
  simulate_once <- function(m, s_start) {
    beta0_m   <- post$beta0[m, ]
    gamma0_m  <- post$gamma0[m, ]
    beta_X_m  <- post$beta_X[m, ]
    beta_Y_m  <- post$beta_Y[m]
    beta_A_m  <- post$beta_A[m]
    gamma_X_m <- post$gamma_X[m, ]
    gamma_Y_m <- post$gamma_Y[m]
    gamma_A_m <- post$gamma_A[m]

    # baseline contribution
    if (p_baseline == 0 || length(beta_X_m) != p_baseline || length(gamma_X_m) != p_baseline) {
      x_beta  <- numeric(n_pat)
      x_gamma <- numeric(n_pat)
    } else {
      x_beta  <- as.numeric(X_cov_baseline %*% beta_X_m)   # n_pat
      x_gamma <- as.numeric(X_cov_baseline %*% gamma_X_m)
    }

    # replicate across bootstrap replicates
    x_beta_mat  <- matrix(x_beta,  n_pat, B)
    x_gamma_mat <- matrix(x_gamma, n_pat, B)

    # state matrices
    alive   <- matrix(TRUE,  n_pat, B)
    Y_prev  <- matrix(0L,    n_pat, B)
    num_Y   <- matrix(0L,    n_pat, B)
    num_atR <- matrix(0L,    n_pat, B)

    for (k in seq_len(K)) {
      if (!any(alive)) break

      A_k <- as.integer(k >= s_start)

      # (a) terminal event
      p_death <- stats::plogis(beta0_m[k] + x_beta_mat + beta_Y_m * Y_prev + beta_A_m * A_k)
      death_draw <- matrix(stats::rbinom(length(p_death), 1, p_death), n_pat, B)
      alive <- alive & (death_draw == 0L)

      # (b) recurrent events
      if (any(alive)) {
        mu_k <- exp(gamma0_m[k] + x_gamma_mat + gamma_Y_m * Y_prev + gamma_A_m * A_k)
        mu_k <- pmin(pmax(mu_k, 1e-6), 1e4)
        Y_k  <- matrix(stats::rpois(length(mu_k), mu_k), n_pat, B)
        Y_k[!alive] <- 0L

        num_atR[alive] <- num_atR[alive] + 1L
        num_Y   <- num_Y + Y_k
        Y_prev  <- Y_k
      }
    }

    # per-subject mean over B replicates
    rate_rep <- ifelse(num_atR > 0, num_Y / num_atR, 0)
    rowMeans(rate_rep)
  }

  S        <- length(s_vec)
  R_mat    <- matrix(NA_real_, nrow = M, ncol = S + 1)
  colnames(R_mat) <- c(paste0("s=", s_vec), paste0("s=", K + 1))

  # Parallel over posterior draws
  if (requireNamespace("future.apply", quietly = TRUE)) {
    R_rows <- future.apply::future_lapply(seq_len(M), function(m) {
      out <- numeric(S + 1)
      out[S + 1] <- stats::weighted.mean(simulate_once(m, K + 1), w = pi_mat[m, ])
      for (j in seq_len(S))
        out[j] <- stats::weighted.mean(simulate_once(m, s_vec[j]), w = pi_mat[m, ])
      out
    }, future.seed = TRUE)
    R_mat <- do.call(rbind, R_rows)
  } else {
    warning("future.apply not installed – falling back to sequential loop; performance will be slower.")
    for (m in seq_len(M)) {
      R_mat[m, S + 1] <- stats::weighted.mean(simulate_once(m, K + 1), w = pi_mat[m, ])
      for (j in seq_len(S))
        R_mat[m, j] <- stats::weighted.mean(simulate_once(m, s_vec[j]), w = pi_mat[m, ])
    }
  }

  # Contrasts
  delta <- setNames(vector("list", S), paste0("s=", s_vec))
  for (j in seq_len(S)) {
    diff <- R_mat[, j] - R_mat[, S + 1]
    good <- diff[is.finite(diff)]
    if (length(good)) {
      qs <- stats::quantile(good, c(.025, .975))
      delta[[j]] <- list(draws = good, mean = mean(good), CI_lower = qs[1], CI_upper = qs[2])
    } else {
      delta[[j]] <- list(draws = numeric(0), mean = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_)
    }
  }

  structure(list(R_mat = R_mat, delta = delta), class = "gcomp_out")
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
