# g_computation

#' Bayesian g‑Computation for Recurrent‑Event Rate Contrasts
#'
#' Perform Bayesian g‑computation to estimate average recurrent‑event
#' rates under different treatment–initiation times versus never treating,
#' using the **teacher‑style** model parameters returned by
#' \code{\link{fit_causal_recur}}.
#'
#' @param fit_out Output list from \code{\link{fit_causal_recur}}.
#' @param s_vec Integer vector of treatment‑start intervals.
#' @param B     Number of Monte‑Carlo replicates per posterior draw (default
#'   \code{50}).
#'
#' @return An object of class \code{gcomp_out}, a list with components:
#'   \describe{
#'     \item{R_mat}{Numeric matrix with one row per posterior draw and one
#'       column per treatment strategy (the last column is “never treat”,
#'       \code{s = K + 1}).  Each entry is the Dirichlet‑weighted mean event
#'       rate under that strategy.}
#'     \item{delta}{Named list of length \code{length(s_vec)}.  Each element
#'       contains the posterior draws of \eqn{\Delta(s)=R(s)-R(K+1)} together
#'       with its mean and 95% equal‑tail interval.}
#'   }
#'
#' @details
#' For each posterior draw \eqn{m=1,\dots,M} the algorithm:
#' \enumerate{
#'   \item Generates Dirichlet weights across subjects (to mimic the Bayesian
#'         bootstrap).
#'   \item For every treatment start time \eqn{s\in s\_vec} \emph{and} for the
#'         “never‑treat” strategy (\eqn{s=K+1}) it calls an internal simulator
#'         that propagates the teacher‑style parameters
#'         \eqn{\{\beta_0(k),\beta_1,\beta_L,\theta_0(k),\theta_1,\theta_L,
#'         \theta_{\mathrm{lag}}\}} over \eqn{B} replicate paths of
#'         death/recurrent events.
#'   \item Computes the subject‑level mean event rate in each replicate, then
#'         Dirichlet‑averages those rates to obtain \eqn{R_m(s)}.
#'   \item Stores all \eqn{R_m(s)} and their contrasts
#'         \eqn{\Delta_m(s)=R_m(s)-R_m(K+1)}.
#' }
#'
#' S3 methods \code{print}, \code{summary}, and \code{plot} are available; e.g.
#' \code{print(gcomp_out)}, \code{plot(gcomp_out,interactive=TRUE)}.
#'
#' @examples
#' \dontrun{
#' gcomp_out <- g_computation(fit_out = fit, s_vec = 1:K, B = 20)
#' print(gcomp_out)
#' plot(gcomp_out, ref_line = 0)
#' }
#'
#' @importFrom rstan extract
#' @importFrom stats rgamma rbinom rpois weighted.mean quantile
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar geom_line geom_point
#'   geom_text labs scale_color_manual scale_shape_manual theme position_dodge
#'   ggsave
#' @importFrom dplyr bind_rows
#' @importFrom plotly ggplotly
#' @keywords internal
#' @export


g_computation <- function(fit_out, s_vec, B = 50) {

  if (!inherits(fit_out$stan_fit, "stanfit"))
    stop("fit_out$stan_fit must be a 'stanfit' object")

  df      <- fit_out$data_preprocessed
  n_pat   <- fit_out$n_pat
  K       <- fit_out$K

  X_full  <- model.matrix(fit_out$design_info$formula_T, data = df)
  Lmat    <- X_full[df$k_idx == 1, -1, drop = FALSE]
  P_data  <- ncol(Lmat)
  if (P_data == 0) Lmat <- matrix(0, n_pat, 0)

  post <- rstan::extract(
    fit_out$stan_fit,
    pars = c("beta0","beta1","betaL",
             "theta0","theta1","thetaL","theta_lag"),
    permuted = TRUE
  )
  ndraws <- length(post$beta1)

  P_stan <- if (is.null(dim(post$betaL))) 0 else dim(post$betaL)[2]

  take_betaL  <- function(m) {
    if (P_stan == 0L || P_data == 0L) numeric(0)
    else as.numeric(post$betaL[m, 1:min(P_data, P_stan), drop = TRUE])
  }
  take_thetaL <- function(m) {
    if (P_stan == 0L || P_data == 0L) numeric(0)
    else as.numeric(post$thetaL[m, 1:min(P_data, P_stan), drop = TRUE])
  }

  ## helpers
  invlogit <- function(x) .Call(stats:::C_logit_linkinv, x)
  rbern    <- function(n, p) stats::rbinom(n, 1, p)

  sim_interv <- function(Lmat, n,
                         beta0, beta1, betaL,
                         theta0, theta1, thetaL, theta_lag,
                         s, K, B) {

    eta_beta  <- if (length(betaL))  as.numeric(Lmat %*% betaL)  else rep(0, n)
    eta_theta <- if (length(thetaL)) as.numeric(Lmat %*% thetaL) else rep(0, n)

    ir_shell <- matrix(NA_real_, n, B)

    for (b in 1:B) {
      ybar <- tbar <- abar <- matrix(0, n, K)

      ## k = 1
      abar[, 1] <- as.integer(1 >= s)
      ybar[, 1] <- rpois(
        n,
        exp(theta0[1] + abar[, 1] * theta1 + eta_theta)
      )

      ## k = 2..K
      for (k in 2:K) {
        abar[, k] <- as.integer(k >= s)

        ## death process
        eta_k <- beta0[k] + abar[, k] * beta1 + eta_beta     # length n
        tbar[, k] <- rbern(n, plogis(eta_k))                 # use plogis for safety

        ## recurrent count
        mu_k <- theta0[k] + abar[, k] * theta1 + eta_theta +
          theta_lag * (ybar[, k - 1] > 0)
        ybar[, k] <- rpois(n, exp(mu_k))

        ## propagate death forward
        dead_prev <- tbar[, k - 1] == 1
        tbar[dead_prev, k] <- 1
        ybar[dead_prev, k] <- 0
      }

      ir_shell[, b] <- rowSums(ybar) / (K - rowSums(tbar))
    }

    rowMeans(ir_shell)
  }

  one_draw <- function(m, s) {
    sim_interv(Lmat, n_pat,
               beta0      = post$beta0[m, ],
               beta1      = post$beta1[m],
               betaL      = take_betaL(m),
               theta0     = post$theta0[m, ],
               theta1     = post$theta1[m],
               thetaL     = take_thetaL(m),
               theta_lag  = post$theta_lag[m],
               s = s, K = K, B = B)
  }

  ## Dirichlet weighting
  W      <- matrix(stats::rgamma(ndraws * n_pat, 1, 1), ndraws, n_pat)
  pi_mat <- W / rowSums(W)

  calc_row <- function(m) {
    base <- one_draw(m, K + 1)  # never treat
    c(vapply(s_vec, function(s) one_draw(m, s), numeric(n_pat)),
      base)
  }

  R_mat <- if (requireNamespace("future.apply", quietly = TRUE)) {
    do.call(rbind,
            future.apply::future_lapply(seq_len(ndraws), calc_row,
                                        future.seed = TRUE))
  } else {
    t(vapply(seq_len(ndraws), calc_row,
             numeric(length(s_vec) + 1)))
  }

  R_mat <- R_mat * pi_mat
  R_mat <- t(apply(R_mat, 1, colSums))

  colnames(R_mat) <- c(paste0("s=", s_vec), paste0("s=", K + 1))

  delta <- lapply(seq_along(s_vec), function(j) {
    d <- R_mat[, j] - R_mat[, length(s_vec) + 1]
    qs <- stats::quantile(d, c(.025, .975), na.rm = TRUE)
    list(draws = d,
         mean  = mean(d, na.rm = TRUE),
         CI_lower = qs[1], CI_upper = qs[2])
  })
  names(delta) <- paste0("s=", s_vec)

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
