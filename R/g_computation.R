#' Bayesian g-computation for recurrent-event rate contrasts
#'
#' @description
#' Estimate average recurrent-event rates under different treatment-start times
#' versus never treating, using teacher-style parameters from [fit_causal_recur()].
#'
#' @param fit_out Output list from [fit_causal_recur()].
#' @param s_vec Integer vector of treatment-start intervals.
#' @param B Integer. Number of Monte Carlo replicates per posterior draw (default 50).
#' @param Lmat Optional baseline design matrix for level/lag effects; if NULL,
#'   it is built from the first row per subject and the right-hand side of
#'   `fit_out$design_info$formula_T` (no intercept).
#' @param cores Integer. If > 1, simulate in parallel via `parallel::parLapply()`.
#'
#' @return A `gcomp_out` list with two components:
#'   `R_mat` (rows = posterior draws, columns = strategies `s = 1..K, K+1` for
#'   never treat) and `delta` (named list with draws and summaries of
#'   `Delta(s) = R(s) - R(K+1)`).
#'
#' @details
#' For each posterior draw, the algorithm draws Dirichlet weights over subjects,
#' simulates B replicated paths under each strategy using the teacher-style
#' parameters (logistic death hazard; Poisson log-mean for recurrent counts with
#' an optional lag term), averages subject-level rates with the Dirichlet weights,
#' and stores the contrasts `Delta(s)`.
#'
#' @examples
#' \donttest{
#' # g <- g_computation(fit_out = fit, s_vec = 1:K, B = 20)
#' # print(g)
#' # plot(g, ref_line = 0)
#' }
#'
#' @importFrom rstan extract
#' @importFrom stats rgamma rbinom rpois weighted.mean quantile
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar geom_line geom_point
#'   geom_text labs scale_color_manual scale_shape_manual theme position_dodge ggsave
#' @importFrom dplyr bind_rows
#' @importFrom plotly ggplotly
#' @keywords internal
#' @export



g_computation <- function(Lmat = NULL,
                          fit_out,
                          s_vec,
                          B     = 50,
                          cores = 1) {

  if (!inherits(fit_out$stan_fit, "stanfit"))
    stop("fit_out$stan_fit must be a 'stanfit' object")

  if (is.null(Lmat)) {
    df_all <- fit_out$data_preprocessed

    id_col <- if (!is.null(fit_out$design_info$id_col)) {
      fit_out$design_info$id_col
    } else if ("pat_id" %in% names(df_all)) {
      "pat_id"
    } else if ("id" %in% names(df_all)) {
      "id"
    } else {
      stop("Cannot infer id column; please supply Lmat explicitly.")
    }

    time_col <- if (!is.null(fit_out$design_info$time_col)) {
      fit_out$design_info$time_col
    } else if ("k_idx" %in% names(df_all)) {
      "k_idx"
    } else if ("k" %in% names(df_all)) {
      "k"
    } else {
      stop("Cannot infer time column; please supply Lmat explicitly.")
    }

    baseline_df <- df_all |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::slice_min(order_by = .data[[time_col]], n = 1, with_ties = FALSE) |>
      dplyr::ungroup()

    rhs_terms <- stats::delete.response(stats::terms(fit_out$design_info$formula_T))
    Lmat <- stats::model.matrix(rhs_terms, data = baseline_df)
    if ("(Intercept)" %in% colnames(Lmat)) {
      Lmat <- Lmat[, colnames(Lmat) != "(Intercept)", drop = FALSE]
    }
  }

  df      <- fit_out$data_preprocessed
  n_pat   <- fit_out$n_pat
  K       <- fit_out$K
  P_data  <- ncol(Lmat)

  post <- rstan::extract(
    fit_out$stan_fit,
    pars = c("beta0", "beta1", "betaL",
             "theta0", "theta1", "thetaL", "thetaLag"),
    permuted = TRUE
  )
  ndraws <- length(post$beta1)
  P_stan <- if (is.null(dim(post$betaL))) 0 else dim(post$betaL)[2]

  take_betaL <- function(m) {
    if (P_stan == 0L || P_data == 0L) numeric(0)
    else as.numeric(post$betaL[m, 1:min(P_data, P_stan), drop = TRUE])
  }
  take_thetaL <- function(m) {
    if (P_stan == 0L || P_data == 0L) numeric(0)
    else as.numeric(post$thetaL[m, 1:min(P_data, P_stan), drop = TRUE])
  }
  take_thetaLag <- function(m) {
    if (is.null(dim(post$thetaLag))) post$thetaLag[m]        # QlagY == 1
    else if (length(post$thetaLag) == 0) 0                   # QlagY == 0
    else post$thetaLag[m, 1]
  }

  sim_interv <- function(Lmat, n,
                         beta0, beta1, betaL,
                         theta0, theta1, thetaL, theta_lag_scalar,
                         s, K, B) {

    eta_beta  <- if (length(betaL))  as.numeric(Lmat %*% betaL)  else rep(0, n)
    eta_theta <- if (length(thetaL)) as.numeric(Lmat %*% thetaL) else rep(0, n)

    abar_vec <- as.integer(seq_len(K) >= s)

    ybar <- array(0, dim = c(n, K, B))
    tbar <- array(0, dim = c(n, K, B))

    ## k = 1 with no lag
    lambda1 <- exp(theta0[1] + abar_vec[1] * theta1 + eta_theta)
    ybar[, 1, ] <- matrix(stats::rpois(n * B, lambda1), n, B)

    ## k = 2 ... K
    for (k in 2:K) {
      abar_k <- abar_vec[k]

      ## death
      p_death <- plogis(beta0[k] + abar_k * beta1 + eta_beta)
      tbar[, k, ] <- matrix(stats::rbinom(n * B, 1, p_death), n, B)

      ## recurrent events
      mu_k <- theta0[k] + abar_k * theta1 + eta_theta +
        theta_lag_scalar * log1p(ybar[, k - 1, ])
      mu_k <- pmin(mu_k, 700)
      ybar[, k, ] <- matrix(stats::rpois(n * B, exp(mu_k)), n, B)

      ## propagate death
      dead_prev <- tbar[, k - 1, ] == 1
      tbar[dead_prev] <- 1
      ybar[dead_prev] <- 0
    }

    denom <- K - apply(tbar, c(1, 3), sum)
    num   <- apply(ybar, c(1, 3), sum)
    ir    <- ifelse(denom == 0, NA_real_, num / denom)
    rowMeans(ir, na.rm = TRUE)
  }

  ## Dirichlet Weights
  W      <- matrix(stats::rgamma(ndraws * n_pat, 1, 1), ndraws, n_pat)
  pi_mat <- W / rowSums(W)

  cl <- NULL
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit({ parallel::stopCluster(cl) }, add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("sim_interv", "Lmat", "n_pat", "K", "B",
                  "take_betaL", "take_thetaL", "take_thetaLag", "post"),
      envir = environment()
    )
  }
  apply_fun <- if (is.null(cl)) lapply else function(X, FUN)
    parallel::parLapply(cl, X, FUN)

  s_vec_base <- c(s_vec, K + 1)
  R_mat_list <- lapply(s_vec_base, function(s) {
    out_mat <- do.call(rbind,
                       apply_fun(seq_len(ndraws), function(m) {
                         sim_interv(
                           Lmat, n_pat,
                           beta0  = post$beta0[m, ],
                           beta1  = post$beta1[m],
                           betaL  = take_betaL(m),
                           theta0 = post$theta0[m, ],
                           theta1 = post$theta1[m],
                           thetaL = take_thetaL(m),
                           theta_lag_scalar = take_thetaLag(m),
                           s = s, K = K, B = B
                         )
                       })
    )
    out_mat * pi_mat     # Dirichlet re-weight
  })

  R_mat <- do.call(rbind, lapply(R_mat_list, rowSums))
  R_mat <- do.call(cbind, lapply(R_mat_list, rowSums))
  colnames(R_mat) <- c(paste0("s=", s_vec), paste0("s=", K + 1))

  delta <- lapply(seq_along(s_vec), function(j) {
    d  <- R_mat[, j] - R_mat[, length(s_vec) + 1]
    qs <- stats::quantile(d, c(.025, .975), na.rm = TRUE)
    list(draws = d,
         mean       = mean(d, na.rm = TRUE),
         CI_lower   = qs[1],
         CI_upper   = qs[2])
  })
  names(delta) <- paste0("s=", s_vec)

  structure(list(R_mat = R_mat, delta = delta),
            class = "gcomp_out")
}


# print / summary / plot methods

#' Print method for g_computation output
#'
#' @description
#' Print a summary table of causal contrasts \code{delta(s, K+1)} for a
#' \code{gcomp_out} object.
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
#' Plot the causal contrast \code{delta(s, K+1)} versus treatment-start
#' interval \code{s} for a \code{gcomp_out} object. By default, returns a static
#' \pkg{ggplot2} graphic. To obtain an interactive Plotly plot or save to file,
#' set \code{interactive = TRUE} or provide \code{save_file}.
#'
#' @param x An object of class \code{gcomp_out}.
#' @param s_vec Integer vector of intervals to plot. If \code{NULL}, all
#'   available \code{s} are used.
#' @param ref_line Numeric or \code{NULL}. y-coordinate for a horizontal
#'   reference line (e.g., \code{0}). Default \code{0}.
#' @param theme_fn A \pkg{ggplot2} theme function (default:
#'   \code{ggplot2::theme_minimal}).
#' @param interactive Logical. If \code{TRUE}, returns a Plotly object via
#'   \code{plotly::ggplotly()}. Default \code{FALSE}.
#' @param save_file Optional character. File path to save the static ggplot
#'   (e.g., \code{"plot.png"}). Default \code{NULL}.
#' @param width Numeric. Width in inches for saving (default: \code{8}).
#' @param height Numeric. Height in inches for saving (default: \code{5}).
#' @param dpi Numeric. Resolution in dpi for saving (default: \code{300}).
#' @param ... Additional arguments forwarded to the underlying static or
#'   interactive plotting helpers (e.g., \code{line_size}, \code{ribbon_alpha},
#'   \code{show_points}, \code{label_points}).
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

  if (inherits(x, "gcomp_out"))    contrast_list <- list(default = x)
  else                              contrast_list <- x

  if (interactive || !is.null(save_file)) {
    p <- plot_posterior_causal_contrast_interactive(
      contrast_list = contrast_list,
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
    p <- plot_posterior_causal_contrast_static(
      contrast_list = contrast_list,
      s_vec         = s_vec,
      ref_line      = ref_line,
      theme_fn      = theme_fn,
      ...
    )
  }
  print(p); invisible(p)
}
