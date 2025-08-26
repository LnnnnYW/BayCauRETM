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

  `%||%` <- function(x, y) if (is.null(x)) y else x

  df_all   <- fit_out$data_preprocessed
  n_pat    <- fit_out$n_pat
  K        <- fit_out$K
  info     <- fit_out$design_info
  labels   <- fit_out$param_labels
  lag_col  <- info$lag_col %||% "lagYk"

  if (is.null(Lmat)) {
    id_col   <- info$id_col %||% "pat_id"
    time_col <- info$time_col %||% "k_idx"

    baseline_df <- df_all |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::slice_min(order_by = .data[[time_col]], n = 1, with_ties = FALSE) |>
      dplyr::ungroup()

    build_mm <- function(terms, d, env) {
      if (!length(terms)) return(matrix(0, nrow(d), 0))
      f <- stats::reformulate(terms, intercept = FALSE); environment(f) <- env
      mm <- stats::model.matrix(f, d, na.action = stats::na.pass)
      mm[is.na(mm)] <- 0
      storage.mode(mm) <- "double"
      mm
    }
    env_T <- attr(stats::terms(info$formula_T), ".Environment")
    env_Y <- attr(stats::terms(info$formula_Y), ".Environment")

    cov_terms <- labels$cov_terms %||% character(0)
    mm_base   <- build_mm(cov_terms, baseline_df, env_T %||% env_Y)

    cov_cols  <- labels$cov_cols %||% colnames(mm_base) %||% character(0)
    if (length(cov_cols)) {
      miss <- setdiff(cov_cols, colnames(mm_base))
      if (length(miss))
        mm_base <- cbind(mm_base, matrix(0, nrow(mm_base), length(miss), dimnames = list(NULL, miss)))
      Lmat <- mm_base[, cov_cols, drop = FALSE]
    } else {
      Lmat <- matrix(0, nrow(baseline_df), 0)
    }
  }
  P_data <- ncol(Lmat)

  post <- rstan::extract(
    fit_out$stan_fit,
    pars = c("beta0","beta1","betaL","theta0","theta1","thetaL","thetaLag"),
    permuted = TRUE
  )
  ndraws <- length(post$beta1)
  P_stan <- if (is.null(dim(post$betaL))) 0L else dim(post$betaL)[2]

  take_betaL <- function(m) {
    if (P_stan == 0L || P_data == 0L) numeric(0) else as.numeric(post$betaL[m, 1:min(P_data, P_stan), drop = TRUE])
  }
  take_thetaL <- function(m) {
    if (P_stan == 0L || P_data == 0L) numeric(0) else as.numeric(post$thetaL[m, 1:min(P_data, P_stan), drop = TRUE])
  }

  lag_terms_used <- labels$Y_lag %||% character(0)
  env_Y <- attr(stats::terms(info$formula_Y), ".Environment")
  build_lag_mm_from_y <- function(y_prev) {
    n <- length(y_prev)
    if (!length(lag_terms_used)) return(matrix(0, n, 0))
    newdat <- setNames(data.frame(y_prev), lag_col)
    f <- stats::reformulate(unique(gsub("^I\\((.*)\\)$", "\\1", lag_terms_used)), intercept = FALSE)
    environment(f) <- env_Y
    mm_full <- stats::model.matrix(f, newdat, na.action = stats::na.pass)

    miss <- setdiff(lag_terms_used, colnames(mm_full))
    mm   <- if (length(miss)) cbind(mm_full, matrix(0, n, length(miss), dimnames = list(NULL, miss))) else mm_full
    mm[, lag_terms_used, drop = FALSE]
  }

  sim_interv <- function(L_TY, n,
                         beta0, beta1, betaL,
                         theta0, theta1, thetaL, thetaLag_vec,
                         s, K, B) {
    eta_beta  <- if (length(betaL))  as.numeric(L_TY %*% betaL)  else rep(0, n)
    eta_theta <- if (length(thetaL)) as.numeric(L_TY %*% thetaL) else rep(0, n)

    abar_vec <- as.integer(seq_len(K) >= s)

    ybar <- array(0, dim = c(n, K, B))
    tbar <- array(0, dim = c(n, K, B))

    # k = 1
    lambda1 <- exp(theta0[1] + abar_vec[1] * theta1 + eta_theta)
    ybar[, 1, ] <- matrix(stats::rpois(n * B, lambda1), n, B)

    # k = 2...K
    for (k in 2:K) {
      abar_k  <- abar_vec[k]
      p_death <- plogis(beta0[k] + abar_k * beta1 + eta_beta)
      tbar[, k, ] <- matrix(stats::rbinom(n * B, 1, p_death), n, B)

      lag_term <- matrix(0, n, B)
      if (length(thetaLag_vec)) {
        for (b in 1:B) {
          F_b <- build_lag_mm_from_y(ybar[, k - 1, b])
          if (ncol(F_b)) lag_term[, b] <- F_b %*% thetaLag_vec[seq_len(ncol(F_b))]
        }
      }
      mu_k <- theta0[k] + abar_k * theta1 + eta_theta + lag_term
      mu_k <- pmin(mu_k, 700)
      ybar[, k, ] <- matrix(stats::rpois(n * B, exp(mu_k)), n, B)
    }

    alive <- array(1L, dim = c(n, K, B))
    if (K >= 2) for (k in 2:K) alive[, k, ] <- (alive[, k - 1, ] == 1L) & (tbar[, k - 1, ] == 0L)

    num_alive <- apply(ybar * alive, c(1, 3), sum)
    den_alive <- apply(alive,        c(1, 3), sum)
    ir <- ifelse(den_alive > 0, num_alive / den_alive, NA_real_)
    rowMeans(ir, na.rm = TRUE)
  }

  # Dirichlet weighting
  W      <- matrix(stats::rgamma(ndraws * n_pat, 1, 1), ndraws, n_pat)
  pi_mat <- W / rowSums(W)

  cl <- NULL
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("sim_interv", "build_lag_mm_from_y", "Lmat",
                  "n_pat", "K", "B", "post"),
      envir = environment()
    )
  }
  par_apply <- if (is.null(cl)) lapply else function(X, FUN) parallel::parLapply(cl, X, FUN)

  s_all <- c(s_vec, K + 1)  # ref = K+1
  R_draws <- sapply(s_all, function(s) {
    out_mat <- do.call(rbind, par_apply(seq_len(ndraws), function(m) {
      sim_interv(
        L_TY = Lmat, n = n_pat,
        beta0  = post$beta0[m, ],
        beta1  = post$beta1[m],
        betaL  = take_betaL(m),
        theta0 = post$theta0[m, ],
        theta1 = post$theta1[m],
        thetaL = take_thetaL(m),
        thetaLag_vec = {
          th <- post$thetaLag[m, , drop = FALSE]
          as.numeric(th)
        },
        s = s, K = K, B = B
      )
    }))
    rowSums(out_mat * pi_mat)
  })

  colnames(R_draws) <- paste0("s=", s_all)

  # delta(s, K+1)
  ref_col <- ncol(R_draws)
  delta <- lapply(seq_along(s_vec), function(j) {
    d  <- R_draws[, j] - R_draws[, ref_col]
    qs <- stats::quantile(d, c(.025, .975), na.rm = TRUE)
    list(draws = d,
         mean     = mean(d, na.rm = TRUE),
         CI_lower = qs[1],
         CI_upper = qs[2])
  })
  names(delta) <- paste0("s=", s_vec)

  structure(list(R_mat = R_draws, delta = delta), class = "gcomp_out")
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
