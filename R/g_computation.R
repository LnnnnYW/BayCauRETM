#' Bayesian g-Computation for Recurrent-Event Rate Contrasts
#'
#' Perform Bayesian g-computation to estimate average recurrent-event rates
#' under different treatment initiation times versus never treating.
#'
#' @param fit_out Output list from \code{fit_causal_recur()}, containing
#'   at least \code{stan_fit}, \code{data_preprocessed}, \code{n_pat}, \code{K}.
#' @param s_vec Integer vector of treatment initiation intervals (from 1 to K).
#' @param B Integer. Number of Monte Carlo replicates per posterior draw (default: 50).
#'
#' @return A list with:
#'   \describe{
#'     \item{R_mat}{Matrix of dimension M × (length(s_vec)+1), where M is the number of posterior draws.
#'     Columns are mean rates under each treatment strategy.}
#'     \item{delta}{Named list of length \code{length(s_vec)}, each an element list with
#'     \code{draws}, \code{mean}, \code{CI_lower}, \code{CI_upper}.}
#'   }
#'
#' @details
#' For each posterior draw:
#' \enumerate{
#'   \item Sample Dirichlet weights for subjects.
#'   \item Simulate B replicate paths of (T_k, Y_k) conditional on treatment at s.
#'   \item Compute subject-level event rates and average with weights.
#'   \item Compute contrasts \eqn{\Delta(s,K+1)=R(s)-R(K+1)}.
#' }
#'
#' @examples
#' \dontrun{
#' # Minimal reproducible workflow
#' df <- data.frame(
#'   patient_id = rep(1:2, each = 2),
#'   k_idx      = rep(1:2, 2),
#'   Y          = rpois(4, 1),
#'   T          = rbinom(4, 1, 0.2),
#'   C          = rbinom(4, 1, 0.05),
#'   A          = rbinom(4, 1, 0.5)
#' )
#' prior <- list(
#'   eta_beta    = 0, sigma_beta  = 1, rho_beta   = 0.5,
#'   eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = 0.5
#' )
#' fit <- fit_causal_recur(
#'   data       = df,
#'   K          = 2,
#'   id_col     = "patient_id",
#'   k_col      = "k_idx",
#'   y_col      = "Y",
#'   t_col      = "T",
#'   c_col      = "C",
#'   a_col      = "A",
#'   x_cols     = NULL,
#'   formula_T  = T_obs ~ Y_prev + A + k_idx,
#'   formula_Y  = Y_obs ~ Y_prev + A + k_idx,
#'   prior      = prior,
#'   num_chains = 1,
#'   iter       = 100
#' )
#' gcomp <- g_computation(
#'   fit_out = fit,
#'   s_vec   = 1:2,
#'   B       = 10
#' )
#' print(gcomp)
#' }
#'
#' @importFrom rstan extract
#' @importFrom stats rgamma rbinom rpois weighted.mean setNames
#' @importFrom utils capture.output
#' @import methods
#' @import cmdstanr
#' @import coda
#' @import survival
#' @importFrom stats lag na.omit
#' @export




g_computation <- function(fit_out, s_vec, B = 50) {
  if (is.null(fit_out$stan_fit)) {
    stop("stan_fit is missing in fit_out")
  }
  if (!inherits(fit_out$stan_fit, "stanfit")) {
    stop("stan_fit must be a stanfit object")
  }

  # Extract inputs
  stan_fit      <- fit_out$stan_fit
  df            <- fit_out$data_preprocessed
  n_pat         <- fit_out$n_pat
  K             <- fit_out$K
  design_info   <- fit_out$design_info

  # 1) Extract posterior draws
  post <- rstan::extract(
    stan_fit,
    pars = c(
      "beta0", "gamma0",
      "beta_X", "beta_Y", "beta_A",
      "gamma_X", "gamma_Y", "gamma_A"
    ),
    permuted = TRUE
  )
  M <- length(post$beta_Y)

  # 2) Dirichlet weights for each draw
  W      <- matrix(rgamma(M * n_pat, shape = 1, rate = 1), nrow = M, ncol = n_pat)
  pi_mat <- W / rowSums(W)

  # 3) Baseline covariate matrix (we use the same covariates for all k)
  mat_cov           <- model.matrix(design_info$formula_T, data = df)
  X_cov             <- mat_cov[, -1, drop = FALSE]  # drop intercept
  idx_baseline      <- which(df$k_idx == 1)
  X_cov_baseline    <- X_cov[idx_baseline, , drop = FALSE]

  # 4) Simulation function for one posterior draw and one treatment time s
  simulate_once <- function(m, s) {
    # unpack parameters for draw m
    beta0_m   <- post$beta0[m, ]
    beta_X_m  <- post$beta_X[m, ]
    beta_Y_m  <- post$beta_Y[m]
    beta_A_m  <- post$beta_A[m]

    gamma0_m  <- post$gamma0[m, ]
    gamma_X_m <- post$gamma_X[m, ]
    gamma_Y_m <- post$gamma_Y[m]
    gamma_A_m <- post$gamma_A[m]

    R_i <- numeric(n_pat)

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

          # define current treatment indicator
          A_k <- as.integer(k >= s)

          # 1) simulate death
          logit_h <- beta0_m[k] +
            sum(Xi_cov * beta_X_m) +
            beta_Y_m * Y_prev_sim +
            beta_A_m * A_k

          if (!is.finite(logit_h)) {
            p_death <- 0
          } else {
            p_death <- 1 / (1 + exp(-logit_h))
          }

          p_death <- pmin(pmax(p_death, 0), 1L)

          T_k_sim <- rbinom(1, 1, p_death)
          if (is.na(T_k_sim)) {
            T_k_sim <- 0L
          }

          if (T_k_sim == 1L) {
            T_prev_sim <- 1L
            next
          }


          # 2) simulate recurrent event
          num_risk    <- num_risk + 1L
          log_mu      <- gamma0_m[k] +
            sum(Xi_cov * gamma_X_m) +
            gamma_Y_m * Y_prev_sim +
            gamma_A_m * A_k
          mu_k        <- exp(log_mu)
          Y_k_sim     <- rpois(1, mu_k)
          Y_prev_sim  <- Y_k_sim
          num_Y       <- num_Y + Y_k_sim
        }

        rate_vec[b] <- if (num_risk > 0L) num_Y / num_risk else 0
      }

      R_i[i] <- mean(rate_vec)
    }

    R_i
  }

  # 5) Build R_mat: rows = posterior draws; cols = s_vec and never-treat (K+1)
  S             <- length(s_vec)
  R_mat         <- matrix(NA_real_, nrow = M, ncol = S + 1)
  colnames(R_mat) <- c(paste0("s=", s_vec), paste0("s=", K + 1))

  for (m in seq_len(M)) {
    # never treat
    R_never           <- simulate_once(m, K + 1)
    R_mat[m, S + 1]   <- weighted.mean(R_never, w = pi_mat[m, ])

    # each treatment-initiation time
    for (j in seq_len(S)) {
      R_s           <- simulate_once(m, s_vec[j])
      R_mat[m, j]   <- weighted.mean(R_s, w = pi_mat[m, ])
    }
  }

  # 6) Compute contrasts delta(s, K+1)
  delta <- setNames(vector("list", S), paste0("s=", s_vec))
  for (j in seq_len(S)) {
    d_vec <- R_mat[, j] - R_mat[, S + 1]
    delta[[j]] <- list(
      draws    = d_vec,
      mean     = mean(d_vec),
      CI_lower = as.numeric(quantile(d_vec, 0.025)),
      CI_upper = as.numeric(quantile(d_vec, 0.975))
    )
  }

  list(
    R_mat = R_mat,
    delta = delta
  )
}

