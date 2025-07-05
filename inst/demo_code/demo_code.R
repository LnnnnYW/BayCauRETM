# inst/demo/demo_code.R
library(BayCauRETM)
set.seed(42)

# 1) Data-generating
N <- 150          # subjects
K <- 4            # intervals
beta_A  <- -0.5   # A -> death   (logit)
gamma_A <- -0.4   # A -> events  (log)

sim_df <- do.call(
  rbind,
  lapply(seq_len(N), function(i) {

    Y_prev <- 0
    alive  <- TRUE
    rows   <- vector("list", K)

    ## subject-specific treatment start s  (K+1 == never treat)
    s_start <- sample(1:(K + 1), 1, prob = rep(1 / (K + 1), K + 1))

    for (k in 1:K) {

      A_k <- as.integer(k >= s_start)

      if (alive) {
        ## death
        logit_h <- -3 + 0.4 * Y_prev + beta_A * A_k
        T_k     <- rbinom(1, 1, plogis(logit_h))

        ## recurrent events
        mu_k <- exp(log(1) + 0.3 * Y_prev + gamma_A * A_k)
        Y_k  <- rpois(1, mu_k)
      } else {          # after death: structural values
        T_k <- 1
        Y_k <- 0        # use 0 (not NA) to satisfy preprocess_data checks
      }

      rows[[k]] <- data.frame(
        pat_id = i,
        k_idx  = k,
        Y_obs  = Y_k,
        T_obs  = T_k,
        A      = A_k
      )

      alive  <- (T_k == 0)
      Y_prev <- Y_k
    }

    do.call(rbind, rows)
  })
)

# 2) Model Fitting
stan_file <- system.file("stan", "causal_recur_model.stan",
                         package = "BayCauRETM")

prior_default <- list(
  eta_beta  = 0, sigma_beta  = 1, rho_beta   = .5,
  eta_gamma = 0, sigma_gamma = 1, rho_gamma  = .5
)

iter_demo   <- 1000
chains_demo <- 1

fit <- fit_causal_recur(
  data            = sim_df,
  K               = K,
  formula_T       = T_obs ~ Y_prev + A + k_idx,
  formula_Y       = Y_obs ~ Y_prev + A + k_idx,
  prior           = prior_default,
  stan_model_file = stan_file,
  num_chains      = chains_demo,
  iter            = iter_demo,
  verbose         = FALSE
)

# 3) MCMC Diagnosis
diag <- mcmc_diagnosis(fit)
print(diag)     # R-hat / ESS table
plot(diag)      # trace plots

# 4) g-computation
gcomp <- g_computation(fit, s_vec = 1:K, B = 50)
print(gcomp)
plot(gcomp, ref_line = 0)

if (requireNamespace("plotly", quietly = TRUE)) {
  plot(gcomp, interactive = TRUE)
}

# 5) Propensity-score diagnostics
ps_diag <- propensity_score_diagnostics(
  data       = fit$data_preprocessed,
  treat_col  = "A",
  covariates = c("Y_prev", "k_idx")
)
plot(ps_diag)

# 6) Switching-probability diagnostics
sw_diag <- switching_probability_summary(fit$data_preprocessed)
plot(sw_diag)

# 7) Summary table
sum_tbl <- result_summary_table(
  fit_out   = fit,
  gcomp_out = gcomp,
  s_vec     = 1:K,
  format    = "kable"
)
print(sum_tbl)

