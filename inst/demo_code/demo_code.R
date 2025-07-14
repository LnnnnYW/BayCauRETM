# inst/demo/demo_code.R

library(BayCauRETM)
library(dplyr)
library(future)
library(future.apply)
plan(multisession)

set.seed(42)

# Generate a long-format data with non-standard names
N       <- 150   # subjects
K       <- 4     # intervals
beta_A  <- -0.5  # effect on death (logit)
gamma_A <- -0.4  # effect on events (log)

raw_df <- do.call(rbind, lapply(seq_len(N), function(i) {
  prev_events <- 0
  rows <- list()
  start_trt <- sample(1:(K + 1), 1)

  for (t in 1:K) {
    trt_flag <- as.integer(t >= start_trt)
    # death
    logit_p <- -3 + 0.4 * prev_events + beta_A * trt_flag
    died    <- rbinom(1, 1, plogis(logit_p))
    # events
    mu      <- exp(log(1) + 0.3 * prev_events + gamma_A * trt_flag)
    cnt     <- rpois(1, mu)

    rows[[t]] <- data.frame(
      sid          = i,       # id
      period       = t,       # time
      death_flag   = died,    # death
      event_count  = cnt,     # event count
      trt_arm      = trt_flag # treatment
    )

    if (died == 1) break
    prev_events <- cnt
  }
  do.call(rbind, rows)
}))

# Introduce exactly the three canonical names we need to bypass A/k_idx/pat_id errors
demo_df <- raw_df %>%
  mutate(
    A     = trt_arm,
    k_idx = period,
    pat_id= sid
  )
# demo_df now has:
# [ sid, period, death_flag, event_count, trt_arm, A, k_idx, pat_id ]

# Fit the model and watch the auto-mapping in action
fit <- fit_causal_recur(
  data            = demo_df,
  K               = K,
  x_cols          = NULL,
  formula_T       = death_flag  ~ event_count + A + k_idx,
  formula_Y       = event_count ~ event_count + A + k_idx,
  prior           = list(
    eta_beta    = 0, sigma_beta  = 1, rho_beta   = .5,
    eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = .5
  ),
  stan_model_file = system.file("stan","causal_recur_model.stan",
                                package="BayCauRETM"),
  num_chains      = 2,
  iter            = 500,
  verbose         = FALSE        # prints the mapping steps
)

# MCMC diagnostics
diag <- mcmc_diagnosis(fit)
print(diag); plot(diag)

# g-computation
gcomp <- g_computation(fit, s_vec = 1:K, B = 50)
print(gcomp); plot(gcomp, ref_line = 0)

# Propensity-score diagnostics
ps_diag <- propensity_score_diagnostics(
  data       = fit$data_preprocessed,
  treat_col  = "A",
  covariates = c("Y_prev", "k_idx")
)
plot(ps_diag)

# Switching-probability diagnostics
sw_diag <- switching_probability_summary(fit$data_preprocessed)
plot(sw_diag)

# Summary table
sum_tbl <- result_summary_table(
  fit_out   = fit,
  gcomp_out = gcomp,
  s_vec     = 1:K,
  format    = "kable"
)
print(sum_tbl)
