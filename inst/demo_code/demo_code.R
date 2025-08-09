# inst/demo/demo_code.R

library(BayCauRETM)
library(dplyr)
library(tidyr)
library(rstan)
library(parallel)

set.seed(432)

use_multicore <- TRUE   # FALSE for single core
if (use_multicore) {
  cores <- max(1, detectCores() - 2)
} else {
  cores <- 1
}
options(mc.cores = cores)
rstan_options(auto_write = TRUE)

warmup <- 500
M      <- 500
iter   <- warmup + M
B      <- 50
s_vec  <- c(3, 6, 9)

df_fit <- df %>%
  filter(id %in% 1:100) %>%
  arrange(id, k) %>%
  mutate(k_fac = as.integer(factor(k, levels = sort(unique(k))))) %>%
  group_by(id) %>%
  mutate(
    lagYk = if ("lagYk" %in% names(.)) replace_na(lagYk, 0) else lag(Yk, default = 0)
  ) %>%
  ungroup() %>%
  drop_na(Tk, Yk, Ak, L.1, L.2) %>%
  mutate(
    L.1 = as.numeric(scale(L.1)),
    L.2 = as.numeric(scale(L.2))
  )
K <- length(unique(df_fit$k_fac))

# 1) Stan fit
fit <- fit_causal_recur(
  data      = df_fit,
  K         = K,
  id_col    = "id",
  time_col  = "k_fac",
  treat_col = "Ak",
  lag_col   = "lagYk",
  formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
  formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
  cores     = cores,
  verbose   = TRUE
)

# 2) MCMC Diagnosis
message("Checking convergence...")
rstan::check_hmc_diagnostics(fit$stan_fit)
diag <- mcmc_diagnosis(fit, pars_to_check = c("beta0","beta1","theta0","theta1","theta_lag"))
plot(diag)

baseline_df <- fit$data_preprocessed %>%
  group_by(pat_id) %>%
  slice_min(order_by = k_idx, n = 1) %>%
  arrange(pat_id) %>%
  ungroup()

# set it into gcomp
rhs_terms <- stats::delete.response(stats::terms(fit$design_info$formula_T))
Lmat <- model.matrix(rhs_terms, data = baseline_df)
Lmat <- Lmat[, colnames(Lmat) != "(Intercept)", drop = FALSE]

# 4) g-computation
message("Running g-computation...")
gcomp <- g_computation(
  Lmat, fit,
  s_vec = s_vec,
  B     = B,
  cores = cores
)
print(gcomp)
plot(gcomp, ref_line = 0)
plot(gcomp, interactive = TRUE, ref_line = 0)

# 5) PS Diagnosis
ps_diag <- propensity_score_diagnostics(
  fit$data_preprocessed,
  treat_col = "A",
  covariates = c("lagYk","k_idx")
)
plot(ps_diag, type = "histogram")
plot(ps_diag, type = "density")

# 6) Switching Probability
sw_diag <- switching_probability_summary(fit$data_preprocessed)
plot(sw_diag, type = "boxplot")

# 7) Results
sum_tbl <- result_summary_table(
  fit_out   = fit,
  gcomp_out = gcomp,
  s_vec     = s_vec,
  format    = "kable",
  pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag")
)
print(sum_tbl)

