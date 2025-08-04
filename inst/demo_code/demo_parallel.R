# inst/demo/demo_code.R

library(dplyr)
library(tidyr)
library(rstan)
library(future)
library(future.apply)
library(parallel)
library(doParallel)

cores = detectCores() - 2

warmup <- 500
M      <- 500
B      <- 100
s_vec  <- c(3, 6, 9)

load("data.Rdata")

df_clean <- df %>%
  dplyr::filter(id %in% 1:100) %>%
  dplyr::filter(k > 1) %>%
  tidyr::drop_na(Tk, Yk, lagYk, Ak, L.1, L.2)

df_bsl = df %>%
  dplyr::filter(id %in% 1:100) %>%
  filter(k==1) %>%
  select(starts_with("L.")) %>%
  as.matrix()


K <- max(df_clean$k)

#LagYK can be nothing, a vector, column name or the formula likesummary "~variable 1 + variable 2"
message("Fitting model with cached compiled Stan...")
fit <- fit_causal_recur(
  data       = df_clean,
  K          = K,
  id_col     = "id",
  time_col   = "k",
  treat_col  = "Ak",
  x_cols     = c("L.1","L.2"),
  lagYK = NA,
  formula_T  = Tk ~ Ak + k + `L.1` + `L.2`,
  formula_Y  = Yk ~ Ak + k + `L.1` + `L.2`,
  prior      = list(eta_beta = 0, sigma_beta = 1, rho_beta = 0.5,
                    eta_gamma = 0, sigma_gamma = 1, rho_gamma = 0.5),
  num_chains = 4,
  iter       = 10,
  cores = cores,
  #compiled_model_cache = "compiled_model.rds",
  stan_model_file = "inst/stan/causal_recur_model.stan",
  verbose    = TRUE
)

message("Checking convergence...")
diag <- mcmc_diagnosis(fit, pars_to_check = c("beta0","beta1","theta0","theta1"))
print(diag)
plot(diag)

message("Running g-computation...")
gcomp <- g_computation(df_bsl, fit, s_vec = s_vec, B = B, cores = cores)
print(gcomp)
plot(gcomp, ref_line = 0)

ps_diag <- propensity_score_diagnostics(
  fit$data_preprocessed,
  treat_col = "A",
  covariates = c("lagYk","k_idx")
)
plot(ps_diag)

sw_diag <- switching_probability_summary(fit$data_preprocessed)
plot(sw_diag)

sum_tbl <- result_summary_table(
  fit_out   = fit,
  gcomp_out = gcomp,
  s_vec     = s_vec,
  format    = "kable",
  pars_to_report = c("beta0","beta1","theta0","theta1")
)
print(sum_tbl)
