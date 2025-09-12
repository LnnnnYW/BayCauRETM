library(testthat)
library(tidyverse)

df <- readRDS(testthat::test_path("data","data.rds"))

df_clean <- df %>%
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
K <- length(unique(df_clean$k_fac))

formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2
formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2

id_col    = "id"
time_col  = "k_fac"
treat_col = "Ak"
lag_col   = "lagYk"
event_col <- all.vars(formula_T)[1]   # Tk
count_col <- all.vars(formula_Y)[1]   # Yk

df <- as.data.frame(df_clean)
df$T_obs  <- df[[event_col]]
df$Y_obs  <- df[[count_col]]
df$pat_id <- df[[id_col]]
df$k_idx  <- df[[time_col]]
df$A      <- df[[treat_col]]

fit <- fit_causal_recur(
  data      = df,
  K         = K,
  id_col    = "id",
  time_col  = "k_fac",
  treat_col = "Ak",
  lag_col   = "lagYk",
  formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
  formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
  cores     = 1,
  iter      = 1000,
  num_chains = 1,
  verbose   = TRUE
)

#expected input
test_that("mcmc_diagnosis expected input", {
  diagnosis <- mcmc_diagnosis(fit)
  expect_type(diagnosis, "list")
  expect_type(print(diagnosis), "list")
})


#extra params in pars_to_check
test_that("mcmc_diagnosis no params", {
  expect_warning(mcmc_diagnosis(fit, pars_to_check = c("beta0", "theta0", "beta1", "theta1", "thetaLag", "theta2")),
                 "The following 'pars_to_check' were not found in the fitted Stan object and will be skipped: theta2")
})

#no params in pars_to_check
test_that("mcmc_diagnosis no params", {
  expect_error(mcmc_diagnosis(fit, pars_to_check = c()),
               "None of the specified 'pars_to_check' exist in the fitted Stan object.")
})

#wrong type of input
test_that("mcmc_diagnosis wrong type of input", {
  expect_error(mcmc_diagnosis(as.matrix(fit)), "input must contain a valid 'stan_fit' rstan::stanfit object")
  expect_error(mcmc_diagnosis(list(a=1, b=2)), "input must contain a valid 'stan_fit' rstan::stanfit object")
})


###test plot
diagnosis <- mcmc_diagnosis(fit)
test_that("mcmc_diagnosis plot works", {
  p <- plot(diagnosis)
  expect_type(p, "list")
})

###not found parameter
test_that("mcmc_diagnosis plot not found parameter", {
  expect_warning(plot(diagnosis, pars = c("L1")),
                 "Some 'pars' were not found in the available plots and will be skipped.")
})
