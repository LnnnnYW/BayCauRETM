library(testthat)
library(tidyverse)

load("data.Rdata")

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
  iter      = 2000,
  verbose   = TRUE
)

B      <- 10
s_vec  <- c(3, 6, 9)
cores <- 1

gcomp <- g_computation(fit_out = fit,
                       s_vec   = s_vec,
                       B       = B,
                       cores   = cores)

# expected input
test_that("result_summary_table expected input", {
  result_summary <- result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL)
  expect_type(result_summary, "list")
})

# wrong type of input
test_that("result_summary_table wrong type of input", {
  expect_error(result_summary_table(fit_out = as.matrix(fit),
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL), "fit_out must contain be a 'stanfit' object")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = as.matrix(gcomp),
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL), "gcomp_out must be the output of g_computation")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = NULL,
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL), "'pars_to_report' must be a non-empty character vector")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = c(3,-6,9),
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL), "'s_vec' must be NULL or a numeric vector of positive integers")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = 5,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL), "'filter_pars' must be NULL, a character vector, or a dplyr filter expression")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = 5,
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = NULL), "'sort_by' must be NULL or a character string")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = c(TRUE, FALSE),
                                         format        = "data.frame",
                                         export_file   = NULL), "'sort_desc' must be a single logical value")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = c("data.frame", "list"),
                                         export_file   = NULL), "'format' must be a single character string")
  expect_error(result_summary_table(fit_out = fit,
                                         gcomp_out = gcomp,
                                         pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                         s_vec         = NULL,
                                         filter_pars   = NULL,
                                         sort_by       = "Mean",
                                         sort_desc     = TRUE,
                                         format        = "data.frame",
                                         export_file   = c("file1", "file2")), "'export_file' must be NULL or a single character string")
})

