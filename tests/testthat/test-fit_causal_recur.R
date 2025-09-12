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

#expected input
test_that("fit_causal_recur expected input", {
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
    iter      = 500,
    verbose   = TRUE
  )
  expect_s3_class(fit, "causal_recur_fit")
})

#wrong input type
test_that("fit_causal_recur wrong type of input", {
  expect_error(fit_causal_recur(
    data      = as.matrix(df),
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "data must be a data.frame")
  expect_error(fit_causal_recur(
    data      = df,
    K         = "3",
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "K must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = c(3,4),
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "K must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = 0,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "K must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = 3.5,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "K must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = 3,
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "id_col must be a single character string")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = 3,
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "time_col must be a single character string")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = 3,
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "treat_col must be a single character string")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = 3,
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "lag_col must be NULL or a single character string")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = "Tk ~ Ak + I(lagYk^2) + L.1 + L.2",
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "formula_T must be a formula")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = "Yk ~ Ak + I(lagYk^2) + L.1 + L.2",
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "formula_Y must be a formula")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = "1",
    iter      = 500,
    verbose   = TRUE
  ), "cores must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = c(1,2),
    iter      = 500,
    verbose   = TRUE
  ), "cores must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 0,
    iter      = 500,
    verbose   = TRUE
  ), "cores must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1.5,
    iter      = 500,
    verbose   = TRUE
  ), "cores must be a positive integer")
  expect_error(fit_causal_recur(
    data      = df,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = "500",
    verbose   = TRUE
  ), "iter must be a positive integer")
})


#columns missing
test_that("fit_causal_recur missing required columns", {
  df2 <- df %>% select(-Tk)
  expect_error(fit_causal_recur(
    data      = df2,
    K         = K,
    id_col    = "id",
    time_col  = "k_fac",
    treat_col = "Ak",
    lag_col   = "lagYk",
    formula_T = Tk ~ Ak + I(lagYk^2) + L.1 + L.2,
    formula_Y = Yk ~ Ak + I(lagYk^2) + L.1 + L.2,
    cores     = 1,
    iter      = 500,
    verbose   = TRUE
  ), "Event column 'Tk' not found")
})

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
  iter      = 500,
  verbose   = TRUE
)

#test summary
test_that("summary.causal_recur_fit works", {
  summ <- summary(fit)
  expect_type(summ, "list")
})
#not found parameter
test_that("summary.causal_recur_fit no parameter", {
  expect_error(summary(fit, pars_to_report = "L1"), "Error extracting summary for parameters L1")
})

#test print
test_that("print.causal_recur_fit works", {
  expect_type(print(fit), "list")
})
#test plot
test_that("plot.causal_recur_fit works", {
  expect_null(plot(fit))
})


