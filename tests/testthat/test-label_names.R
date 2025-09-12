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
  iter      = 500,
  verbose   = TRUE
)

test_that("result_summary_table works", {
  expect_type(Summary_lable(fit), "list")
})
