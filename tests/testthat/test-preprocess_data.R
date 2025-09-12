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
test_that("preprocess_data expected input", {
  expect_type(preprocess_data(df, K)$processed_df, "list")
  expect_type(preprocess_data(df, K)$n_pat, "integer")
})
test_that("preprocess_data expected input", {
  expect_type(preprocess_data(df, K, lag_col = "lag", need_lag = TRUE)$processed_df, "list")
  expect_type(preprocess_data(df, K, lag_col = "lag", need_lag = TRUE)$n_pat, "integer")
})

#columns missing
test_that("preprocess_data missing required columns", {
  df2 <- df %>% select(-T_obs)
  expect_error(preprocess_data(df2, K), "df is missing required columns: T_obs")
})

#wrong type of input
test_that("preprocess_data wrong type of input", {
  expect_error(preprocess_data(as.matrix(df), K), "is.data.frame")
  expect_error(preprocess_data(df, "3"), "is.numeric")
  expect_error(preprocess_data(df, c(3,4)), "length")
  expect_error(preprocess_data(df, 0), "K >= 1")
})

#no NA
test_that("preprocess_data conflict input", {
  df2 <- df %>% mutate(lagYk = ifelse(runif(n()) < 0.1, NA, lagYk))
  expect_error(preprocess_data(df2, K), "df must not contain any missing values")
})


