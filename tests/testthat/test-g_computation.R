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
  num_chains = 1,
  verbose   = TRUE
)

B      <- 10
s_vec  <- c(3, 6, 9)
cores <- 1

#expected input
test_that("g_computation expected input", {
  gcomp <- g_computation(fit_out = fit,
                         s_vec   = s_vec,
                         B       = B,
                         cores   = cores)
  expect_type(gcomp, "list")
  expect_type(print(gcomp), "list")
  p <- plot(gcomp)
  expect_s3_class(p, "ggplot")
  p <- plot(gcomp, interactive = TRUE)
  expect_s3_class(p, "plotly")
})
test_that("g_computation expected input: parallel", {
  gcomp <- g_computation(fit_out = fit,
                         s_vec   = s_vec,
                         B       = B,
                         cores   = 4)
  expect_type(gcomp, "list")
})

# wrong type of input
test_that("g_computation wrong type of input", {
  expect_error(g_computation(fit_out = as.matrix(fit),
                             s_vec   = s_vec,
                             B       = B,
                             cores   = cores), "fit_out must contain be a 'stanfit' object")
  expect_error(g_computation(fit_out = fit,
                             s_vec   = "3,6,9",
                             B       = B,
                             cores   = cores), "s_vec must be a vector of positive integers")
  expect_error(g_computation(fit_out = fit,
                             s_vec   = c(3,6,9),
                             B       = "10",
                             cores   = cores), "B must be a single positive integer")
  expect_error(g_computation(fit_out = fit,
                             s_vec   = c(3,6,9),
                             B       = c(10,20),
                             cores   = cores), "B must be a single positive integer")
  expect_error(g_computation(fit_out = fit,
                             s_vec   = c(3,6,9),
                             B       = 10,
                             cores   = "1"), "cores must be a single positive integer")
  expect_error(g_computation(fit_out = fit,
                             s_vec   = c(3,6,9),
                             B       = 10,
                             cores   = c(1,2)), "cores must be a single positive integer")
})

#test ci
test_that("g_computation: ci validation works", {
  gcomp <- g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                         ci = c(0.025, 0.975))
  expect_type(gcomp, "list")
  expect_true(all(c("mean", "CI_lower", "CI_upper") %in% names(gcomp$delta[[1]])))

  # Not numeric
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c("0.025", "0.975")),
    "'ci' must be a numeric length-2 vector"
  )

  # Wrong length
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = 0.95),
    "'ci' must be a numeric length-2 vector"
  )
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(0.025, 0.5, 0.975)),
    "'ci' must be a numeric length-2 vector"
  )

  # NA present
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(NA_real_, 0.975)),
    "'ci' must be a numeric length-2 vector"
  )

  # Out of (0,1)
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(0, 0.975)),
    "'ci' must be within \\(0,1\\)"
  )
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(0.025, 1)),
    "'ci' must be within \\(0,1\\)"
  )
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(-1, 0.975)),
    "'ci' must be within \\(0,1\\)"
  )

  # Bad order / equal
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(0.975, 0.025)),
    "ci\\[1\\] < ci\\[2\\]"
  )
  expect_error(
    g_computation(fit_out = fit, s_vec = s_vec, B = B, cores = cores,
                  ci = c(0.5, 0.5)),
    "ci\\[1\\] < ci\\[2\\]"
  )
})
