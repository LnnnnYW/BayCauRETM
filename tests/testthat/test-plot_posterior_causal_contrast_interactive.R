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

gcomp <- g_computation(fit_out = fit,
                       s_vec   = s_vec,
                       B       = B,
                       cores   = cores)

#test expected input
test_that("plot_posterior_causal_contrast_interactive works", {
  p <- plot_posterior_causal_contrast_interactive(gcomp, s_vec = 3)
  expect_s3_class(p, "ggplot")
})

#test wrong type of input
test_that("plot_posterior_causal_contrast_interactive wrong type of input", {
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, s_vec = -3), "s_vec must be a vector of positive integers")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, theme_fn = "a"), "theme_fn must be a ggplot2 theme function")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, point_size = -1), "point_size must be a single positive number")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, error_width = 2), "error_width must be a single number in \\(0, 1\\]")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, ref_line = "a"), "ref_line must be NULL or a single number")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, interactive = "a"), "interactive must be a single logical value")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, save_file = 1), "save_file must be NULL or a character string")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, width = -1), "width must be a single positive number")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, height = -1), "height must be a single positive number")
  expect_error(plot_posterior_causal_contrast_interactive(gcomp, dpi = -1), "dpi must be a single positive number")
})


