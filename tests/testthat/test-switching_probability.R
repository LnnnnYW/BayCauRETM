# tests/testthat/test-switching_probability.R

library(testthat)
library(BayCauRETM)

test_that("switching_probability_diagnostics returns summary, plot, and model", {
  raw_df <- data.frame(
    patient_id = rep(1:2, each = 3),
    k_idx      = rep(1:3, times = 2),
    Y          = rep(0:2, 2),
    T          = rep(c(0,1,0), 2),
    C          = rep(0, 6),
    A          = rep(c(0,1,0), 2),
    X1         = seq_len(6)
  )
  prep <- preprocess_data(
    df      = raw_df,
    id_col  = "patient_id",
    k_col   = "k_idx",
    y_col   = "Y",
    t_col   = "T",
    c_col   = "C",
    a_col   = "A",
    x_cols  = "X1",
    K       = 3
  )
  df <- prep$processed_df

  expect_silent({
    res <- suppressMessages(
      suppressWarnings(
        switching_probability_diagnostics(df, covariates = NULL)
      )
    )
  })

  expect_named(res, c("summary", "plot", "model"))
  expect_s3_class(res$summary, "data.frame")
  expect_s3_class(res$plot,    "ggplot")
})
