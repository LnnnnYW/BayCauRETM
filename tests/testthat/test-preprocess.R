# tests/testthat/test-preprocess.R

library(testthat)
library(BayCauRETM)

test_that("preprocess_data works with minimal data and x_cols=NULL", {

  df <- data.frame(
    patient_id = rep(1:2, each = 3),
    k_idx      = rep(1:3, 2),
    Y          = sample(0:2, 6, TRUE),
    T          = rbinom(6, 1, 0.2),
    C          = rbinom(6, 1, 0.1),
    A          = rbinom(6, 1, 0.5)
  )

  ## 抑制可能出现的提示性 warning（例如“无变异”或“缺失区间”）
  expect_silent(
    out <- suppressWarnings(
      preprocess_data(
        df      = df,
        id_col  = "patient_id",
        k_col   = "k_idx",
        y_col   = "Y",
        t_col   = "T",
        c_col   = "C",
        a_col   = "A",
        x_cols  = NULL,
        K       = 3
      )
    )
  )

  expect_s3_class(out, "list")
  expect_named(out, c("processed_df", "n_pat"))
  expect_equal(nrow(out$processed_df), 6)
  expect_equal(out$n_pat, 2)
  expect_true(all(c("Y_prev", "T_prev", "C_prev") %in% names(out$processed_df)))
})

