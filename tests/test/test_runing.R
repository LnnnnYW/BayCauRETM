#' simulate_data()
#'
#' 按照论文补充材料 S3 描述的流程，生成示例的 long-format 数据集
#' （N 个患者，K 个区间），包含以下列：
#'   patient_id, k_idx, Y, T, C, A, L1, L2, L3, L4, L5
#'
#' 具体流程：
#'   1) 生成 N 名患者的基线协变量 L1 至 L5，其中：
#'        L1, L2, L3 ~ N(0, 1)；L4, L5 ~ Bernoulli(0.5)
#'
#'   2) 在 k = 1 时：
#'        T1 = 0；C1 = 0
#'        A1 ~ Bernoulli( expit(-1 + sum(L * gamma_A)) )
#'        Y1 ~ Poisson( exp(theta0 + A1 + sum(L * theta_L)) )
#'
#'   3) 在 k = 2 时：
#'        如果 A1 == 1，则 A2 = 1
#'        否则 A2 ~ Bernoulli( expit(-1 + sum(L * gamma_A)) )
#'
#'        T2 ~ Bernoulli( expit(-beta0 + sum(L * beta_L) + A2) )
#'        C2 ~ Bernoulli( expit(censor0 + sum(L * theta_L) + gamma_lag * I(A1 > 0)) )
#'
#'        若 T2 == 1 或 C2 == 1，则该患者后续不再模拟；
#'        否则 Y2 ~ Poisson( exp(theta0 + A2 + sum(L * theta_L) + theta_lag * I(A1 > 0)) )
#'
#'   4) 对于 k = 3 到 K，重复同样逻辑：
#'        如果 A_{k-1} == 1，则 A_k = 1
#'        否则 A_k ~ Bernoulli( expit(-1 + sum(L * gamma_A)) )
#'
#'        T_k ~ Bernoulli( expit(-beta0 + sum(L * beta_L) + A_k) )
#'        C_k ~ Bernoulli( expit(censor0 + sum(L * theta_L) + gamma_lag * I(A_{k-1} > 0)) )
#'
#'        若 T_k == 1 或 C_k == 1，则停止模拟；
#'        否则 Y_k ~ Poisson( exp(theta0 + A_k + sum(L * theta_L) + theta_lag * I(A_{k-1} > 0)) )
#'
#' 参数：
#'   N    : 患者数量（默认 500）
#'   K    : 区间数量（默认 12）
#'   seed : 随机种子（默认 12345）
#'
#' 返回值：
#'   一个长格式数据框 df_long，包含以下列：
#'   patient_id, k_idx, Y, T, C, A, L1, L2, L3, L4, L5



simulate_data <- function(N = 500, K = 12, seed = 42) {
  set.seed(seed)

  # baseline covariates L1~L5
  L1 <- rnorm(N, mean = 0, sd = 1)
  L2 <- rnorm(N, mean = 0, sd = 1)
  L3 <- rnorm(N, mean = 0, sd = 1)
  L4 <- rbinom(N, size = 1, prob = 0.5)
  L5 <- rbinom(N, size = 1, prob = 0.5)
  L_mat <- cbind(L1, L2, L3, L4, L5)

  # set params in model
  gamma_A_vec <- c(-0.5, 0.5, -0.2, 0.2, -0.2)  # A_k 的 L 线性效应
  theta_L_vec <- c(-1, -1, -1, 1, -0.5)         # Y_k 的 L 线性效应
  theta_lag    <- 1.5                           # 上一期 gamma>0 对 Y_k 的影响
  gamma_lag    <- 1.5                           # 上一期 gamma>0 对 C_k 的影响
  beta_L_vec   <- c(-1, -1, -1, 0.5, -0.5)      # T_k 的 L 线性效应

  theta0_vec <- rep(0, K)   # Y_k intercept
  beta0      <- 0           # T_k intercept
  censor0    <- -2.5        # C_k intercept

  # initialization
  records <- vector("list", N * K)
  idx     <- 1

  # simulate for every patient i
  for (i in 1:N) {
    Li      <- L_mat[i, ]
    A_prev  <- 0
    T_prev  <- 0
    C_prev  <- 0
    Y_prev  <- 0

    for (k in 1:K) {
      if (T_prev == 1 || C_prev == 1) {
        # if die or censoring, break here
        break
      }

      # A_k
      if (k == 1) {
        # k=1: A1 ~ Bernoulli(expit(-1 + sum L*gamma_A))
        p_A <- 1 / (1 + exp(-(-1 + sum(Li * gamma_A_vec))))
        A_k <- rbinom(1, size = 1, prob = p_A)
      } else {
        # k>=2: if A_prev==1，then A_k=1；or else same as k=1
        if (A_prev == 1) {
          A_k <- 1
        } else {
          p_A <- 1 / (1 + exp(-(-1 + sum(Li * gamma_A_vec))))
          A_k <- rbinom(1, size = 1, prob = p_A)
        }
      }

      # T_k
      p_T <- 1 / (1 + exp(-(-beta0 + sum(Li * beta_L_vec) + A_k)))
      T_k <- rbinom(1, size = 1, prob = p_T)

      # C_k
      p_C <- 1 / (1 + exp(-(censor0 + sum(Li * theta_L_vec) + gamma_lag * as.numeric(Y_prev > 0))))
      C_k <- rbinom(1, size = 1, prob = p_C)

      # Simulate Y_k
      if (T_k == 1) {
        Y_k <- 0
      } else if (C_k == 1) {
        # if censoring，Y_k set to be NA
        Y_k <- NA
      } else {
        #  Y_k ~ Poisson(exp(θ0 + A_k + sum L*theta_L + thata_lag I(Y_prev>0)))
        log_mu <- theta0_vec[k] + A_k + sum(Li * theta_L_vec) + theta_lag * as.numeric(Y_prev > 0)
        mu_k   <- exp(log_mu)
        Y_k    <- rpois(1, lambda = mu_k)
      }

      records[[idx]] <- data.frame(
        patient_id = i,
        k_idx      = k,
        Y          = ifelse(is.na(Y_k), NA_integer_, as.integer(Y_k)),
        T          = as.integer(T_k),
        C          = as.integer(C_k),
        A          = as.integer(A_k),
        L1         = Li[1],
        L2         = Li[2],
        L3         = Li[3],
        L4         = Li[4],
        L5         = Li[5],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1

      A_prev <- A_k
      T_prev <- T_k
      C_prev <- C_k
      Y_prev <- ifelse(is.na(Y_k), 0, Y_k)

      if (T_k == 1 || C_k == 1) {
        break
      }
    } # end for k
  } # end for i
  df_long <- do.call(rbind, records[!sapply(records, is.null)])
  return(df_long)
}








df_long <- simulate_data(N = 500, K = 12, seed = 42)

head(df_long)


formula_T <- as.formula("T_obs ~ Y_prev + A + L1 + L2 + L3 + L4 + L5 + k_idx")
formula_Y <- as.formula("Y_obs ~ Y_prev + A + L1 + L2 + L3 + L4 + L5 + k_idx")

prior_list <- list(
  eta_beta    = 0,
  sigma_beta  = 1,
  rho_beta    = 0.9,
  eta_gamma   = 0,
  sigma_gamma = 1,
  rho_gamma   = 0.9
)

# 3 拟合模型
fit_out <- fit_causal_recur(
  data            = df_long,
  K               = 12,
  id_col          = "patient_id",
  k_col           = "k_idx",
  y_col           = "Y",
  t_col           = "T",
  c_col           = "C",
  a_col           = "A",
  x_cols          = c("L1", "L2", "L3", "L4", "L5"),
  formula_T       = formula_T,
  formula_Y       = formula_Y,
  prior           = prior_list,
  num_chains      = 4,
  iter            = 2000,
  stan_model_file = "inst/stan/recur_model.stan"
)

# 4 MCMC 收敛诊断
mcmc_diagnosis(fit_out, save_plots = FALSE)

# 5 贝叶斯 g-computation：计算 s = 1..12 vs s = 13（永不转药）
gcomp_out <- g_computation(fit_out, s_vec = 1:12, B = 50)

# 6 汇总结果
res_summary <- summarize_results(fit_out, gcomp_out)
res_summary$delta_summary

# 7 绘制 delta(s,13) 曲线
plot_delta_vs_s(gcomp_out, s_vec = 1:12)

