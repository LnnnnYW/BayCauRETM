// Model：
// If in the kth interval, T_k-1=0 & C_k-1=0: T_obs ~ Bernoulli_logit( beta0[k] + X_T·beta_X + beta_Y*Y_prev + beta_A*A )
// If T_obs == 0，then Y_obs ~ Poisson_log( gamma0[k] + X_Y·gamma_X + gamma_Y*Y_prev + gamma_A*A )
//    where beta0[1:K], gamma0[1:K] follows gAR(1) smoothing prior
//    X_T, X_Y: covariates matrix


data {
  int<lower=0> NY1;
  int<lower=0> NYk;
  int<lower=0> NTk;
  int<lower=1> K;
  int<lower=1> P;

  int<lower=1,upper=K>         kvecT[NTk];
  matrix[NTk,P]                L_Tk;
  vector<lower=0,upper=1>[NTk] A_Tk;
  int<lower=0,upper=1>         Tk[NTk];

  int<lower=1,upper=K>         kvecY[NYk];
  matrix[NYk,P]                L_Yk;
  vector<lower=0,upper=1>[NYk] lagYk;
  vector<lower=0,upper=1>[NYk] A_Yk;
  int<lower=0>                 Yk[NYk];

  matrix[NY1,P]                L_Y1;
  vector<lower=0,upper=1>[NY1] A_Y1;
  int<lower=0>                 Y1[NY1];
}

parameters {
  real                beta1;
  vector[K]           beta_eps;
  real                beta0_star;
  vector[P]           betaL;
  vector<lower=0>[K]  sigma_beta;
  real<lower=0, upper=1> rho_beta_star;

  real                theta1;
  vector[K]           theta_eps;
  real                theta0_star;
  vector[P]           thetaL;
  real                theta_lag;
  vector<lower=0>[K]  sigma_theta;
  real<lower=0, upper=1> rho_theta_star;
}

transformed parameters {
  real<lower=-1, upper=1> rho_beta  = 2 * (rho_beta_star  - 0.5);
  real<lower=-1, upper=1> rho_theta = 2 * (rho_theta_star - 0.5);

  vector[K] beta0;
  vector[K] theta0;

  beta0[1]  = beta0_star  + sigma_beta[1]  * beta_eps[1];
  theta0[1] = theta0_star + sigma_theta[1] * theta_eps[1];

  for (k in 2:K) {
    beta0[k]  = beta0_star  * (1 - rho_beta)  + rho_beta  * beta0[k-1]
                + sigma_beta[k]  * beta_eps[k];
    theta0[k] = theta0_star * (1 - rho_theta) + rho_theta * theta0[k-1]
                + sigma_theta[k] * theta_eps[k];
  }
}

model {
  beta_eps  ~ normal(0, 1);
  theta_eps ~ normal(0, 1);

  beta0_star  ~ normal(0, 1);
  theta0_star ~ normal(0, 1);

  rho_beta_star  ~ beta(1, 1);
  rho_theta_star ~ beta(1, 1);

  beta1  ~ normal(0, 1);
  theta1 ~ normal(0, 1);

  betaL  ~ normal(0, 1);
  thetaL ~ normal(0, 1);
  theta_lag ~ normal(0, 1);

  // ---- Vectorized likelihood ----------------------------------------------
  vector[NTk] eta_T;
  for (i in 1:NTk)
    eta_T[i] = beta0[kvecT[i]] + A_Tk[i] * beta1 + dot_product(L_Tk[i], betaL);
  Tk ~ bernoulli_logit(eta_T);

  if (NY1 > 0) {
    vector[NY1] eta_Y1;
    for (i in 1:NY1)
      eta_Y1[i] = theta0[1] + A_Y1[i] * theta1 + dot_product(L_Y1[i], thetaL);
    Y1 ~ poisson_log(eta_Y1);
  }

  if (NYk > 0) {
    vector[NYk] eta_Yk;
    for (i in 1:NYk)
      eta_Yk[i] = theta0[kvecY[i]] + A_Yk[i] * theta1 +
                  dot_product(L_Yk[i], thetaL) + theta_lag * lagYk[i];
    Yk ~ poisson_log(eta_Yk);
  }
}
