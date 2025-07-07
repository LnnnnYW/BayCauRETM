// Model：
// If in the kth interval, T_k-1=0 & C_k-1=0: T_obs ~ Bernoulli_logit( beta0[k] + X_T·beta_X + beta_Y*Y_prev + beta_A*A )
// If T_obs == 0，then Y_obs ~ Poisson_log( gamma0[k] + X_Y·gamma_X + gamma_Y*Y_prev + gamma_A*A )
//    where beta0[1:K], gamma0[1:K] follows gAR(1) smoothing prior
//    X_T, X_Y: covariates matrix


// Discrete-time Bayesian model for recurrent events and terminal event
// - Terminal event hazard: Bernoulli-logit
// - Recurrent-event count: Poisson-log, conditional on survival and non-censoring
// - Time-varying intercepts follow generalized AR(1) priors

data {
  int<lower=1> N;           // number of rows
  int<lower=1> n_pat;       // number of subjects
  int<lower=1> K;           // number of intervals
  int<lower=1> p_T;         // covariate dimension for hazard
  int<lower=1> p_Y;         // covariate dimension for count

  int<lower=1,upper=n_pat> pat_id[N];
  int<lower=1,upper=K>   k_idx[N];

  int<lower=0,upper=1> T_prev[N];
  int<lower=0,upper=1> A[N];

  matrix[N,p_T] X_T;
  int<lower=0> Y_prev[N];

  matrix[N,p_Y] X_Y;
  int<lower=0,upper=1> T_obs[N];
  int<lower=0>         Y_obs[N];

  real eta_beta;
  real<lower=0> sigma_beta;
  real<lower=-1,upper=1> rho_beta;
  real eta_gamma;
  real<lower=0> sigma_gamma;
  real<lower=-1,upper=1> rho_gamma;
}

parameters {
  vector[K] beta0;
  vector[K] gamma0;
  vector[p_T] beta_X;
  real beta_Y;
  real beta_A;
  vector[p_Y] gamma_X;
  real gamma_Y;
  real gamma_A;
}

model {
  // gAR(1) priors on time-varying intercepts
  beta0[1] ~ normal(eta_beta, sigma_beta);
  for (k in 2:K)
    beta0[k] ~ normal(eta_beta + rho_beta * (beta0[k-1] - eta_beta), sigma_beta);

  gamma0[1] ~ normal(eta_gamma, sigma_gamma);
  for (k in 2:K)
    gamma0[k] ~ normal(eta_gamma + rho_gamma * (gamma0[k-1] - eta_gamma), sigma_gamma);

  // weakly informative priors
  beta_X ~ normal(0, 1);
  beta_Y ~ normal(0, 1);
  beta_A ~ normal(0, 1);
  gamma_X ~ normal(0, 1);
  gamma_Y ~ normal(0, 1);
  gamma_A ~ normal(0, 1);

  // likelihood
  for (n in 1:N) {
    real lp = beta0[k_idx[n]]
              + X_T[n] * beta_X
              + beta_Y * Y_prev[n]
              + beta_A * A[n];

    if (T_prev[n] == 0) {                 // subject still alive at start of interval
      T_obs[n] ~ bernoulli_logit(lp);

      if (T_obs[n] == 0) {                // if not dead, model recurrent events
        Y_obs[n] ~ poisson_log(
          gamma0[k_idx[n]]
          + X_Y[n] * gamma_X
          + gamma_Y * Y_prev[n]
          + gamma_A * A[n]);
      }
    }
  }
}

