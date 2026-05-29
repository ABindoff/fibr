// Connection-corrected (horizontal) reparameterisation of the centred GLMM.
//
// The fiber variable is split into:
//   alpha_vert[j] = alpha[j] - A_mu[j]*mu - A_sigma[j]*sigma
//
// where A_mu[j] and A_sigma[j] are the posterior-mean connection coefficients
// computed from the Milestone 2 connection analysis.  This removes the
// dominant linear coupling between the fiber alpha and the base (mu, sigma),
// reducing the non-trivial holonomy and improving NUTS mixing.
//
// For A_mu = 0, A_sigma = 0: recovers the centred parameterisation.
// For A_mu = 1/sigma, A_sigma = 0: approximates the non-centred parameterisation.
// For A_mu, A_sigma from the Fisher metric connection: the horizontal parameterisation.

data {
  int<lower=0> N;
  int<lower=1> J;
  array[N] int<lower=1,upper=J> group;
  matrix[N, 2] X;
  array[N] int<lower=0,upper=1> y;
  // Connection coefficients from compute_connection() — fixed at posterior mean
  vector[J] A_mu;     // A[j, mu]   = 1 / (sigma_bar^2 * G_FF_bar[j])
  vector[J] A_sigma;  // A[j, sigma] = 2*(alpha_bar[j]-mu_bar) / (sigma_bar^3 * G_FF_bar[j])
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector[J] alpha_vert;   // vertical (connection-corrected) fiber component
  vector[2] beta;
}
transformed parameters {
  // Reconstruct alpha from vertical component + horizontal base tracking
  vector[J] alpha = alpha_vert + A_mu * mu + A_sigma * sigma;
}
model {
  mu ~ normal(0, 5);
  sigma ~ exponential(1);

  // Prior on alpha_vert derived from alpha ~ normal(mu, sigma):
  //   alpha_vert[j] = alpha[j] - A_mu[j]*mu - A_sigma[j]*sigma
  //   => alpha_vert[j] ~ normal(mu*(1-A_mu[j]) - A_sigma[j]*sigma, sigma)
  alpha_vert ~ normal(mu * (1 - A_mu) - A_sigma * sigma, sigma);

  beta ~ normal(0, 2);
  y ~ bernoulli_logit(alpha[group] + X * beta);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[group[n]] + X[n] * beta);
}
