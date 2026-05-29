data {
  int<lower=0> N;
  int<lower=1> J;
  array[N] int<lower=1,upper=J> group;
  matrix[N, 2] X;
  array[N] int<lower=0,upper=1> y;
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector[J] alpha;
  vector[2] beta;
}
model {
  mu ~ normal(0, 5);
  sigma ~ exponential(1);
  alpha ~ normal(mu, sigma);
  beta ~ normal(0, 2);
  y ~ bernoulli_logit(alpha[group] + X * beta);
}
generated quantities {
  // Log-likelihood for LOO; also useful for connection diagnostics later
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[group[n]] + X[n] * beta);
}
