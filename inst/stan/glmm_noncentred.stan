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
  vector[J] alpha_tilde;   // standardised group effects; independent of (mu, sigma)
  vector[2] beta;
}
transformed parameters {
  // alpha lives in the same space as the centred model for direct comparison
  vector[J] alpha = mu + sigma * alpha_tilde;
}
model {
  mu ~ normal(0, 5);
  sigma ~ exponential(1);
  alpha_tilde ~ normal(0, 1);
  beta ~ normal(0, 2);
  y ~ bernoulli_logit(alpha[group] + X * beta);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[group[n]] + X[n] * beta);
}
