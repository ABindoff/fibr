// Partial non-centred parameterisation of the centred GLMM.
//
// Each group j uses a blending weight w[j] in [0, 1], supplied as data
// (typically w[j] = pi_j = prior fraction from a pilot run).
//
// Parameterisation (Betancourt & Girolami 2013):
//   psi[j] ~ N(0, sigma^(1 - w[j]))
//   alpha[j] = mu + sigma^w[j] * psi[j]
//
// Special cases:
//   w[j] = 0 -> psi[j] ~ N(0, sigma),   alpha[j] = mu + psi[j]   [centred]
//   w[j] = 1 -> psi[j] ~ N(0, 1),       alpha[j] = mu + sigma*psi[j]  [non-centred]
//
// The marginal alpha[j] ~ N(mu, sigma^2) for all w[j].
// w[j] is kept as real-valued data to allow fractional blending.

data {
  int<lower=0> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> group;
  matrix[N, 2] X;
  array[N] int<lower=0, upper=1> y;
  vector<lower=0, upper=1>[J] w;   // per-group blending weights
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector[J] psi;    // partially whitened group effects
  vector[2] beta;
}
transformed parameters {
  vector[J] alpha;
  for (j in 1:J)
    alpha[j] = mu + pow(sigma, w[j]) * psi[j];
}
model {
  mu    ~ normal(0, 5);
  sigma ~ exponential(1);
  for (j in 1:J)
    psi[j] ~ normal(0, pow(sigma, 1.0 - w[j]));
  beta  ~ normal(0, 2);
  y     ~ bernoulli_logit(alpha[group] + X * beta);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha[group[n]] + X[n] * beta);
}
