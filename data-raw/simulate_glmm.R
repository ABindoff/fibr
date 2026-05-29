## Simulate data from the toy GLMM and fit the centred Stan model.
##
## Produces:  data-raw/glmm_data.rds   -- simulated data list (Stan-ready)
##            data-raw/glmm_draws.rds  -- posterior draws as a draws_array
##            data-raw/glmm_fit.rds    -- full CmdStanMCMC object
##
## Run from the package root, e.g.:
##   source("data-raw/simulate_glmm.R")

library(cmdstanr)
library(posterior)

# ── 1. Simulate ───────────────────────────────────────────────────────────────

set.seed(42)

J      <- 8L        # groups
n_j    <- 20L       # observations per group
N      <- J * n_j

mu_true    <- 0
sigma_true <- 1.5
beta_true  <- c(0.8, -0.5)

alpha_true <- rnorm(J, mu_true, sigma_true)
group_id   <- rep(seq_len(J), each = n_j)
X          <- matrix(rnorm(N * 2L), ncol = 2L)
eta        <- alpha_true[group_id] + X %*% beta_true
y          <- rbinom(N, 1L, plogis(eta))

cat(sprintf(
  "Simulated: N=%d, J=%d, mean(y)=%.3f\n",
  N, J, mean(y)
))
cat(sprintf(
  "True alpha range: [%.2f, %.2f]\n",
  min(alpha_true), max(alpha_true)
))

stan_data <- list(
  N     = N,
  J     = J,
  group = group_id,
  X     = X,
  y     = y
)

saveRDS(
  list(stan_data = stan_data, alpha_true = alpha_true,
       mu_true = mu_true, sigma_true = sigma_true, beta_true = beta_true),
  "data-raw/glmm_data.rds"
)
cat("Saved: data-raw/glmm_data.rds\n")

# ── 2. Compile and fit ────────────────────────────────────────────────────────

stan_file <- system.file("stan", "glmm_centred.stan", package = "fibr")
if (!nzchar(stan_file)) {
  # Fall back to relative path when running from package root before install
  stan_file <- "inst/stan/glmm_centred.stan"
}

mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data            = stan_data,
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 2000L,
  seed            = 42L,
  refresh         = 500L
)

# ── 3. Basic diagnostics ──────────────────────────────────────────────────────

cat("\n── Sampler diagnostics ──────────────────────────────────────────────────\n")
fit$diagnostic_summary()

draws <- fit$draws()   # draws_array [iter, chain, variable]

rhat_vals  <- summarise_draws(draws, "rhat")
ess_vals   <- summarise_draws(draws, "ess_bulk", "ess_tail")

# Flag any parameter with R-hat > 1.01 or ESS < 400
rhat_bad <- rhat_vals[!is.na(rhat_vals$rhat) & rhat_vals$rhat > 1.01, ]
ess_bad  <- ess_vals[
  (!is.na(ess_vals$ess_bulk) & ess_vals$ess_bulk < 400) |
  (!is.na(ess_vals$ess_tail)  & ess_vals$ess_tail  < 400), ]

if (nrow(rhat_bad) > 0L) {
  cat("\nWARNING: R-hat > 1.01 for:\n")
  print(rhat_bad)
} else {
  cat("\nAll R-hat <= 1.01\n")
}

if (nrow(ess_bad) > 0L) {
  cat("\nWARNING: ESS < 400 for:\n")
  print(ess_bad)
} else {
  cat("All ESS >= 400\n")
}

# Quick recovery check vs truth
key_params <- c("mu", "sigma", "beta[1]", "beta[2]")
cat("\n── Parameter recovery ───────────────────────────────────────────────────\n")
summ <- summarise_draws(
  subset_draws(draws, variable = key_params),
  "mean", "sd", ~quantile(.x, c(0.05, 0.95))
)
truth_row <- c(mu_true, sigma_true, beta_true)
summ$truth <- c(truth_row)
print(summ)

# ── 4. Save outputs ───────────────────────────────────────────────────────────

saveRDS(draws, "data-raw/glmm_draws.rds")
cat("\nSaved: data-raw/glmm_draws.rds\n")

saveRDS(fit,   "data-raw/glmm_fit.rds")
cat("Saved: data-raw/glmm_fit.rds\n")

cat("\nDone. draws_array dimensions: ", paste(dim(draws), collapse = " x "), "\n")
cat("  [iterations x chains x variables]\n")
cat("  Variables include: mu, sigma, alpha[1..J], beta[1..2], log_lik[1..N]\n")
