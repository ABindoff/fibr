## Sparse-data GLMM simulation designed to stress the centred parameterisation.
##
## n_j = 3 obs/group + sigma_true = 3 creates a strong funnel: the posterior
## of alpha | (mu, sigma) is dominated by the prior, so the centred NUTS chain
## mixes slowly in the fiber.  The holonomy diagnostic should now see non-trivial
## structure because loop endpoints are no longer approximately independent.
##
## Produces:
##   data-raw/glmm_sparse_data.rds
##   data-raw/glmm_sparse_draws.rds
##   data-raw/glmm_sparse_fit.rds
##
## Run from package root: Rscript data-raw/simulate_glmm_sparse.R

library(cmdstanr)
library(posterior)

# ── 1. Simulate ───────────────────────────────────────────────────────────────

set.seed(123)

J          <- 8L
n_j        <- 3L          # was 20; now sparse — prior dominates
N          <- J * n_j

mu_true    <- 0
sigma_true <- 3           # was 1.5; now wide — strong funnel in centred space
beta_true  <- c(0.8, -0.5)

alpha_true <- rnorm(J, mu_true, sigma_true)
group_id   <- rep(seq_len(J), each = n_j)
X          <- matrix(rnorm(N * 2L), ncol = 2L)
eta        <- alpha_true[group_id] + X %*% beta_true
y          <- rbinom(N, 1L, plogis(eta))

cat(sprintf("Sparse model: N=%d, J=%d, n_j=%d, sigma_true=%.1f\n",
            N, J, n_j, sigma_true))
cat(sprintf("mean(y) = %.3f\n", mean(y)))
cat(sprintf("alpha range: [%.2f, %.2f]\n", min(alpha_true), max(alpha_true)))

stan_data <- list(N = N, J = J, group = group_id, X = X, y = y)

saveRDS(
  list(stan_data = stan_data, alpha_true = alpha_true,
       mu_true = mu_true, sigma_true = sigma_true, beta_true = beta_true),
  "data-raw/glmm_sparse_data.rds"
)

# ── 2. Fit ────────────────────────────────────────────────────────────────────

stan_file <- "inst/stan/glmm_centred.stan"
mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data            = stan_data,
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 2000L,
  seed            = 123L,
  refresh         = 500L
)

# ── 3. Diagnostics ────────────────────────────────────────────────────────────

cat("\n── Sampler diagnostics ──────────────────────────────────────────────────\n")
fit$diagnostic_summary()

draws <- fit$draws()
summ  <- summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail")

alpha_summ <- summ[grep("^alpha", summ$variable), ]
cat("\n── Alpha parameter summary ──────────────────────────────────────────────\n")
print(alpha_summ, n = Inf)

cat(sprintf("\nMin ESS_bulk (alpha): %.0f  (%.1f%% of nominal)\n",
            min(alpha_summ$ess_bulk, na.rm = TRUE),
            100 * min(alpha_summ$ess_bulk, na.rm = TRUE) / (4 * 2000)))

cat(sprintf("Max R-hat    (alpha): %.4f\n",
            max(alpha_summ$rhat, na.rm = TRUE)))

# ── 4. Save ───────────────────────────────────────────────────────────────────

saveRDS(draws, "data-raw/glmm_sparse_draws.rds")
saveRDS(fit,   "data-raw/glmm_sparse_fit.rds")
cat("\nSaved sparse draws and fit.\n")
