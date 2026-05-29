## Milestone 4: four-way comparison on the sparse GLMM.
##
## Centred (Stan NUTS, Euclidean)
## Non-centred (Stan NUTS, Euclidean)
## H-corrected (Stan NUTS, fixed connection reparameterisation)
## RMHMC (explicit Riemannian leapfrog, per-trajectory Fisher metric)
##
## Run from package root: Rscript data-raw/run_rmhmc.R

library(cmdstanr)
library(posterior)
library(ggplot2)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir   <- file.path(pkg_root, "data-raw")
saved     <- readRDS(file.path(out_dir, "glmm_sparse_data.rds"))
stan_data <- saved$stan_data

base_vars  <- c("mu", "sigma")
fiber_vars <- paste0("alpha[", 1:8, "]")
key_vars   <- c("mu", "sigma", fiber_vars, "beta[1]", "beta[2]")

# ── 1. Sanity-check the metric ────────────────────────────────────────────────

cat("── 1. Metric sanity check ────────────────────────────────────────────────\n")

draws_c    <- readRDS(file.path(out_dir, "glmm_sparse_draws.rds"))
full_mat_c <- do.call(rbind, lapply(seq_len(dim(draws_c)[2L]), function(ci) {
  m   <- as_draws_matrix(subset_draws(draws_c, chain=ci))
  mat <- as.matrix(m)
  colnames(mat) <- variables(m)
  mat
}))

set.seed(1L)
check_idx <- sample.int(nrow(full_mat_c), 5L)
cat("Diagonal metric entries (always > 0) at 5 chain draws:\n")
for (i in check_idx) {
  row  <- full_mat_c[i, ]
  mu_i <- row[["mu"]]
  ls_i <- log(row[["sigma"]])
  al_i <- as.vector(row[paste0("alpha[", 1:8, "]")])
  be_i <- as.vector(row[c("beta[1]", "beta[2]")])
  G_d  <- .glmm_diag_metric(mu_i, ls_i, al_i, be_i, stan_data)
  cat(sprintf("  Draw %d: min G_diag = %.4f  (G_ls = %.4f, sigma = %.3f)\n",
              i, min(G_d), G_d[2L], row[["sigma"]]))
}
cat("\nFull G_BB off-diagonal causes non-PD outside mode (centred GLMM funnel).\n")
cat("Diagonal metric used for RMHMC: always PD, captures key adaptation.\n\n")

# ── 2. RMHMC chains ───────────────────────────────────────────────────────────

cat("\n── 2. RMHMC samplers ────────────────────────────────────────────────────\n")

# Initialise near the Stan posterior mean so we measure mixing, not burn-in
summ_c   <- summarise_draws(draws_c, "mean")
get_mean <- function(pat) summ_c$mean[grep(pat, summ_c$variable)]
init_pm  <- list(
  mu    = get_mean("^mu$"),
  sigma = get_mean("^sigma$"),
  alpha = get_mean("^alpha\\["),
  beta  = get_mean("^beta\\[")
)

RMHMC_ARGS <- list(
  stan_data   = stan_data,
  n_iter      = 2000L,
  n_warmup    = 2000L,
  n_chains    = 4L,
  L           = 1L,        # MALA-style: single step, stable in funnels
  epsilon     = 0.10,
  target_rate = 0.57,      # MALA optimum
  init        = init_pm,
  seed        = 42L,
  verbose     = TRUE
)

cat("\n2a. Diagonal metric (per-parameter step sizes, no off-diagonal):\n")
draws_diag <- do.call(riemannian_mcmc, c(RMHMC_ARGS, list(method = "diagonal")))

cat("\n2b. SoftAbs metric (full metric with off-diagonal coupling, always PD):\n")
draws_sa <- do.call(riemannian_mcmc, c(RMHMC_ARGS, list(method = "softabs",
                                                          softabs_alpha = 1.0)))

cat("\nDiagonal RMHMC summary:\n")
print(summarise_draws(subset_draws(draws_diag, variable = key_vars),
                      "ess_bulk", "rhat"), n = Inf)

cat("\nSoftAbs RMHMC summary:\n")
print(summarise_draws(subset_draws(draws_sa, variable = key_vars),
                      "ess_bulk", "rhat"), n = Inf)

saveRDS(draws_diag, file.path(out_dir, "glmm_sparse_rmhmc_diag.rds"))
saveRDS(draws_sa,   file.path(out_dir, "glmm_sparse_rmhmc_sa.rds"))

# ── 3. Reload other chains ────────────────────────────────────────────────────

cat("\n── 3. Four-way ESS comparison ───────────────────────────────────────────\n")

draws_nc <- readRDS(file.path(out_dir, "glmm_sparse_nc_draws.rds"))

# Re-run H-corrected (quick Stan run)
conn <- compute_connection(draws_c, base_vars, fiber_vars,
                           method="analytic_glmm", stan_data=stan_data,
                           n_subsample=500L)
A_mu    <- colMeans(conn$A[, , 1L, drop=TRUE])
A_sigma <- colMeans(conn$A[, , 2L, drop=TRUE])
mod_hc  <- cmdstan_model(file.path(pkg_root, "inst", "stan", "glmm_hconnected.stan"))
fit_hc  <- mod_hc$sample(
  data=c(stan_data, list(A_mu=A_mu, A_sigma=A_sigma)),
  chains=4L, parallel_chains=4L,
  iter_warmup=1000L, iter_sampling=2000L,
  seed=42L, refresh=0L
)
draws_hc <- fit_hc$draws()

.get_ess <- function(draws, label) {
  s <- summarise_draws(subset_draws(draws, variable=intersect(key_vars, variables(draws))),
                       "ess_bulk")
  data.frame(variable=s$variable, ess=s$ess_bulk, model=label)
}

nominal <- 4L * 2000L
df_ess <- rbind(
  .get_ess(draws_c,    "Centred"),
  .get_ess(draws_hc,   "H-corrected"),
  .get_ess(draws_diag, "RMHMC-diag"),
  .get_ess(draws_sa,   "RMHMC-SoftAbs"),
  .get_ess(draws_nc,   "Non-centred")
)
df_ess$pct      <- df_ess$ess / nominal * 100
df_ess$model    <- factor(df_ess$model,
                           levels=c("Centred","H-corrected",
                                    "RMHMC-diag","RMHMC-SoftAbs","Non-centred"))
df_ess$variable <- factor(df_ess$variable, levels=key_vars)

cat("\nMin ESS_bulk per model:\n")
agg <- aggregate(pct ~ model, data=df_ess, FUN=min)
print(agg, row.names=FALSE)

# ── 4. Bar chart ──────────────────────────────────────────────────────────────

p_ess <- ggplot(df_ess, aes(x=variable, y=pct, fill=model)) +
  geom_col(position="dodge", width=0.72) +
  geom_hline(yintercept=100, linetype="dotted", colour="grey40") +
  scale_fill_manual(
    values=c("Centred"       = "#4682B4",
             "H-corrected"   = "#B22222",
             "RMHMC-diag"    = "#FF8C00",
             "RMHMC-SoftAbs" = "#9932CC",
             "Non-centred"   = "#2E8B57")
  ) +
  labs(
    title="Five-way ESS comparison — sparse GLMM (J=8, n_j=3, sigma_true=3)",
    subtitle=paste("SoftAbs: full G with off-diagonal coupling, always PD;",
                   "Diagonal: per-parameter scaling only"),
    x=NULL, y="ESS bulk (% of nominal)", fill=NULL
  ) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=40, hjust=1),
        legend.position="bottom")

ggsave(file.path(out_dir, "milestone4_ess.png"),
       plot=p_ess, width=10, height=5, dpi=150)
cat("Saved: data-raw/milestone4_ess.png\n")

cat("\nDone.\n")
