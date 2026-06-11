## ESS comparison: centred vs horizontally-corrected vs non-centred,
## on the sparse GLMM benchmark (J = 8, n_j = 3, sigma_true = 3).
##
## Fills Table "ess_comparison" in manuscript/fibr_paper.tex.
## Writes:
##   data-raw/ess_comparison.rds   full ESS table, all parameters x 3 models
##   data-raw/ess_table.tex        LaTeX rows for copy-paste into the paper
##
## Requires data-raw/glmm_sparse_data.rds and glmm_sparse_draws.rds
## (produced by data-raw/simulate_glmm_sparse.R and run_diagnostic_sparse.R).
##
## Run from package root:  Rscript data-raw/run_ess_comparison.R

library(cmdstanr)
library(posterior)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir <- file.path(pkg_root, "data-raw")

SEED        <- 123L
N_WARMUP    <- 1000L
N_SAMPLING  <- 2000L
N_CHAINS    <- 4L
NOMINAL_ESS <- N_CHAINS * N_SAMPLING

# ── 1. Load sparse benchmark data and existing centred draws ─────────────────

saved     <- readRDS(file.path(out_dir, "glmm_sparse_data.rds"))
stan_data <- saved$stan_data
J         <- stan_data$J

draws_c <- readRDS(file.path(out_dir, "glmm_sparse_draws.rds"))

# ── 2. Connection coefficients at the posterior mean (for h-connected) ───────

conn <- compute_connection(
  chain       = draws_c,
  base_vars   = c("mu", "sigma"),
  fiber_vars  = paste0("alpha[", seq_len(J), "]"),
  method      = "analytic_glmm",
  stan_data   = stan_data,
  n_subsample = 500L
)

# conn$A is [n_sub x J x 2]; average over draws
A_bar   <- apply(conn$A, c(2L, 3L), mean)
A_mu    <- A_bar[, 1L]
A_sigma <- A_bar[, 2L]

cat("Posterior-mean connection coefficients:\n")
print(round(cbind(A_mu = A_mu, A_sigma = A_sigma), 4))

# ── 3. Fit the three models ───────────────────────────────────────────────────

fit_model <- function(stan_file, data, label) {
  cat(sprintf("\n── Fitting %s ───────────────────────────────────────\n", label))
  mod <- cmdstan_model(file.path(pkg_root, "inst", "stan", stan_file))
  fit <- mod$sample(
    data            = data,
    chains          = N_CHAINS,
    parallel_chains = N_CHAINS,
    iter_warmup     = N_WARMUP,
    iter_sampling   = N_SAMPLING,
    seed            = SEED,
    refresh         = 0L
  )
  print(fit$diagnostic_summary())
  fit
}

fit_c  <- fit_model("glmm_centred.stan",    stan_data, "centred")
fit_h  <- fit_model("glmm_hconnected.stan",
                    c(stan_data, list(A_mu = A_mu, A_sigma = A_sigma)),
                    "h-connected")
fit_nc <- fit_model("glmm_noncentred.stan", stan_data, "non-centred")

# ── 4. ESS table ──────────────────────────────────────────────────────────────

# Parameters reported in the manuscript table.  alpha[1] is defined in all
# three models (the non-centred and h-connected models reconstruct it in
# transformed parameters / generated quantities).
PARS <- c("mu", "sigma", "alpha[1]", "beta[1]")

ess_pct <- function(fit) {
  s <- summarise_draws(fit$draws(PARS), "ess_bulk")
  setNames(100 * s$ess_bulk / NOMINAL_ESS, s$variable)[PARS]
}

tab <- data.frame(
  parameter   = PARS,
  centred     = round(ess_pct(fit_c),  1),
  h_corrected = round(ess_pct(fit_h),  1),
  non_centred = round(ess_pct(fit_nc), 1),
  row.names   = NULL
)

cat("\nESS bulk (% of nominal", NOMINAL_ESS, "):\n")
print(tab, row.names = FALSE)

saveRDS(tab, file.path(out_dir, "ess_comparison.rds"))

# ── 5. LaTeX rows for the manuscript ──────────────────────────────────────────

tex_par <- c("$\\mu$", "$\\sigma$", "$\\alpha[1]$", "$\\beta[1]$")
tex <- paste0(
  sprintf("    %s & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\",
          tex_par, tab$centred, tab$h_corrected, tab$non_centred),
  collapse = "\n"
)
writeLines(tex, file.path(out_dir, "ess_table.tex"))
cat("\nWrote LaTeX rows to data-raw/ess_table.tex:\n", tex, "\n")
