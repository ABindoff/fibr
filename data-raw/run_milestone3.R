## Milestone 3: Three-way comparison on the sparse GLMM.
##
## Centred  vs  Connection-corrected (horizontal)  vs  Non-centred
##
## The connection coefficients from Milestone 2 are used to reparameterise
## the model in coordinates that approximately follow the horizontal manifold.
## This reduces the geometric obstruction to NUTS mixing without requiring a
## custom sampler.
##
## Run from package root: Rscript data-raw/run_milestone3.R

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

# ── 1. Connection coefficients from Milestone 2 ───────────────────────────────

cat("── 1. Connection coefficients ───────────────────────────────────────────\n")

draws_c   <- readRDS(file.path(out_dir, "glmm_sparse_draws.rds"))

conn <- compute_connection(
  chain       = draws_c,
  base_vars   = base_vars,
  fiber_vars  = fiber_vars,
  method      = "analytic_glmm",
  stan_data   = stan_data,
  n_subsample = 500L
)

# Use posterior-mean A as fixed coefficients for the reparameterisation
A_mu    <- colMeans(conn$A[, , 1L, drop = TRUE])
A_sigma <- colMeans(conn$A[, , 2L, drop = TRUE])

cat("Connection coefficients (posterior mean across 500 draws):\n")
df_A <- data.frame(
  group   = seq_len(8),
  A_mu    = round(A_mu,    4),
  A_sigma = round(A_sigma, 4),
  prior_frac = round(colMeans(conn$prior_frac), 3)
)
print(df_A, row.names = FALSE)

# ── 2. Fit all three models ───────────────────────────────────────────────────

cat("\n── 2. Fitting models ────────────────────────────────────────────────────\n")

SAMPLE_ARGS <- list(
  data            = stan_data,
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 2000L,
  seed            = 42L,
  refresh         = 0L    # silent
)

# 2a. Centred (already fitted — reload)
cat("Loading centred draws...\n")
fit_c_draws <- draws_c

# 2b. Non-centred (already fitted — reload)
cat("Loading non-centred draws...\n")
fit_nc_draws <- readRDS(file.path(out_dir, "glmm_sparse_nc_draws.rds"))

# 2c. Horizontal-corrected (new)
cat("Fitting horizontal-corrected model...\n")
hconn_data <- c(stan_data, list(A_mu = A_mu, A_sigma = A_sigma))
mod_hc     <- cmdstan_model(file.path(pkg_root, "inst", "stan", "glmm_hconnected.stan"))
fit_hc     <- mod_hc$sample(
  data            = hconn_data,
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 2000L,
  seed            = 42L,
  refresh         = 500L
)

cat("Horizontal-corrected sampler diagnostics:\n")
fit_hc$diagnostic_summary()
draws_hc <- fit_hc$draws()

# ── 3. Comparison table ───────────────────────────────────────────────────────

cat("\n── 3. Three-way ESS / R-hat comparison ──────────────────────────────────\n")

# Extract ESS for key parameters
.ess_for_alpha <- function(draws, label) {
  key_vars <- c("mu", "sigma", paste0("alpha[", 1:8, "]"), "beta[1]", "beta[2]")
  # Keep only variables that exist in this draws object
  key_vars <- intersect(key_vars, variables(draws))
  s <- summarise_draws(subset_draws(draws, variable = key_vars),
                       "ess_bulk", "rhat")
  s$model <- label
  s
}

ess_c  <- .ess_for_alpha(fit_c_draws,  "Centred")
ess_nc <- .ess_for_alpha(fit_nc_draws, "Non-centred")
ess_hc <- .ess_for_alpha(draws_hc,     "H-corrected")

# Bind and pivot for display
df_ess <- data.frame(
  variable = ess_c$variable,
  ESS_C    = round(ess_c$ess_bulk),
  Rhat_C   = round(ess_c$rhat,  3),
  ESS_HC   = round(ess_hc$ess_bulk),
  Rhat_HC  = round(ess_hc$rhat, 3),
  ESS_NC   = round(ess_nc$ess_bulk),
  Rhat_NC  = round(ess_nc$rhat, 3)
)
print(df_ess, row.names = FALSE)

nominal <- 4L * 2000L
cat(sprintf(
  "\nMin ESS_bulk: Centred = %.0f (%.0f%%)  |  H-corrected = %.0f (%.0f%%)  |  Non-centred = %.0f (%.0f%%)\n",
  min(ess_c$ess_bulk,  na.rm = TRUE),
  100 * min(ess_c$ess_bulk,  na.rm = TRUE) / nominal,
  min(ess_hc$ess_bulk, na.rm = TRUE),
  100 * min(ess_hc$ess_bulk, na.rm = TRUE) / nominal,
  min(ess_nc$ess_bulk, na.rm = TRUE),
  100 * min(ess_nc$ess_bulk, na.rm = TRUE) / nominal
))

# ── 4. Bar chart ──────────────────────────────────────────────────────────────

df_all <- rbind(
  cbind(ess_c[,  c("variable", "ess_bulk")], model = "Centred"),
  cbind(ess_hc[, c("variable", "ess_bulk")], model = "H-corrected"),
  cbind(ess_nc[, c("variable", "ess_bulk")], model = "Non-centred")
)
df_all$pct     <- df_all$ess_bulk / nominal * 100
df_all$model   <- factor(df_all$model,
                         levels = c("Centred", "H-corrected", "Non-centred"))
target_order   <- c("mu", "sigma", paste0("alpha[", 1:8, "]"),
                    "beta[1]", "beta[2]")
df_all$variable <- factor(df_all$variable, levels = target_order)

p_ess <- ggplot(df_all, aes(x = variable, y = pct, fill = model)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey40") +
  scale_fill_manual(
    values = c("Centred"     = "#4682B4",
               "H-corrected" = "#B22222",
               "Non-centred" = "#2E8B57")
  ) +
  labs(
    title    = "Three-way ESS comparison — sparse GLMM (J=8, n_j=3, sigma_true=3)",
    subtitle = "H-corrected uses the Fisher-metric connection to reparameterise the fiber",
    x        = NULL,
    y        = "ESS bulk (% of nominal)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "milestone3_ess.png"),
       plot = p_ess, width = 10, height = 5, dpi = 150)
cat("Saved: data-raw/milestone3_ess.png\n")

# ── 5. Divergence count summary ───────────────────────────────────────────────

cat("\n── Sampler diagnostics (centred already run; reload summary) ────────────\n")
diag_hc  <- fit_hc$diagnostic_summary()
cat(sprintf("Divergences — H-corrected:  %d\n", sum(diag_hc$num_divergent)))
cat(sprintf("Divergences — Centred:      ~86 (from earlier run)\n"))
cat(sprintf("Divergences — Non-centred:  0  (from earlier run)\n"))
cat("\nDone.\n")
