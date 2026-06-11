## Milestone 2 end-to-end: connection form, curvature, and parallel transport.
##
## Loads the sparse GLMM chain (centred parameterisation), computes the
## analytic Fisher metric connection, integrates it along detected loops,
## and compares the geometric prediction against the Milestone 1 empirical
## transport.
##
## Run from package root: Rscript data-raw/run_connection.R

library(posterior)
library(FNN)
library(ggplot2)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir    <- file.path(pkg_root, "data-raw")
draws_file <- file.path(out_dir, "glmm_sparse_draws.rds")
data_file  <- file.path(out_dir, "glmm_sparse_data.rds")

draws     <- readRDS(draws_file)
saved     <- readRDS(data_file)
stan_data <- saved$stan_data

base_vars  <- c("mu", "sigma")
fiber_vars <- paste0("alpha[", 1:8, "]")
beta_vars  <- c("beta[1]", "beta[2]")

# ── 1. Compute connection ─────────────────────────────────────────────────────

cat("── 1. Connection form ───────────────────────────────────────────────────\n")

conn <- compute_connection(
  chain       = draws,
  base_vars   = base_vars,
  fiber_vars  = fiber_vars,
  method      = "analytic_glmm",
  stan_data   = stan_data,
  beta_vars   = beta_vars,
  n_subsample = 500L
)

print(conn)

p_conn <- plot(conn)
ggsave(file.path(out_dir, "connection_form.png"),
       plot = p_conn, width = 10, height = 5, dpi = 150)
cat("Saved: data-raw/connection_form.png\n\n")

# ── 2. Milestone 1 diagnostic (for loop pairs) ────────────────────────────────

cat("── 2. Holonomy diagnostic (Milestone 1, min_gap=3) ──────────────────────\n")

hd <- holonomy_diagnostic(
  chain             = draws,
  base_vars         = base_vars,
  fiber_vars        = fiber_vars,
  min_gap           = 3L,
  n_bootstrap       = 100L,
  residualize_fiber = TRUE
)

cat(sprintf("Loops available for transport integration: %d\n\n", hd$n_loops))

# ── 3. Synthetic loop validation (Stokes vs numerical integration) ────────────
#
# MCMC loops are too short to enclose meaningful area.  Instead, we trace a
# synthetic circular path in (mu, sigma) space and compare two predictions:
#   (a) numerical: integrate A along the circle discretely
#   (b) Stokes:    exp(F_mean * pi * r^2)
# Agreement between (a) and (b) validates both the connection formula and
# the curvature formula independently of MCMC dynamics.

cat("── 3. Synthetic loop validation ─────────────────────────────────────────\n")

# Small-loop regime only (r <= 0.2): Stokes is accurate here.
radii <- c(0.05, 0.10, 0.15, 0.20)

for (r in radii) {
  res <- synthetic_holonomy_loop(conn, radius = r, n_steps = 400L)
  cat(sprintf("\nRadius = %.2f  (area = %.5f)\n", r, pi * r^2))
  print(res[, c("j", "H_numerical", "H_stokes", "F_j", "delta_alpha")])
}

# ── 4. Multi-radius holonomy curve ────────────────────────────────────────────

cat("\n── 4. Holonomy vs loop radius ────────────────────────────────────────────\n")

# Use small to moderate radii (beyond ~0.3 the linearised Stokes diverges)
r_grid <- seq(0.02, 0.30, by = 0.02)
df_curve <- do.call(rbind, lapply(r_grid, function(r) {
  res <- synthetic_holonomy_loop(conn, radius = r, n_steps = 400L)
  data.frame(
    radius      = r,
    group       = factor(res$j),
    H_numerical = res$H_numerical,
    H_stokes    = res$H_stokes
  )
}))

cat("Holonomy vs radius (per group):\n")
print(df_curve, row.names = FALSE, digits = 4)

p_curve <- ggplot(df_curve, aes(x = radius, colour = group)) +
  geom_vline(xintercept = 0.20, linetype = "dashed", colour = "grey70") +
  geom_line(aes(y = H_stokes), linewidth = 0.7) +
  geom_point(aes(y = H_numerical), size = 1.8, shape = 16) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  annotate("text", x = 0.20, y = Inf, label = "r = 0.20",
           hjust = -0.1, vjust = 1.5, size = 3.2, colour = "grey50") +
  scale_colour_viridis_d(option = "turbo", name = "Group j") +
  labs(
    title    = "Holonomy vs loop radius — GLMM centred parameterisation",
    subtitle = "Lines: Stokes  H = 1 + Fⱼ·πr²/α₀ⱼ  |  Points: numerical integration  |  centre = posterior mean",
    x        = expression(r ~ "(loop radius, posterior SD units)"),
    y        = expression(H[j] ~ "(holonomy per group)")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "holonomy_vs_radius.png"),
       plot = p_curve, width = 7, height = 5, dpi = 150)
cat("Saved: data-raw/holonomy_vs_radius.png\n")

cat(sprintf(
  "\nMilestone 1  lambda_1 (min_gap=3) = %.4f\n",
  Mod(hd$eigenvalues[1L])
))
cat("The Stokes curve above shows the GEOMETRIC holonomy prediction;\n")
cat("Milestone 1 measures STATISTICAL fiber memory (autocorrelation + geometry).\n")
cat("Their difference quantifies the MCMC baseline contamination.\n")
cat("\nDone.\n")
