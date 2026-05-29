## Centred vs non-centred holonomy comparison on the sparse GLMM.
##
## Fits the non-centred model to the existing sparse data, runs the holonomy
## diagnostic on both chains at min_gap=3 (inside the autocorrelation time),
## and produces a side-by-side eigenspectrum plot.
##
## Run from package root: Rscript data-raw/run_comparison.R

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

MIN_GAP     <- 3L
N_BOOTSTRAP <- 200L

# ── 1. Load sparse data and centred draws ─────────────────────────────────────

saved   <- readRDS(file.path(pkg_root, "data-raw", "glmm_sparse_data.rds"))
draws_c <- readRDS(file.path(pkg_root, "data-raw", "glmm_sparse_draws.rds"))

cat("── Non-centred fit ───────────────────────────────────────────────────────\n")

# ── 2. Fit non-centred model ──────────────────────────────────────────────────

mod_nc <- cmdstan_model(file.path(pkg_root, "inst", "stan", "glmm_noncentred.stan"))

fit_nc <- mod_nc$sample(
  data            = saved$stan_data,
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 2000L,
  seed            = 123L,
  refresh         = 500L
)

cat("\nSampler diagnostics (non-centred):\n")
fit_nc$diagnostic_summary()

draws_nc <- fit_nc$draws()

# Quick ESS comparison
summ_c  <- summarise_draws(draws_c,  "ess_bulk")[grep("^alpha\\[", summarise_draws(draws_c,  "ess_bulk")$variable), ]
summ_nc <- summarise_draws(draws_nc, "ess_bulk")[grep("^alpha_tilde\\[", summarise_draws(draws_nc, "ess_bulk")$variable), ]

cat(sprintf("\nMin ESS_bulk — centred alpha:       %.0f  (%.1f%%)\n",
            min(summ_c$ess_bulk,  na.rm = TRUE),
            100 * min(summ_c$ess_bulk,  na.rm = TRUE) / (4 * 2000)))
cat(sprintf("Min ESS_bulk — non-centred tilde:   %.0f  (%.1f%%)\n",
            min(summ_nc$ess_bulk, na.rm = TRUE),
            100 * min(summ_nc$ess_bulk, na.rm = TRUE) / (4 * 2000)))

saveRDS(draws_nc, file.path(pkg_root, "data-raw", "glmm_sparse_nc_draws.rds"))

# ── 3. Run holonomy diagnostic on both ───────────────────────────────────────

base_vars  <- c("mu", "sigma")

# Centred: fiber = alpha[1..8] (raw intercepts; residualised against base)
fiber_c    <- paste0("alpha[",       1:8, "]")

# Non-centred: fiber = alpha_tilde[1..8] (standardised intercepts)
# alpha_tilde is independent of (mu,sigma) by construction;
# residualisation has no effect but we apply it for consistency.
fiber_nc   <- paste0("alpha_tilde[", 1:8, "]")

cat("\n── Holonomy diagnostic: centred (min_gap = ", MIN_GAP, ") ────────────────\n", sep = "")
hd_c <- holonomy_diagnostic(
  chain              = draws_c,
  base_vars          = base_vars,
  fiber_vars         = fiber_c,
  epsilon            = NULL,
  n_bootstrap        = N_BOOTSTRAP,
  min_gap            = MIN_GAP,
  residualize_fiber  = TRUE
)
print(hd_c)

cat("\n── Holonomy diagnostic: non-centred (min_gap = ", MIN_GAP, ") ──────────────\n", sep = "")
hd_nc <- holonomy_diagnostic(
  chain              = draws_nc,
  base_vars          = base_vars,
  fiber_vars         = fiber_nc,
  epsilon            = NULL,
  n_bootstrap        = N_BOOTSTRAP,
  min_gap            = MIN_GAP,
  residualize_fiber  = TRUE
)
print(hd_nc)

# ── 4. Side-by-side eigenspectrum plot ────────────────────────────────────────

.evals_panel <- function(hd, label) {
  evals <- hd$eigenvalues
  boot  <- hd$boot_eigenvalues

  boot_df <- if (!is.null(boot)) {
    df <- data.frame(x = as.vector(Re(boot)), y = as.vector(Im(boot)))
    df[complete.cases(df), ]
  } else {
    data.frame(x = numeric(0), y = numeric(0))
  }
  boot_df$panel <- label

  pts_df <- data.frame(
    x     = Re(evals),
    y     = Im(evals),
    idx   = seq_along(evals),
    panel = label
  )

  list(
    boot = boot_df,
    pts  = pts_df,
    sub  = sprintf("||H-I||_F = %.3f  |  %d loops", hd$frobenius_dev, hd$n_loops)
  )
}

p_c  <- .evals_panel(hd_c,  sprintf("Centred  (min_gap=%d)", MIN_GAP))
p_nc <- .evals_panel(hd_nc, sprintf("Non-centred  (min_gap=%d)", MIN_GAP))

boot_all <- rbind(p_c$boot, p_nc$boot)
pts_all  <- rbind(p_c$pts,  p_nc$pts)

# Subtitle per panel via labeller
sub_map <- c(p_c$sub, p_nc$sub)
names(sub_map) <- c(p_c$pts$panel[1], p_nc$pts$panel[1])
panel_labeller <- labeller(panel = sub_map)

theta  <- seq(0, 2 * pi, length.out = 300)
circle <- data.frame(x = cos(theta), y = sin(theta),
                     panel = p_c$pts$panel[1])  # drawn in both via facet

comparison_plot <-
  ggplot() +
  # Unit circle (one data frame per panel via dummy facet)
  geom_path(
    data = rbind(
      transform(circle, panel = p_c$pts$panel[1]),
      transform(circle, panel = p_nc$pts$panel[1])
    ),
    aes(x = x, y = y),
    colour = "grey60", linetype = "dashed", linewidth = 0.5
  ) +
  # Bootstrap clouds
  geom_point(
    data = boot_all,
    aes(x = x, y = y),
    colour = "steelblue", alpha = 0.06, size = 0.6
  ) +
  # Point estimates
  geom_point(
    data = pts_all,
    aes(x = x, y = y),
    colour = "firebrick", size = 3
  ) +
  geom_text(
    data = pts_all,
    aes(x = x, y = y, label = idx),
    nudge_y = 0.05, size = 2.8, colour = "firebrick"
  ) +
  # Identity point
  geom_point(
    data = rbind(
      data.frame(x = 1, y = 0, panel = p_c$pts$panel[1]),
      data.frame(x = 1, y = 0, panel = p_nc$pts$panel[1])
    ),
    aes(x = x, y = y),
    shape = 3, size = 4, colour = "black", stroke = 1
  ) +
  facet_wrap(~ panel, labeller = panel_labeller) +
  coord_equal() +
  labs(
    title    = "Centred vs non-centred: holonomy eigenspectrum",
    subtitle = sprintf("Sparse GLMM  |  J=8, n_j=3, sigma_true=3  |  min_gap=%d", MIN_GAP),
    x        = expression(Re(lambda)),
    y        = expression(Im(lambda))
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

out_file <- file.path(pkg_root, "data-raw", "holonomy_comparison.png")
ggsave(out_file, plot = comparison_plot, width = 12, height = 6, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_file))
cat("Done.\n")
