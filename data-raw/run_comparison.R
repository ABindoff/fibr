## Centred vs non-centred holonomy comparison on the sparse GLMM.
##
## Runs the diagonal holonomy diagnostic at gaps {3,10,25,50} for both
## the centred and non-centred parameterisations and produces a two-facet
## gap-profile figure with per-group h_j lines and bootstrap 90% ribbons.
##
## Non-centred draws are loaded from cache (glmm_sparse_nc_draws.rds);
## if the cache is absent the non-centred model is refit from scratch.
##
## Writes:
##   data-raw/holonomy_comparison.png
##   data-raw/holonomy_comparison.rds   (comparison data frame for paper numbers)
##
## Run from package root:  Rscript data-raw/run_comparison.R

library(posterior)
library(ggplot2)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir     <- file.path(pkg_root, "data-raw")
GAP_GRID    <- c(3L, 10L, 25L, 50L)
N_BOOTSTRAP <- 200L
J           <- 8L
BASE_VARS   <- c("mu", "sigma")
FIBER_C     <- paste0("alpha[",       seq_len(J), "]")
FIBER_NC    <- paste0("alpha_tilde[", seq_len(J), "]")

# ── 1. Load draws ─────────────────────────────────────────────────────────────

draws_c  <- readRDS(file.path(out_dir, "glmm_sparse_draws.rds"))
nc_cache <- file.path(out_dir, "glmm_sparse_nc_draws.rds")

if (file.exists(nc_cache)) {
  cat("Loading cached non-centred draws.\n")
  draws_nc <- readRDS(nc_cache)
} else {
  cat("Cache absent — fitting non-centred model.\n")
  library(cmdstanr)
  saved  <- readRDS(file.path(out_dir, "glmm_sparse_data.rds"))
  mod_nc <- cmdstan_model(file.path(pkg_root, "inst", "stan",
                                    "glmm_noncentred.stan"))
  fit_nc <- mod_nc$sample(
    data            = saved$stan_data,
    chains          = 4L,
    parallel_chains = 4L,
    iter_warmup     = 1000L,
    iter_sampling   = 2000L,
    seed            = 123L,
    refresh         = 500L
  )
  draws_nc <- fit_nc$draws()
  saveRDS(draws_nc, nc_cache)
}

# ── 2. Run holonomy_diagnostic at each gap for both parameterisations ─────────

run_gaps <- function(draws, fiber_vars, label) {
  lapply(GAP_GRID, function(g) {
    cat(sprintf("  %s  gap = %d\n", label, g))
    hd <- tryCatch(
      suppressMessages(holonomy_diagnostic(
        chain             = draws,
        base_vars         = BASE_VARS,
        fiber_vars        = fiber_vars,
        min_gap           = g,
        n_bootstrap       = N_BOOTSTRAP,
        structure         = "diagonal",
        weights           = "distance",
        residualize_fiber = TRUE
      )),
      error = function(e) {
        cat(sprintf("    Error: %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(hd)) return(NULL)

    h_point <- Re(hd$eigenvalues)          # J-vector

    # 90% bootstrap CI per group from boot_eigenvalues [n_boot x J]
    boot_re <- apply(Re(hd$boot_eigenvalues), 2L, function(x)
      quantile(x, c(0.05, 0.95), na.rm = TRUE))

    data.frame(
      param   = label,
      gap     = g,
      group   = seq_len(J),
      h       = h_point,
      h_lo    = boot_re[1L, ],
      h_hi    = boot_re[2L, ],
      n_loops = hd$n_loops
    )
  })
}

cat("── Centred ──────────────────────────────────────────────────────────────\n")
rows_c  <- run_gaps(draws_c,  FIBER_C,  "Centred")

cat("── Non-centred ──────────────────────────────────────────────────────────\n")
rows_nc <- run_gaps(draws_nc, FIBER_NC, "Non-centred")

comp_df <- do.call(rbind, c(rows_c, rows_nc))
comp_df$gap_f   <- factor(comp_df$gap, levels = GAP_GRID)
comp_df$group_f <- factor(comp_df$group)
comp_df$param   <- factor(comp_df$param, levels = c("Centred", "Non-centred"))

saveRDS(comp_df, file.path(out_dir, "holonomy_comparison.rds"))
cat("\nSaved: data-raw/holonomy_comparison.rds\n")

# ── 3. Print summary numbers for paper paragraph ──────────────────────────────

cat("\nMean |h_j| by parameterisation and gap:\n")
agg_h <- tapply(abs(comp_df[["h"]]), paste(comp_df[["param"]], comp_df[["gap"]]), mean)
print(round(sort(agg_h), 4))

# Quick check of "substantially larger" claim
for (g in GAP_GRID) {
  h_c  <- comp_df$h[comp_df$param == "Centred"     & comp_df$gap == g]
  h_nc <- comp_df$h[comp_df$param == "Non-centred" & comp_df$gap == g]
  cat(sprintf("gap=%2d  mean h_c=%.4f  mean h_nc=%.4f  ratio=%.1f\n",
              g, mean(h_c), mean(h_nc),
              mean(abs(h_c)) / max(mean(abs(h_nc)), 1e-6)))
}

# ── 4. Figure: h_j vs gap, lines per group, ribbons, two facets ───────────────

pal <- scales::hue_pal()(J)

fig <- ggplot(comp_df, aes(x = gap_f, y = h,
                            colour = group_f, group = group_f)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_ribbon(aes(ymin = h_lo, ymax = h_hi, fill = group_f),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ param, nrow = 1L, scales = "fixed") +
  scale_colour_manual(values = pal, name = "Group") +
  scale_fill_manual(  values = pal, name = "Group") +
  scale_x_discrete(labels = as.character(GAP_GRID)) +
  labs(
    title    = "Centred vs non-centred: per-group transport factors across gap grid",
    subtitle = sprintf(
      "Sparse GLMM  |  J=%d, n_j=3, σ_true=3  |  diagonal estimator, bootstrap 90%% CI",
      J
    ),
    x = "Gap  g",
    y = expression(hat(h)[j] ~ "(transport factor)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "right",
    strip.text       = element_text(face = "bold", size = 12),
    panel.spacing    = unit(1.5, "lines")
  )

ggsave(file.path(out_dir, "holonomy_comparison.png"),
       plot = fig, width = 11, height = 5, dpi = 150)
cat("Saved: data-raw/holonomy_comparison.png\n")
cat("\nDone.\n")
