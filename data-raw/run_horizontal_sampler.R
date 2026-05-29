## Milestone 3: compare the horizontal-corrected sampler against the standard
## centred sampler on the sparse GLMM.
##
## Run from package root: Rscript data-raw/run_horizontal_sampler.R

library(posterior)
library(ggplot2)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir    <- file.path(pkg_root, "data-raw")
saved      <- readRDS(file.path(out_dir, "glmm_sparse_data.rds"))
stan_data  <- saved$stan_data

base_vars  <- c("mu", "sigma")
fiber_vars <- paste0("alpha[", 1:8, "]")

# ── 1. Run both samplers ───────────────────────────────────────────────────────

N_WARMUP <- 1000L
N_ITER   <- 2000L
N_CHAINS <- 4L

cat("── Standard centred sampler (no horizontal correction) ──────────────────\n")
draws_std <- horizontal_mcmc(
  stan_data      = stan_data,
  n_iter         = N_ITER,
  n_warmup       = N_WARMUP,
  n_chains       = N_CHAINS,
  use_correction = FALSE,
  seed           = 42L,
  verbose        = TRUE
)

cat("\n── Horizontal-corrected sampler ─────────────────────────────────────────\n")
draws_hor <- horizontal_mcmc(
  stan_data      = stan_data,
  n_iter         = N_ITER,
  n_warmup       = N_WARMUP,
  n_chains       = N_CHAINS,
  use_correction = TRUE,
  seed           = 42L,
  verbose        = TRUE
)

# ── 2. ESS comparison ─────────────────────────────────────────────────────────

cat("\n── ESS comparison ───────────────────────────────────────────────────────\n")

ess_std <- summarise_draws(draws_std, "ess_bulk", "ess_tail", "rhat")
ess_hor <- summarise_draws(draws_hor, "ess_bulk", "ess_tail", "rhat")

target_vars <- c("mu", "sigma", fiber_vars, "beta[1]", "beta[2]")

ess_std_sub <- ess_std[ess_std$variable %in% target_vars, ]
ess_hor_sub <- ess_hor[ess_hor$variable %in% target_vars, ]

nominal <- N_ITER * N_CHAINS

df_ess <- data.frame(
  variable  = ess_std_sub$variable,
  ESS_std   = round(ess_std_sub$ess_bulk),
  pct_std   = round(100 * ess_std_sub$ess_bulk / nominal, 1),
  ESS_hor   = round(ess_hor_sub$ess_bulk),
  pct_hor   = round(100 * ess_hor_sub$ess_bulk / nominal, 1),
  gain      = round(ess_hor_sub$ess_bulk / ess_std_sub$ess_bulk, 2)
)
print(df_ess, row.names = FALSE)

# Save draws for downstream analysis
saveRDS(draws_std, file.path(out_dir, "glmm_sparse_rwm_std.rds"))
saveRDS(draws_hor, file.path(out_dir, "glmm_sparse_rwm_hor.rds"))

# ── 3. ESS bar chart ──────────────────────────────────────────────────────────

df_plot <- rbind(
  data.frame(variable = ess_std_sub$variable,
             ess      = ess_std_sub$ess_bulk / nominal * 100,
             sampler  = "Standard (no correction)"),
  data.frame(variable = ess_hor_sub$variable,
             ess      = ess_hor_sub$ess_bulk / nominal * 100,
             sampler  = "Horizontal correction")
)
df_plot$variable <- factor(df_plot$variable, levels = target_vars)

p_ess <- ggplot(df_plot, aes(x = variable, y = ess, fill = sampler)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey40") +
  scale_fill_manual(values = c("Standard (no correction)" = "steelblue",
                               "Horizontal correction"    = "firebrick")) +
  labs(
    title    = "ESS (% of nominal) — standard vs horizontal-corrected sampler",
    subtitle = sprintf("Sparse GLMM  |  J=8, n_j=3, sigma_true=3  |  %d iterations x %d chains",
                       N_ITER, N_CHAINS),
    x        = NULL,
    y        = "ESS bulk (% of nominal)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "ess_comparison.png"),
       plot = p_ess, width = 9, height = 5, dpi = 150)
cat("Saved: data-raw/ess_comparison.png\n")

# ── 4. Trace plots for worst alpha group ──────────────────────────────────────

# Pick the alpha variable with lowest ESS in the standard sampler
worst_var <- ess_std_sub$variable[which.min(ess_std_sub$ess_bulk)]
cat(sprintf("\nTrace for %s (lowest ESS in standard sampler):\n", worst_var))

# Flatten to data frame for ggplot (first 2 chains only for clarity)
.draws_to_df <- function(draws, var, n_chain = 2L, label) {
  arr <- as_draws_matrix(subset_draws(draws, variable = var))
  do.call(rbind, lapply(seq_len(min(n_chain, dim(draws)[2L])), function(c) {
    data.frame(
      iteration = seq_len(nrow(arr)),
      value     = as.vector(
        as_draws_matrix(subset_draws(draws, variable = var, chain = c))
      ),
      chain   = sprintf("chain %d", c),
      sampler = label
    )
  }))
}

trace_std <- .draws_to_df(draws_std, worst_var, label = "Standard")
trace_hor <- .draws_to_df(draws_hor, worst_var, label = "Horizontal")
trace_all <- rbind(trace_std, trace_hor)

p_trace <- ggplot(trace_all, aes(x = iteration, y = value, colour = chain)) +
  geom_line(alpha = 0.8, linewidth = 0.4) +
  facet_wrap(~ sampler, ncol = 1L) +
  labs(
    title    = sprintf("Trace plot — %s", worst_var),
    subtitle = "Horizontal correction improves autocorrelation structure",
    x        = "Iteration (post-warmup)",
    y        = worst_var,
    colour   = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "trace_comparison.png"),
       plot = p_trace, width = 8, height = 5, dpi = 150)
cat("Saved: data-raw/trace_comparison.png\n")

cat("\nDone.\n")
