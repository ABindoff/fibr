## Attribution analysis: is the transport signal geometric or just
## autocorrelation?  (Manuscript Section "Attribution".)
##
## PRIMARY ANALYSIS (restricted): short loops (length ≤ 2*MIN_GAP) in the
## tightest distance quartile; Stokes prediction (exp(F_mean*area)-1)*alpha_start;
## displacements in units of per-group posterior SD.
##
## SECONDARY ANALYSIS (full): all loops, full path-integral prediction, raw units.
## Kept as a check; null reported in Section 4.7 is from the full analysis.
##
## Decision rule (HANDOFF_revision2.md):
##   Restricted/Stokes slopes > 0 -> signal detectable for tight short loops only;
##                                    report as consistent with footprint framing.
##   Still null -> add caveat (test harshness) + cite control-pair result.
##
## Writes:
##   data-raw/attribution_results.rds
##   data-raw/attribution_scatter.png   (primary: Stokes, scaled, restricted)
##   data-raw/attribution_area.png      (primary: area vs emp displacement)
##   data-raw/attribution_scatter_full.png   (secondary: original, all loops)
##
## Requires data-raw/glmm_sparse_data.rds and glmm_sparse_draws.rds.
## Run from package root:  Rscript data-raw/run_attribution.R

library(posterior)
library(ggplot2)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir <- file.path(pkg_root, "data-raw")

MIN_GAP     <- 25L    # past the base IACT so generic autocorrelation is weak
N_LOOPS_INT <- 2000L  # integrate more loops so the restricted filter has enough

# ── 1. Load sparse benchmark ──────────────────────────────────────────────────

saved     <- readRDS(file.path(out_dir, "glmm_sparse_data.rds"))
stan_data <- saved$stan_data
J         <- stan_data$J

draws <- readRDS(file.path(out_dir, "glmm_sparse_draws.rds"))

BASE_VARS  <- c("mu", "sigma")
FIBER_VARS <- paste0("alpha[", seq_len(J), "]")

# ── 2. Per-group posterior SD (for units scaling) ─────────────────────────────

chains_list <- .split_chains(draws)
full_mat    <- do.call(rbind, chains_list)
post_sd     <- apply(full_mat[, FIBER_VARS, drop = FALSE], 2L, sd)
names(post_sd) <- FIBER_VARS

cat("Per-group posterior SD:\n"); print(round(post_sd, 3))

# ── 3. Detect loops and estimate per-group transport factors ─────────────────

hd <- holonomy_diagnostic(
  chain             = draws,
  base_vars         = BASE_VARS,
  fiber_vars        = FIBER_VARS,
  min_gap           = MIN_GAP,
  n_bootstrap       = 200L,
  structure         = "diagonal",
  residualize_fiber = TRUE
)
print(hd)

# ── 4. Analytic connection and integrated transport along each loop ──────────

conn <- compute_connection(
  chain       = draws,
  base_vars   = BASE_VARS,
  fiber_vars  = FIBER_VARS,
  method      = "analytic_glmm",
  stan_data   = stan_data,
  n_subsample = 500L
)

# Integrate along tightest N_LOOPS_INT loops (sorted by distance inside fn)
tr <- integrate_transport(conn, hd$loops, draws, n_loops = N_LOOPS_INT)
cat(sprintf("Integrated %d loops (min_gap = %d).\n", nrow(tr), MIN_GAP))

# ── 5. Build long format: one row per loop x group ───────────────────────────

safe_lab <- gsub("[^[:alnum:]]", "_", FIBER_VARS)

long_rows <- lapply(seq_len(J), function(j) {
  jl <- safe_lab[j]
  data.frame(
    group         = j,
    loop          = seq_len(nrow(tr)),
    area          = tr$area,
    length        = tr$length,
    distance      = tr$distance,
    # Full path-integral prediction (secondary analysis)
    pred_full     = tr[[paste0("transport_error_", jl)]],
    # Stokes prediction: (exp(F_mean*area) - 1) * alpha_start
    pred_stokes   = (tr[[paste0("holonomy_stokes_", jl)]] - 1) *
                      tr[[paste0("alpha_start_", jl)]],
    emp           = tr[[paste0("empirical_disp_",  jl)]],
    post_sd       = post_sd[j]
  )
})
long_df <- do.call(rbind, long_rows)
long_df$group_f <- factor(long_df$group)

# Scaled versions (units of per-group posterior SD)
long_df$emp_sc         <- long_df$emp         / long_df$post_sd
long_df$pred_stokes_sc <- long_df$pred_stokes / long_df$post_sd
long_df$pred_full_sc   <- long_df$pred_full   / long_df$post_sd

# ── 6. PRIMARY ANALYSIS: restricted loops ─────────────────────────────────────
# (a) path length <= 2 * MIN_GAP
# (b) distance in tightest quartile of the INTEGRATED loops

dist_q25    <- quantile(tr$distance, 0.25)
max_length  <- 2L * MIN_GAP
restr_df <- long_df[
  long_df$length   <= max_length &
  long_df$distance <= dist_q25,
]
cat(sprintf(
  "\nRestricted filter: length <= %d AND distance <= %.4f  ->  %d loops.\n",
  max_length, dist_q25, length(unique(restr_df$loop))
))

if (nrow(restr_df) < J * 5L) {
  warning("Very few loops after restriction — primary test unreliable.")
}

cat("\n── PRIMARY: Stokes prediction test (restricted loops, scaled units) ──\n")
primary_test <- do.call(rbind, lapply(seq_len(J), function(j) {
  d <- restr_df[restr_df$group == j, ]
  if (nrow(d) < 10L) {
    return(data.frame(group=j, n=nrow(d), slope=NA, se=NA, t=NA, p=NA, cor=NA))
  }
  fit <- lm(emp_sc ~ pred_stokes_sc, data = d)
  cf  <- summary(fit)$coefficients
  data.frame(
    group = j,
    n     = nrow(d),
    slope = cf["pred_stokes_sc", "Estimate"],
    se    = cf["pred_stokes_sc", "Std. Error"],
    t     = cf["pred_stokes_sc", "t value"],
    p     = cf["pred_stokes_sc", "Pr(>|t|)"],
    cor   = cor(d$pred_stokes_sc, d$emp_sc, use = "complete.obs"),
    h_j   = Re(hd$eigenvalues)[j]
  )
}))
print(round(primary_test, 4), row.names = FALSE)

fit_pooled_primary <- lm(emp_sc ~ pred_stokes_sc, data = restr_df)
cat(sprintf("\nPooled slope (primary): %.3f (SE %.3f, p = %.4f)\n",
            coef(fit_pooled_primary)[2L],
            summary(fit_pooled_primary)$coefficients["pred_stokes_sc", "Std. Error"],
            summary(fit_pooled_primary)$coefficients["pred_stokes_sc", "Pr(>|t|)"]))

# ── 7. SECONDARY ANALYSIS: full loops, original prediction (original result) ──

cat("\n── SECONDARY: full path-integral prediction (all loops, raw units) ──\n")
secondary_test <- do.call(rbind, lapply(seq_len(J), function(j) {
  d   <- long_df[long_df$group == j, ]
  fit <- lm(emp ~ pred_full, data = d)
  cf  <- summary(fit)$coefficients
  data.frame(
    group = j,
    n     = nrow(d),
    slope = cf["pred_full", "Estimate"],
    se    = cf["pred_full", "Std. Error"],
    t     = cf["pred_full", "t value"],
    p     = cf["pred_full", "Pr(>|t|)"],
    cor   = cor(d$pred_full, d$emp, use = "complete.obs"),
    h_j   = Re(hd$eigenvalues)[j]
  )
}))
print(round(secondary_test, 4), row.names = FALSE)

fit_pooled_sec <- lm(emp ~ pred_full, data = long_df)
cat(sprintf("\nPooled slope (secondary, full): %.5f (SE %.5f)\n",
            coef(fit_pooled_sec)[2L],
            summary(fit_pooled_sec)$coefficients["pred_full", "Std. Error"]))

# ── 8. Orientation test: emp vs signed area (all loops) ───────────────────────

cat("\n── Orientation test: emp ~ area, per group ────────────────────────\n")
area_test <- do.call(rbind, lapply(seq_len(J), function(j) {
  d   <- long_df[long_df$group == j, ]
  fit <- lm(emp_sc ~ area, data = d)
  cf  <- summary(fit)$coefficients
  data.frame(
    group      = j,
    slope_area = cf["area", "Estimate"],
    se         = cf["area", "Std. Error"],
    p          = cf["area", "Pr(>|t|)"],
    F_j_mean   = mean(conn$curvature[, j], na.rm = TRUE)
  )
}))
print(round(area_test, 4), row.names = FALSE)

saveRDS(
  list(
    primary   = primary_test,
    secondary = secondary_test,
    orientation = area_test,
    min_gap     = MIN_GAP,
    n_all       = length(unique(long_df$loop)),
    n_restricted = length(unique(restr_df$loop)),
    max_length  = max_length,
    dist_q25    = dist_q25
  ),
  file.path(out_dir, "attribution_results.rds")
)

# ── 9. Figures ────────────────────────────────────────────────────────────────

# Primary: Stokes prediction vs empirical (restricted, scaled)
fig_primary <- ggplot(restr_df, aes(x = pred_stokes_sc, y = emp_sc)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(alpha = 0.35, size = 0.9, colour = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "firebrick",
              linewidth = 0.7, se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50") +
  facet_wrap(~ group_f, ncol = 4L, scales = "free",
             labeller = label_both) +
  labs(
    title    = "Primary attribution: Stokes prediction vs empirical displacement",
    subtitle = sprintf(
      paste0("Restricted loops (length ≤ %d, tightest distance quartile):  ",
             "%d loops  |  units = posterior SD"),
      max_length, length(unique(restr_df$loop))
    ),
    x = expression(Delta * alpha[j]^{Stokes} / SD[j] ~ "(Stokes prediction, scaled)"),
    y = expression(Delta * alpha[j]^{emp} / SD[j])
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "attribution_scatter.png"),
       plot = fig_primary, width = 10, height = 6, dpi = 150)
cat("Saved: data-raw/attribution_scatter.png  (primary/Stokes/restricted)\n")

# Secondary: full path-integral prediction, all loops, raw units
fig_secondary <- ggplot(long_df, aes(x = pred_full, y = emp)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(alpha = 0.2, size = 0.7, colour = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "firebrick",
              linewidth = 0.7, se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50") +
  facet_wrap(~ group_f, ncol = 4L, scales = "free",
             labeller = label_both) +
  labs(
    title    = "Secondary: full path-integral prediction vs empirical displacement",
    subtitle = sprintf(
      "All %d loops (null result: slope ~0)", length(unique(long_df$loop))
    ),
    x = expression(Delta * alpha[j]^{pred} ~ "(integrated connection)"),
    y = expression(Delta * alpha[j]^{emp})
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "attribution_scatter_full.png"),
       plot = fig_secondary, width = 10, height = 6, dpi = 150)
cat("Saved: data-raw/attribution_scatter_full.png  (secondary/full)\n")

# Orientation test: emp_sc vs signed area, all loops
fig_area <- ggplot(long_df, aes(x = area, y = emp_sc)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(alpha = 0.25, size = 0.9, colour = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "firebrick",
              linewidth = 0.7, se = TRUE) +
  facet_wrap(~ group_f, ncol = 4L, scales = "free",
             labeller = label_both) +
  labs(
    title    = "Orientation test: displacement vs signed loop area",
    subtitle = "Holonomy is orientation-aware; autocorrelation is not  |  units = posterior SD",
    x = "Signed enclosed area (shoelace, base space)",
    y = expression(Delta * alpha[j]^{emp} / SD[j])
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "attribution_area.png"),
       plot = fig_area, width = 10, height = 6, dpi = 150)
cat("Saved: data-raw/attribution_area.png\n")

cat("\nDone.\n")
