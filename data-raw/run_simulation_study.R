## Simulation study: holonomy diagnostic sensitivity across identifiability regimes.
##
## Tests the prediction that ||H - I||_F increases with prior dominance of the
## fiber parameters, as quantified by the analytic prior fraction
##   pi_j = (1/sigma^2) / G_FF_j.
##
## Design:
##   n_j        in {3, 10, 30, 100}       (observations per group)
##   sigma_true in {0.5, 1.0, 2.0, 3.0}  (true hierarchical SD)
##   J = 8 groups (fixed), N_REP = 10 replicates per cell
##
## For each cell-replicate:
##   1. Simulate GLMM data with beta_true = (0.8, -0.5)
##   2. Fit the centred Stan model (4 chains x 1000 warmup + 2000 sampling)
##   3. holonomy_diagnostic() evaluated at a grid of min_gap values (see
##      GAP_GRID below). Small gaps sit inside the autocorrelation time in
##      prior-dominated cells; reporting the full gap profile separates the
##      autocorrelation contribution from the geometric signal.
##   4. compute_connection() for the analytic prior fraction per group
##
## Results are checkpointed to disk after every replicate. Re-running the
## script resumes from the last completed cell.
##
## Outputs (written to data-raw/):
##   simstud_results.rds              scalar summary, one row per cell-replicate
##   simstud_pergroup.rds             per-group metrics, one row per cell x rep x group
##   simstud_hd_corners.rds           fibr_holonomy objects for four extremal cells (rep 1)
##   simstud_heatmap.png              Figure A: ||H-I||_F heatmap over design grid
##   simstud_theory_vs_empirical.png  Figure B: prior_frac vs ||H-I||_F scatter
##   simstud_eigenspectra.png         Figure C: eigenspectrum panels at four corners
##
## Estimated wall time: ~4 hours on a modern laptop (160 Stan fits,
## parallel_chains = 4 each). Set N_REP = 3 for a quick smoke test.
##
## Run from package root:  Rscript data-raw/run_simulation_study.R

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

out_dir <- file.path(pkg_root, "data-raw")

# ── Study parameters ──────────────────────────────────────────────────────────

NJ_GRID    <- c(3L, 10L, 30L, 100L)
SIGMA_GRID <- c(0.5, 1.0, 2.0, 3.0)
J          <- 8L
N_REP      <- 10L
MU_TRUE    <- 0
BETA_TRUE  <- c(0.8, -0.5)

BASE_VARS  <- c("mu", "sigma")
FIBER_VARS <- paste0("alpha[", seq_len(J), "]")
BETA_VARS  <- c("beta[1]", "beta[2]")

GAP_GRID    <- c(3L, 10L, 25L, 50L)  # loop gaps; profile separates
                                     # autocorrelation (decays with gap)
                                     # from geometric signal (persists)
N_BOOTSTRAP <- 100L  # transport-factor bootstrap resamples

# Draws are saved per cell-replicate so the diagnostic can be re-run
# post hoc (new gaps, estimators) without refitting Stan.
draws_dir <- file.path(out_dir, "simstud_draws")
dir.create(draws_dir, showWarnings = FALSE)

# Corner cells for the eigenspectrum figure (Figure C).
# Use replicate 1 only; listed in reading order (left→right, top→bottom).
CORNERS <- list(
  list(n_j = 3L,   sigma = 3.0, label = "High holonomy\n(n[j]==3, sigma==3)"),
  list(n_j = 3L,   sigma = 0.5, label = "Low prior dominance\n(n[j]==3, sigma==0.5)"),
  list(n_j = 100L, sigma = 3.0, label = "Wide prior, dense data\n(n[j]==100, sigma==3)"),
  list(n_j = 100L, sigma = 0.5, label = "Near-zero holonomy\n(n[j]==100, sigma==0.5)")
)

.is_corner <- function(n_j, sigma_t) {
  any(vapply(CORNERS,
             function(co) co$n_j == n_j && abs(co$sigma - sigma_t) < 1e-9,
             logical(1L)))
}

# ── Compile Stan model once ───────────────────────────────────────────────────

stan_file <- file.path(pkg_root, "inst", "stan", "glmm_centred.stan")
cat("Compiling Stan model...\n")
mod <- cmdstan_model(stan_file)
cat("Done.\n\n")

# ── Checkpoint: load previously completed results ─────────────────────────────

result_file <- file.path(out_dir, "simstud_results.rds")
pergrp_file <- file.path(out_dir, "simstud_pergroup.rds")
corner_file <- file.path(out_dir, "simstud_hd_corners.rds")

results_df <- if (file.exists(result_file)) readRDS(result_file) else NULL
pergrp_df  <- if (file.exists(pergrp_file)) readRDS(pergrp_file) else NULL
corner_hd  <- if (file.exists(corner_file)) readRDS(corner_file) else list()

# Schema guard: results from the old single-gap (g = 3), full-estimator run
# lack the `gap` column. Archive them and start fresh.
if (!is.null(results_df) && !"gap" %in% names(results_df)) {
  cat("Old-schema results found; archiving as *_v1_gap3.rds and starting fresh.\n")
  file.rename(result_file, sub("\\.rds$", "_v1_gap3.rds", result_file))
  if (file.exists(pergrp_file))
    file.rename(pergrp_file, sub("\\.rds$", "_v1_gap3.rds", pergrp_file))
  if (file.exists(corner_file))
    file.rename(corner_file, sub("\\.rds$", "_v1_gap3.rds", corner_file))
  results_df <- NULL
  pergrp_df  <- NULL
  corner_hd  <- list()
}

.already_done <- function(n_j, sigma_t, r_idx) {
  if (is.null(results_df) || nrow(results_df) == 0L) return(FALSE)
  any(results_df$n_j == n_j &
      abs(results_df$sigma_true - sigma_t) < 1e-9 &
      results_df$replicate == r_idx)
}

total_cells <- length(NJ_GRID) * length(SIGMA_GRID) * N_REP
done_count  <- if (!is.null(results_df)) nrow(results_df) else 0L
cat(sprintf("Simulation study: %d cell-replicates total, %d already complete.\n\n",
            total_cells, done_count))

# ── Main grid loop ────────────────────────────────────────────────────────────

for (n_j in NJ_GRID) {
  for (sigma_true in SIGMA_GRID) {
    for (r_idx in seq_len(N_REP)) {

      if (.already_done(n_j, sigma_true, r_idx)) {
        cat(sprintf("  skip  n_j=%3d  sigma=%.1f  rep=%02d\n",
                    n_j, sigma_true, r_idx))
        next
      }

      cat(sprintf(
        "\n── n_j=%3d  sigma_true=%.1f  rep=%02d ──────────────────────────────\n",
        n_j, sigma_true, r_idx
      ))

      # Deterministic unique seed: encodes cell position and replicate.
      nj_idx    <- which(NJ_GRID    == n_j)
      sig_idx   <- which(abs(SIGMA_GRID - sigma_true) < 1e-9)
      cell_seed <- nj_idx * 10000L + sig_idx * 1000L + r_idx

      # ── 1. Simulate GLMM data ──────────────────────────────────────────────
      set.seed(cell_seed)
      N          <- J * n_j
      alpha_true <- rnorm(J,  MU_TRUE,   sigma_true)
      group_id   <- rep(seq_len(J), each = n_j)
      X          <- matrix(rnorm(N * 2L), ncol = 2L)
      eta        <- alpha_true[group_id] + X %*% BETA_TRUE
      y          <- rbinom(N, 1L, plogis(eta))
      stan_data  <- list(N = N, J = J, group = group_id, X = X, y = y)

      # ── 2. Fit centred Stan model ──────────────────────────────────────────
      fit <- tryCatch(
        mod$sample(
          data            = stan_data,
          chains          = 4L,
          parallel_chains = 4L,
          iter_warmup     = 1000L,
          iter_sampling   = 2000L,
          seed            = cell_seed,
          refresh         = 0L
        ),
        error = function(e) {
          cat(sprintf("    Stan error: %s\n", conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit)) next

      invisible(suppressMessages(capture.output(
        diag_summ <- fit$diagnostic_summary()
      )))
      n_divergent <- sum(diag_summ$num_divergent)
      draws       <- fit$draws()
      ess_tbl     <- summarise_draws(draws, "ess_bulk")
      min_ess     <- min(ess_tbl$ess_bulk, na.rm = TRUE)

      # Integrated autocorrelation time proxy per variable block:
      #   IACT ~= n_draws_total / ess_bulk
      n_draws_tot <- prod(dim(draws)[1:2])
      iact_of <- function(vars) {
        e <- ess_tbl$ess_bulk[ess_tbl$variable %in% vars]
        max(n_draws_tot / e, na.rm = TRUE)
      }
      iact_base  <- iact_of(BASE_VARS)
      iact_fiber <- iact_of(FIBER_VARS)

      cat(sprintf("    divergences: %d  |  min ESS: %.0f  |  IACT base/fiber: %.1f / %.1f\n",
                  n_divergent, min_ess, iact_base, iact_fiber))

      # Save draws for post-hoc re-analysis (new gaps/estimators, no refit)
      saveRDS(draws, file.path(draws_dir,
        sprintf("nj%d_s%.1f_r%02d.rds", n_j, sigma_true, r_idx)))

      # ── 3. Holonomy diagnostic over the gap grid ──────────────────────────
      # Diagonal estimator: per-group transport factors h_j (group-aligned).
      # Primary scalar: mean h = mean_j h_j, per gap.
      #   ~0  when loop endpoints are independent (mixed at loop time-scale)
      #   ->1 with loop-conditional persistence
      hd_by_gap <- list()
      for (gap in GAP_GRID) {
        hd_by_gap[[as.character(gap)]] <- tryCatch(
          suppressMessages(holonomy_diagnostic(
            chain             = draws,
            base_vars         = BASE_VARS,
            fiber_vars        = FIBER_VARS,
            epsilon           = NULL,
            n_bootstrap       = N_BOOTSTRAP,
            min_gap           = gap,
            k                 = 100L,
            max_loops         = 5000L,
            structure         = "diagonal",
            residualize_fiber = TRUE
          )),
          error = function(e) {
            cat(sprintf("    holonomy_diagnostic (gap %d): %s\n",
                        gap, conditionMessage(e)))
            NULL
          }
        )
      }

      # Save the full hd objects for Figure C (corner cells, rep 1 only).
      if (r_idx == 1L && .is_corner(n_j, sigma_true)) {
        corner_key <- sprintf("nj%d_s%.1f", n_j, sigma_true)
        corner_hd[[corner_key]] <- hd_by_gap
        saveRDS(corner_hd, corner_file)
      }

      # ── 4. Analytic connection: prior fraction ─────────────────────────────
      conn <- tryCatch(
        suppressMessages(compute_connection(
          chain       = draws,
          base_vars   = BASE_VARS,
          fiber_vars  = FIBER_VARS,
          method      = "analytic_glmm",
          stan_data   = stan_data,
          beta_vars   = BETA_VARS,
          n_subsample = 200L
        )),
        error = function(e) NULL
      )

      mean_prior_frac <- if (!is.null(conn)) {
        mean(colMeans(conn$prior_frac))
      } else NA_real_

      # ── 5. Record and checkpoint (one row per gap) ─────────────────────────
      for (gap in GAP_GRID) {
        hd <- hd_by_gap[[as.character(gap)]]
        h_j      <- if (!is.null(hd)) Re(hd$eigenvalues) else rep(NA_real_, J)
        mean_h   <- if (!is.null(hd)) mean(h_j)           else NA_real_
        h_max    <- if (!is.null(hd)) max(h_j)            else NA_real_
        frob_dev <- if (!is.null(hd)) hd$frobenius_dev    else NA_real_
        n_loops  <- if (!is.null(hd)) hd$n_loops          else NA_integer_

        cat(sprintf("    gap %2d: loops %5s  |  mean h: %s\n", gap,
                    if (is.na(n_loops)) "NA" else as.character(n_loops),
                    if (is.na(mean_h)) "NA" else sprintf("%.4f", mean_h)))

        new_row <- data.frame(
          n_j             = n_j,
          sigma_true      = sigma_true,
          replicate       = r_idx,
          gap             = gap,
          seed            = cell_seed,
          n_loops         = n_loops,
          mean_h          = mean_h,
          h_max           = h_max,
          frobenius_dev   = frob_dev,
          mean_prior_frac = mean_prior_frac,
          iact_base       = iact_base,
          iact_fiber      = iact_fiber,
          n_divergent     = n_divergent,
          min_ess_bulk    = min_ess
        )
        results_df <- if (is.null(results_df)) new_row else rbind(results_df, new_row)

        if (!is.null(conn) && !is.null(hd)) {
          pf_j   <- colMeans(conn$prior_frac)
          grp_df <- data.frame(
            n_j        = n_j,
            sigma_true = sigma_true,
            replicate  = r_idx,
            gap        = gap,
            group      = seq_len(J),
            prior_frac = pf_j,
            h_j        = h_j
          )
          pergrp_df <- if (is.null(pergrp_df)) grp_df else rbind(pergrp_df, grp_df)
        }
      }
      saveRDS(results_df, result_file)
      if (!is.null(pergrp_df)) saveRDS(pergrp_df, pergrp_file)

    }  # r_idx
  }    # sigma_true
}      # n_j

cat("\n══════════════════════════════════════════════════════════════\n")
cat("All cell-replicates complete. Generating figures.\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Reload from disk (supports use as a standalone figure-regeneration script)
results_df <- readRDS(result_file)
pergrp_df  <- if (file.exists(pergrp_file)) readRDS(pergrp_file) else NULL
corner_hd  <- if (file.exists(corner_file)) readRDS(corner_file) else list()

# ── Summary table ─────────────────────────────────────────────────────────────

cat("Cell-level summary (median over replicates):\n\n")

summary_rows <- vector("list", length(NJ_GRID) * length(SIGMA_GRID))
k_row <- 0L
for (nj in NJ_GRID) {
  for (sig in SIGMA_GRID) {
    k_row <- k_row + 1L
    sub   <- results_df[results_df$n_j == nj &
                        abs(results_df$sigma_true - sig) < 1e-9, ]
    sub50 <- sub[sub$gap == 50L, ]   # headline gap: past the IACT in most cells
    summary_rows[[k_row]] <- data.frame(
      n_j              = nj,
      sigma_true       = sig,
      n_complete       = sum(!is.na(sub50$mean_h)),
      median_mean_h    = round(median(sub50$mean_h, na.rm = TRUE), 4),
      median_h_max     = round(median(sub50$h_max,  na.rm = TRUE), 4),
      median_pf        = round(median(sub50$mean_prior_frac, na.rm = TRUE), 3),
      median_iact_fib  = round(median(sub50$iact_fiber, na.rm = TRUE), 1),
      pct_div          = round(100 * mean(sub50$n_divergent > 0L, na.rm = TRUE))
    )
  }
}
print(do.call(rbind, summary_rows), row.names = FALSE)

# ── Figure A: mean transport factor heatmap, faceted by gap ──────────────────

agg_rows <- vector("list",
                   length(NJ_GRID) * length(SIGMA_GRID) * length(GAP_GRID))
k_row <- 0L
for (g in GAP_GRID) {
  for (nj in NJ_GRID) {
    for (sig in SIGMA_GRID) {
      k_row <- k_row + 1L
      sub   <- results_df[results_df$n_j == nj &
                          abs(results_df$sigma_true - sig) < 1e-9 &
                          results_df$gap == g, ]
      mh <- sub$mean_h
      agg_rows[[k_row]] <- data.frame(
        gap_f      = factor(sprintf("gap = %d", g),
                            levels = sprintf("gap = %d", GAP_GRID)),
        n_j_f      = factor(nj,  levels = rev(NJ_GRID)),
        sigma_f    = factor(sig, levels = SIGMA_GRID),
        med_h      = median(mh, na.rm = TRUE),
        q25_h      = quantile(mh, 0.25, na.rm = TRUE),
        q75_h      = quantile(mh, 0.75, na.rm = TRUE),
        n_complete = sum(!is.na(mh))
      )
    }
  }
}
agg_df <- do.call(rbind, agg_rows)

agg_df$tile_label <- sprintf(
  "%.3f\n[%.3f, %.3f]",
  agg_df$med_h, agg_df$q25_h, agg_df$q75_h
)

fig_a <- ggplot(agg_df, aes(x = sigma_f, y = n_j_f)) +
  geom_tile(aes(fill = med_h), colour = "white", linewidth = 0.9) +
  geom_text(aes(label = tile_label), size = 2.4, lineheight = 1.1,
            colour = "black") +
  facet_wrap(~ gap_f, ncol = 2L) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    name      = expression(bar(h)),
    na.value  = "grey85",
    limits    = c(0, max(agg_df$med_h, na.rm = TRUE))
  ) +
  scale_x_discrete(
    name   = expression(sigma[true] ~ "(hierarchical SD)"),
    labels = paste0("σ = ", SIGMA_GRID)
  ) +
  scale_y_discrete(
    name   = expression(n[j] ~ "(obs per group)"),
    labels = paste0("n[j] = ", rev(NJ_GRID))
  ) +
  labs(
    title    = expression("Mean transport factor" ~ bar(h) ~
                          "by prior dominance and loop gap"),
    subtitle = sprintf(
      "Centred GLMM  |  J = %d groups  |  %d replicates per cell  |  diagonal estimator",
      J, N_REP
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right",
    axis.text       = element_text(size = 9)
  )

ggsave(file.path(out_dir, "simstud_heatmap.png"),
       plot = fig_a, width = 10, height = 8, dpi = 150)
cat("Saved: data-raw/simstud_heatmap.png\n")

# ── Figure B: analytic prior fraction vs empirical mean transport factor ─────

fig_b_df <- results_df[
  !is.na(results_df$mean_h) & !is.na(results_df$mean_prior_frac), ]
fig_b_df$sigma_f <- factor(
  sprintf("sigma = %.1f", fig_b_df$sigma_true),
  levels = sprintf("sigma = %.1f", SIGMA_GRID)
)
fig_b_df$nj_f <- factor(
  as.character(fig_b_df$n_j),
  levels = as.character(NJ_GRID)
)
fig_b_df$gap_f <- factor(
  sprintf("gap = %d", fig_b_df$gap),
  levels = sprintf("gap = %d", GAP_GRID)
)

fig_b <- ggplot(fig_b_df,
                aes(x = mean_prior_frac, y = mean_h)) +
  geom_point(aes(colour = sigma_f, shape = nj_f), alpha = 0.75, size = 2.0) +
  geom_smooth(
    method = "loess", formula = y ~ x,
    colour = "grey30", fill = "grey80",
    linewidth = 0.8, alpha = 0.25, se = TRUE
  ) +
  facet_wrap(~ gap_f, ncol = 2L) +
  scale_colour_brewer(
    palette = "OrRd",
    name    = expression(sigma[true])
  ) +
  scale_shape_manual(
    values = c(`3` = 16L, `10` = 17L, `30` = 15L, `100` = 18L),
    name   = expression(n[j])
  ) +
  labs(
    title    = "Analytic prior fraction vs empirical transport factor, by gap",
    subtitle = "mean pi_j = (1/sigma^2) / G_FF_j  (analytic)  vs  mean h_j  (diagnostic)",
    x = expression(bar(pi) ~ "(mean prior fraction, analytic)"),
    y = expression(bar(h) ~ "(mean transport factor)")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "simstud_theory_vs_empirical.png"),
       plot = fig_b, width = 9, height = 7, dpi = 150)
cat("Saved: data-raw/simstud_theory_vs_empirical.png\n")

# ── Figure C: per-group transport profiles over gaps, four corner cells ──────
# Each panel: h_j (with bootstrap 90% intervals) against gap, one line per
# group. Autocorrelation decays with gap; a persistent component is the
# candidate geometric signal.

if (length(corner_hd) == 0L) {
  cat("No corner fibr_holonomy objects found; skipping Figure C.\n")
} else {
  prof_rows <- list()
  for (co in CORNERS) {
    key       <- sprintf("nj%d_s%.1f", co$n_j, co$sigma)
    hd_by_gap <- corner_hd[[key]]
    if (is.null(hd_by_gap)) next
    cell_lab <- sprintf("n[j] = %d,  sigma = %.1f", co$n_j, co$sigma)
    for (g in names(hd_by_gap)) {
      hd <- hd_by_gap[[g]]
      if (is.null(hd)) next
      h_j <- Re(hd$eigenvalues)
      ci  <- apply(Re(hd$boot_eigenvalues), 2L, quantile,
                   probs = c(0.05, 0.95), na.rm = TRUE)
      prof_rows[[length(prof_rows) + 1L]] <- data.frame(
        cell  = cell_lab,
        gap   = as.integer(g),
        group = factor(seq_along(h_j)),
        h     = h_j,
        lo    = ci[1L, ],
        hi    = ci[2L, ]
      )
    }
  }

  if (length(prof_rows) > 0L) {
    prof_df <- do.call(rbind, prof_rows)
    fig_c <- ggplot(prof_df, aes(x = gap, y = h, colour = group)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = group),
                  alpha = 0.12, colour = NA) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.8) +
      facet_wrap(~ cell, ncol = 2L) +
      scale_x_continuous(breaks = GAP_GRID) +
      scale_colour_viridis_d(name = "Group") +
      scale_fill_viridis_d(guide = "none") +
      labs(
        title    = "Per-group transport factors across loop gaps",
        subtitle = sprintf(
          "Centred GLMM  |  J = %d groups  |  diagonal estimator  |  replicate 1",
          J
        ),
        x = "Minimum loop gap g (iterations)",
        y = expression(h[j])
      ) +
      theme_minimal(base_size = 12)

    ggsave(file.path(out_dir, "simstud_gap_profiles.png"),
           plot = fig_c, width = 9, height = 7, dpi = 150)
    cat("Saved: data-raw/simstud_gap_profiles.png\n")
  }
}

cat("\nDone.\n")

