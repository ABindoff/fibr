## Negative control: base-closure conditioning vs matched-gap non-loop pairs.
##
## Purpose: determine whether the h_j signal in true loops is explained by
## lag autocorrelation + residualisation leakage alone, or whether
## base-closure conditioning (requiring theta_e ~ theta_s) adds signal.
##
## For each dataset and gap:
##   TRUE LOOPS:    detect_loops() as usual; estimate diagonal h_j
##   CONTROL PAIRS: same K pairs, gap distribution matched to true loops,
##                  but base distance LARGE (above the median of all pairs at
##                  that gap). Uniform weights (distance-weighting is
##                  meaningless when we are deliberately choosing far pairs).
##
## Decision rule:
##   h_true >> h_control (prior-dominated cell) -> base-closure carries signal
##   h_true ~ h_control               -> lag autocorrelation / leakage only
##
## Writes:
##   data-raw/control_pairs.rds    data frame (cell, gap, group, h_true, h_control)
##   data-raw/control_pairs.png    scatter h_true vs h_control, faceted by cell x gap
##
## Run from package root:  Rscript data-raw/run_control_pairs.R

library(posterior)
library(ggplot2)

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

out_dir  <- file.path(pkg_root, "data-raw")
draw_dir <- file.path(out_dir, "simstud_draws")

GAP_GRID    <- c(3L, 10L, 25L, 50L)
N_BOOTSTRAP <- 100L
SEED        <- 42L

BASE_VARS  <- c("mu", "sigma")
J          <- 8L
FIBER_VARS <- paste0("alpha[", seq_len(J), "]")

# ── Datasets ------------------------------------------------------------------

datasets <- list(
  list(
    label = "sparse_benchmark",
    file  = file.path(out_dir, "glmm_sparse_draws.rds"),
    fiber_vars = FIBER_VARS
  ),
  list(
    label = "nj3_s0.5",
    file  = file.path(draw_dir, "nj3_s0.5_r01.rds"),
    fiber_vars = FIBER_VARS
  ),
  list(
    label = "nj3_s3.0",
    file  = file.path(draw_dir, "nj3_s3.0_r01.rds"),
    fiber_vars = FIBER_VARS
  ),
  list(
    label = "nj100_s0.5",
    file  = file.path(draw_dir, "nj100_s0.5_r01.rds"),
    fiber_vars = FIBER_VARS
  )
)

# ── Helper: build control pairs -----------------------------------------------
# For each unique gap value in sampled_gaps, collect all within-chain pairs
# at that lag, keep those with base distance above dist_threshold, then
# sample n_needed from the far set.
#
# Arguments:
#   base_std     : standardised base matrix (n_total x K_base)
#   sampled_gaps : integer vector of K gap values drawn from true-loop gaps
#   chain_starts, chain_ends : integer vectors of per-chain boundaries
#                              (1-indexed, inclusive)
#   seed         : optional RNG seed

.build_control_pairs <- function(base_std, sampled_gaps,
                                  chain_starts, chain_ends,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  K    <- length(sampled_gaps)
  ugs  <- sort(unique(sampled_gaps))

  # Estimate median pairwise distance (pooled over gaps/chains, subsampled)
  samp_dists <- numeric(0)
  for (g in ugs) {
    for (ci in seq_along(chain_starts)) {
      cs <- chain_starts[ci]; ce <- chain_ends[ci]
      if (ce - cs < g) next
      idx_s <- seq(cs, ce - g)
      n_samp <- min(200L, length(idx_s))
      idx_s2 <- idx_s[sample.int(length(idx_s), n_samp)]
      d <- sqrt(rowSums(
        (base_std[idx_s2 + g, , drop = FALSE] -
           base_std[idx_s2,   , drop = FALSE])^2
      ))
      samp_dists <- c(samp_dists, d)
    }
  }
  dist_threshold <- median(samp_dists, na.rm = TRUE)

  # Build the control data frame gap by gap
  rows <- vector("list", length(ugs))
  for (gi in seq_along(ugs)) {
    g        <- ugs[gi]
    n_needed <- sum(sampled_gaps == g)

    # Collect all far pairs across all chains at this gap
    far_s <- integer(0); far_d <- numeric(0)
    for (ci in seq_along(chain_starts)) {
      cs <- chain_starts[ci]; ce <- chain_ends[ci]
      if (ce - cs < g) next
      idx_s <- seq(cs, ce - g)
      d <- sqrt(rowSums(
        (base_std[idx_s + g, , drop = FALSE] -
           base_std[idx_s,   , drop = FALSE])^2
      ))
      mask <- d > dist_threshold
      far_s <- c(far_s, idx_s[mask])
      far_d <- c(far_d, d[mask])
    }

    # Fall back to top 50% if threshold too strict
    if (length(far_s) == 0L) {
      all_s <- integer(0); all_d <- numeric(0)
      for (ci in seq_along(chain_starts)) {
        cs <- chain_starts[ci]; ce <- chain_ends[ci]
        if (ce - cs < g) next
        idx_s <- seq(cs, ce - g)
        d <- sqrt(rowSums(
          (base_std[idx_s + g, , drop = FALSE] -
             base_std[idx_s,   , drop = FALSE])^2
        ))
        all_s <- c(all_s, idx_s); all_d <- c(all_d, d)
      }
      half <- order(all_d, decreasing = TRUE)[seq_len(max(1L, length(all_d) %/% 2L))]
      far_s <- all_s[half]; far_d <- all_d[half]
    }

    chosen <- sample.int(length(far_s), n_needed, replace = length(far_s) < n_needed)
    rows[[gi]] <- data.frame(
      start    = far_s[chosen],
      end      = far_s[chosen] + g,
      distance = far_d[chosen]
    )
  }
  do.call(rbind, rows)
}

# ── Main loop -----------------------------------------------------------------

results <- vector("list", length(datasets) * length(GAP_GRID))
k_row <- 0L

for (ds in datasets) {
  cat(sprintf("\n══ %s ══\n", ds$label))

  if (!file.exists(ds$file)) {
    cat(sprintf("  File not found: %s  — skipping.\n", ds$file))
    next
  }

  draws <- readRDS(ds$file)
  chains_list <- .split_chains(draws)
  chain_sizes <- vapply(chains_list, nrow, integer(1L))

  full_mat   <- do.call(rbind, chains_list)
  full_base  <- full_mat[, BASE_VARS, drop = FALSE]
  full_fiber <- full_mat[, ds$fiber_vars, drop = FALSE]

  # Residualise once: shared by both true-loop and control estimation
  fiber_resid <- .residualize(full_fiber, full_base)

  # Standardise base for distance computation
  base_std <- scale(full_base)

  # Chain boundary indices (1-indexed, inclusive)
  chain_ends   <- cumsum(chain_sizes)
  chain_starts <- c(1L, chain_ends[-length(chain_ends)] + 1L)

  for (gap in GAP_GRID) {
    cat(sprintf("  gap = %d\n", gap))

    # ── True loops ------------------------------------------------------------
    hd <- tryCatch(
      suppressMessages(holonomy_diagnostic(
        chain             = draws,
        base_vars         = BASE_VARS,
        fiber_vars        = ds$fiber_vars,
        min_gap           = gap,
        n_bootstrap       = N_BOOTSTRAP,
        k                 = 100L,
        max_loops         = 5000L,
        structure         = "diagonal",
        weights           = "distance",
        residualize_fiber = TRUE
      )),
      error = function(e) {
        cat(sprintf("    holonomy_diagnostic error: %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(hd)) next

    h_true  <- Re(hd$eigenvalues)   # J-vector, group-aligned
    K_loops <- nrow(hd$loops)
    true_gaps_vec <- hd$loops$end - hd$loops$start

    cat(sprintf("    true loops: %d  |  mean h: %.4f\n",
                K_loops, mean(h_true)))

    # ── Control pairs --------------------------------------------------------
    set.seed(SEED + gap)
    sampled_gaps <- sample(true_gaps_vec, K_loops, replace = TRUE)

    ctrl_loops <- tryCatch(
      .build_control_pairs(base_std, sampled_gaps,
                            chain_starts, chain_ends,
                            seed = SEED + gap + 1L),
      error = function(e) {
        cat(sprintf("    control-pair error: %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(ctrl_loops) || nrow(ctrl_loops) == 0L) next

    tm_ctrl <- tryCatch(
      estimate_transport_map(
        fiber_draws = fiber_resid,
        loops       = ctrl_loops,
        n_bootstrap = 0L,        # no bootstrap needed for controls
        structure   = "diagonal",
        weights     = "uniform"
      ),
      error = function(e) {
        cat(sprintf("    estimate_transport_map (control) error: %s\n",
                    conditionMessage(e)))
        NULL
      }
    )
    if (is.null(tm_ctrl)) next

    h_ctrl <- Re(tm_ctrl$eigenvalues)

    cat(sprintf("    control pairs: %d  |  mean h: %.4f\n",
                nrow(ctrl_loops), mean(h_ctrl)))

    k_row <- k_row + 1L
    results[[k_row]] <- data.frame(
      cell      = ds$label,
      gap       = gap,
      group     = seq_len(J),
      h_true    = h_true,
      h_control = h_ctrl
    )
  }
}

results <- do.call(rbind, results[seq_len(k_row)])
saveRDS(results, file.path(out_dir, "control_pairs.rds"))
cat(sprintf("\nSaved: data-raw/control_pairs.rds  (%d rows)\n", nrow(results)))

# ── Figure --------------------------------------------------------------------

results$gap_f  <- factor(sprintf("gap = %d",  results$gap),
                          levels = sprintf("gap = %d", GAP_GRID))
results$cell_f <- factor(results$cell, levels = vapply(datasets, `[[`, "", "label"))

fig <- ggplot(results, aes(x = h_control, y = h_true, colour = factor(group))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 2.5, alpha = 0.85) +
  facet_grid(cell_f ~ gap_f, scales = "free") +
  scale_colour_viridis_d(name = "Group") +
  labs(
    title    = "Negative control: true loops vs matched-gap far pairs",
    subtitle = paste0(
      "Points above the dashed y = x line indicate base-closure carries signal",
      " beyond lag autocorrelation"
    ),
    x = expression(h[j] ~ "(control: large base distance)"),
    y = expression(h[j] ~ "(true loops: small base distance)")
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", strip.text = element_text(size = 8))

ggsave(file.path(out_dir, "control_pairs.png"),
       plot = fig, width = 12, height = 9, dpi = 150)
cat("Saved: data-raw/control_pairs.png\n")
cat("\nDone.\n")
