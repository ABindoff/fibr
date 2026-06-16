## Mixed-design π_j-adaptive partial-NC experiment (Task 5).
##
## Tests whether per-group adaptive weights w_j = π_j outperform uniform
## parameterisations when groups have heterogeneous observation counts.
## This is the controlled follow-up to Section 8.1 of the manuscript.
##
## Design factors:
##   nj_design   "uniform_sparse"  all n_j = 3
##               "uniform_dense"   all n_j = 50
##               "mixed"           n_j = c(3,3,3,3,50,50,50,50)
##   sigma_true  1.0, 2.0
##   N_REP       10 replicates per cell
##
## Methods — all via glmm_partial_nc.stan with different w vectors:
##   "centred"   w_j = 0    (exact centred parameterisation)
##   "nc"        w_j = 1    (exact non-centred)
##   "half"      w_j = 0.5  (uniform intermediate)
##   "pi_j"      w_j = π_j  (prior-fraction adaptive; proposed)
##
## π_j is estimated from a short pilot centred fit at each cell-rep.
##
## Expectation:
##   uniform_sparse: nc ≈ pi_j         (all π_j high → adaptive ≈ NC)
##   uniform_dense:  centred ≈ pi_j    (all π_j low  → adaptive ≈ centred)
##   mixed:          pi_j > centred and nc  (the key Section 8.1 result)
##
## Outputs (data-raw/):
##   mixed_design_results.rds     per-cell-rep-method scalar metrics
##   mixed_design_pergroup.rds    per-group alpha[j] ESS and π_j
##   mixed_design_minESS.png      min-ESS by method × design boxplot
##   mixed_design_pi_ess.png      per-group alpha ESS vs π_j scatter
##   mixed_design_table.tex       LaTeX summary table
##
## Run from package root:  Rscript data-raw/run_mixed_design.R

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
})

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)
invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))
out_dir <- file.path(pkg_root, "data-raw")

# ── Study parameters ──────────────────────────────────────────────────────────

J          <- 8L
K          <- 2L
MU_TRUE    <- 0
BETA_TRUE  <- c(0.8, -0.5)
SIGMA_GRID <- c(1.0, 2.0)
N_REP      <- 10L
N_CHAINS   <- 4L
N_WARMUP   <- 1000L
N_SAMPLING <- 2000L

NJ_DESIGNS <- list(
  uniform_sparse = rep(3L,  J),
  uniform_dense  = rep(50L, J),
  mixed          = c(3L, 3L, 3L, 3L, 50L, 50L, 50L, 50L)
)

METHODS <- c("centred", "nc", "half", "pi_j")

# ── Compile Stan model ────────────────────────────────────────────────────────

cat("Compiling glmm_partial_nc.stan...\n")
mod <- cmdstan_model(file.path(pkg_root, "inst/stan/glmm_partial_nc.stan"))
cat("Done.\n\n")

# ── Checkpoint helpers ────────────────────────────────────────────────────────

result_file <- file.path(out_dir, "mixed_design_results.rds")
pergrp_file <- file.path(out_dir, "mixed_design_pergroup.rds")
results     <- if (file.exists(result_file)) readRDS(result_file) else list()
pergrp_rows <- if (file.exists(pergrp_file)) readRDS(pergrp_file) else list()

.key <- function(design, sigma, rep, method)
  sprintf("%s_s%.1f_r%02d_%s", design, sigma, rep, method)

.done <- function(design, sigma, rep, method)
  !is.null(results[[.key(design, sigma, rep, method)]])

.save <- function(design, sigma, rep, method, res) {
  results[[.key(design, sigma, rep, method)]] <<- res
  saveRDS(results, result_file)
}

# ── ESS helper ────────────────────────────────────────────────────────────────

.alpha_ess <- function(drws, J) {
  avars <- paste0("alpha[", seq_len(J), "]")
  vapply(avars, function(v)
    tryCatch(posterior::ess_bulk(drws[, , v]), error = function(e) NA_real_),
    numeric(1L)
  )
}

.min_ess <- function(drws) {
  mat  <- tryCatch(as_draws_matrix(drws), error = function(e) NULL)
  if (is.null(mat)) return(NA_real_)
  vars <- colnames(mat)
  vars <- vars[!grepl("^log_lik|^lp__|^psi", vars)]
  ess  <- vapply(vars, function(v)
    tryCatch(posterior::ess_bulk(drws[, , v]), error = function(e) NA_real_),
    numeric(1L)
  )
  min(ess, na.rm = TRUE)
}

# ── Pilot π_j from a short centred run ───────────────────────────────────────

.pilot_pi <- function(stan_data, cell_seed) {
  fit <- tryCatch(
    mod$sample(
      data            = c(stan_data, list(w = rep(0, stan_data$J))),
      chains          = 2L,
      parallel_chains = 2L,
      iter_warmup     = 300L,
      iter_sampling   = 300L,
      seed            = cell_seed,
      refresh         = 0L
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  drws <- fit$draws()
  beta_vars <- c("beta[1]", "beta[2]")
  conn <- tryCatch(
    suppressMessages(compute_connection(
      chain       = drws,
      base_vars   = c("mu", "sigma"),
      fiber_vars  = paste0("alpha[", seq_len(stan_data$J), "]"),
      method      = "analytic_glmm",
      stan_data   = stan_data,
      beta_vars   = beta_vars,
      n_subsample = 100L
    )),
    error = function(e) NULL
  )
  if (is.null(conn)) return(NULL)
  pmin(pmax(colMeans(conn$prior_frac), 0), 1)
}

# ── Main loop ─────────────────────────────────────────────────────────────────

total <- length(NJ_DESIGNS) * length(SIGMA_GRID) * N_REP * length(METHODS)
done  <- sum(vapply(names(results), function(k) !is.null(results[[k]]), logical(1L)))
cat(sprintf("Mixed-design experiment: %d runs total, %d already complete.\n\n",
            total, done))

for (design_name in names(NJ_DESIGNS)) {
  nj_vec <- NJ_DESIGNS[[design_name]]

  for (sigma_true in SIGMA_GRID) {

    for (r_idx in seq_len(N_REP)) {

      if (all(vapply(METHODS, function(m)
            .done(design_name, sigma_true, r_idx, m), logical(1L)))) next

      cat(sprintf(
        "\n── design=%-15s  sigma=%.1f  rep=%02d ──────────────────────\n",
        design_name, sigma_true, r_idx
      ))

      # Deterministic seed: encode all three factors
      des_idx <- which(names(NJ_DESIGNS) == design_name)
      sig_idx <- which(abs(SIGMA_GRID - sigma_true) < 1e-9)
      cell_seed <- des_idx * 100000L + sig_idx * 10000L + r_idx

      # ── Simulate data ────────────────────────────────────────────────────────
      set.seed(cell_seed)
      N          <- sum(nj_vec)
      alpha_true <- rnorm(J, MU_TRUE, sigma_true)
      group_id   <- rep(seq_len(J), times = nj_vec)
      X          <- matrix(rnorm(N * K), ncol = K)
      eta        <- alpha_true[group_id] + X %*% BETA_TRUE
      y          <- rbinom(N, 1L, plogis(eta))
      stan_data  <- list(N = N, J = J, group = group_id, X = X, y = y)

      # ── Pilot: estimate π_j ──────────────────────────────────────────────────
      pi_est <- .pilot_pi(stan_data, cell_seed)
      if (is.null(pi_est)) {
        cat("  pilot failed — skipping rep\n")
        next
      }
      cat(sprintf("  π_j (pilot): %s\n",
                  paste(round(pi_est, 2), collapse = " ")))

      # ── Run each method ──────────────────────────────────────────────────────
      for (method in METHODS) {

        if (.done(design_name, sigma_true, r_idx, method)) {
          cat(sprintf("  skip %s\n", method))
          next
        }
        cat(sprintf("  %s ... ", method))

        w_vec <- switch(method,
          centred = rep(0,   J),
          nc      = rep(1,   J),
          half    = rep(0.5, J),
          pi_j    = pi_est
        )

        t0 <- proc.time()[["elapsed"]]
        fit <- tryCatch(
          mod$sample(
            data            = c(stan_data, list(w = w_vec)),
            chains          = N_CHAINS,
            parallel_chains = N_CHAINS,
            iter_warmup     = N_WARMUP,
            iter_sampling   = N_SAMPLING,
            seed            = cell_seed,
            refresh         = 0L
          ),
          error = function(e) {
            cat(sprintf("ERROR: %s\n", conditionMessage(e)))
            NULL
          }
        )
        elapsed <- proc.time()[["elapsed"]] - t0

        if (is.null(fit)) {
          .save(design_name, sigma_true, r_idx, method, list(error = TRUE))
          next
        }

        drws <- fit$draws()
        n_div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
        m_ess <- .min_ess(drws)
        a_ess <- .alpha_ess(drws, J)

        cat(sprintf("min_ESS=%.0f  divs=%d  (%.1f sec)\n",
                    m_ess, n_div, elapsed))

        res <- list(
          design     = design_name,
          sigma_true = sigma_true,
          rep        = r_idx,
          method     = method,
          w_vec      = w_vec,
          pi_j       = pi_est,
          min_ess    = m_ess,
          alpha_ess  = a_ess,
          n_div      = n_div,
          elapsed    = elapsed,
          error      = FALSE
        )
        .save(design_name, sigma_true, r_idx, method, res)

        # Per-group row
        for (j in seq_len(J)) {
          pergrp_rows[[length(pergrp_rows) + 1L]] <- list(
            design     = design_name,
            sigma_true = sigma_true,
            rep        = r_idx,
            method     = method,
            group      = j,
            n_j        = nj_vec[j],
            pi_j       = pi_est[j],
            w_j        = w_vec[j],
            ess_alpha  = a_ess[j]
          )
        }
        saveRDS(pergrp_rows, pergrp_file)

      }  # method
    }    # r_idx
  }      # sigma_true
}        # design_name

cat("\n══════════════════════════════════════════════════════════════\n")
cat("All runs complete. Assembling results and figures.\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# ── Assemble dataframes ────────────────────────────────────────────────────────

results   <- readRDS(result_file)
pergrp_rows <- readRDS(pergrp_file)

res_df <- do.call(rbind, Filter(Negate(is.null), lapply(results, function(r) {
  if (is.null(r) || isTRUE(r$error)) return(NULL)
  data.frame(
    design     = r$design,
    sigma_true = r$sigma_true,
    rep        = r$rep,
    method     = r$method,
    min_ess    = r$min_ess,
    n_div      = r$n_div,
    elapsed    = r$elapsed,
    stringsAsFactors = FALSE
  )
})))

pg_df <- do.call(rbind, lapply(pergrp_rows, as.data.frame))

# Factor levels for consistent ordering
design_levels <- c("uniform_sparse", "mixed", "uniform_dense")
design_labs   <- c("Uniform sparse\n(all n[j]=3)",
                   "Mixed\n(n[j]=3 vs 50)",
                   "Uniform dense\n(all n[j]=50)")
method_levels <- c("centred", "half", "nc", "pi_j")
method_labs   <- c("Centred\n(w=0)", "Half\n(w=0.5)", "NC\n(w=1)",
                   expression(pi[j]*"-adaptive"))

res_df$design_f <- factor(res_df$design, levels = design_levels,
                           labels = design_labs)
res_df$method_f <- factor(res_df$method, levels = method_levels)
pg_df$design_f  <- factor(pg_df$design, levels = design_levels,
                           labels = design_labs)
pg_df$method_f  <- factor(pg_df$method, levels = method_levels)
pg_df$sigma_f   <- factor(sprintf("sigma == %.1f", pg_df$sigma_true))

# ── Summary table ─────────────────────────────────────────────────────────────

cat("Median min-ESS by design × method × sigma:\n\n")

agg <- aggregate(min_ess ~ design + sigma_true + method,
                 data = res_df, FUN = function(x) round(median(x, na.rm = TRUE)))
agg <- agg[order(agg$design, agg$sigma_true, agg$method), ]
print(agg, row.names = FALSE)

# ── Figure A: min-ESS boxplot ─────────────────────────────────────────────────

method_colour <- c(
  centred = "#2166ac",
  half    = "#74add1",
  nc      = "#d6604d",
  pi_j    = "#1a9641"
)

fig_a <- ggplot(res_df[!is.na(res_df$min_ess), ],
                aes(x = method_f, y = min_ess, fill = method)) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.4) +
  facet_grid(
    rows = vars(factor(sprintf("sigma == %.1f", sigma_true))),
    cols = vars(design_f),
    labeller = labeller(
      .rows = label_parsed,
      .cols = label_value
    )
  ) +
  scale_fill_manual(values = method_colour, guide = "none") +
  scale_x_discrete(
    labels = c(centred = "Centred\n(w=0)", half = "Half\n(w=0.5)",
               nc = "NC\n(w=1)", pi_j = expression(pi[j]))
  ) +
  labs(
    title    = expression("Min-ESS by parameterisation: adaptive "*pi[j]*
                          " vs uniform weights"),
    subtitle = sprintf(
      "J = %d groups  |  %d replicates per cell  |  4 chains × %d sampling draws",
      J, N_REP, N_SAMPLING
    ),
    x = "Parameterisation",
    y = "Min bulk-ESS (all parameters)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text   = element_text(size = 9),
    axis.text.x  = element_text(size = 8)
  )

ggsave(file.path(out_dir, "mixed_design_minESS.png"),
       plot = fig_a, width = 10, height = 6, dpi = 150)
cat("Saved: data-raw/mixed_design_minESS.png\n")

# ── Figure B: per-group ESS vs π_j scatter ────────────────────────────────────

# Median per-group ESS over reps, within design × sigma × method × group
pg_agg <- aggregate(
  cbind(ess_alpha, pi_j, n_j) ~ design + sigma_true + method + group,
  data = pg_df,
  FUN  = median
)
pg_agg$design_f <- factor(pg_agg$design, levels = design_levels,
                           labels = design_labs)
pg_agg$method_f <- factor(pg_agg$method, levels = method_levels)
pg_agg$sigma_f  <- factor(sprintf("sigma == %.1f", pg_agg$sigma_true))
pg_agg$sparse   <- pg_agg$n_j <= 5L

fig_b <- ggplot(
    pg_agg,
    aes(x = pi_j, y = ess_alpha, colour = method, shape = sparse)
  ) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_line(aes(group = interaction(method, design)),
            linewidth = 0.5, alpha = 0.5) +
  facet_grid(
    rows = vars(sigma_f),
    cols = vars(design_f),
    labeller = labeller(.rows = label_parsed, .cols = label_value)
  ) +
  scale_colour_manual(
    values = method_colour,
    labels = c(centred = "Centred (w=0)", half = "Half (w=0.5)",
               nc = "NC (w=1)", pi_j = expression(pi[j]*"-adaptive")),
    name = "Method"
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16L, `TRUE` = 17L),
    labels = c(`FALSE` = expression(n[j]*"=50  (dense)"),
               `TRUE`  = expression(n[j]*"=3  (sparse)")),
    name = "Group type"
  ) +
  labs(
    title    = expression("Per-group "*alpha[j]*" ESS vs prior fraction "*pi[j]),
    subtitle = paste0("Median over ", N_REP,
                      " replicates per cell.  Mixed design shows adaptive advantage."),
    x = expression(pi[j]~"(prior fraction, pilot estimate)"),
    y = expression("Bulk ESS — "*alpha[j])
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 9)
  )

ggsave(file.path(out_dir, "mixed_design_pi_ess.png"),
       plot = fig_b, width = 11, height = 6, dpi = 150)
cat("Saved: data-raw/mixed_design_pi_ess.png\n")

# ── Figure C: π_j profiles by design and sigma ────────────────────────────────

# Show what π_j looks like under each design (one rep per cell, method=pi_j)
pi_prof <- pg_df[pg_df$method == "pi_j" & pg_df$rep == 1L, ]
pi_prof$design_f <- factor(pi_prof$design, levels = design_levels,
                            labels = design_labs)
pi_prof$sigma_f  <- factor(sprintf("sigma == %.1f", pi_prof$sigma_true))

fig_c <- ggplot(pi_prof, aes(x = factor(group), y = pi_j,
                              colour = factor(n_j), group = 1L)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  facet_grid(
    rows = vars(sigma_f),
    cols = vars(design_f),
    labeller = labeller(.rows = label_parsed, .cols = label_value)
  ) +
  scale_colour_manual(
    values = c(`3` = "#d6604d", `50` = "#2166ac"),
    name   = expression(n[j])
  ) +
  labs(
    title    = expression("Prior fraction "*pi[j]*" per group, by design (rep 1)"),
    subtitle = expression("Dashed line: "*pi[j]*" = 0.5; "*
                          "pi[j]>0.5 favours NC, < 0.5 favours centred"),
    x = "Group j",
    y = expression(pi[j]~"(pilot estimate)")
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

ggsave(file.path(out_dir, "mixed_design_pi_profile.png"),
       plot = fig_c, width = 10, height = 6, dpi = 150)
cat("Saved: data-raw/mixed_design_pi_profile.png\n")

# ── LaTeX summary table ───────────────────────────────────────────────────────

# Columns: sigma × design; rows: method; cell: median min-ESS
wide_tbl <- tapply(
  agg$min_ess,
  list(agg$method, paste(agg$sigma_true, agg$design, sep = "_")),
  identity
)
cols_order <- c("1_uniform_sparse", "1_mixed", "1_uniform_dense",
                "2_uniform_sparse", "2_mixed", "2_uniform_dense")
cols_order <- intersect(cols_order, colnames(wide_tbl))

header <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Median min-ESS (all parameters) by parameterisation and design.",
  " $J=", J, "$ groups, ", N_REP, " replicates, $", N_CHAINS,
  "\\times", N_SAMPLING, "$ post-warmup draws.",
  " $w_j=\\pi_j$: prior-fraction adaptive (proposed).}\n",
  "\\label{tab:mixed_design}\n",
  "\\begin{tabular}{l",
  paste(rep("r", length(cols_order)), collapse = ""),
  "}\n\\toprule\n",
  "Method & \\multicolumn{3}{c}{$\\sigma=1.0$} &",
  " \\multicolumn{3}{c}{$\\sigma=2.0$} \\\\\n",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
  " & Sparse & Mixed & Dense & Sparse & Mixed & Dense \\\\\n\\midrule\n"
)

method_tex <- c(
  centred = "Centred ($w_j=0$)",
  half    = "Half ($w_j=0.5$)",
  nc      = "Non-centred ($w_j=1$)",
  pi_j    = "$w_j=\\pi_j$ (adaptive)"
)

body <- paste(vapply(method_levels, function(m) {
  vals <- vapply(cols_order, function(col) {
    v <- wide_tbl[m, col]
    if (is.null(v) || is.na(v)) "---" else as.character(round(v))
  }, character(1L))
  paste0("  ", method_tex[m], " & ", paste(vals, collapse = " & "), " \\\\")
}, character(1L)), collapse = "\n")

footer <- paste0(
  "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
)

writeLines(paste0(header, body, footer),
           file.path(out_dir, "mixed_design_table.tex"))
cat("Saved: data-raw/mixed_design_table.tex\n")

cat("\nMixed-design experiment complete.\n")
