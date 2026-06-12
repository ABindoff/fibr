## M4a: sanity benchmarks — relative efficiency across the 4×4 simulation grid.
##
## Reuses the grid and cell seeds from run_simulation_study.R.
## Comparators:
##   centred_nuts      NUTS on glmm_centred.stan       (cached draws reused)
##   nc_nuts           NUTS on glmm_noncentred.stan
##   partial_nc        NUTS on glmm_partial_nc.stan     w_j = pi_j from pilot
##   asis              asis_mcmc()                      MwG + NC interweaving
##   m2a_conditional   horizontal_mcmc()                conditional transport MwG
##   m2b_reparam       horizontal_hmc()                 reparameterised HMC
##   riemannian        riemannian_mcmc()                SoftAbs metric HMC
##   m2c_marginal      marginal_mcmc()                  reference (GH + exact alpha)
##
## Metrics: ESS/sec and ESS/n_grad_equiv (1 per log-post eval).
##
## Headline figure: relative efficiency of m2a/m2b vs partial_nc vs asis
## as a function of mean(pi_j), mean(|kappa_j|), mean(|F_j|) at each cell.
## Check collinearity of |kappa_j| and |F_j| first; if r > 0.95 flag and
## stop for Aidan before fitting separate-predictor models.
##
## Expectation (M4a is the control layer): on this grid M2c should dominate;
## the real question lives in M4b.
##
## Run from package root:  Rscript data-raw/run_m4a.R
## Resume:  re-running the script skips completed (cell, rep, method) triples.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)
invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))
out_dir    <- file.path(pkg_root, "data-raw")
draws_dir  <- file.path(out_dir, "simstud_draws")

# ── Study parameters (match run_simulation_study.R) ──────────────────────────

NJ_GRID    <- c(3L, 10L, 30L, 100L)
SIGMA_GRID <- c(0.5, 1.0, 2.0, 3.0)
J          <- 8L
K          <- 2L
N_REP      <- 10L
MU_TRUE    <- 0
BETA_TRUE  <- c(0.8, -0.5)

METHODS <- c("centred_nuts", "nc_nuts", "partial_nc",
             "asis", "m2a_conditional", "m2b_reparam",
             "riemannian", "m2c_marginal")

# Post-warmup iterations per method (n_chains = 4 for all).
# marginal uses more raw draws so ESS-based thinning still yields ~100 IID draws.
METHOD_ITERS <- list(
  centred_nuts   = list(n_warmup = 1000L, n_iter = 2000L, L = NA_integer_),
  nc_nuts        = list(n_warmup = 1000L, n_iter = 2000L, L = NA_integer_),
  partial_nc     = list(n_warmup = 1000L, n_iter = 2000L, L = NA_integer_),
  asis           = list(n_warmup = 1000L, n_iter = 2000L, L = NA_integer_),
  m2a_conditional= list(n_warmup = 1000L, n_iter = 2000L, L = NA_integer_),
  m2b_reparam    = list(n_warmup =  200L, n_iter =  500L, L = 5L),
  riemannian     = list(n_warmup = 1000L, n_iter = 2000L, L = 10L),
  m2c_marginal   = list(n_warmup = 2000L, n_iter = 4000L, L = NA_integer_)
)
N_CHAINS <- 4L
GH_NODES <- 15L   # Gauss-Hermite nodes used in marginal_mcmc

# ── Compile Stan models ───────────────────────────────────────────────────────

cat("Compiling Stan models...\n")
mod_centred    <- cmdstan_model(file.path(pkg_root, "inst/stan/glmm_centred.stan"))
mod_nc         <- cmdstan_model(file.path(pkg_root, "inst/stan/glmm_noncentred.stan"))
mod_partial_nc <- cmdstan_model(file.path(pkg_root, "inst/stan/glmm_partial_nc.stan"))
cat("Done.\n\n")

# ── Checkpoint helpers ────────────────────────────────────────────────────────

result_file <- file.path(out_dir, "m4a_results.rds")
results     <- if (file.exists(result_file)) readRDS(result_file) else list()

.result_key <- function(n_j, sigma, rep, method)
  sprintf("nj%d_s%.1f_r%02d_%s", n_j, sigma, rep, method)

.already_done <- function(n_j, sigma, rep, method)
  !is.null(results[[.result_key(n_j, sigma, rep, method)]])

.save_result <- function(n_j, sigma, rep, method, res) {
  results[[.result_key(n_j, sigma, rep, method)]] <<- res
  saveRDS(results, result_file)
}

# ── ESS computation helper ────────────────────────────────────────────────────

.summarise_draws <- function(drws, n_grad_equiv, elapsed_sec) {
  if (is.null(drws)) return(NULL)
  mat   <- tryCatch(as_draws_matrix(drws), error = function(e) NULL)
  if (is.null(mat)) return(NULL)
  vars  <- dimnames(mat)[[2]]
  vars  <- vars[!grepl("^log_lik|^lp__|^alpha_tilde", vars)]
  alpha_vars <- vars[grepl("^alpha\\[", vars)]
  base_vars  <- c("mu", "sigma")
  beta_vars  <- vars[grepl("^beta\\[", vars)]

  .ess <- function(v) tryCatch(
    max(0, posterior::ess_bulk(drws[, , v])),
    error = function(e) NA_real_
  )
  ess_all   <- vapply(vars, .ess, numeric(1L))
  ess_alpha <- vapply(alpha_vars, .ess, numeric(1L))

  list(
    min_ess       = min(ess_all, na.rm = TRUE),
    alpha_min_ess = if (length(ess_alpha) > 0) min(ess_alpha, na.rm = TRUE) else NA_real_,
    mu_ess        = .ess("mu"),
    sigma_ess     = .ess("sigma"),
    n_grad_equiv  = n_grad_equiv,
    elapsed_sec   = elapsed_sec,
    ess_per_sec   = min(ess_all, na.rm = TRUE) / max(elapsed_sec, 1e-6),
    ess_per_grad  = min(ess_all, na.rm = TRUE) / max(n_grad_equiv, 1L)
  )
}

# ── Cell-level geometry helper ────────────────────────────────────────────────

.cell_geometry <- function(centred_draws, stan_data) {
  # Posterior mean of base parameters from centred NUTS pilot
  mat   <- as.matrix(as_draws_matrix(centred_draws[, ,
    c("mu", "sigma", "beta[1]", "beta[2]")]))
  mu_m  <- mean(mat[, "mu"])
  sig_m <- mean(mat[, "sigma"])
  be_m  <- c(mean(mat[, "beta[1]"]), mean(mat[, "beta[2]"]))
  ls_m  <- log(sig_m)

  lap <- tryCatch(
    .glmm_laplace(c(mu_m, ls_m), be_m, stan_data),
    error = function(e) NULL
  )
  if (is.null(lap)) return(NULL)

  # G_FF_j = 1/s_j^2  (conditional precision at the Laplace mode)
  G_FF <- 1 / lap$s^2

  # pi_j (prior fraction)
  pi_j <- .glmm_prior_fraction(G_FF, sig_m)

  # |F_j| (curvature magnitude)
  F_j  <- abs(.glmm_curvature(G_FF, sig_m))

  # kappa_j = -t_j(m_j) * s_j^3
  # t_j(a) = sum_{i in j} p_i(1-p_i)(1-2p_i)
  eta   <- lap$m[stan_data$group] + as.vector(stan_data$X %*% be_m)
  pv    <- plogis(eta)
  t_j   <- as.vector(tapply(pv * (1 - pv) * (1 - 2 * pv), stan_data$group, sum))
  kappa_j <- -t_j * lap$s^3

  # Per-group pi_j as the Stan partial-NC weights (clamp to [0,1])
  w_j <- pmin(pmax(pi_j, 0), 1)

  list(
    mu = mu_m, sigma = sig_m, beta = be_m,
    pi_j = pi_j, F_j = F_j, kappa_j = kappa_j,
    w_j = w_j,
    mean_pi    = mean(pi_j),
    mean_abs_F = mean(F_j),
    mean_abs_kappa = mean(abs(kappa_j))
  )
}

# ── Main grid loop ─────────────────────────────────────────────────────────────

total_runs <- length(NJ_GRID) * length(SIGMA_GRID) * N_REP * length(METHODS)
done_runs  <- sum(vapply(names(results), function(k) !is.null(results[[k]]), logical(1L)))
cat(sprintf("M4a: %d total (method × cell × rep), %d already complete.\n\n",
            total_runs, done_runs))

for (n_j in NJ_GRID) {
  for (sigma_true in SIGMA_GRID) {
    for (r_idx in seq_len(N_REP)) {

      # Skip if all methods done for this cell-rep
      if (all(vapply(METHODS, function(m) .already_done(n_j, sigma_true, r_idx, m),
                     logical(1L)))) {
        next
      }

      cat(sprintf("\n── n_j=%3d  sigma=%.1f  rep=%02d ──────────────────────────────────\n",
                  n_j, sigma_true, r_idx))

      # ── Reproduce data from simulation study seed ──────────────────────────
      nj_idx    <- which(NJ_GRID    == n_j)
      sig_idx   <- which(abs(SIGMA_GRID - sigma_true) < 1e-9)
      cell_seed <- nj_idx * 10000L + sig_idx * 1000L + r_idx

      set.seed(cell_seed)
      N          <- J * n_j
      alpha_true <- rnorm(J,  MU_TRUE,   sigma_true)
      group_id   <- rep(seq_len(J), each = n_j)
      X          <- matrix(rnorm(N * K), ncol = K)
      eta        <- alpha_true[group_id] + X %*% BETA_TRUE
      y          <- rbinom(N, 1L, plogis(eta))
      stan_data  <- list(N = N, J = J, group = group_id, X = X, y = y)

      # ── Load centred NUTS (cached) ─────────────────────────────────────────
      draw_file <- file.path(draws_dir,
                             sprintf("nj%d_s%.1f_r%02d.rds", n_j, sigma_true, r_idx))
      centred_draws <- if (file.exists(draw_file)) readRDS(draw_file) else NULL

      # ── Cell geometry (pi_j, kappa_j, |F_j|) from centred NUTS ────────────
      geom <- if (!is.null(centred_draws)) {
        tryCatch(.cell_geometry(centred_draws, stan_data), error = function(e) NULL)
      } else NULL

      # ── Method runners ─────────────────────────────────────────────────────

      run_method <- function(method) {
        if (.already_done(n_j, sigma_true, r_idx, method)) {
          cat(sprintf("  skip %s\n", method))
          return(invisible(NULL))
        }
        cat(sprintf("  %s ...", method))
        spec <- METHOD_ITERS[[method]]
        t0   <- proc.time()[["elapsed"]]
        drws <- NULL
        n_grad <- NA_integer_

        if (method == "centred_nuts") {
          drws    <- centred_draws
          n_grad  <- if (!is.null(drws)) N_CHAINS * spec$n_iter * 2L else NA_integer_
          elapsed <- if (!is.null(drws)) NA_real_ else NA_real_
        } else {
          elapsed <- tryCatch({
            if (method == "nc_nuts") {
              fit <- mod_nc$sample(
                data            = stan_data,
                chains          = N_CHAINS,
                parallel_chains = min(N_CHAINS, 4L),
                iter_warmup     = spec$n_warmup,
                iter_sampling   = spec$n_iter,
                seed            = cell_seed,
                refresh         = 0L
              )
              drws   <<- fit$draws()
              n_grad <<- sum(fit$sampler_diagnostics()[, , "n_leapfrog__"])
            } else if (method == "partial_nc") {
              if (is.null(geom)) stop("no geometry")
              fit <- mod_partial_nc$sample(
                data            = c(stan_data, list(w = geom$w_j)),
                chains          = N_CHAINS,
                parallel_chains = min(N_CHAINS, 4L),
                iter_warmup     = spec$n_warmup,
                iter_sampling   = spec$n_iter,
                seed            = cell_seed,
                refresh         = 0L
              )
              drws   <<- fit$draws()
              n_grad <<- sum(fit$sampler_diagnostics()[, , "n_leapfrog__"])
            } else if (method == "asis") {
              drws <<- asis_mcmc(stan_data,
                n_iter = spec$n_iter, n_warmup = spec$n_warmup,
                n_chains = N_CHAINS, seed = cell_seed, verbose = FALSE)
              n_grad <<- N_CHAINS * (spec$n_warmup + spec$n_iter) * 3L
            } else if (method == "m2a_conditional") {
              drws <<- horizontal_mcmc(stan_data,
                transport = "conditional",
                n_iter = spec$n_iter, n_warmup = spec$n_warmup,
                n_chains = N_CHAINS, seed = cell_seed, verbose = FALSE)
              n_grad <<- N_CHAINS * (spec$n_warmup + spec$n_iter) * 1L
            } else if (method == "m2b_reparam") {
              drws <<- horizontal_hmc(stan_data,
                transport = "reparam", L = spec$L,
                n_iter = spec$n_iter, n_warmup = spec$n_warmup,
                n_chains = N_CHAINS, seed = cell_seed, verbose = FALSE)
              n_grad <<- N_CHAINS * (spec$n_warmup + spec$n_iter) * spec$L
            } else if (method == "riemannian") {
              drws <<- riemannian_mcmc(stan_data,
                n_iter = spec$n_iter, n_warmup = spec$n_warmup,
                n_chains = N_CHAINS, L = spec$L,
                seed = cell_seed, verbose = FALSE)
              n_grad <<- N_CHAINS * (spec$n_warmup + spec$n_iter) * spec$L
            } else if (method == "m2c_marginal") {
              drws <<- marginal_mcmc(stan_data,
                n_iter = spec$n_iter, n_warmup = spec$n_warmup,
                n_chains = N_CHAINS, seed = cell_seed, verbose = FALSE)
              n_grad <<- N_CHAINS * (spec$n_warmup + spec$n_iter) * J * GH_NODES
            }
            proc.time()[["elapsed"]] - t0
          }, error = function(e) {
            cat(sprintf(" ERROR: %s\n", conditionMessage(e)))
            NA_real_
          })
        }

        elapsed_final <- if (method == "centred_nuts") NA_real_ else elapsed

        smry <- if (!is.null(drws)) {
          .summarise_draws(drws, n_grad, elapsed_final)
        } else NULL

        res <- list(
          n_j = n_j, sigma_true = sigma_true, rep = r_idx, method = method,
          cell_seed = cell_seed,
          n_warmup = spec$n_warmup, n_iter = spec$n_iter,
          elapsed_sec  = elapsed_final,
          n_grad_equiv = n_grad,
          min_ess      = if (!is.null(smry)) smry$min_ess       else NA_real_,
          alpha_min_ess= if (!is.null(smry)) smry$alpha_min_ess else NA_real_,
          mu_ess       = if (!is.null(smry)) smry$mu_ess        else NA_real_,
          sigma_ess    = if (!is.null(smry)) smry$sigma_ess     else NA_real_,
          ess_per_sec  = if (!is.null(smry)) smry$ess_per_sec   else NA_real_,
          ess_per_grad = if (!is.null(smry)) smry$ess_per_grad  else NA_real_,
          mean_pi      = if (!is.null(geom)) geom$mean_pi            else NA_real_,
          mean_abs_F   = if (!is.null(geom)) geom$mean_abs_F         else NA_real_,
          mean_abs_kappa = if (!is.null(geom)) geom$mean_abs_kappa   else NA_real_,
          error = is.null(smry)
        )

        .save_result(n_j, sigma_true, r_idx, method, res)
        if (!is.null(smry))
          cat(sprintf(" min_ESS=%.0f  ESS/sec=%.1f\n",
                      smry$min_ess, smry$ess_per_sec))
        invisible(NULL)
      }

      for (method in METHODS) run_method(method)
    }
  }
}

# ── Assemble results dataframe ─────────────────────────────────────────────────

cat("\n── Assembling results ────────────────────────────────────────────────────\n")

res_df <- dplyr::bind_rows(lapply(results, function(r) {
  if (is.null(r) || r$error) return(NULL)
  data.frame(
    n_j          = r$n_j,
    sigma_true   = r$sigma_true,
    rep          = r$rep,
    method       = r$method,
    n_iter       = r$n_iter,
    elapsed_sec  = r$elapsed_sec,
    n_grad_equiv = r$n_grad_equiv,
    min_ess      = r$min_ess,
    alpha_ess    = r$alpha_min_ess,
    mu_ess       = r$mu_ess,
    sigma_ess    = r$sigma_ess,
    ess_per_sec  = r$ess_per_sec,
    ess_per_grad = r$ess_per_grad,
    mean_pi      = r$mean_pi,
    mean_abs_F   = r$mean_abs_F,
    mean_abs_kappa = r$mean_abs_kappa,
    stringsAsFactors = FALSE
  )
}))

if (is.null(res_df) || nrow(res_df) == 0L) {
  cat("No complete results yet.\n")
  quit(save = "no")
}

write.csv(res_df, file.path(out_dir, "m4a_summary.csv"), row.names = FALSE)
cat("Saved: data-raw/m4a_summary.csv\n")

# ── Collinearity check: |kappa_j| vs |F_j| ────────────────────────────────────

cell_geom <- res_df |>
  dplyr::filter(method == "centred_nuts") |>
  dplyr::select(n_j, sigma_true, mean_pi, mean_abs_F, mean_abs_kappa) |>
  dplyr::distinct()

if (nrow(cell_geom) >= 4L) {
  r_kappa_F <- cor(cell_geom$mean_abs_kappa, cell_geom$mean_abs_F, use = "complete.obs")
  cat(sprintf("\n|kappa_j| vs |F_j| correlation across cells: r = %.3f\n", r_kappa_F))
  if (abs(r_kappa_F) > 0.95)
    cat("WARNING: |kappa_j| and |F_j| are too collinear (|r| > 0.95).\n")
  cat("         Stop and ask Aidan whether to extend the grid.\n")
} else {
  cat("\nToo few cells for collinearity check.\n")
}

# ── Efficiency table ──────────────────────────────────────────────────────────

cat("\n── ESS/sec by method (median over reps) ─────────────────────────────────\n")

eff_tbl <- res_df |>
  dplyr::group_by(n_j, sigma_true, method) |>
  dplyr::summarise(
    n_rep      = dplyr::n(),
    med_ess_sec  = median(ess_per_sec, na.rm = TRUE),
    med_min_ess  = median(min_ess, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(n_j, sigma_true, dplyr::desc(med_ess_sec))

print(eff_tbl, n = 80)

# ── Relative efficiency vs partial_nc ─────────────────────────────────────────

cat("\n── Relative efficiency (vs partial_nc, median over reps) ────────────────\n")

ref_eff <- res_df |>
  dplyr::filter(method == "partial_nc") |>
  dplyr::select(n_j, sigma_true, rep, ref_ess_sec = ess_per_sec)

rel_eff <- res_df |>
  dplyr::left_join(ref_eff, by = c("n_j", "sigma_true", "rep")) |>
  dplyr::mutate(rel_eff = ess_per_sec / ref_ess_sec) |>
  dplyr::group_by(n_j, sigma_true, method) |>
  dplyr::summarise(med_rel_eff = median(rel_eff, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = method, values_from = med_rel_eff)

print(rel_eff, width = 100)

# ── Figure: relative efficiency vs mean_pi ────────────────────────────────────

plot_df <- res_df |>
  dplyr::left_join(ref_eff, by = c("n_j", "sigma_true", "rep")) |>
  dplyr::mutate(rel_eff = ess_per_sec / ref_ess_sec) |>
  dplyr::filter(method %in% c("m2a_conditional", "m2b_reparam", "asis",
                               "m2c_marginal", "riemannian")) |>
  dplyr::group_by(n_j, sigma_true, method) |>
  dplyr::summarise(
    med_rel = median(rel_eff, na.rm = TRUE),
    mean_pi = mean(mean_pi, na.rm = TRUE),
    mean_abs_kappa = mean(mean_abs_kappa, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(plot_df) > 0L) {
  p_pi <- ggplot(plot_df, aes(x = mean_pi, y = med_rel, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    scale_y_log10(name = "Relative ESS/sec (vs partial-NC, log scale)") +
    scale_x_continuous(name = expression(bar(pi)[j] ~ "(mean prior fraction)")) +
    scale_colour_brewer(palette = "Dark2", name = "Method") +
    labs(title = "M4a: Relative efficiency vs prior fraction",
         subtitle = "Points = median over 10 reps per cell; line = loess") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "m4a_rel_eff_pi.png"), p_pi,
         width = 8, height = 5, dpi = 150)
  cat("Saved: data-raw/m4a_rel_eff_pi.png\n")

  p_kappa <- ggplot(plot_df, aes(x = mean_abs_kappa, y = med_rel, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    scale_y_log10(name = "Relative ESS/sec (vs partial-NC, log scale)") +
    scale_x_continuous(name = expression(bar(abs)(kappa)[j] ~ "(mean non-Gaussianity)")) +
    scale_colour_brewer(palette = "Dark2", name = "Method") +
    labs(title = "M4a: Relative efficiency vs non-Gaussianity",
         subtitle = "Triage hypothesis: M2a/M2b advantage predicted by |kappa_j|, not |F_j|") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "m4a_rel_eff_kappa.png"), p_kappa,
         width = 8, height = 5, dpi = 150)
  cat("Saved: data-raw/m4a_rel_eff_kappa.png\n")
}

cat("\nM4a complete.\n")
