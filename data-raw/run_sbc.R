## M3: SBC validation battery for the conditional-transport samplers.
##
## Tests four samplers on J=4, n_j=5 GLMMs drawn from the prior:
##   M2a-conditional  : horizontal_mcmc(transport = "conditional")  [must pass]
##   M2a-fisher_legacy: horizontal_mcmc(transport = "fisher_legacy") [expected fail; paper fig]
##   M2b-reparam      : horizontal_hmc()                            [must pass]
##   M2c-marginal     : marginal_mcmc()                             [must pass]
##
## Algorithm: for each of N_SBC prior draws,
##   1. Sample (mu, sigma, beta, alpha) from the joint prior.
##   2. Simulate y from the likelihood.
##   3. Run sampler; thin chain to ESS(mu) draws; compute fractional rank of
##      each true parameter in the thinned posterior draws.
##   4. Fractional ranks should be Uniform[0,1] under correct targeting.
##
## Thinning rationale: SBC rank uniformity requires approximately IID posterior
## draws (Talts et al. 2018). With thin_factor = floor(n_iter / ESS_mu) the
## thinned draws are approximately independent and the chi-sq GOF is valid.
## Fractional ranks (rank / L_thinned) are stored so the chi-sq test is
## invariant to the per-dataset L_thinned.
##
## Warmup note: MwG (M2a) and RWM (M2c) need more warmup than HMC because the
## simple multiplicative step-size tuner requires ~10 windows (1000 iters) to
## grow the initial step 0.10 up to ~0.4-0.8.
##
## Outputs (written to data-raw/):
##   sbc_{method}.rds          per-dataset fractional rank matrices (checkpoint)
##   sbc_ranks_{method}.png    rank histogram panels (one per method)
##   sbc_legacy_fail.png       paper figure: fisher_legacy failure
##   sbc_summary.csv           chi-squared GOF results per (method, parameter)
##
## Use the cell_seed scheme: seed_s = SBC_BASE + s (reproducible, resumable).
## SBC_BASE = 271828 (fresh seed; previous runs used 314159).
##
## Run from the package root:  Rscript data-raw/run_sbc.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(posterior)
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

N_SBC      <- 400L        # prior draws
N_CHAINS   <- 1L          # single chain per dataset (SBC is valid)
J          <- 4L
N_J        <- 5L
K          <- 2L
N          <- J * N_J
N_BINS     <- 20L         # bins for rank histograms
SBC_BASE   <- 271828L     # fresh seed (e); prior run used 314159

# Per-method settings. Marginal uses 4000 post-warmup iterations to provide
# enough raw draws for ESS-based thinning to yield ~100 approximately-IID draws.
METHOD_SPECS <- list(
  conditional   = list(n_iter = 4000L, n_warmup = 2000L),
  fisher_legacy = list(n_iter = 1000L, n_warmup =  500L),
  reparam_hmc   = list(n_iter =  500L, n_warmup =  200L),
  marginal      = list(n_iter = 4000L, n_warmup = 2000L)
)

# ── Fixed design matrix ───────────────────────────────────────────────────────

set.seed(SBC_BASE)
X_sbc     <- matrix(rnorm(N * K), N, K)
group_sbc <- rep(seq_len(J), each = N_J)

# ── Parameter names ───────────────────────────────────────────────────────────

par_names_full <- c(
  "mu", "sigma",
  paste0("alpha[", seq_len(J), "]"),
  paste0("beta[",  seq_len(K), "]")
)

# ── Prior simulator ───────────────────────────────────────────────────────────

.sim_glmm_prior <- function() {
  mu    <- rnorm(1L, 0, 5)
  sigma <- rexp(1L, rate = 1)
  beta  <- rnorm(K, 0, 2)
  alpha <- rnorm(J, mu, sigma)
  list(mu = mu, sigma = sigma, beta = beta, alpha = alpha)
}

# ── Likelihood simulator ──────────────────────────────────────────────────────

.sim_glmm_y <- function(theta) {
  eta <- theta$alpha[group_sbc] + as.vector(X_sbc %*% theta$beta)
  as.integer(rbinom(N, 1L, plogis(eta)))
}

# ── Rank computation: thin to ESS, return fractional rank in [0,1] ───────────
#
# Thinning to ESS(mu) gives approximately IID draws, making the chi-sq GOF
# test for rank uniformity valid. Fractional rank is scale-invariant.

.compute_ranks <- function(draws_arr, theta_true) {
  mat <- as.matrix(as_draws_matrix(draws_arr))
  n   <- nrow(mat)

  # Thin factor = floor(n / ESS_mu); fall back to no thinning on error
  ess_mu <- tryCatch(
    max(1L, floor(posterior::ess_bulk(draws_arr[, , "mu"]))),
    error = function(e) n
  )
  thin_by <- max(1L, floor(n / ess_mu))
  idx     <- seq(1L, n, by = thin_by)
  mat_t   <- mat[idx, , drop = FALSE]
  L_t     <- nrow(mat_t)

  ranks <- setNames(numeric(length(par_names_full)), par_names_full)
  for (v in par_names_full)
    if (v %in% colnames(mat_t))
      ranks[v] <- sum(mat_t[, v] < theta_true[[v]]) / L_t
  ranks   # fractional, in [0, 1]
}

# ── True-value list from a prior draw ────────────────────────────────────────

.theta_to_named <- function(theta) {
  nms <- list(mu = theta$mu, sigma = theta$sigma)
  for (j in seq_len(J)) nms[[sprintf("alpha[%d]", j)]] <- theta$alpha[j]
  for (k in seq_len(K)) nms[[sprintf("beta[%d]",  k)]] <- theta$beta[k]
  nms
}

# ── Generic SBC runner ────────────────────────────────────────────────────────

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_sbc <- function(method_name, sampler_fn, ...) {
  spec     <- METHOD_SPECS[[method_name]]
  n_iter_m <- spec$n_iter
  n_warm_m <- spec$n_warmup

  out_file <- file.path(out_dir, sprintf("sbc_%s.rds", method_name))
  results  <- if (file.exists(out_file)) readRDS(out_file) else vector("list", N_SBC)

  if (length(results) < N_SBC) length(results) <- N_SBC

  n_done <- sum(!vapply(results, is.null, logical(1L)))
  if (n_done > 0L) cat(sprintf("  resuming %s from dataset %d\n", method_name, n_done + 1L))

  for (s in seq_len(N_SBC)) {
    if (!is.null(results[[s]])) next

    set.seed(SBC_BASE + s)
    theta_s <- .sim_glmm_prior()
    y_s     <- .sim_glmm_y(theta_s)
    data_s  <- list(N = N, J = J, group = group_sbc, X = X_sbc, y = y_s)

    draws_s <- tryCatch(
      sampler_fn(data_s, n_iter = n_iter_m, n_warmup = n_warm_m,
                 n_chains = N_CHAINS, verbose = FALSE, seed = SBC_BASE + s, ...),
      error = function(e) NULL
    )

    if (!is.null(draws_s) && all(is.finite(as.array(draws_s)))) {
      results[[s]] <- list(
        ranks  = .compute_ranks(draws_s, .theta_to_named(theta_s)),
        n_iter = n_iter_m,
        error  = FALSE
      )
    } else {
      results[[s]] <- list(ranks = NULL, n_iter = n_iter_m, error = TRUE)
    }

    saveRDS(results, out_file)

    if (s %% 20L == 0L || s == N_SBC)
      cat(sprintf("  %s: %d/%d done\n", method_name, s, N_SBC))
  }

  results
}

# ── Assemble fractional-rank matrix from checkpoint list ──────────────────────

.ranks_matrix <- function(results) {
  ok <- vapply(results, function(r) !is.null(r) && !r$error && !is.null(r$ranks),
               logical(1L))
  if (!any(ok)) return(NULL)
  mat    <- do.call(rbind, lapply(results[ok], `[[`, "ranks"))
  n_iter <- results[ok][[1L]]$n_iter %||% NA_integer_
  list(mat = mat, n_iter = n_iter)   # n_iter kept for reporting only
}

# ── Chi-squared GOF: Uniform[0,1] on fractional ranks ────────────────────────

.chisq_sbc <- function(ranks_vec) {
  breaks <- seq(0, 1, length.out = N_BINS + 1L)
  obs    <- as.vector(table(cut(ranks_vec, breaks = breaks, include.lowest = TRUE)))
  exp_c  <- length(ranks_vec) / N_BINS
  stat   <- sum((obs - exp_c)^2 / exp_c)
  list(stat = stat, df = N_BINS - 1L, pvalue = pchisq(stat, df = N_BINS - 1L, lower.tail = FALSE))
}

# ── Rank histogram (fractional ranks) ─────────────────────────────────────────

.plot_ranks <- function(rm_obj, method_label, fill_col = "#4477AA") {
  rank_mat <- rm_obj$mat
  n_iter   <- rm_obj$n_iter
  breaks   <- seq(0, 1, length.out = N_BINS + 1L)
  expected <- nrow(rank_mat) / N_BINS

  df <- as.data.frame(rank_mat) |>
    tidyr::pivot_longer(everything(), names_to = "parameter", values_to = "rank") |>
    dplyr::mutate(
      parameter = factor(parameter, levels = par_names_full),
      bin = cut(rank, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    ) |>
    dplyr::count(parameter, bin)

  subtitle <- if (!is.na(n_iter))
    sprintf("N = %d draws, %d post-warmup iter, ESS-thinned; dashed = expected",
            nrow(rank_mat), n_iter)
  else
    sprintf("N = %d draws, ESS-thinned; dashed = expected", nrow(rank_mat))

  ggplot(df, aes(x = bin, y = n)) +
    geom_col(fill = fill_col, colour = "white", width = 0.9) +
    geom_hline(yintercept = expected, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ parameter, ncol = 4L) +
    scale_x_continuous(
      name   = sprintf("Rank fractile (%d bins)", N_BINS),
      breaks = c(1, N_BINS),
      labels = c("0", "1")
    ) +
    scale_y_continuous(name = "Count") +
    labs(title = sprintf("SBC rank histograms: %s", method_label),
         subtitle = subtitle) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 7))
}

# ── Run all four methods ──────────────────────────────────────────────────────

cat("\n── M2a: transport = 'conditional' ──────────────────────────────────────────\n")
res_cond <- run_sbc(
  "conditional",
  function(...) horizontal_mcmc(..., transport = "conditional")
)

cat("\n── M2a: transport = 'fisher_legacy' ────────────────────────────────────────\n")
res_legacy <- run_sbc(
  "fisher_legacy",
  function(...) horizontal_mcmc(..., transport = "fisher_legacy")
)

cat("\n── M2b: horizontal_hmc (reparam) ────────────────────────────────────────────\n")
res_hmc <- run_sbc(
  "reparam_hmc",
  function(...) horizontal_hmc(..., transport = "reparam", L = 5L)
)

cat("\n── M2c: marginal_mcmc ───────────────────────────────────────────────────────\n")
res_marg <- run_sbc(
  "marginal",
  function(...) marginal_mcmc(...)
)

# ── Build rank matrices ───────────────────────────────────────────────────────

rm_cond   <- .ranks_matrix(res_cond)
rm_legacy <- .ranks_matrix(res_legacy)
rm_hmc    <- .ranks_matrix(res_hmc)
rm_marg   <- .ranks_matrix(res_marg)

# ── Chi-squared summary ───────────────────────────────────────────────────────

.chisq_table <- function(rm_obj, method_name) {
  if (is.null(rm_obj)) return(NULL)
  rank_mat <- rm_obj$mat
  rows <- lapply(par_names_full, function(v) {
    cq <- .chisq_sbc(rank_mat[, v])
    data.frame(
      method    = method_name,
      parameter = v,
      n_iter    = rm_obj$n_iter,
      chisq     = round(cq$stat, 2),
      df        = cq$df,
      pvalue    = signif(cq$pvalue, 3),
      pass      = cq$pvalue > 0.01,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

chisq_summary <- dplyr::bind_rows(
  .chisq_table(rm_cond,   "M2a-conditional"),
  .chisq_table(rm_legacy, "M2a-fisher_legacy"),
  .chisq_table(rm_hmc,    "M2b-reparam_hmc"),
  .chisq_table(rm_marg,   "M2c-marginal")
)

cat("\n── Chi-squared GOF summary (pass = p > 0.01) ────────────────────────────────\n")
print(chisq_summary[, setdiff(names(chisq_summary), "n_iter")], row.names = FALSE)

write.csv(chisq_summary, file.path(out_dir, "sbc_summary.csv"), row.names = FALSE)
cat("Saved: data-raw/sbc_summary.csv\n")

# ── Rank histogram plots ──────────────────────────────────────────────────────

.save_plot <- function(rm_obj, method_name, label, filename, fill_col = "#4477AA") {
  if (is.null(rm_obj)) {
    cat(sprintf("  skipping plot for %s (no valid results)\n", method_name))
    return(invisible(NULL))
  }
  p    <- .plot_ranks(rm_obj, label, fill_col)
  path <- file.path(out_dir, filename)
  ggsave(path, p, width = 10, height = 6, dpi = 150)
  cat(sprintf("Saved: data-raw/%s\n", filename))
}

.save_plot(rm_cond,   "conditional",   "M2a — conditional transport",
           "sbc_ranks_conditional.png")
.save_plot(rm_legacy, "fisher_legacy",
           "M2a — Fisher legacy (biased; no Jacobian correction)",
           "sbc_ranks_fisher_legacy.png", fill_col = "#CC3333")
.save_plot(rm_hmc,    "reparam_hmc",   "M2b — reparameterised HMC",
           "sbc_ranks_reparam_hmc.png")
.save_plot(rm_marg,   "marginal",      "M2c — marginal MCMC (gold standard)",
           "sbc_ranks_marginal.png")

# ── Paper figure: fisher_legacy failure ──────────────────────────────────────

if (!is.null(rm_legacy)) {
  p_fail <- .plot_ranks(rm_legacy,
                        "Fisher-legacy transport (no Jacobian correction) — SBC failure",
                        fill_col = "#CC3333") +
    labs(title    = "SBC failure: Fisher-legacy transport (no Jacobian correction)",
         subtitle = sprintf(
           "Rank non-uniformity indicates incorrect targeting (%d prior draws)",
           nrow(rm_legacy$mat)))
  ggsave(file.path(out_dir, "sbc_legacy_fail.png"), p_fail,
         width = 10, height = 6, dpi = 150)
  cat("Saved: data-raw/sbc_legacy_fail.png\n")
}

# ── Method-level pass/fail summary ───────────────────────────────────────────

cat("\n── Method-level summary (pass = ALL parameters have p > 0.01) ──────────────\n")
method_summary <- chisq_summary |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    n_params   = dplyr::n(),
    n_pass     = sum(pass),
    all_pass   = all(pass),
    min_pvalue = signif(min(pvalue), 3),
    .groups    = "drop"
  )
print(method_summary, width = 80)

cat("\nM3 SBC complete.\n")
cat(sprintf("  Pass: %s\n",
            paste(method_summary$method[method_summary$all_pass], collapse = ", ")))
cat(sprintf("  Fail: %s\n",
            paste(method_summary$method[!method_summary$all_pass], collapse = ", ")))
