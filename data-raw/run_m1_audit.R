# M1.3 audit: is horizontal_mcmc() biased, and by how much?
#
# Parts:
#   A. Finite-difference verification of the analytic dA_j/dalpha_j
#      (manuscript2/M1_exactness_memo.md, Section 2).
#   B. Proposal-level diagnostics at a representative point: log|det DT|
#      and round-trip defect ||alpha - T_rev(T(alpha))|| for typical steps.
#   C. Long-run bias audit on the sparse benchmark: legacy corrected
#      (biased?), uncorrected (exact baseline), and the exact stochastic
#      transport-guided fix (memo 4a), all vs NUTS-non-centred ground truth.
#
# Checkpointed: writes data-raw/m1_audit_<variant>.rds per variant and
# data-raw/m1_audit_summary.rds + data-raw/m1_audit_bias.png at the end.
# Safe to re-run; finished variants are skipped.

suppressPackageStartupMessages({
  has_tidy <- requireNamespace("dplyr", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
root <- "."
if (!file.exists("R/connection_glmm.R")) {
  if (file.exists("../R/connection_glmm.R")) root <- ".."
  else stop("run from the package root or data-raw/")
}
source(file.path(root, "R", "connection_glmm.R"))
source(file.path(root, "R", "log_posterior_glmm.R"))

stan_data <- readRDS(file.path(root, "data-raw", "glmm_sparse_data.rds"))
if (!is.null(stan_data$stan_data)) stan_data <- stan_data$stan_data
J <- max(stan_data$group)

cat("Sparse benchmark: N =", length(stan_data$y), " J =", J, "\n")

# ── Analytic dA/dalpha (memo Section 2) ──────────────────────────────────────

# Returns list(A = J x 2 [mu, sigma], dA = J x 2 d/dalpha_j of each column)
.connection_and_deriv <- function(sigma, mu, alpha, X, group, beta) {
  eta <- as.vector(alpha[group] + X %*% beta)
  p   <- plogis(eta)
  s   <- as.vector(tapply(p * (1 - p), group, sum))
  sp  <- as.vector(tapply(p * (1 - p) * (1 - 2 * p), group, sum))
  G   <- 1 / sigma^2 + s
  d   <- alpha - mu
  A_mu  <- 1 / (sigma^2 * G)
  A_sig <- 2 * d / (sigma^3 * G)
  dA_mu  <- -sp / (sigma^2 * G^2)
  dA_sig <- (2 / (sigma^3 * G)) * (1 - d * sp / G)
  list(A = cbind(A_mu, A_sig), dA = cbind(dA_mu, dA_sig))
}

# ── Part A: finite-difference check ──────────────────────────────────────────

set.seed(101)
max_rel_err <- 0
for (rep in 1:20) {
  mu    <- rnorm(1); sigma <- exp(rnorm(1, 0, 0.5))
  alpha <- rnorm(J, mu, sigma); beta <- rnorm(2, 0, 0.5)
  cd <- .connection_and_deriv(sigma, mu, alpha, stan_data$X,
                              stan_data$group, beta)
  h <- 1e-6
  for (j in seq_len(J)) {
    ap <- alpha; ap[j] <- ap[j] + h
    am <- alpha; am[j] <- am[j] - h
    Gp <- .glmm_G_FF(sigma, ap, stan_data$X, stan_data$group, beta)
    Gm <- .glmm_G_FF(sigma, am, stan_data$X, stan_data$group, beta)
    Ap <- .glmm_connection(Gp, .glmm_G_BF(sigma, mu, ap))
    Am <- .glmm_connection(Gm, .glmm_G_BF(sigma, mu, am))
    fd <- (Ap[j, ] - Am[j, ]) / (2 * h)
    rel <- abs(fd - cd$dA[j, ]) / pmax(abs(fd), 1e-8)
    max_rel_err <- max(max_rel_err, rel)
  }
}
cat(sprintf("Part A: max relative error analytic vs FD dA/dalpha: %.2e\n",
            max_rel_err))
stopifnot(max_rel_err < 1e-5)

# ── Shared transport map (exactly as in horizontal_mcmc) ─────────────────────

.transport <- function(alpha, mu, sigma, dmu, dsig, beta, sd_clamp = TRUE) {
  G_FF <- .glmm_G_FF(sigma, alpha, stan_data$X, stan_data$group, beta)
  G_BF <- .glmm_G_BF(sigma, mu, alpha)
  A    <- .glmm_connection(G_FF, G_BF)
  da   <- as.vector(A %*% c(dmu, dsig))
  if (sd_clamp) {
    mc <- 3 * sigma
    da <- pmax(pmin(da, mc), -mc)
  }
  alpha + da
}

# ── Part B: proposal-level defect sizes at a representative point ───────────

set.seed(202)
mu0 <- 0; sig0 <- 1.5
al0 <- rnorm(J, mu0, sig0); be0 <- c(0, 0)
s_base <- 0.10
defect <- logdet <- numeric(200)
for (k in 1:200) {
  dmu  <- rnorm(1, 0, s_base); dls <- rnorm(1, 0, s_base)
  sig1 <- sig0 * exp(dls); dsig <- sig1 - sig0
  a1   <- .transport(al0, mu0, sig0, dmu, dsig, be0)
  aRT  <- .transport(a1, mu0 + dmu, sig1, -dmu, sig0 - sig1, be0)
  defect[k] <- sqrt(sum((aRT - al0)^2))
  cd <- .connection_and_deriv(sig0, mu0, al0, stan_data$X,
                              stan_data$group, be0)
  logdet[k] <- sum(log(abs(1 + dmu * cd$dA[, 1] + dsig * cd$dA[, 2])))
}
cat(sprintf("Part B (step %.2f): median |log det DT| = %.4f,  median round-trip defect = %.4g\n",
            s_base, median(abs(logdet)), median(defect)))

# ── Part C: long-run bias audit ──────────────────────────────────────────────

# One MwG implementation, three variants of the base block:
#   "plain"   : alpha untouched (exact)
#   "legacy"  : deterministic clamped Euler transport, no correction (current)
#   "exact4a" : transport + N(0, s_noise^2) fiber noise, MH-corrected (exact)
run_variant <- function(variant, n_iter = 40000L, n_warmup = 5000L,
                        n_chains = 4L, step_base = 0.10, step_alpha = 0.30,
                        step_beta = 0.15, s_noise = 0.20,
                        target_rate = 0.30, seed = 1) {
  set.seed(seed)
  P <- 2L + J + 2L
  out <- array(NA_real_, dim = c(n_iter, n_chains, P))
  acc_rates <- matrix(NA_real_, n_chains, 3)

  for (ch in seq_len(n_chains)) {
    mu_c <- rnorm(1, 0, 1); sigma_c <- rexp(1, 2)
    alpha_c <- rnorm(J, mu_c, sigma_c); beta_c <- rnorm(2, 0, 0.5)
    ls_c <- log(sigma_c)
    lp_c <- .glmm_log_post(mu_c, ls_c, alpha_c, beta_c, stan_data)
    s_b <- step_base; s_a <- step_alpha; s_be <- step_beta
    ab <- nb <- aa <- na_ <- abe <- nbe <- 0L

    for (it in seq_len(n_warmup + n_iter)) {
      # base block
      dmu <- rnorm(1, 0, s_b); dls <- rnorm(1, 0, s_b)
      pmu <- mu_c + dmu; pls <- ls_c + dls; psig <- exp(pls)
      dsig <- psig - sigma_c
      log_q_corr <- 0
      if (variant == "plain") {
        a_prop <- alpha_c
      } else if (variant == "legacy") {
        a_prop <- .transport(alpha_c, mu_c, sigma_c, dmu, dsig, beta_c)
      } else { # exact4a
        Tf <- .transport(alpha_c, mu_c, sigma_c, dmu, dsig, beta_c)
        a_prop <- Tf + rnorm(J, 0, s_noise)
        Tr <- .transport(a_prop, pmu, psig, -dmu, sigma_c - psig, beta_c)
        log_q_corr <- sum(dnorm(alpha_c, Tr, s_noise, log = TRUE)) -
                      sum(dnorm(a_prop, Tf, s_noise, log = TRUE))
      }
      lp_p <- tryCatch(.glmm_log_post(pmu, pls, a_prop, beta_c, stan_data),
                       error = function(e) -Inf)
      if (!is.finite(lp_p)) lp_p <- -Inf
      lr <- lp_p - lp_c + log_q_corr
      if (is.finite(lr) && log(runif(1)) < lr) {
        mu_c <- pmu; ls_c <- pls; sigma_c <- psig
        alpha_c <- a_prop; lp_c <- lp_p; ab <- ab + 1L
      }
      nb <- nb + 1L

      # fiber block (per-group RW)
      for (j in seq_len(J)) {
        atry <- alpha_c; atry[j] <- atry[j] + rnorm(1, 0, s_a)
        lpt <- tryCatch(.glmm_log_post(mu_c, ls_c, atry, beta_c, stan_data),
                        error = function(e) -Inf)
        if (!is.finite(lpt)) lpt <- -Inf
        if (log(runif(1)) < lpt - lp_c) {
          alpha_c <- atry; lp_c <- lpt; aa <- aa + 1L
        }
        na_ <- na_ + 1L
      }

      # beta block
      btry <- beta_c + rnorm(2, 0, s_be)
      lpb <- tryCatch(.glmm_log_post(mu_c, ls_c, alpha_c, btry, stan_data),
                      error = function(e) -Inf)
      if (!is.finite(lpb)) lpb <- -Inf
      if (log(runif(1)) < lpb - lp_c) {
        beta_c <- btry; lp_c <- lpb; abe <- abe + 1L
      }
      nbe <- nbe + 1L

      # warmup tuning (as in horizontal_mcmc)
      if (it <= n_warmup && it %% 100L == 0L) {
        s_b  <- min(if (ab / nb   > target_rate) s_b  * 1.3 else s_b  / 1.3, 0.5)
        s_a  <- min(if (aa / na_  > target_rate) s_a  * 1.3 else s_a  / 1.3, 0.5)
        s_be <- min(if (abe / nbe > target_rate) s_be * 1.3 else s_be / 1.3, 0.5)
        ab <- nb <- aa <- na_ <- abe <- nbe <- 0L
      }
      if (it > n_warmup) out[it - n_warmup, ch, ] <-
        c(mu_c, exp(ls_c), alpha_c, beta_c)
    }
    acc_rates[ch, ] <- c(ab / max(nb, 1), aa / max(na_, 1), abe / max(nbe, 1))
  }
  dimnames(out) <- list(NULL, NULL,
                        c("mu", "sigma", paste0("alpha[", 1:J, "]"),
                          paste0("beta[", 1:2, "]")))
  list(draws = out, acc = acc_rates)
}

variants <- c("plain", "legacy", "exact4a")
for (v in variants) {
  f <- file.path(root, "data-raw", paste0("m1_audit_", v, ".rds"))
  if (file.exists(f)) { cat("skip (cached):", v, "\n"); next }
  cat("running:", v, "...\n"); t0 <- Sys.time()
  res <- run_variant(v, seed = match(v, variants))
  cat(sprintf("  done in %.1f min; base acc %.2f\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$acc[, 1])))
  saveRDS(res, f)
}

# ── Ground truth: NUTS non-centred ───────────────────────────────────────────

nc <- readRDS(file.path(root, "data-raw", "glmm_sparse_nc_draws.rds"))
nc_mat <- if (requireNamespace("posterior", quietly = TRUE)) {
  as.matrix(posterior::as_draws_matrix(nc))
} else {
  apply(unclass(nc), 3, c)  # iterations*chains x variables
}
gt_names <- colnames(nc_mat)
get_gt <- function(var) {
  hit <- gt_names[grepl(paste0("^", gsub("\\[", "\\\\[", gsub("\\]", "\\\\]", var)), "$"),
                        gt_names)]
  if (length(hit) == 1) return(nc_mat[, hit])
  NULL
}
# reconstruct alpha from alpha_tilde if necessary
gt_draw <- function(var) {
  x <- get_gt(var)
  if (!is.null(x)) return(x)
  if (grepl("^alpha\\[", var)) {
    j  <- sub("alpha\\[(\\d+)\\]", "\\1", var)
    at <- get_gt(paste0("alpha_tilde[", j, "]"))
    if (!is.null(at)) return(get_gt("mu") + get_gt("sigma") * at)
  }
  stop("ground-truth variable not found: ", var)
}

# batch-means MCSE
mcse <- function(x, nb = 50) {
  m <- floor(length(x) / nb)
  bm <- colMeans(matrix(x[1:(m * nb)], nrow = m))
  sd(bm) / sqrt(nb)
}

vars <- c("mu", "sigma", paste0("alpha[", 1:J, "]"))
rows <- list()
for (v in variants) {
  res <- readRDS(file.path(root, "data-raw", paste0("m1_audit_", v, ".rds")))
  for (vn in vars) {
    x  <- as.vector(res$draws[, , vn])
    gt <- gt_draw(vn)
    rows[[length(rows) + 1]] <- data.frame(
      variant = v, var = vn,
      mean = mean(x), sd = sd(x),
      gt_mean = mean(gt), gt_sd = sd(gt),
      mcse_mean = mcse(as.vector(res$draws[, , vn]))
    )
  }
}
summ <- do.call(rbind, rows)
summ$bias_mean <- summ$mean - summ$gt_mean
summ$bias_sd   <- summ$sd - summ$gt_sd
summ$z_mean    <- summ$bias_mean / pmax(summ$mcse_mean, 1e-12)

# pi_j at ground-truth posterior means for stratification
pm <- vapply(vars, function(v) mean(gt_draw(v)), numeric(1))
G_pm <- .glmm_G_FF(pm["sigma"], pm[paste0("alpha[", 1:J, "]")],
                   stan_data$X, stan_data$group,
                   c(mean(gt_draw("beta[1]")), mean(gt_draw("beta[2]"))))
pi_j <- .glmm_prior_fraction(G_pm, pm["sigma"])
summ$pi_j <- NA_real_
for (j in 1:J) summ$pi_j[summ$var == paste0("alpha[", j, "]")] <- pi_j[j]

saveRDS(list(summary = summ, pi_j = pi_j,
             partA_max_rel_err = max_rel_err,
             partB = list(logdet = logdet, defect = defect)),
        file.path(root, "data-raw", "m1_audit_summary.rds"))
print(summ, digits = 3)

if (has_tidy) {
  library(ggplot2)
  d <- subset(summ, grepl("alpha", var))
  p <- ggplot(d, aes(pi_j, bias_sd, colour = variant)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2) +
    labs(x = expression(pi[j]), y = "bias in SD(alpha_j) vs NUTS-NC",
         title = "M1 audit: per-group SD bias by prior fraction") +
    theme_minimal()
  ggsave(file.path(root, "data-raw", "m1_audit_bias.png"), p,
         width = 7, height = 4.5, dpi = 150)
}
cat("M1 audit complete.\n")
