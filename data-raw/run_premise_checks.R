## M0 premise checks — reproduces manuscript2/M1_exactness_memo.md Section 7.
##
## Checks:
##  1. Pure funnel: A^sigma/(d_j/sigma) = 2 analytically.
##  2. Mean-shift ratio = 1.000 for all groups at all tested sigma values.
##  3. Deviation-scaling ratio table (sigma in {0.3, 1, 1.73, 4, 8}).
##     Expected from memo: ~2.00 at sigma=0.3; 1.25-1.88 at sigma=1.73;
##     drops below 1 at large sigma.
##
## STOP and report to Aidan if any assertion fails or if numbers disagree
## with memo Section 7 values.
##
## Run from package root:  Rscript data-raw/run_premise_checks.R

root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)
source(file.path(root, "R", "connection_glmm.R"))
source(file.path(root, "R", "log_posterior_glmm.R"))

# ── Load sparse benchmark ──────────────────────────────────────────────────────

raw       <- readRDS(file.path(root, "data-raw", "glmm_sparse_data.rds"))
stan_data <- if (!is.null(raw$stan_data)) raw$stan_data else raw
J         <- max(stan_data$group)
cat(sprintf("Sparse benchmark: N = %d, J = %d\n\n", length(stan_data$y), J))

# ── Local Newton / IFT helpers (same logic as test-conditional-transport.R) ───

.newton_laplace <- function(mu, sigma, beta, sd, m_init = NULL) {
  group <- sd$group; X <- sd$X; y <- sd$y
  J     <- max(group)
  m     <- if (!is.null(m_init)) m_init else rep(mu, J)

  for (iter in seq_len(50L)) {
    eta  <- m[group] + as.vector(X %*% beta)
    p    <- plogis(eta)
    r_j  <- as.vector(tapply(y - p,       group, sum))
    s_jm <- as.vector(tapply(p * (1-p),   group, sum))
    G_j  <- 1 / sigma^2 + s_jm
    step <- ((mu - m) / sigma^2 + r_j) / G_j
    al   <- 1.0
    for (k in seq_len(10L)) { if (all(is.finite(m + al*step))) break; al <- al/2 }
    m <- m + al * step
    if (max(abs(al * step)) < 1e-12) break
  }

  eta  <- m[group] + as.vector(X %*% beta)
  p    <- plogis(eta)
  s_jm <- as.vector(tapply(p * (1-p),             group, sum))
  t_j  <- as.vector(tapply(p * (1-p) * (1-2*p),  group, sum))
  G_j  <- 1 / sigma^2 + s_jm
  list(m = m, s = 1/sqrt(G_j), G = G_j, p = p, t_j = t_j)
}

.ift_derivs <- function(mu, sigma, beta, sd, lap) {
  group <- sd$group; X <- sd$X
  J     <- max(group)
  m <- lap$m; G_j <- lap$G; p <- lap$p; t_j <- lap$t_j

  dm_dmu  <- (1 / sigma^2) / G_j
  dm_dsig <- 2 * (m - mu) / (sigma^3 * G_j)
  dm_db   <- matrix(0, J, length(beta))
  for (k in seq_along(beta)) {
    wt <- p * (1-p) * X[,k]
    dm_db[,k] <- -as.vector(tapply(wt, group, sum)) / G_j
  }
  dG_dmu  <- t_j * dm_dmu
  dG_dsig <- -2/sigma^3 + t_j * dm_dsig
  dG_db   <- matrix(0, J, length(beta))
  for (k in seq_along(beta)) {
    wt <- p * (1-p) * (1-2*p) * X[,k]
    dG_db[,k] <- as.vector(tapply(wt, group, sum)) + t_j * dm_db[,k]
  }
  list(dm_dmu  = dm_dmu,   dm_dsig  = dm_dsig,
       dls_dmu = -0.5 * dG_dmu  / G_j,
       dls_dsig= -0.5 * dG_dsig / G_j)
}

# ── CHECK 1: Pure funnel (direct formula) ─────────────────────────────────────

cat("── CHECK 1: Pure funnel (s_j = 0) ─────────────────────────────────────────\n")
cat("Expected: A^sigma / (d_j/sigma) = 2.000, F_j = -2/sigma\n")
sigma_pf <- 1.5; mu_pf <- 0.3; alpha_pf <- c(0.7, -1.2, 0.4)
d_pf <- alpha_pf - mu_pf
G_pf <- rep(1/sigma_pf^2, 3L)
ratio_pf <- (2*d_pf/(sigma_pf^3*G_pf)) / (d_pf/sigma_pf)
F_pf     <- -2/(sigma_pf^5 * G_pf^2)
cat(sprintf("  A^sigma/(d_j/sigma): %s  (all should be 2.000)\n",
            paste(round(ratio_pf, 6), collapse = ", ")))
cat(sprintf("  F_j / (-2/sigma):    %s  (all should be 1.000)\n",
            paste(round(F_pf / (-2/sigma_pf), 6), collapse = ", ")))
stopifnot(all(abs(ratio_pf - 2) < 1e-10))
stopifnot(all(abs(F_pf / (-2/sigma_pf) - 1) < 1e-10))
cat("  PASS\n\n")

# ── SHARED PARAMETER POINT for checks 2 and 3 ─────────────────────────────────

mu0   <- 0.0
beta0 <- c(0.0, 0.0)   # zero fixed effects for clean decomposition

# ── CHECK 2 & 3: Deviation-scaling and mean-shift tables ──────────────────────

sigma_grid <- c(0.3, 1.0, 1.73, 4.0, 8.0)

cat("── CHECK 2: Mean-shift ratio = 1.000 (Fisher A at m_j vs dm_j/dtheta) ─────\n")
cat(sprintf("  %-8s  %-30s  %-30s\n",
            "sigma", "ratio(mu): range [min, max]", "ratio(sigma): range [min, max]"))

cat("── CHECK 3: Deviation-scaling ratio table ───────────────────────────────────\n")
cat("  Fisher dev coeff: 2/(sigma^3 G_j)\n")
cat("  Correct dev coeff: dlog s_j/dsigma\n")
cat(sprintf("  %-8s  %-40s\n",
            "sigma", "ratio = Fisher_dev / correct_dev  (per-group range)"))

for (sigma in sigma_grid) {
  lap <- .newton_laplace(mu0, sigma, beta0, stan_data)
  d   <- .ift_derivs(mu0, sigma, beta0, stan_data, lap)

  # Mean-shift ratios (Fisher A at m_j vs IFT dm/dtheta)
  A_mu  <- 1 / (sigma^2 * lap$G)
  A_sig <- 2 * (lap$m - mu0) / (sigma^3 * lap$G)
  ratio_mu  <- A_mu  / d$dm_dmu
  ratio_sig <- ifelse(abs(d$dm_dsig) < 1e-10, 1,
                      A_sig / d$dm_dsig)   # 0/0 -> 1 when d_j = 0
  cat(sprintf("  sigma=%-4.2f  mu-dir [%.4f, %.4f]  sigma-dir [%.4f, %.4f]\n",
              sigma,
              min(ratio_mu), max(ratio_mu),
              min(ratio_sig[is.finite(ratio_sig)]),
              max(ratio_sig[is.finite(ratio_sig)])))

  stopifnot(all(abs(ratio_mu - 1) < 1e-4))
  stopifnot(all(abs(ratio_sig[is.finite(ratio_sig)] - 1) < 1e-4))

  # Deviation-scaling ratio
  fisher_dev   <- 2 / (sigma^3 * lap$G)          # Fisher coeff for d_j in sigma direction
  correct_dev  <- d$dls_dsig                      # d(log s_j)/dsigma
  dev_ratio    <- fisher_dev / correct_dev        # should -> 2 at sigma=0.3
  cat(sprintf("  sigma=%-4.2f  dev-scaling ratio [%.3f, %.3f]\n",
              sigma, min(dev_ratio), max(dev_ratio)))
}
cat("\n")

cat("── Expected from memo Section 7 ────────────────────────────────────────────\n")
cat("  Mean-shift ratio: 1.000 in both directions at all sigma  (algebraically exact)\n")
cat("  Dev-scaling ratio: ~2.00 at sigma=0.3; 1.25-1.88 at sigma=1.73;\n")
cat("                     drops below 1 at sigma=4 or 8\n\n")

# ── CHECK 3b: FD verification of IFT derivative formulas ──────────────────────

cat("── CHECK 3b: IFT derivatives vs central FD (tolerance 1e-5 relative) ───────\n")
sigma_test <- 1.73; h <- 1e-5
lap0  <- .newton_laplace(mu0, sigma_test, beta0, stan_data)
d0    <- .ift_derivs(mu0, sigma_test, beta0, stan_data, lap0)

# mu direction
lap_pm <- .newton_laplace(mu0 + h, sigma_test, beta0, stan_data, m_init = lap0$m)
lap_mm <- .newton_laplace(mu0 - h, sigma_test, beta0, stan_data, m_init = lap0$m)
fd_dm_dmu  <- (lap_pm$m - lap_mm$m)           / (2*h)
fd_dls_dmu <- (log(lap_pm$s) - log(lap_mm$s)) / (2*h)
err_dm_mu  <- max(abs(fd_dm_dmu  - d0$dm_dmu)  / pmax(abs(fd_dm_dmu),  1e-8))
err_dls_mu <- max(abs(fd_dls_dmu - d0$dls_dmu) / pmax(abs(fd_dls_dmu), 1e-8))
cat(sprintf("  dm/dmu    max rel err: %.2e  %s\n", err_dm_mu,
            if (err_dm_mu  < 1e-5) "PASS" else "*** FAIL ***"))
cat(sprintf("  dls/dmu   max rel err: %.2e  %s\n", err_dls_mu,
            if (err_dls_mu < 1e-5) "PASS" else "*** FAIL ***"))

# sigma direction
lap_ps <- .newton_laplace(mu0, sigma_test + h, beta0, stan_data, m_init = lap0$m)
lap_ms <- .newton_laplace(mu0, sigma_test - h, beta0, stan_data, m_init = lap0$m)
fd_dm_dsig  <- (lap_ps$m - lap_ms$m)           / (2*h)
fd_dls_dsig <- (log(lap_ps$s) - log(lap_ms$s)) / (2*h)
err_dm_sig  <- max(abs(fd_dm_dsig  - d0$dm_dsig)  / pmax(abs(fd_dm_dsig),  1e-8))
err_dls_sig <- max(abs(fd_dls_dsig - d0$dls_dsig) / pmax(abs(fd_dls_dsig), 1e-8))
cat(sprintf("  dm/dsigma max rel err: %.2e  %s\n", err_dm_sig,
            if (err_dm_sig  < 1e-5) "PASS" else "*** FAIL ***"))
cat(sprintf("  dls/dsig  max rel err: %.2e  %s\n", err_dls_sig,
            if (err_dls_sig < 1e-5) "PASS" else "*** FAIL ***"))

stopifnot(err_dm_mu  < 1e-5, err_dls_mu < 1e-5,
          err_dm_sig < 1e-5, err_dls_sig < 1e-5)

cat("\nAll M0 premise checks PASSED. Proceed to M1.\n")
