# M0 premise checks: conditional transport for the centred GLMM.
#
# Three checks from HANDOFF_conditional_transport.md M0:
#  1. Pure funnel (s_j = 0): A^sigma / (d_j/sigma) = 2 exactly; F_j = -2/sigma != 0.
#  2. Mean-shift ratio: Fisher A evaluated at the Laplace mode m_j agrees with
#     dm_j/dtheta (IFT) to 1e-4.
#  3. All implicit-function derivatives match central FD at 1e-5 relative.
#
# Helper functions are defined locally here; they move to R/conditional_transport.R
# in M1. These do NOT export; no roxygen tags in this file.

# ── Local helpers ─────────────────────────────────────────────────────────────

# Newton solver for the Laplace mode m_j and scale s_j.
# Returns list(m, s, G, p, t_j) — all J-vectors except p (N-vector).
.local_laplace <- function(mu, sigma, beta, sd, m_init = NULL) {
  J     <- max(sd$group); group <- sd$group; X <- sd$X; y <- sd$y
  m     <- if (!is.null(m_init)) m_init else rep(mu, J)

  for (iter in seq_len(50L)) {
    eta  <- m[group] + as.vector(X %*% beta)
    p    <- plogis(eta)
    r_j  <- as.vector(tapply(y - p,       group, sum))
    s_jm <- as.vector(tapply(p * (1 - p), group, sum))
    G_j  <- 1 / sigma^2 + s_jm
    step <- ((mu - m) / sigma^2 + r_j) / G_j
    al   <- 1.0
    for (k in seq_len(10L)) {                # step-halving guard
      if (all(is.finite(m + al * step))) break
      al <- al / 2
    }
    m <- m + al * step
    if (max(abs(al * step)) < 1e-12) break
  }

  eta  <- m[group] + as.vector(X %*% beta)
  p    <- plogis(eta)
  s_jm <- as.vector(tapply(p * (1 - p),              group, sum))
  t_j  <- as.vector(tapply(p * (1 - p) * (1 - 2*p),  group, sum))
  G_j  <- 1 / sigma^2 + s_jm

  list(m = m, s = 1 / sqrt(G_j), G = G_j, p = p, t_j = t_j)
}

# Implicit-function derivatives of m_j and log s_j w.r.t. (mu, sigma, beta).
# Returns a named list; all quantities are at a = m_j.
.local_ift_derivs <- function(mu, sigma, beta, sd, lap) {
  J <- max(sd$group); group <- sd$group; X <- sd$X
  m <- lap$m; G_j <- lap$G; p <- lap$p; t_j <- lap$t_j

  dm_dmu  <- (1 / sigma^2) / G_j
  dm_dsig <- 2 * (m - mu) / (sigma^3 * G_j)
  dm_db   <- matrix(0, J, length(beta))
  for (k in seq_along(beta)) {
    wt <- p * (1 - p) * X[, k]
    dm_db[, k] <- -as.vector(tapply(wt, group, sum)) / G_j
  }

  dG_dmu  <- t_j * dm_dmu
  dG_dsig <- -2 / sigma^3 + t_j * dm_dsig
  dG_db   <- matrix(0, J, length(beta))
  for (k in seq_along(beta)) {
    wt <- p * (1 - p) * (1 - 2*p) * X[, k]
    dG_db[, k] <- as.vector(tapply(wt, group, sum)) + t_j * dm_db[, k]
  }

  dls_dmu  <- -0.5 * dG_dmu  / G_j
  dls_dsig <- -0.5 * dG_dsig / G_j
  dls_db   <- -0.5 * dG_db   / G_j

  list(dm_dmu  = dm_dmu,   dm_dsig  = dm_dsig,   dm_db  = dm_db,
       dls_dmu = dls_dmu,  dls_dsig = dls_dsig,  dls_db = dls_db,
       t_j = t_j)
}

# ── Tiny test data (J=4, n_j=4 each) ─────────────────────────────────────────

tiny <- list(
  N     = 16L,
  J     = 4L,
  group = rep(1:4, each = 4L),
  X     = matrix(c(rep(c(1,-1), 8L), rep(c(-1,1), 8L)), 16L, 2L),
  y     = c(1,0,1,0, 0,1,1,0, 1,1,0,0, 0,0,1,1)
)
mu_t <- 0; sig_t <- 1.5; beta_t <- c(0.3, -0.2)
alpha_t <- c(0.5, -0.5, 1.0, -1.0)

# ── Check 1: Pure funnel ─────────────────────────────────────────────────────

test_that("Pure funnel (direct formulas): A^sigma/(d_j/sigma) = 2, F_j = -2/sigma", {
  # Package functions don't support empty-group data (tapply drops missing levels);
  # we implement the pure-funnel limit (s_j = 0) directly.
  sigma <- 1.5; mu <- 0.3; alpha <- c(0.7, -1.2, 0.4)
  d_j   <- alpha - mu
  G_j   <- rep(1 / sigma^2, 3L)

  A_sigma <- 2 * d_j / (sigma^3 * G_j)          # = 2 d_j / sigma
  A_mu    <- 1 / (sigma^2 * G_j)                 # = 1
  F_j     <- -2 / (sigma^5 * G_j^2)              # = -2 / sigma

  expect_equal(A_sigma / (d_j / sigma), rep(2, 3L), tolerance = 1e-12)
  expect_equal(A_mu,                    rep(1, 3L),  tolerance = 1e-12)
  expect_equal(F_j,                     rep(-2 / sigma, 3L), tolerance = 1e-12)
  expect_true(all(F_j != 0))
})

test_that("Package functions A^sigma/ratio -> 2 at large |alpha| (s_j -> 0)", {
  # With alpha >> 1, p*(1-p) -> 0, G_j -> 1/sigma^2 and ratio -> 2
  sigma <- 1.5; mu <- 0.0
  alpha_ext <- c(12.0, -12.0, 12.0, -12.0)       # extreme: p*(1-p) < 1e-9 per obs
  beta  <- c(0.0, 0.0)
  X_z   <- matrix(0.0, nrow = 4L, ncol = 2L)     # zero design: eta = alpha_j
  g_vec <- 1:4
  tiny_z <- list(X = X_z, y = c(1L,0L,1L,0L), group = g_vec)

  G_j <- .glmm_G_FF(sigma, alpha_ext, tiny_z$X, tiny_z$group, beta)
  A   <- .glmm_connection(G_j, .glmm_G_BF(sigma, mu, alpha_ext))
  F_j <- .glmm_curvature_linearised(G_j, sigma)

  ratio <- A[, 2L] / ((alpha_ext - mu) / sigma)
  expect_equal(ratio, rep(2, 4L), tolerance = 1e-4)
  expect_true(all(F_j < 0))
  expect_equal(F_j, rep(-2 / sigma, 4L), tolerance = 1e-4)
})

# ── Check 2: Mean-shift ratio = 1 ────────────────────────────────────────────

test_that("Mean-shift ratio: Fisher A at m_j equals dm_j/dtheta (1e-4)", {
  lap <- .local_laplace(mu_t, sig_t, beta_t, tiny)
  d   <- .local_ift_derivs(mu_t, sig_t, beta_t, tiny, lap)
  m   <- lap$m; G_j <- lap$G

  A_mu  <- 1 / (sig_t^2 * G_j)
  A_sig <- 2 * (m - mu_t) / (sig_t^3 * G_j)

  # A_j^mu at alpha_j = m_j equals dm_j/dmu
  expect_equal(A_mu,  d$dm_dmu,  tolerance = 1e-4)
  # A_j^sigma at alpha_j = m_j equals dm_j/dsigma
  expect_equal(A_sig, d$dm_dsig, tolerance = 1e-4)
})

# ── Check 3: Implicit-function derivatives vs central FD ──────────────────────

test_that("dm/dmu and dlog_s/dmu match FD at 1e-5 relative tolerance", {
  h <- 1e-5
  lap <- .local_laplace(mu_t, sig_t, beta_t, tiny)
  d   <- .local_ift_derivs(mu_t, sig_t, beta_t, tiny, lap)

  lap_p <- .local_laplace(mu_t + h, sig_t, beta_t, tiny, m_init = lap$m)
  lap_m <- .local_laplace(mu_t - h, sig_t, beta_t, tiny, m_init = lap$m)

  fd_dm  <- (lap_p$m - lap_m$m)   / (2 * h)
  fd_dls <- (log(lap_p$s) - log(lap_m$s)) / (2 * h)

  rel_dm  <- abs(fd_dm  - d$dm_dmu)  / pmax(abs(fd_dm),  abs(d$dm_dmu),  1e-8)
  rel_dls <- abs(fd_dls - d$dls_dmu) / pmax(abs(fd_dls), abs(d$dls_dmu), 1e-6)
  expect_true(max(rel_dm)  < 1e-5, label = "dm/dmu max rel err")
  expect_true(max(rel_dls) < 1e-5, label = "dlog_s/dmu max rel err")
})

test_that("dm/dsigma and dlog_s/dsigma match FD at 1e-5 relative tolerance", {
  h <- 1e-5
  lap <- .local_laplace(mu_t, sig_t, beta_t, tiny)
  d   <- .local_ift_derivs(mu_t, sig_t, beta_t, tiny, lap)

  lap_p <- .local_laplace(mu_t, sig_t + h, beta_t, tiny, m_init = lap$m)
  lap_m <- .local_laplace(mu_t, sig_t - h, beta_t, tiny, m_init = lap$m)

  fd_dm  <- (lap_p$m - lap_m$m)           / (2 * h)
  fd_dls <- (log(lap_p$s) - log(lap_m$s)) / (2 * h)

  rel_dm  <- abs(fd_dm  - d$dm_dsig)  / pmax(abs(fd_dm),  1e-8)
  rel_dls <- abs(fd_dls - d$dls_dsig) / pmax(abs(fd_dls), 1e-8)
  expect_true(max(rel_dm)  < 1e-5, label = "dm/dsigma max rel err")
  expect_true(max(rel_dls) < 1e-5, label = "dlog_s/dsigma max rel err")
})

test_that("dm/dbeta and dlog_s/dbeta match FD at 1e-5 relative tolerance", {
  h <- 1e-5
  lap <- .local_laplace(mu_t, sig_t, beta_t, tiny)
  d   <- .local_ift_derivs(mu_t, sig_t, beta_t, tiny, lap)

  for (k in 1:2) {
    bp <- beta_t; bp[k] <- bp[k] + h
    bm <- beta_t; bm[k] <- bm[k] - h
    lap_p <- .local_laplace(mu_t, sig_t, bp, tiny, m_init = lap$m)
    lap_m <- .local_laplace(mu_t, sig_t, bm, tiny, m_init = lap$m)

    fd_dm  <- (lap_p$m - lap_m$m)           / (2 * h)
    fd_dls <- (log(lap_p$s) - log(lap_m$s)) / (2 * h)

    rel_dm  <- abs(fd_dm  - d$dm_db[, k])  / pmax(abs(fd_dm),  1e-8)
    rel_dls <- abs(fd_dls - d$dls_db[, k]) / pmax(abs(fd_dls), 1e-8)
    expect_true(max(rel_dm)  < 1e-5,
                label = sprintf("dm/dbeta[%d] max rel err", k))
    expect_true(max(rel_dls) < 1e-5,
                label = sprintf("dlog_s/dbeta[%d] max rel err", k))
  }
})

test_that("Newton solver converges from a distant start", {
  # m_init = 0 (far from mode for extreme sigma)
  lap_cold <- .local_laplace(mu_t, 0.5, beta_t, tiny, m_init = rep(0, 4L))
  lap_warm <- .local_laplace(mu_t, 0.5, beta_t, tiny, m_init = alpha_t)
  expect_equal(lap_cold$m, lap_warm$m, tolerance = 1e-8)
})
