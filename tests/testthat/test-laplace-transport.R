# M1 tests: .glmm_laplace, .glmm_cond_transport, .glmm_lap_derivs, .glmm_kappa
#
# Checks (per HANDOFF_conditional_transport.md M1):
#  1. Round-trip exactness: T_{theta->theta'} then T_{theta'->theta} recovers alpha to 1e-10.
#  2. Log-Jacobian equals log|det| computed by FD of the map.
#  3. sigma -> 0 limit: transport approaches non-centering (NC) map.
#  4. Derivative FD checks: dm and dlog_s in all directions at 1e-5 relative.
#  5. Newton convergence from a bad start.
#  6. .glmm_kappa returns a finite J-vector with expected sign.

# ── Tiny test data (J = 4, n_j = 4 each) ─────────────────────────────────────

tiny <- list(
  N     = 16L,
  J     = 4L,
  group = rep(1:4, each = 4L),
  X     = matrix(c(rep(c(1, -1), 8L), rep(c(-1, 1), 8L)), 16L, 2L),
  y     = c(1, 0, 1, 0,  0, 1, 1, 0,  1, 1, 0, 0,  0, 0, 1, 1)
)

theta0 <- c(0.0, log(1.5))          # (mu, log_sigma)
beta0  <- c(0.3, -0.2)
alpha0 <- c(0.5, -0.5, 1.0, -1.0)

theta1 <- c(0.4, log(0.9))          # a second distinct base point

# ── Check 1: Round-trip exactness ─────────────────────────────────────────────

test_that("round-trip T_{theta->theta'} then T_{theta'->theta} is exact (1e-10)", {
  lap0 <- .glmm_laplace(theta0, beta0, tiny)
  lap1 <- .glmm_laplace(theta1, beta0, tiny)

  fwd  <- .glmm_cond_transport(alpha0, lap0, lap1)
  back <- .glmm_cond_transport(fwd$alpha_new, lap1, lap0)

  expect_equal(back$alpha_new, alpha0, tolerance = 1e-10)
})

test_that("round-trip log-Jacobians cancel to zero (1e-12)", {
  lap0 <- .glmm_laplace(theta0, beta0, tiny)
  lap1 <- .glmm_laplace(theta1, beta0, tiny)

  fwd  <- .glmm_cond_transport(alpha0, lap0, lap1)
  back <- .glmm_cond_transport(fwd$alpha_new, lap1, lap0)

  expect_equal(fwd$log_jac + back$log_jac, 0, tolerance = 1e-12)
})

# ── Check 2: Log-Jacobian equals log|det DT| from FD ─────────────────────────

test_that("log_jac equals log|det DT| (FD, 1e-4)", {
  J    <- max(tiny$group)
  lap0 <- .glmm_laplace(theta0, beta0, tiny)
  lap1 <- .glmm_laplace(theta1, beta0, tiny)

  h <- 1e-5
  # Build Jacobian matrix by perturbing each alpha component
  DT <- matrix(0, J, J)
  for (k in seq_len(J)) {
    ep <- rep(0, J); ep[k] <- h
    T_p <- .glmm_cond_transport(alpha0 + ep, lap0, lap1)$alpha_new
    T_m <- .glmm_cond_transport(alpha0 - ep, lap0, lap1)$alpha_new
    DT[, k] <- (T_p - T_m) / (2 * h)
  }

  log_det_fd <- determinant(DT, logarithm = TRUE)$modulus
  expect_equal(
    .glmm_cond_transport(alpha0, lap0, lap1)$log_jac,
    as.numeric(log_det_fd),
    tolerance = 1e-4
  )
})

# ── Check 3: sigma -> 0 approaches the NC map ─────────────────────────────────

test_that("transport approaches NC map as sigma -> 0 (no data)", {
  # NC map: alpha'_j = mu' + (alpha_j - mu) [pure funnel: m_j = mu, s_j = sigma]
  # Use zero design matrix so likelihood plays no role
  tiny_nc <- list(
    N = 4L, J = 4L,
    group = 1:4,
    X = matrix(0, 4L, 2L),
    y = c(0L, 0L, 0L, 0L)
  )
  theta_from <- c(0.5, log(0.001))    # sigma = 0.001; O(sigma^2) NC deviation < 1e-5
  theta_to   <- c(1.2, log(0.001))
  beta_nc    <- c(0, 0)

  lap_from <- .glmm_laplace(theta_from, beta_nc, tiny_nc)
  lap_to   <- .glmm_laplace(theta_to,   beta_nc, tiny_nc)

  T_cond <- .glmm_cond_transport(alpha0, lap_from, lap_to)$alpha_new

  # NC map: alpha' = mu' + (alpha - mu)
  mu_from  <- theta_from[1L];  mu_to <- theta_to[1L]
  T_nc     <- mu_to + (alpha0 - mu_from)

  expect_equal(T_cond, T_nc, tolerance = 1e-6)
})

# ── Check 4: IFT derivatives vs central FD ────────────────────────────────────

test_that("dm/dmu and dlog_s/dmu match FD at 1e-5 relative", {
  h    <- 1e-5
  lap0 <- .glmm_laplace(theta0, beta0, tiny)
  d0   <- .glmm_lap_derivs(theta0, beta0, tiny, lap0)

  theta_p <- theta0; theta_p[1L] <- theta0[1L] + h
  theta_m <- theta0; theta_m[1L] <- theta0[1L] - h
  lap_p <- .glmm_laplace(theta_p, beta0, tiny, m_init = lap0$m)
  lap_m <- .glmm_laplace(theta_m, beta0, tiny, m_init = lap0$m)

  fd_dm  <- (lap_p$m - lap_m$m) / (2 * h)
  fd_dls <- (log(lap_p$s) - log(lap_m$s)) / (2 * h)

  rel_dm  <- abs(fd_dm  - d0$dm_dmu)  / pmax(abs(fd_dm),  abs(d0$dm_dmu),  1e-8)
  rel_dls <- abs(fd_dls - d0$dls_dmu) / pmax(abs(fd_dls), abs(d0$dls_dmu), 1e-6)
  expect_true(max(rel_dm)  < 1e-5, label = "dm/dmu max rel err")
  expect_true(max(rel_dls) < 1e-5, label = "dlog_s/dmu max rel err")
})

test_that("dm/dsigma and dlog_s/dsigma match FD at 1e-5 relative", {
  # Perturb sigma directly (theta[2] = log_sigma -> theta[2] +/- h changes sigma)
  h    <- 1e-5
  sigma0 <- exp(theta0[2L])
  lap0 <- .glmm_laplace(theta0, beta0, tiny)
  d0   <- .glmm_lap_derivs(theta0, beta0, tiny, lap0)

  # Perturb sigma (not log_sigma) to match the dsigma derivative
  theta_ps <- c(theta0[1L], log(sigma0 + h))
  theta_ms <- c(theta0[1L], log(sigma0 - h))
  lap_ps <- .glmm_laplace(theta_ps, beta0, tiny, m_init = lap0$m)
  lap_ms <- .glmm_laplace(theta_ms, beta0, tiny, m_init = lap0$m)

  fd_dm  <- (lap_ps$m - lap_ms$m) / (2 * h)
  fd_dls <- (log(lap_ps$s) - log(lap_ms$s)) / (2 * h)

  rel_dm  <- abs(fd_dm  - d0$dm_dsig)  / pmax(abs(fd_dm),  abs(d0$dm_dsig),  1e-8)
  rel_dls <- abs(fd_dls - d0$dls_dsig) / pmax(abs(fd_dls), abs(d0$dls_dsig), 1e-8)
  expect_true(max(rel_dm)  < 1e-5, label = "dm/dsigma max rel err")
  expect_true(max(rel_dls) < 1e-5, label = "dlog_s/dsigma max rel err")
})

test_that("dm/dbeta and dlog_s/dbeta match FD at 1e-5 relative", {
  h    <- 1e-5
  lap0 <- .glmm_laplace(theta0, beta0, tiny)
  d0   <- .glmm_lap_derivs(theta0, beta0, tiny, lap0)

  for (k in 1:2) {
    bp <- beta0; bp[k] <- bp[k] + h
    bm <- beta0; bm[k] <- bm[k] - h
    lap_p <- .glmm_laplace(theta0, bp, tiny, m_init = lap0$m)
    lap_m <- .glmm_laplace(theta0, bm, tiny, m_init = lap0$m)

    fd_dm  <- (lap_p$m - lap_m$m) / (2 * h)
    fd_dls <- (log(lap_p$s) - log(lap_m$s)) / (2 * h)

    rel_dm  <- abs(fd_dm  - d0$dm_db[, k])  / pmax(abs(fd_dm),  abs(d0$dm_db[, k]),  1e-8)
    rel_dls <- abs(fd_dls - d0$dls_db[, k]) / pmax(abs(fd_dls), abs(d0$dls_db[, k]), 1e-6)
    expect_true(max(rel_dm)  < 1e-5, label = sprintf("dm/dbeta[%d] max rel err",  k))
    expect_true(max(rel_dls) < 1e-5, label = sprintf("dls/dbeta[%d] max rel err", k))
  }
})

# ── Check 5: Newton convergence from bad start ────────────────────────────────

test_that("Newton converges from a distant start (cold vs warm agree to 1e-8)", {
  lap_cold <- .glmm_laplace(theta0, beta0, tiny, m_init = rep(5, 4L))
  lap_warm <- .glmm_laplace(theta0, beta0, tiny)    # default warm start = mu
  expect_equal(lap_cold$m, lap_warm$m, tolerance = 1e-8)
  expect_equal(lap_cold$s, lap_warm$s, tolerance = 1e-8)
})

test_that("Newton converges at extreme sigma (funnel neck)", {
  theta_funnel <- c(0.0, log(0.05))   # sigma = 0.05
  lap_cold <- .glmm_laplace(theta_funnel, beta0, tiny, m_init = rep(3, 4L))
  lap_warm <- .glmm_laplace(theta_funnel, beta0, tiny)
  expect_equal(lap_cold$m, lap_warm$m, tolerance = 1e-6)
})

# ── Check 6: .glmm_kappa finite and correct sign ──────────────────────────────

test_that(".glmm_kappa returns finite J-vector; sign consistent with t_j", {
  lap  <- .glmm_laplace(theta0, beta0, tiny)
  kap  <- .glmm_kappa(lap, beta0, tiny)

  expect_length(kap, max(tiny$group))
  expect_true(all(is.finite(kap)))

  # kappa_j = -t_j * s_j^3; s_j^3 > 0 always, so sign(kappa) = -sign(t_j)
  eta  <- lap$m[tiny$group] + as.vector(tiny$X %*% beta0)
  p    <- plogis(eta)
  t_j  <- as.vector(tapply(p * (1 - p) * (1 - 2 * p), tiny$group, sum))
  expect_equal(sign(kap), -sign(t_j))
})

test_that(".glmm_kappa is near zero for prior-dominated groups (sigma -> 0)", {
  # With very small sigma, all data is irrelevant, m_j -> mu, p -> constant,
  # t_j -> 0 because p*(1-p)*(1-2p) is evaluated at near-constant p
  theta_prior <- c(0.0, log(0.01))
  lap <- .glmm_laplace(theta_prior, beta0, tiny)
  kap <- .glmm_kappa(lap, beta0, tiny)
  # All kappa should be small (prior-dominated -> Gaussian conditional)
  expect_true(all(abs(kap) < 0.01))
})
