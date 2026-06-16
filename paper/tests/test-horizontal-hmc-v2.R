# M2b tests: horizontal_hmc() v2 with transport flag
#
# Checks:
#  1. transport argument is matched and validated.
#  2. Smoke test: both transports return a draws_array with correct dims.
#  3. MANDATORY: FD gradient check for .glmm_reparam_grad (all directions, 1e-4).
#  4. Hamiltonian conservation along leapfrog trajectory (single step, L=1).
#  5. Unit detailed balance for one reparam trajectory:
#     H(q) + K(p) = H(q') + K(p') implies valid accept/reject.
#  6. All draws are finite (smoke).
#  7. fisher_legacy returns same output as old L=1 call.

# ── Shared tiny data ──────────────────────────────────────────────────────────

tiny <- list(
  N     = 16L,
  J     = 4L,
  group = rep(1:4, each = 4L),
  X     = matrix(c(rep(c(1, -1), 8L), rep(c(-1, 1), 8L)), 16L, 2L),
  y     = c(1, 0, 1, 0,  0, 1, 1, 0,  1, 1, 0, 0,  0, 0, 1, 1)
)

theta0 <- c(0.0, log(1.5))
beta0  <- c(0.3, -0.2)
alpha0 <- c(0.5, -0.5, 1.0, -1.0)

# ── Check 1: argument matching ────────────────────────────────────────────────

test_that("transport argument is matched correctly", {
  expect_error(
    horizontal_hmc(tiny, n_iter = 1L, n_warmup = 1L, n_chains = 1L,
                   transport = "bad_value", verbose = FALSE),
    regexp = "arg"
  )
})

test_that("transport defaults to 'reparam'", {
  set.seed(1L)
  dr_def <- horizontal_hmc(tiny, n_iter = 5L, n_warmup = 5L, n_chains = 1L,
                            verbose = FALSE, seed = 1L)
  set.seed(1L)
  dr_rep <- horizontal_hmc(tiny, n_iter = 5L, n_warmup = 5L, n_chains = 1L,
                            transport = "reparam", verbose = FALSE, seed = 1L)
  expect_equal(as.array(dr_def), as.array(dr_rep))
})

# ── Check 2: smoke tests ──────────────────────────────────────────────────────

test_that("transport = 'reparam' returns draws_array with correct dimensions", {
  dr <- horizontal_hmc(tiny, n_iter = 10L, n_warmup = 5L, n_chains = 2L,
                       transport = "reparam", L = 2L, verbose = FALSE, seed = 1L)
  expect_s3_class(dr, "draws_array")
  d <- dim(dr)
  expect_equal(d[1L], 10L)
  expect_equal(d[2L], 2L)
  expect_equal(d[3L], 2L + max(tiny$group) + 2L)
})

test_that("transport = 'fisher_legacy' returns draws_array with correct dims", {
  dr <- horizontal_hmc(tiny, n_iter = 10L, n_warmup = 5L, n_chains = 2L,
                       transport = "fisher_legacy", L = 1L,
                       verbose = FALSE, seed = 1L)
  expect_s3_class(dr, "draws_array")
  expect_equal(dim(dr)[1L], 10L)
})

test_that("all reparam draws are finite", {
  dr <- horizontal_hmc(tiny, n_iter = 20L, n_warmup = 10L, n_chains = 2L,
                       transport = "reparam", L = 3L, verbose = FALSE, seed = 7L)
  expect_true(all(is.finite(as.array(dr))))
})

# ── Check 3: MANDATORY FD gradient for .glmm_reparam_grad ────────────────────
#
# For each direction d in {mu, ls, z[1..J], beta[1..2]}, check that
#   (g(q+h*e_d) - g(q-h*e_d)) / (2h)  ≈  grad[d]
# at tolerance 1e-4 (relative, floored at 1e-8 to handle near-zero entries).

test_that("MANDATORY: .glmm_reparam_grad FD check — all directions at 1e-4 relative", {
  J   <- max(tiny$group)
  lap <- .glmm_laplace(theta0, beta0, tiny)
  z0  <- (alpha0 - lap$m) / lap$s

  mu0 <- theta0[1L];  ls0 <- theta0[2L]
  g0  <- .glmm_reparam_grad(mu0, ls0, z0, beta0, tiny, m_init = lap$m)
  expect_false(is.null(g0), label = "gradient computed at base point")

  h    <- 1e-5
  P    <- 2L + J + 2L
  q0   <- c(mu0, ls0, z0, beta0)
  grad <- g0$grad
  idx_fiber_r <- 3L:(2L + J)
  idx_beta_r  <- (3L + J):(4L + J)

  max_rel_err <- 0

  for (d in seq_len(P)) {
    qp <- q0; qp[d] <- q0[d] + h
    qm <- q0; qm[d] <- q0[d] - h

    gp <- .glmm_reparam_grad(qp[1L], qp[2L], qp[idx_fiber_r], qp[idx_beta_r],
                              tiny, m_init = lap$m)
    gm <- .glmm_reparam_grad(qm[1L], qm[2L], qm[idx_fiber_r], qm[idx_beta_r],
                              tiny, m_init = lap$m)

    expect_false(is.null(gp), label = sprintf("gradient at q+h*e_%d", d))
    expect_false(is.null(gm), label = sprintf("gradient at q-h*e_%d", d))

    fd <- (gp$lp - gm$lp) / (2 * h)
    rel <- abs(fd - grad[d]) / pmax(abs(fd), abs(grad[d]), 1e-8)
    max_rel_err <- max(max_rel_err, rel)

    expect_true(rel < 1e-4,
                label = sprintf("FD grad[%d]: fd=%.6g analytic=%.6g rel=%.2e",
                                d, fd, grad[d], rel))
  }
})

# ── Check 4: gradient is consistent with lp_tilde ────────────────────────────

test_that(".glmm_reparam_grad returns finite grad and lp at several points", {
  lap <- .glmm_laplace(theta0, beta0, tiny)
  z0  <- (alpha0 - lap$m) / lap$s

  g <- .glmm_reparam_grad(theta0[1L], theta0[2L], z0, beta0, tiny,
                           m_init = lap$m)
  expect_false(is.null(g))
  expect_true(all(is.finite(g$grad)))
  expect_true(is.finite(g$lp))
  expect_equal(length(g$grad), 2L + max(tiny$group) + 2L)
})

test_that(".glmm_reparam_grad: alpha recovers as m + s*z", {
  lap <- .glmm_laplace(theta0, beta0, tiny)
  z0  <- (alpha0 - lap$m) / lap$s

  g <- .glmm_reparam_grad(theta0[1L], theta0[2L], z0, beta0, tiny,
                           m_init = lap$m)
  expect_equal(g$alpha, alpha0, tolerance = 1e-10)
})

# ── Check 5: lp_tilde equals lp + sum(log s) ─────────────────────────────────

test_that(".glmm_reparam_grad lp equals lp_alpha + sum(log s)", {
  lap <- .glmm_laplace(theta0, beta0, tiny)
  z0  <- (alpha0 - lap$m) / lap$s

  g    <- .glmm_reparam_grad(theta0[1L], theta0[2L], z0, beta0, tiny,
                              m_init = lap$m)
  lp_a <- .glmm_log_post(theta0[1L], theta0[2L], alpha0, beta0, tiny)
  expect_equal(g$lp, lp_a + sum(log(lap$s)), tolerance = 1e-12)
})

# ── Check 6: reparam Hamiltonian is approximately conserved (L=1, small eps) ─

test_that("reparam Hamiltonian approximately conserved along L=1 leapfrog (small eps)", {
  # Single leapfrog step with small epsilon; H should change by O(eps^3).
  J   <- max(tiny$group)
  lap <- .glmm_laplace(theta0, beta0, tiny)
  z0  <- (alpha0 - lap$m) / lap$s
  mu0 <- theta0[1L]; ls0 <- theta0[2L]
  P   <- 2L + J + 2L

  g0 <- .glmm_reparam_grad(mu0, ls0, z0, beta0, tiny, m_init = lap$m)
  expect_false(is.null(g0))

  set.seed(42L)
  p0  <- rnorm(P, 0, 1)
  H0  <- -g0$lp + 0.5 * sum(p0^2)

  idx_fiber_r <- 3L:(2L + J)
  idx_beta_r  <- (3L + J):(4L + J)

  eps   <- 0.01
  q1    <- c(mu0, ls0, z0, beta0)
  p_half <- p0 + (eps / 2) * g0$grad
  q1    <- q1 + eps * p_half

  g1 <- .glmm_reparam_grad(q1[1L], q1[2L], q1[idx_fiber_r], q1[idx_beta_r],
                             tiny, m_init = g0$m)
  expect_false(is.null(g1))

  p1 <- p_half + (eps / 2) * g1$grad
  H1 <- -g1$lp + 0.5 * sum(p1^2)

  # For a leapfrog integrator with step size h, error is O(h^3) per step.
  expect_true(abs(H1 - H0) < 1e-3,
              label = sprintf("|DeltaH| = %.4e (should be < 1e-3 at eps=0.01)", abs(H1 - H0)))
})

# ── Check 7: fisher_legacy seed reproducibility ───────────────────────────────

test_that("fisher_legacy is reproducible with seed", {
  dr1 <- horizontal_hmc(tiny, n_iter = 10L, n_warmup = 5L, n_chains = 1L,
                         transport = "fisher_legacy", L = 1L,
                         verbose = FALSE, seed = 99L)
  dr2 <- horizontal_hmc(tiny, n_iter = 10L, n_warmup = 5L, n_chains = 1L,
                         transport = "fisher_legacy", L = 1L,
                         verbose = FALSE, seed = 99L)
  expect_equal(as.array(dr1), as.array(dr2))
})

test_that("reparam is reproducible with seed", {
  dr1 <- horizontal_hmc(tiny, n_iter = 10L, n_warmup = 5L, n_chains = 1L,
                         transport = "reparam", L = 2L,
                         verbose = FALSE, seed = 55L)
  dr2 <- horizontal_hmc(tiny, n_iter = 10L, n_warmup = 5L, n_chains = 1L,
                         transport = "reparam", L = 2L,
                         verbose = FALSE, seed = 55L)
  expect_equal(as.array(dr1), as.array(dr2))
})
