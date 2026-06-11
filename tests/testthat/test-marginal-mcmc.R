# M2c tests: .gh_nodes_weights, .marginal_log_post, .sample_alpha_cond,
#            marginal_mcmc()
#
# Checks:
#  1. GH nodes/weights: sum(w) = sqrt(pi), symmetry, exactness for even moments.
#  2. .marginal_log_post returns finite lp at a well-identified point.
#  3. In the prior-dominated limit (sigma -> 0), marginal log-lik -> joint
#     log-lik evaluated at alpha = mu.
#  4. .marginal_log_post is <= joint log-posterior (marginalisation reduces to
#     the marginal, which must be <= the max over alpha).
#  5. .sample_alpha_cond returns J values in the expected grid range.
#  6. smoke: marginal_mcmc() returns draws_array with correct dims.
#  7. All draws are finite; sigma draws are positive.
#  8. Reproducibility with seed.

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

gh15 <- .gh_nodes_weights(15L)

# ── Check 1: GH nodes and weights ────────────────────────────────────────────

test_that("GH weights sum to sqrt(pi) (1e-12)", {
  expect_equal(sum(gh15$w), sqrt(pi), tolerance = 1e-12)
})

test_that("GH nodes are symmetric around 0 (1e-12)", {
  n <- length(gh15$x)
  # Nodes should come in ± pairs (n=15 is odd, so one node at 0)
  expect_equal(gh15$x[1L], -gh15$x[n], tolerance = 1e-12)
  expect_equal(abs(gh15$x[ceiling(n / 2)]), 0, tolerance = 1e-12)
})

test_that("GH is exact for even moments: sum(w * x^2) = sqrt(pi)/2", {
  # int x^2 exp(-x^2) dx = Gamma(3/2) = sqrt(pi)/2
  expect_equal(sum(gh15$w * gh15$x^2), sqrt(pi) / 2, tolerance = 1e-12)
})

test_that("GH is exact for even moments: sum(w * x^4) = 3*sqrt(pi)/4", {
  # int x^4 exp(-x^2) dx = Gamma(5/2) = 3 sqrt(pi)/4
  expect_equal(sum(gh15$w * gh15$x^4), 3 * sqrt(pi) / 4, tolerance = 1e-12)
})

test_that("GH is exact for even moments: sum(w * x^6) = 15*sqrt(pi)/8", {
  # int x^6 exp(-x^2) dx = Gamma(7/2) = 15 sqrt(pi)/8
  expect_equal(sum(gh15$w * gh15$x^6), 15 * sqrt(pi) / 8, tolerance = 1e-10)
})

test_that("GH all weights are positive", {
  expect_true(all(gh15$w > 0))
})

# ── Check 2: .marginal_log_post is finite ─────────────────────────────────────

test_that(".marginal_log_post returns finite lp and m at a well-posed point", {
  r <- .marginal_log_post(theta0[1L], theta0[2L], beta0, tiny, gh15)
  expect_true(is.finite(r$lp))
  expect_false(is.null(r$m))
  expect_length(r$m, max(tiny$group))
  expect_true(all(is.finite(r$m)))
})

test_that(".marginal_log_post returns very low lp for extreme sigma", {
  # ls = 20 => sigma = exp(20) ≈ 5e8: dexp(sigma) = -sigma dominates
  r <- .marginal_log_post(0, 20, beta0, tiny, gh15)
  # The exp(-sigma) prior makes the log-posterior extremely negative
  expect_true(is.numeric(r$lp))
  expect_true(r$lp < -1e6)
})

# ── Check 3: prior-dominated limit ───────────────────────────────────────────

test_that("marginal log-lik -> joint log-lik at alpha=mu as sigma -> 0", {
  # With sigma = 0.001, all alpha_j -> mu; marginal lik -> joint(alpha=mu)
  mu0   <- 0.5; sigma0 <- 0.001; ls0 <- log(sigma0)
  beta_ <- c(0.0, 0.0)
  gh    <- .gh_nodes_weights(15L)

  r     <- .marginal_log_post(mu0, ls0, beta_, tiny, gh)

  # Joint log-lik at alpha = mu (likelihood only, no per-alpha priors)
  alpha_mu <- rep(mu0, max(tiny$group))
  eta      <- alpha_mu[tiny$group] + as.vector(tiny$X %*% beta_)
  ll_joint <- sum(tiny$y * eta - log1p(exp(eta)))

  # Base priors (to subtract from r$lp)
  lp_base <- dnorm(mu0, 0, 5, log = TRUE) +
             dexp(sigma0, rate = 1, log = TRUE) + ls0 +
             sum(dnorm(beta_, 0, 2, log = TRUE))

  # Marginal log-lik (sans base priors) should be close to ll_joint
  expect_equal(r$lp - lp_base, ll_joint, tolerance = 0.02)
})

# ── Check 4: marginal log-lik <= 0 (it's a log-probability) ──────────────────

test_that("per-group marginal log-lik is <= 0 (expectation of a probability)", {
  # sum_j log p(y_j | theta, beta) = sum_j log E_{alpha_j~prior}[p(y_j|alpha_j,beta)]
  # Each p(y_j|alpha_j,beta) <= 1, so the expectation <= 1, log <= 0.
  r <- .marginal_log_post(theta0[1L], theta0[2L], beta0, tiny, gh15)

  lp_base <- dnorm(theta0[1L], 0, 5, log = TRUE) +
             dexp(exp(theta0[2L]), rate = 1, log = TRUE) + theta0[2L] +
             sum(dnorm(beta0, 0, 2, log = TRUE))

  log_ml <- r$lp - lp_base   # marginal log-likelihood (sans base priors)
  expect_true(log_ml <= 1e-6,
              label = sprintf("log marginal lik = %.6f > 0", log_ml))
})

# ── Check 5: .sample_alpha_cond returns values in grid range ─────────────────

test_that(".sample_alpha_cond returns values in [m_j - 6*s_j, m_j + 6*s_j]", {
  set.seed(1L)
  lap <- .glmm_laplace(theta0, beta0, tiny)

  # Draw 20 independent samples and check they are all in range
  J <- max(tiny$group)
  for (rep in seq_len(20L)) {
    alpha_draw <- .sample_alpha_cond(theta0, beta0, tiny, lap)
    expect_length(alpha_draw, J)
    for (j in seq_len(J)) {
      lo <- lap$m[j] - 6 * lap$s[j]
      hi <- lap$m[j] + 6 * lap$s[j]
      expect_true(alpha_draw[j] >= lo && alpha_draw[j] <= hi,
                  label = sprintf("alpha[%d] = %.4f not in [%.4f, %.4f]",
                                  j, alpha_draw[j], lo, hi))
    }
  }
})

test_that(".sample_alpha_cond returns finite values", {
  set.seed(2L)
  lap        <- .glmm_laplace(theta0, beta0, tiny)
  alpha_draw <- .sample_alpha_cond(theta0, beta0, tiny, lap)
  expect_true(all(is.finite(alpha_draw)))
})

# ── Check 6: smoke tests ──────────────────────────────────────────────────────

test_that("marginal_mcmc returns draws_array with correct dimensions", {
  dr <- marginal_mcmc(tiny, n_iter = 10L, n_warmup = 10L, n_chains = 2L,
                      verbose = FALSE, seed = 1L)
  expect_s3_class(dr, "draws_array")
  d <- dim(dr)
  expect_equal(d[1L], 10L)
  expect_equal(d[2L], 2L)
  expect_equal(d[3L], 2L + max(tiny$group) + ncol(tiny$X))
})

test_that("marginal_mcmc variable names match convention", {
  dr <- marginal_mcmc(tiny, n_iter = 5L, n_warmup = 5L, n_chains = 1L,
                      verbose = FALSE, seed = 2L)
  vn <- posterior::variables(dr)
  expect_true("mu"       %in% vn)
  expect_true("sigma"    %in% vn)
  expect_true("alpha[1]" %in% vn)
  expect_true("beta[1]"  %in% vn)
})

# ── Check 7: finite draws and positive sigma ──────────────────────────────────

test_that("marginal_mcmc all draws are finite", {
  dr <- marginal_mcmc(tiny, n_iter = 20L, n_warmup = 20L, n_chains = 2L,
                      verbose = FALSE, seed = 3L)
  expect_true(all(is.finite(as.array(dr))))
})

test_that("marginal_mcmc sigma draws are all positive", {
  dr     <- marginal_mcmc(tiny, n_iter = 20L, n_warmup = 20L, n_chains = 1L,
                           verbose = FALSE, seed = 4L)
  sigmas <- as.vector(posterior::as_draws_matrix(
    posterior::subset_draws(dr, "sigma")))
  expect_true(all(sigmas > 0))
})

# ── Check 8: reproducibility ──────────────────────────────────────────────────

test_that("marginal_mcmc is reproducible with seed", {
  dr1 <- marginal_mcmc(tiny, n_iter = 10L, n_warmup = 10L, n_chains = 1L,
                       verbose = FALSE, seed = 99L)
  dr2 <- marginal_mcmc(tiny, n_iter = 10L, n_warmup = 10L, n_chains = 1L,
                       verbose = FALSE, seed = 99L)
  expect_equal(as.array(dr1), as.array(dr2))
})
