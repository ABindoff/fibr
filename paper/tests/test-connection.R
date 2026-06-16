library(posterior)

# Tiny GLMM for fast tests
tiny_stan_data <- list(
  N     = 16L,
  J     = 4L,
  group = rep(1:4, each = 4L),
  X     = matrix(c(rep(c(1,-1), 8L), rep(c(-1,1), 8L)), 16L, 2L),
  y     = c(1,0,1,0, 0,1,1,0, 1,1,0,0, 0,0,1,1)
)

# A parameter point near the posterior mode
mu_t    <- 0; ls_t  <- log(1.5); sig_t <- 1.5
alpha_t <- c(0.5, -0.5, 1.0, -1.0)
beta_t  <- c(0.3, -0.2)

# ── .glmm_G_FF ────────────────────────────────────────────────────────────────

test_that(".glmm_G_FF returns J-vector of positive values", {
  g <- .glmm_G_FF(sig_t, alpha_t, tiny_stan_data$X,
                  tiny_stan_data$group, beta_t)
  expect_length(g, 4L)
  expect_true(all(g > 0))
})

test_that(".glmm_G_FF increases with more data (larger likelihood contribution)", {
  # G_FF[j] = 1/sigma^2 + sum_i p_ij(1-p_ij); prior term dominates for sparse data
  g_sparse <- .glmm_G_FF(sig_t, alpha_t, tiny_stan_data$X[1:8, ],
                          tiny_stan_data$group[1:8], beta_t)
  g_full   <- .glmm_G_FF(sig_t, alpha_t, tiny_stan_data$X,
                          tiny_stan_data$group, beta_t)
  # Full data should give >= sparse for at least some groups
  expect_true(any(g_full >= g_sparse))
})

# ── .glmm_diag_metric ─────────────────────────────────────────────────────────

test_that(".glmm_diag_metric returns P-vector of positive values", {
  P <- 2L + 4L + 2L  # J=4
  g <- .glmm_diag_metric(mu_t, ls_t, alpha_t, beta_t, tiny_stan_data)
  expect_length(g, P)
  expect_true(all(g > 0))
  expect_true(all(is.finite(g)))
})

test_that(".glmm_diag_metric G_ls increases with larger alpha deviations", {
  alpha_far  <- c(5, -5, 5, -5)
  alpha_near <- c(0.1, -0.1, 0.1, -0.1)
  g_far  <- .glmm_diag_metric(mu_t, ls_t, alpha_far,  beta_t, tiny_stan_data)
  g_near <- .glmm_diag_metric(mu_t, ls_t, alpha_near, beta_t, tiny_stan_data)
  # G_ls = sigma + 2*sum(dev^2)/sigma^2; larger deviations -> larger G_ls
  expect_gt(g_far[2L], g_near[2L])
})

# ── .softabs_eval ─────────────────────────────────────────────────────────────

test_that(".softabs_eval is always strictly positive", {
  lambdas <- c(-10, -1, -0.1, 0, 0.1, 1, 10)
  vals    <- .softabs_eval(lambdas, alpha = 1.0)
  expect_true(all(vals > 0))
  expect_true(all(is.finite(vals)))
})

test_that(".softabs_eval(0, alpha) = 1/alpha", {
  alpha <- 2.0
  expect_equal(.softabs_eval(0, alpha), 1 / alpha, tolerance = 1e-8)
})

test_that(".softabs_eval(lambda, large_alpha) ≈ |lambda| for large |lambda|", {
  alpha <- 100
  lambda <- c(-5, 5)
  vals  <- .softabs_eval(lambda, alpha)
  expect_equal(vals, abs(lambda), tolerance = 1e-3)
})

# ── .softabs_decomp ───────────────────────────────────────────────────────────

test_that(".softabs_decomp produces PD metric", {
  blocks <- .glmm_full_metric(mu_t, ls_t, alpha_t, beta_t, tiny_stan_data)
  sa     <- .softabs_decomp(blocks, J = 4L, alpha = 1.0)

  expect_true(all(sa$lambda_sa > 0))
  expect_true(is.finite(sa$log_det))

  # Reconstruct G = U diag(lambda_sa) U^T and check PD
  G_sa <- sa$U %*% diag(sa$lambda_sa) %*% t(sa$U)
  ev   <- eigen(G_sa, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that(".softabs_decomp log_det equals sum(log(lambda_sa))", {
  blocks <- .glmm_full_metric(mu_t, ls_t, alpha_t, beta_t, tiny_stan_data)
  sa     <- .softabs_decomp(blocks, J = 4L, alpha = 1.0)
  expect_equal(sa$log_det, sum(log(sa$lambda_sa)), tolerance = 1e-10)
})

# ── compute_connection ────────────────────────────────────────────────────────

test_that("compute_connection returns fibr_connection with correct dims", {
  # Needs a draws_array — build one synthetically
  set.seed(42L)
  n_draws <- 30L
  J       <- tiny_stan_data$J

  # Simulate posterior-like draws
  draws_mat <- cbind(
    mu      = rnorm(n_draws, 0, 0.5),
    sigma   = rexp(n_draws, 1),
    matrix(rnorm(n_draws * J, 0, 1), n_draws, J,
           dimnames = list(NULL, paste0("alpha[", seq_len(J), "]"))),
    `beta[1]` = rnorm(n_draws, 0.3, 0.2),
    `beta[2]` = rnorm(n_draws, -0.2, 0.2)
  )
  # Wrap as a minimal draws_array (n_draws × 1 chain × P variables)
  arr <- posterior::as_draws_array(
    array(draws_mat, dim = c(n_draws, 1L, ncol(draws_mat)),
          dimnames = list(NULL, NULL, colnames(draws_mat)))
  )

  conn <- compute_connection(
    chain      = arr,
    base_vars  = c("mu", "sigma"),
    fiber_vars = paste0("alpha[", seq_len(J), "]"),
    method     = "analytic_glmm",
    stan_data  = tiny_stan_data,
    beta_vars  = c("beta[1]", "beta[2]"),
    n_subsample = 20L
  )

  expect_s3_class(conn, "fibr_connection")
  expect_equal(dim(conn$A), c(20L, J, 2L))
  expect_equal(ncol(conn$G_FF),       J)
  expect_equal(dim(conn$curvature),   c(20L, J))
  expect_true(all(conn$G_FF > 0))
  expect_true(all(conn$curvature < 0))   # linearised curvature always negative for GLMM
  expect_true(all(conn$prior_frac > 0 & conn$prior_frac < 1))
})

# ── Flatness: true connection has zero holonomy ───────────────────────────────
#
# These tests port the RK4 ODE from data-raw/verify_flat_connection.R.
# Single-group logistic model (n=3, random intercept only):
#   p_i = plogis(alpha + x_i),  alpha ~ N(mu, sigma^2)
#   G_FF(a,s) = 1/s^2 + sum_i p_i(1-p_i)
#   A_mu(a,s)  =  1 / (s^2 G_FF)        [frozen: same with G_0]
#   A_sig(a,s) =  2(a-mu) / (s^3 G_FF)  [frozen: same with G_0]
# RK4 integrates dalpha/dt = A_mu*dmu/dt + A_sig*dsig/dt around a circle
# of radius r in (mu,sigma) space.

local({
  set.seed(1L)
  n <- 3L; xb <- rnorm(n, 0, 1)

  S     <- function(a)     { p <- plogis(a + xb); sum(p * (1 - p)) }
  G_FF  <- function(a, s)  1 / s^2 + S(a)
  A_mu  <- function(a, s)  1 / (s^2 * G_FF(a, s))
  A_sig <- function(a, mu, s) 2 * (a - mu) / (s^3 * G_FF(a, s))

  holonomy <- function(mu0, sigma0, alpha0, r, frozen, nsteps = 20000L) {
    G0 <- G_FF(alpha0, sigma0)
    f  <- function(tk, ak) {
      m  <- mu0    + r * cos(tk)
      s  <- sigma0 + r * sin(tk)
      dm <- -r * sin(tk)
      ds <-  r * cos(tk)
      if (frozen) {
        am  <- 1 / (s^2 * G0)
        asg <- 2 * (ak - m) / (s^3 * G0)
      } else {
        am  <- A_mu(ak, s)
        asg <- A_sig(ak, m, s)
      }
      am * dm + asg * ds
    }
    th <- seq(0, 2 * pi, length.out = nsteps + 1L)
    a  <- alpha0
    for (k in seq_len(nsteps)) {
      dt <- th[k + 1L] - th[k]; t <- th[k]
      k1 <- f(t,          a)
      k2 <- f(t + dt / 2, a + dt / 2 * k1)
      k3 <- f(t + dt / 2, a + dt / 2 * k2)
      k4 <- f(t + dt,     a + dt * k3)
      a  <- a + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    }
    a - alpha0
  }

  test_that("linearised (frozen-G) holonomy is nonzero at r=0.1", {
    hol <- holonomy(mu0 = 0, sigma0 = 1.0, alpha0 = 0.0, r = 0.1,
                    frozen = TRUE)
    expect_gt(abs(hol), 0.01)
  })

  test_that("full connection is flat: true holonomy < 1e-8 at r=0.1", {
    hol <- holonomy(mu0 = 0, sigma0 = 1.0, alpha0 = 0.0, r = 0.1,
                    frozen = FALSE)
    expect_lt(abs(hol), 1e-8)
  })
})
