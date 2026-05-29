library(posterior)

# Tiny dataset for fast sampler tests
tiny_data <- list(
  N     = 16L,
  J     = 2L,
  group = rep(1:2, each = 8L),
  X     = matrix(rnorm(32L, sd = 0.5), 16L, 2L),
  y     = c(1,0,1,1,0,1,0,1, 0,1,0,0,1,0,1,0)
)

# ── .glmm_log_post ────────────────────────────────────────────────────────────

test_that(".glmm_log_post returns a finite scalar", {
  lp <- .glmm_log_post(0, log(1), c(0.5, -0.5), c(0.3, -0.2), tiny_data)
  expect_length(lp, 1L)
  expect_true(is.finite(lp))
})

test_that(".glmm_log_post decreases at extreme parameter values", {
  lp_good    <- .glmm_log_post(0,  log(1), c(0, 0), c(0, 0), tiny_data)
  lp_extreme <- .glmm_log_post(100, log(1), c(0, 0), c(0, 0), tiny_data)
  expect_lt(lp_extreme, lp_good)
})

# ── .glmm_grad_log_post ───────────────────────────────────────────────────────

test_that(".glmm_grad_log_post returns named list of correct lengths", {
  g <- .glmm_grad_log_post(0, log(1.5), c(0.5, -0.5), c(0.3, -0.2), tiny_data)
  expect_named(g, c("grad_mu", "grad_log_sigma", "grad_alpha", "grad_beta"))
  expect_length(g$grad_mu,        1L)
  expect_length(g$grad_log_sigma, 1L)
  expect_length(g$grad_alpha,     tiny_data$J)
  expect_length(g$grad_beta,      2L)
  expect_true(all(is.finite(unlist(g))))
})

test_that(".glmm_grad_vec matches finite-difference gradient", {
  mu0 <- 0.1; ls0 <- log(1.2); al0 <- c(0.4, -0.3); be0 <- c(0.2, -0.1)
  h   <- 1e-5
  J   <- tiny_data$J

  grad_an <- .glmm_grad_vec(mu0, ls0, al0, be0, tiny_data)
  lp0     <- .glmm_log_post(mu0, ls0, al0, be0, tiny_data)

  # Finite difference for mu
  fd_mu <- (.glmm_log_post(mu0 + h, ls0, al0, be0, tiny_data) - lp0) / h
  expect_equal(grad_an[1L], fd_mu, tolerance = 1e-3)

  # Finite difference for alpha[1]
  al_p      <- al0; al_p[1L] <- al_p[1L] + h
  fd_alpha1 <- (.glmm_log_post(mu0, ls0, al_p, be0, tiny_data) - lp0) / h
  expect_equal(grad_an[3L], fd_alpha1, tolerance = 1e-3)
})

# ── riemannian_mcmc ───────────────────────────────────────────────────────────

test_that("riemannian_mcmc returns a draws_array with correct dims (diagonal)", {
  draws <- riemannian_mcmc(
    tiny_data, n_iter = 10L, n_warmup = 10L, n_chains = 2L,
    L = 1L, method = "diagonal", seed = 1L, verbose = FALSE
  )
  expect_true(inherits(draws, "draws_array"))
  expect_equal(dim(draws)[1L], 10L)  # iterations
  expect_equal(dim(draws)[2L], 2L)   # chains
  expect_equal(dim(draws)[3L], 2L + tiny_data$J + 2L)  # parameters
})

test_that("riemannian_mcmc returns a draws_array with correct dims (softabs)", {
  draws <- riemannian_mcmc(
    tiny_data, n_iter = 10L, n_warmup = 10L, n_chains = 2L,
    L = 1L, method = "softabs", seed = 2L, verbose = FALSE
  )
  expect_true(inherits(draws, "draws_array"))
  expect_equal(dim(draws)[1L], 10L)
})

test_that("draws are all finite for both methods", {
  for (meth in c("diagonal", "softabs")) {
    draws <- riemannian_mcmc(
      tiny_data, n_iter = 20L, n_warmup = 20L, n_chains = 1L,
      L = 1L, method = meth, seed = 3L, verbose = FALSE
    )
    mat <- as.matrix(as_draws_matrix(draws))
    expect_true(all(is.finite(mat)),
                info = paste("method =", meth, "produced non-finite draws"))
  }
})

test_that("sigma draws are positive", {
  draws <- riemannian_mcmc(
    tiny_data, n_iter = 30L, n_warmup = 20L, n_chains = 2L,
    L = 1L, method = "diagonal", seed = 4L, verbose = FALSE
  )
  sigma_draws <- as.vector(as_draws_matrix(subset_draws(draws, "sigma")))
  expect_true(all(sigma_draws > 0))
})

test_that("variable names match Stan centred model", {
  draws <- riemannian_mcmc(
    tiny_data, n_iter = 5L, n_warmup = 5L, n_chains = 1L,
    L = 1L, seed = 5L, verbose = FALSE
  )
  vnames <- variables(draws)
  expect_true("mu"      %in% vnames)
  expect_true("sigma"   %in% vnames)
  expect_true("alpha[1]" %in% vnames)
  expect_true("beta[1]"  %in% vnames)
})

test_that("riemannian_mcmc is reproducible with the same seed", {
  args <- list(tiny_data, n_iter=20L, n_warmup=10L, n_chains=1L,
               L=1L, method="diagonal", seed=99L, verbose=FALSE)
  d1 <- do.call(riemannian_mcmc, args)
  d2 <- do.call(riemannian_mcmc, args)
  expect_equal(as.array(d1), as.array(d2))
})
