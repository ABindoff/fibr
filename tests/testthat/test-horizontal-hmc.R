library(posterior)

tiny_data <- list(
  N     = 16L,
  J     = 2L,
  group = rep(1:2, each = 8L),
  X     = matrix(rnorm(32L, sd = 0.5), 16L, 2L),
  y     = c(1,0,1,1,0,1,0,1, 0,1,0,0,1,0,1,0)
)

test_that("horizontal_hmc returns a draws_array with correct dims", {
  draws <- horizontal_hmc(
    tiny_data, n_iter = 10L, n_warmup = 10L, n_chains = 2L,
    L = 1L, seed = 1L, verbose = FALSE
  )
  expect_true(inherits(draws, "draws_array"))
  expect_equal(dim(draws)[1L], 10L)
  expect_equal(dim(draws)[2L], 2L)
  expect_equal(dim(draws)[3L], 2L + tiny_data$J + 2L)
})

test_that("horizontal_hmc draws are all finite", {
  draws <- horizontal_hmc(
    tiny_data, n_iter = 20L, n_warmup = 20L, n_chains = 1L,
    L = 1L, seed = 2L, verbose = FALSE
  )
  expect_true(all(is.finite(as.matrix(as_draws_matrix(draws)))))
})

test_that("horizontal_hmc sigma draws are positive", {
  draws <- horizontal_hmc(
    tiny_data, n_iter = 30L, n_warmup = 20L, n_chains = 2L,
    L = 1L, seed = 3L, verbose = FALSE
  )
  sigma_draws <- as.vector(as_draws_matrix(subset_draws(draws, "sigma")))
  expect_true(all(sigma_draws > 0))
})

test_that("horizontal_hmc variable names match Stan centred model", {
  draws <- horizontal_hmc(
    tiny_data, n_iter = 5L, n_warmup = 5L, n_chains = 1L,
    L = 1L, seed = 4L, verbose = FALSE
  )
  vnames <- variables(draws)
  expect_true("mu"       %in% vnames)
  expect_true("sigma"    %in% vnames)
  expect_true("alpha[1]" %in% vnames)
  expect_true("beta[1]"  %in% vnames)
})

test_that("horizontal_hmc is reproducible with the same seed", {
  args <- list(tiny_data, n_iter = 20L, n_warmup = 10L, n_chains = 1L,
               L = 1L, seed = 99L, verbose = FALSE)
  d1 <- do.call(horizontal_hmc, args)
  d2 <- do.call(horizontal_hmc, args)
  expect_equal(as.array(d1), as.array(d2))
})

test_that("horizontal_hmc L=1 and L=3 both run without error", {
  for (l in c(1L, 3L)) {
    draws <- horizontal_hmc(
      tiny_data, n_iter = 10L, n_warmup = 10L, n_chains = 1L,
      L = l, seed = 5L, verbose = FALSE
    )
    expect_true(inherits(draws, "draws_array"),
                info = paste("L =", l))
  }
})

test_that("horizontal_hmc and riemannian_mcmc produce different draws (connection has effect)", {
  set.seed(7L)
  d_horiz <- horizontal_hmc(
    tiny_data, n_iter = 50L, n_warmup = 50L, n_chains = 1L,
    L = 1L, seed = 7L, verbose = FALSE
  )
  d_riem <- riemannian_mcmc(
    tiny_data, n_iter = 50L, n_warmup = 50L, n_chains = 1L,
    L = 1L, method = "diagonal", seed = 7L, verbose = FALSE
  )
  # They use the same seed and same metric but different position updates;
  # draws should differ (connection correction has a real effect).
  expect_false(
    isTRUE(all.equal(as.array(d_horiz), as.array(d_riem))),
    info = "horizontal correction should produce different trajectories"
  )
})
