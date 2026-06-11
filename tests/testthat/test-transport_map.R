# ── helpers ───────────────────────────────────────────────────────────────────

# Build a loops data frame from paired index ranges
make_loops <- function(n_loops, start_offset = 0L, end_offset = 50L,
                       distance = 0.05) {
  starts <- seq_len(n_loops)
  data.frame(start    = starts + start_offset,
             end      = starts + end_offset,
             distance = distance)
}

# ── estimate_transport_map ────────────────────────────────────────────────────

test_that("output has correct structure", {
  set.seed(1L)
  J <- 4L; N <- 200L; K <- 50L
  fiber <- matrix(rnorm(N * J), N, J)
  loops <- make_loops(K, end_offset = 100L)
  tm    <- estimate_transport_map(fiber, loops, n_bootstrap = 10L)

  expect_equal(dim(tm$H), c(J, J))
  expect_length(tm$eigenvalues, J)
  expect_equal(dim(tm$boot_eigenvalues), c(10L, J))
  expect_true(is.numeric(tm$frobenius_dev) && tm$frobenius_dev >= 0)
  expect_equal(tm$n_loops, K)
})

test_that("full-structure eigenvalues are sorted by decreasing modulus", {
  set.seed(2L)
  J <- 5L; N <- 200L
  fiber <- matrix(rnorm(N * J), N, J)
  loops <- make_loops(60L, end_offset = 80L)
  tm    <- estimate_transport_map(fiber, loops, n_bootstrap = 5L,
                                  structure = "full")
  mods  <- Mod(tm$eigenvalues)
  expect_true(all(diff(mods) <= 1e-10))  # non-increasing
})

test_that("diagonal structure returns a diagonal H with group-aligned factors", {
  set.seed(21L)
  J <- 4L; K <- 300L
  h_true <- c(0.2, 0.5, 0.8, 1.0)   # heterogeneous per-group contraction
  alpha  <- matrix(rnorm(K * J), K, J)
  fiber  <- rbind(alpha,
                  sweep(alpha, 2L, h_true, "*") +
                    matrix(rnorm(K * J, 0, 0.01), K, J))
  loops  <- make_loops(K, start_offset = 0L, end_offset = K)
  tm     <- estimate_transport_map(fiber, loops, n_bootstrap = 20L,
                                   weights = "uniform", structure = "diagonal")

  # H is diagonal
  off_diag <- tm$H; diag(off_diag) <- 0
  expect_equal(max(abs(off_diag)), 0)
  expect_identical(tm$structure, "diagonal")

  # Factors recovered in GROUP order (not modulus-sorted)
  expect_equal(Re(tm$eigenvalues), h_true, tolerance = 0.05)
  expect_true(all(abs(Im(tm$eigenvalues)) < 1e-12))

  # Bootstrap columns stay aligned to groups
  boot_means <- colMeans(Re(tm$boot_eigenvalues), na.rm = TRUE)
  expect_equal(boot_means, h_true, tolerance = 0.1)
})

test_that("diagonal estimator matches full estimator when truth is diagonal", {
  set.seed(22L)
  J <- 3L; K <- 500L
  alpha <- matrix(rnorm(K * J), K, J)
  fiber <- rbind(alpha, 0.6 * alpha + matrix(rnorm(K * J, 0, 0.01), K, J))
  loops <- make_loops(K, start_offset = 0L, end_offset = K)
  tm_d  <- estimate_transport_map(fiber, loops, n_bootstrap = 0L,
                                  weights = "uniform", structure = "diagonal")
  tm_f  <- estimate_transport_map(fiber, loops, n_bootstrap = 0L,
                                  weights = "uniform", structure = "full")
  expect_lt(norm(tm_d$H - tm_f$H, "F"), 0.1)
})

test_that("identity transform gives H near I", {
  # alpha_end = alpha_start + tiny noise => H should be near I
  set.seed(3L)
  J <- 3L; K <- 200L; N <- K * 2L
  alpha <- matrix(rnorm(K * J), K, J)
  fiber <- rbind(alpha, alpha + matrix(rnorm(K * J, 0, 0.01), K, J))
  loops <- make_loops(K, start_offset = 0L, end_offset = K)
  tm    <- estimate_transport_map(fiber, loops, n_bootstrap = 0L, weights = "uniform")
  expect_lt(norm(tm$H - diag(J), "F"), 0.3)  # within 0.3 of identity
})

test_that("doubling transform gives H near 2*I", {
  set.seed(4L)
  J <- 2L; K <- 200L; N <- K * 2L
  alpha <- matrix(rnorm(K * J), K, J)
  fiber <- rbind(alpha, 2 * alpha + matrix(rnorm(K * J, 0, 0.01), K, J))
  loops <- make_loops(K, start_offset = 0L, end_offset = K)
  tm    <- estimate_transport_map(fiber, loops, n_bootstrap = 0L, weights = "uniform")
  expect_lt(norm(tm$H - 2 * diag(J), "F"), 0.5)
})

test_that("uniform and distance weighting give similar results for equal distances", {
  set.seed(5L)
  J <- 3L; N <- 300L
  fiber <- matrix(rnorm(N * J), N, J)
  loops <- make_loops(80L, end_offset = 100L, distance = 0.1)  # all same dist
  tm_u  <- estimate_transport_map(fiber, loops, n_bootstrap = 0L, weights = "uniform")
  tm_d  <- estimate_transport_map(fiber, loops, n_bootstrap = 0L, weights = "distance")
  expect_lt(norm(tm_u$H - tm_d$H, "F"), 0.01)  # should be nearly identical
})
