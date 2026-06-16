library(FNN)

# ── helpers ───────────────────────────────────────────────────────────────────

# Circular chain: lots of near-returns
circle_chain <- function(n = 200L) {
  t <- seq(0, 4 * pi, length.out = n)
  cbind(cos(t), sin(t))
}

# Monotone chain: no loops possible
mono_chain <- function(n = 100L) {
  cbind(seq(0, 10, length.out = n), rep(0, n))
}

# ── detect_loops ──────────────────────────────────────────────────────────────

test_that("detect_loops returns a data frame with correct columns", {
  loops <- detect_loops(circle_chain(), epsilon = 0.3, min_gap = 20L)
  expect_s3_class(loops, "data.frame")
  expect_named(loops, c("start", "end", "distance"))
})

test_that("all loop distances are below epsilon", {
  eps <- 0.3
  loops <- detect_loops(circle_chain(), epsilon = eps, min_gap = 20L)
  expect_true(all(loops$distance < eps))
})

test_that("all loop gaps satisfy min_gap", {
  gap <- 30L
  loops <- detect_loops(circle_chain(), epsilon = 0.4, min_gap = gap)
  expect_true(all(loops$end - loops$start >= gap))
})

test_that("detect_loops preserves epsilon as attribute", {
  eps <- 0.25
  loops <- detect_loops(circle_chain(), epsilon = eps, min_gap = 10L)
  expect_equal(attr(loops, "epsilon"), eps)
})

test_that("auto-epsilon produces a numeric value and detects loops", {
  loops <- detect_loops(circle_chain(), epsilon = NULL, min_gap = 20L)
  eps <- attr(loops, "epsilon")
  expect_true(is.numeric(eps) && is.finite(eps) && eps > 0)
  expect_gt(nrow(loops), 0L)
})

test_that("max_loops caps the result", {
  cap <- 15L
  loops <- detect_loops(circle_chain(300L), epsilon = 0.5, min_gap = 5L,
                        max_loops = cap)
  expect_lte(nrow(loops), cap)
})

test_that("detect_loops errors when no loops found", {
  expect_error(detect_loops(mono_chain(), epsilon = 0.001, min_gap = 5L),
               regexp = "No loops found")
})

test_that("scale=FALSE skips standardisation", {
  b <- circle_chain()
  # With scale=FALSE, epsilon is in raw units (already unit-scale here)
  loops <- detect_loops(b, epsilon = 0.4, min_gap = 10L, scale = FALSE)
  expect_gt(nrow(loops), 0L)
})
