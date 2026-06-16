# Tests for prior_fraction().
# The core tests (arithmetic + family information) run everywhere and are fast.
# The brms-adapter tests fit small models, so they are skipped on CRAN and when
# brms is unavailable. Numeric checks recompute pi independently using the SAME
# posterior draws (ndraws = all), so they are deterministic and tight; a failure
# is a real extraction/arithmetic bug, not MCMC noise. The fully-annotated
# standalone version lives in data-raw/test_prior_fraction.R.

test_that("prior_fraction.default: arithmetic and edge cases", {
  pf <- prior_fraction(c(4, 1, 0.25), lik_information = c(0, 1, 0.75))
  expect_equal(pf$pi, c(1, 0.5, 0.25))
  expect_equal(prior_fraction(0, lik_information = 0)$pi, 1)
  expect_error(prior_fraction(c(1, 2), lik_information = 1))
})

test_that(".glm_information matches the GLM working weights", {
  eta <- c(-1, 0, 2); p <- plogis(eta); m <- exp(eta)
  expect_equal(fibr:::.glm_information("gaussian", eta, dispersion = 4), rep(0.25, 3))
  expect_equal(fibr:::.glm_information("bernoulli", eta), p * (1 - p))
  expect_equal(fibr:::.glm_information("binomial", eta, trials = 5), 5 * p * (1 - p))
  expect_equal(fibr:::.glm_information("poisson", eta), m)
  expect_equal(fibr:::.glm_information("negbinomial", eta, dispersion = 3), m / (1 + m / 3))
  expect_error(fibr:::.glm_information("student", eta))
})

# ---- brms adapter (skipped on CRAN / when brms is unavailable) ---------------

.pf_backend <- function() if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
.pf_fit <- function(formula, data, family) {
  brms::brm(formula, data = data, family = family, chains = 2, iter = 1000,
            warmup = 500, refresh = 0, seed = 1, backend = .pf_backend())
}

test_that("brmsfit gaussian (1|g): pi matches the n_g formula", {
  skip_on_cran(); skip_if_not_installed("brms")
  set.seed(1)
  ng <- rep(c(2L, 5L, 20L, 80L), 3); g <- factor(rep(seq_along(ng), ng))
  re <- rnorm(length(ng)); y <- rnorm(sum(ng), re[as.integer(g)], 1)
  dat <- data.frame(y = y, g = g)
  fit <- .pf_fit(y ~ 1 + (1 | g), dat, gaussian())
  nd  <- nrow(posterior::as_draws_df(fit)); dd <- posterior::as_draws_df(fit)
  pf  <- prior_fraction(fit, ndraws = nd); pfi <- pf[pf$coef == "Intercept", ]
  expect_equal(sort(pfi$n_obs), sort(as.integer(table(dat$g))))
  pp <- 1 / mean(dd$sd_g__Intercept)^2
  pe <- pp / (pp + pfi$n_obs / mean(dd$sigma)^2)
  expect_equal(pfi$pi, pe, tolerance = 1e-6)
  expect_true(all(pfi$pi > 0 & pfi$pi <= 1) && cor(pfi$n_obs, pfi$pi) < 0)
})

test_that("brmsfit bernoulli (1|g): pi matches sum p(1-p)", {
  skip_on_cran(); skip_if_not_installed("brms")
  set.seed(2)
  ng <- rep(c(3L, 8L, 25L, 90L), 3); g <- factor(rep(seq_along(ng), ng))
  re <- rnorm(length(ng)); y <- rbinom(sum(ng), 1, plogis(re[as.integer(g)]))
  dat <- data.frame(y = y, g = g)
  fit <- .pf_fit(y ~ 1 + (1 | g), dat, brms::bernoulli())
  nd  <- nrow(posterior::as_draws_df(fit)); dd <- posterior::as_draws_df(fit)
  pf  <- prior_fraction(fit, ndraws = nd); pfi <- pf[pf$coef == "Intercept", ]
  eta <- colMeans(brms::posterior_linpred(fit, ndraws = nd)); pv <- plogis(eta)
  li  <- tapply(pv * (1 - pv), dat$g, sum)
  pp  <- 1 / mean(dd$sd_g__Intercept)^2
  pe  <- pp / (pp + li)
  pg  <- pfi$pi[match(names(pe), pfi$level)]
  expect_false(anyNA(pg))                       # [brms-internal: ranef levels]
  expect_equal(unname(pg), unname(as.numeric(pe)), tolerance = 1e-4)
})

test_that("brmsfit poisson (1|g): pi matches sum mu", {
  skip_on_cran(); skip_if_not_installed("brms")
  set.seed(3)
  ng <- rep(c(4L, 15L, 60L), 3); g <- factor(rep(seq_along(ng), ng))
  re <- rnorm(length(ng), 0, 0.5); y <- rpois(sum(ng), exp(0.2 + re[as.integer(g)]))
  dat <- data.frame(y = y, g = g)
  fit <- .pf_fit(y ~ 1 + (1 | g), dat, poisson())
  nd  <- nrow(posterior::as_draws_df(fit)); dd <- posterior::as_draws_df(fit)
  pf  <- prior_fraction(fit, ndraws = nd); pfi <- pf[pf$coef == "Intercept", ]
  eta <- colMeans(brms::posterior_linpred(fit, ndraws = nd))
  li  <- tapply(exp(eta), dat$g, sum)
  pp  <- 1 / mean(dd$sd_g__Intercept)^2
  pg  <- pfi$pi[match(names(li), pfi$level)]
  expect_false(anyNA(pg))
  expect_equal(unname(pg), unname(as.numeric(pp / (pp + li))), tolerance = 1e-4)
})

test_that("brmsfit correlated (x|g): slope uses x^2 loadings and warns", {
  skip_on_cran(); skip_if_not_installed("brms")
  set.seed(4)
  ng <- rep(c(6L, 30L), 5); g <- factor(rep(seq_along(ng), ng)); N <- sum(ng)
  x <- rnorm(N); b0 <- rnorm(length(ng)); b1 <- rnorm(length(ng), 0, 0.7)
  y <- rnorm(N, b0[as.integer(g)] + b1[as.integer(g)] * x, 1)
  dat <- data.frame(y = y, g = g, x = x)
  fit <- .pf_fit(y ~ x + (x | g), dat, gaussian())
  expect_message(pf <- prior_fraction(fit), "correlated")
  expect_true(all(c("Intercept", "x") %in% pf$coef))
  dd  <- posterior::as_draws_df(fit); pfx <- pf[pf$coef == "x", ]
  sx2 <- tapply(dat$x^2, dat$g, sum)
  pp  <- 1 / mean(dd$sd_g__x)^2
  pg  <- pfx$pi[match(names(sx2), pfx$level)]
  expect_false(anyNA(pg))                       # [brms-internal: Z loadings]
  expect_equal(unname(pg), unname(as.numeric(pp / (pp + sx2 / mean(dd$sigma)^2))),
               tolerance = 1e-5)
})

test_that("print and plot methods work", {
  skip_on_cran(); skip_if_not_installed("brms")
  set.seed(5)
  g <- factor(rep(1:6, each = 10)); y <- rnorm(60, rnorm(6)[as.integer(g)], 1)
  fit <- .pf_fit(y ~ 1 + (1 | g), data.frame(y = y, g = g), gaussian())
  pf <- prior_fraction(fit)
  expect_output(print(pf))
  expect_s3_class(plot(pf), "ggplot")
})
