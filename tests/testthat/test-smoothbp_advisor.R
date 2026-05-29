library(posterior)

# ── Internal math helpers ─────────────────────────────────────────────────────

test_that(".sbp_dmu_domega gives zero for delta=0 and non-first breakpoint", {
  tau      <- c(1, 2, 3, 4)
  omega_ki <- c(2, 2, 2, 2)
  rho_ki   <- c(1, 1, 1, 1)
  delta_ki <- c(0, 0, 0, 0)
  b1_i     <- c(0.5, 0.5, 0.5, 0.5)

  # k > 1: b1 not included
  dmu <- .sbp_dmu_domega(tau, omega_ki, rho_ki, delta_ki, b1_i, is_first = FALSE)
  expect_true(all(dmu == 0))
})

test_that(".sbp_dmu_domega includes b1 only for first breakpoint", {
  tau <- c(0, 1, 2, 3)
  omega_ki <- c(1, 1, 1, 1)
  rho_ki   <- c(1, 1, 1, 1)
  delta_ki <- c(0, 0, 0, 0)   # zero out delta to isolate b1 effect
  b1_i     <- c(1, 1, 1, 1)

  dmu_first    <- .sbp_dmu_domega(tau, omega_ki, rho_ki, delta_ki, b1_i, is_first = TRUE)
  dmu_nonfirst <- .sbp_dmu_domega(tau, omega_ki, rho_ki, delta_ki, b1_i, is_first = FALSE)

  expect_true(all(dmu_first == -b1_i))     # -b1 only
  expect_true(all(dmu_nonfirst == 0))      # zero
})

test_that(".sbp_dmu_domega is finite for extreme inputs", {
  set.seed(1L)
  tau      <- rnorm(20)
  omega_ki <- rnorm(20)
  rho_ki   <- abs(rnorm(20)) + 0.1
  delta_ki <- rnorm(20)
  b1_i     <- rnorm(20)

  dmu <- .sbp_dmu_domega(tau, omega_ki, rho_ki, delta_ki, b1_i, is_first = TRUE)
  expect_true(all(is.finite(dmu)))
})

test_that(".sbp_lp reconstructs correctly from design matrix", {
  X <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 3L,
              dimnames = list(NULL, c("(Intercept)", "re_g_A")))
  row   <- c("pre_(Intercept)" = 2.0, "pre_re_g_A" = 0.5)
  lp    <- .sbp_lp(row, "pre_", X, c("(Intercept)", "re_g_A"))
  # Row 1: 2 + 0 = 2; Row 2: 2 + 0 = 2; Row 3: 2 + 0.5 = 2.5
  expect_equal(lp, c(2.0, 2.0, 2.5))
})

test_that(".sbp_lp applies spike-and-slab zeroing", {
  X <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 3L,
              dimnames = list(NULL, c("(Intercept)", "re_g_A")))
  row <- c("pre_(Intercept)" = 2.0, "pre_re_g_A" = 0.5,
           "g_(Intercept)" = 0,  "g_re_g_A" = 1)   # gamma: 0 for intercept, 1 for RE
  lp <- .sbp_lp(row, "pre_", X, c("(Intercept)", "re_g_A"),
                gamma_prefix = "g_")
  # Intercept zeroed (gamma=0): Row 1=0+0=0; Row 2=0+0=0; Row 3=0+0.5=0.5
  expect_equal(lp, c(0.0, 0.0, 0.5))
})

# ── smoothbp_advisor: synthetic smoothbp_fit ──────────────────────────────────

# Build a minimal smoothbp_fit-like object for testing the advisor logic.
# We don't have a compiled smoothbp here, so we construct the object structure
# that smoothbp_advisor() expects.

.make_fake_fit <- function(n_obs = 30L, n_groups = 3L,
                            sigma_true = 0.5, sigma_re_true = 1.0,
                            n_draws = 100L, seed = 42L) {
  set.seed(seed)
  J       <- n_groups
  tau     <- seq(0, 5, length.out = n_obs)
  group   <- rep(seq_len(J), length.out = n_obs)
  omega_true  <- 2.5
  u_om_true   <- rnorm(J, 0, sigma_re_true)
  rho_true    <- 1.0
  delta_true  <- 1.0
  b1_true     <- 0.5

  # Group-specific changepoints
  omega_i <- omega_true + u_om_true[group]
  d_i     <- tau - omega_i
  s_i     <- plogis(d_i * rho_true)
  mu_i    <- b1_true * (tau - omega_true) + delta_true * d_i * s_i
  y       <- mu_i + rnorm(n_obs, 0, sigma_true)

  # Design matrices
  X_b1     <- matrix(1, n_obs, 1L, dimnames = list(NULL, "(Intercept)"))
  X_delta  <- matrix(1, n_obs, 1L, dimnames = list(NULL, "(Intercept)"))
  X_rho    <- matrix(1, n_obs, 1L, dimnames = list(NULL, "(Intercept)"))

  # X_om: intercept column + one RE column per group
  lvls  <- paste0("G", seq_len(J))
  X_re  <- matrix(0L, n_obs, J, dimnames = list(NULL, paste0("re_group_", lvls)))
  for (g in seq_len(J)) X_re[group == g, g] <- 1L
  X_om  <- cbind(matrix(1, n_obs, 1L, dimnames = list(NULL, "(Intercept)")), X_re)
  attr(X_om, "re_mask") <- c(0L, rep(1L, J))

  # Fake draws_array: n_draws × 1 chain × n_params
  col_b1    <- "(Intercept)"
  col_del   <- "(Intercept)"
  col_om    <- c("(Intercept)", paste0("re_group_", lvls))
  col_rho   <- "(Intercept)"

  make_pnames <- function(prefix, cols) paste0(prefix, cols)

  pnames <- c(
    make_pnames("b1_",    col_b1),
    make_pnames("delta1_", col_del),
    make_pnames("omega1_", col_om),
    make_pnames("rho1_",  col_rho),
    "sigma", "sigma_re_omega1"
  )

  # Simulate posterior draws centred on truth
  n_p   <- length(pnames)
  draws_arr <- array(NA_real_, dim = c(n_draws, 1L, n_p),
                     dimnames = list(NULL, NULL, pnames))
  draws_arr[, 1L, "b1_(Intercept)"]     <- rnorm(n_draws, b1_true,    0.05)
  draws_arr[, 1L, "delta1_(Intercept)"] <- rnorm(n_draws, delta_true, 0.05)
  draws_arr[, 1L, "omega1_(Intercept)"] <- rnorm(n_draws, omega_true, 0.05)
  for (g in seq_len(J)) {
    vn <- paste0("omega1_re_group_G", g)
    draws_arr[, 1L, vn] <- rnorm(n_draws, u_om_true[g], 0.05)
  }
  draws_arr[, 1L, "rho1_(Intercept)"]  <- abs(rnorm(n_draws, rho_true, 0.05))
  draws_arr[, 1L, "sigma"]             <- abs(rnorm(n_draws, sigma_true, 0.02))
  draws_arr[, 1L, "sigma_re_omega1"]   <- abs(rnorm(n_draws, sigma_re_true, 0.05))

  fit <- structure(
    list(
      draws  = posterior::as_draws_array(draws_arr),
      time   = "tau",
      data   = data.frame(tau = tau, group = group, y = y),
      dm     = list(
        X_b1           = X_b1,
        X_deltas       = list(X_delta),
        X_om           = list(X_om),
        X_rho          = list(X_rho),
        col_names_b1      = col_b1,
        col_names_deltas  = list(col_del),
        col_names_om      = list(col_om),
        col_names_rho     = list(col_rho)
      )
    ),
    class = "smoothbp_fit"
  )
  fit
}

test_that("smoothbp_advisor returns fibr_smoothbp_advice for valid fit", {
  fit    <- .make_fake_fit()
  advice <- smoothbp_advisor(fit, n_draws = 30L)
  expect_s3_class(advice, "fibr_smoothbp_advice")
  expect_length(advice$breakpoints, 1L)
})

test_that("smoothbp_advisor has one result per RE group", {
  J   <- 4L
  fit <- .make_fake_fit(n_groups = J)
  adv <- smoothbp_advisor(fit, n_draws = 30L)
  bp  <- adv$breakpoints[[1L]]
  expect_equal(length(bp$re_vars), J)
  expect_equal(length(bp$prior_frac_mean), J)
})

test_that("prior_frac values are in (0, 1)", {
  fit <- .make_fake_fit()
  adv <- smoothbp_advisor(fit, n_draws = 50L)
  pf  <- adv$breakpoints[[1L]]$prior_frac_mean
  expect_true(all(pf > 0 & pf < 1))
})

test_that("prior-dominated model gives higher prior_frac", {
  # prior_frac = G_prior / (G_prior + G_lik) where G_prior = 1/sigma_re^2
  # Small sigma_re → tight prior → HIGH G_prior → HIGH prior_frac (prior-dominated)
  # Large sigma_re → diffuse prior → LOW  G_prior → LOW  prior_frac (data-dominated)
  fit_narrow <- .make_fake_fit(sigma_re_true = 0.1)   # tight prior
  fit_wide   <- .make_fake_fit(sigma_re_true = 5.0)   # diffuse prior

  pf_narrow <- mean(smoothbp_advisor(fit_narrow, n_draws = 50L)$breakpoints[[1L]]$prior_frac_mean)
  pf_wide   <- mean(smoothbp_advisor(fit_wide,   n_draws = 50L)$breakpoints[[1L]]$prior_frac_mean)

  expect_gt(pf_narrow, pf_wide)   # tight prior → more prior-dominated → higher prior_frac
})

test_that("recommendation is one of the three expected strings", {
  fit <- .make_fake_fit()
  adv <- smoothbp_advisor(fit, n_draws = 30L)
  rec <- adv$breakpoints[[1L]]$recommendation
  expect_true(all(rec %in% c("non-centred", "centred (OK)", "borderline")))
})

test_that("print.fibr_smoothbp_advice runs without error", {
  fit <- .make_fake_fit()
  adv <- smoothbp_advisor(fit, n_draws = 20L)
  expect_output(print(adv), "prior_frac")
})

test_that("smoothbp_advisor returns NULL for fit with no RE on omega", {
  fit <- .make_fake_fit()
  # Remove RE columns from X_om to simulate no-RE fit
  X_plain <- matrix(1, 30L, 1L, dimnames = list(NULL, "(Intercept)"))
  attr(X_plain, "re_mask") <- 0L
  fit$dm$X_om           <- list(X_plain)
  fit$dm$col_names_om   <- list("(Intercept)")
  # Remove RE draws
  fit$draws <- posterior::subset_draws(fit$draws,
                 variable = c("b1_(Intercept)", "delta1_(Intercept)",
                              "omega1_(Intercept)", "rho1_(Intercept)",
                              "sigma", "sigma_re_omega1"))
  expect_message(result <- smoothbp_advisor(fit, n_draws = 10L), "no random effects")
  expect_null(result)
})
