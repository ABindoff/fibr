# =============================================================================
# Test script for fibr::prior_fraction()   (run by a non-sandboxed agent)
#
# HOW TO RUN (from the package root):
#   devtools::document()    # registers prior_fraction() + S3 methods in NAMESPACE
#   devtools::load_all()
#   source("data-raw/test_prior_fraction.R")
#
# Needs: brms + a working backend (cmdstanr preferred, else rstan).
# The script never stops on failure; it prints PASS/FAIL per check and a summary
# at the end. Please paste the full output back.
#
# Design note: the numeric checks compare the function's pi against an
# INDEPENDENT recomputation that uses the SAME posterior draws (ndraws = all),
# so they are deterministic and the tolerances are tight (1e-6). A failure means
# a genuine extraction/arithmetic bug, not sampling noise. The two brms-internal
# things most likely to differ across versions are (a) the level labels in
# attr(fit$ranef, "levels") and (b) the standata names J_<id> / Z_<id>_<cn>;
# tests that depend on those are marked [brms-internal] so you can localise.
# =============================================================================

suppressMessages({
  ok_brms <- requireNamespace("brms", quietly = TRUE)
  backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
})
options(mc.cores = 2)
set.seed(1)

# ---- tiny test harness -------------------------------------------------------
.results <- list()
check <- function(name, pass, info = "") {
  pass <- isTRUE(pass)
  .results[[length(.results) + 1L]] <<- list(name = name, pass = pass, info = info)
  cat(sprintf("[%s] %s%s\n", if (pass) "PASS" else "FAIL", name,
              if (nzchar(info)) paste0("   <-- ", info) else ""))
}
run <- function(name, expr) {
  tryCatch(force(expr),
           error = function(e) check(name, FALSE, paste("ERROR:", conditionMessage(e))))
  invisible(NULL)
}
fit_brm <- function(formula, data, family) {
  brms::brm(formula, data = data, family = family,
            chains = 2, iter = 1000, warmup = 500, refresh = 0,
            seed = 1, backend = backend)
}
# capture messages emitted by an expression
with_messages <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr,
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  list(value = val, messages = msgs)
}

# =============================================================================
# 1. Core arithmetic and family information  (no brms needed)
# =============================================================================
cat("\n--- 1. core: prior_fraction.default + .glm_information ---\n")

run("default: pi = pp/(pp+li)", {
  pf <- prior_fraction(c(4, 1, 0.25), lik_information = c(0, 1, 0.75))
  expect <- c(1, 0.5, 0.25)          # 4/4, 1/2, 0.25/1
  check("default: pi = pp/(pp+li)", max(abs(pf$pi - expect)) < 1e-12)
})
run("default: zero total precision -> pi = 1", {
  pf <- prior_fraction(0, lik_information = 0)
  check("default: zero total precision -> pi = 1", isTRUE(pf$pi == 1))
})
run("default: length mismatch errors", {
  e <- tryCatch({ prior_fraction(c(1, 2), lik_information = 1); FALSE },
                error = function(e) TRUE)
  check("default: length mismatch errors", isTRUE(e))
})

gi <- fibr:::.glm_information
run(".glm_information families", {
  eta <- c(-1, 0, 2)
  p   <- plogis(eta)
  check(".glm_information gaussian = 1/disp",
        max(abs(gi("gaussian", eta, dispersion = 4) - 0.25)) < 1e-12)
  check(".glm_information bernoulli = p(1-p)",
        max(abs(gi("bernoulli", eta) - p * (1 - p))) < 1e-12)
  check(".glm_information binomial = n p(1-p)",
        max(abs(gi("binomial", eta, trials = 5) - 5 * p * (1 - p))) < 1e-12)
  check(".glm_information poisson = exp(eta)",
        max(abs(gi("poisson", eta) - exp(eta))) < 1e-12)
  m <- exp(eta)
  check(".glm_information negbinomial = m/(1+m/shape)",
        max(abs(gi("negbinomial", eta, dispersion = 3) - m / (1 + m / 3))) < 1e-12)
  check(".glm_information unsupported family errors",
        isTRUE(tryCatch({ gi("student", eta); FALSE }, error = function(e) TRUE)))
})

if (!ok_brms) {
  cat("\n[skip] brms not installed; only core tests ran.\n")
} else {

# =============================================================================
# 2. brms gaussian random intercept:  y ~ 1 + (1 | g)
#    info is constant (1/sigma_y^2), so pi depends only on n_g -> tightest anchor
# =============================================================================
cat("\n--- 2. brms gaussian (1|g) ---\n")
run("gaussian (1|g)", {
  G  <- 12
  ng <- rep(c(2L, 5L, 20L, 80L), 3)            # wide range of group sizes
  g  <- factor(rep(seq_len(G), ng))
  re <- rnorm(G, 0, 1.0)
  y  <- rnorm(sum(ng), re[as.integer(g)], 1.0)
  dat <- data.frame(y = y, g = g)
  fit <- fit_brm(y ~ 1 + (1 | g), dat, gaussian())
  nd  <- nrow(posterior::as_draws_df(fit))     # all post-warmup draws (deterministic)
  dd  <- posterior::as_draws_df(fit)

  pf  <- prior_fraction(fit, ndraws = nd)
  pfi <- pf[pf$coef == "Intercept", ]

  check("gaussian (1|g): n_obs matches data",
        identical(sort(pfi$n_obs), sort(as.integer(table(dat$g)))))

  sd_re <- mean(dd$sd_g__Intercept); sy <- mean(dd$sigma)
  pi_exp <- (1 / sd_re^2) / (1 / sd_re^2 + pfi$n_obs / sy^2)
  check("gaussian (1|g): pi = (1/sd_re^2)/(1/sd_re^2 + n/sigma_y^2)",
        max(abs(pfi$pi - pi_exp)) < 1e-6)
  check("gaussian (1|g): pi in (0,1] and decreasing in n_obs",
        all(pfi$pi > 0 & pfi$pi <= 1) &&
          cor(pfi$n_obs, pfi$pi) < 0)
})

# =============================================================================
# 3. brms bernoulli random intercept (the paper's GLMM case)
#    pi depends on sum p(1-p) per group -> recompute with the SAME draws
# =============================================================================
cat("\n--- 3. brms bernoulli (1|g) ---\n")
run("bernoulli (1|g)", {
  G  <- 12
  ng <- rep(c(3L, 8L, 25L, 90L), 3)
  g  <- factor(rep(seq_len(G), ng))
  re <- rnorm(G, 0, 1.0)
  pr <- plogis(re[as.integer(g)])
  y  <- rbinom(sum(ng), 1, pr)
  dat <- data.frame(y = y, g = g)
  fit <- fit_brm(y ~ 1 + (1 | g), dat, brms::bernoulli())
  nd  <- nrow(posterior::as_draws_df(fit))
  dd  <- posterior::as_draws_df(fit)

  pf  <- prior_fraction(fit, ndraws = nd)
  pfi <- pf[pf$coef == "Intercept", ]

  # independent: same eta the function uses, then sum p(1-p) per group
  eta  <- colMeans(brms::posterior_linpred(fit, ndraws = nd))   # data-row order
  pvec <- plogis(eta); info <- pvec * (1 - pvec)
  li   <- tapply(info, dat$g, sum)                              # named by level
  sd_re <- mean(dd$sd_g__Intercept)
  pi_exp <- (1 / sd_re^2) / (1 / sd_re^2 + li)                  # named by level
  pi_got <- pfi$pi[match(names(pi_exp), pfi$level)]             # [brms-internal: levels]

  check("bernoulli (1|g): levels align (attr(ranef,'levels'))",
        all(!is.na(pi_got)))
  check("bernoulli (1|g): pi = (1/sd_re^2)/(1/sd_re^2 + sum p(1-p))",
        all(!is.na(pi_got)) && max(abs(pi_got - pi_exp)) < 1e-4)
})

# =============================================================================
# 4. brms poisson random intercept  (info = mu = exp(eta))
# =============================================================================
cat("\n--- 4. brms poisson (1|g) ---\n")
run("poisson (1|g)", {
  G  <- 9
  ng <- rep(c(4L, 15L, 60L), 3)
  g  <- factor(rep(seq_len(G), ng))
  re <- rnorm(G, 0, 0.5)
  mu <- exp(0.2 + re[as.integer(g)])
  y  <- rpois(sum(ng), mu)
  dat <- data.frame(y = y, g = g)
  fit <- fit_brm(y ~ 1 + (1 | g), dat, poisson())
  nd  <- nrow(posterior::as_draws_df(fit))
  dd  <- posterior::as_draws_df(fit)

  pf  <- prior_fraction(fit, ndraws = nd)
  pfi <- pf[pf$coef == "Intercept", ]
  eta <- colMeans(brms::posterior_linpred(fit, ndraws = nd))
  li  <- tapply(exp(eta), dat$g, sum)
  sd_re <- mean(dd$sd_g__Intercept)
  pi_exp <- (1 / sd_re^2) / (1 / sd_re^2 + li)
  pi_got <- pfi$pi[match(names(pi_exp), pfi$level)]
  check("poisson (1|g): pi = (1/sd_re^2)/(1/sd_re^2 + sum mu)",
        all(!is.na(pi_got)) && max(abs(pi_got - pi_exp)) < 1e-4)
})

# =============================================================================
# 5. correlated random slope:  y ~ x + (x | g)
#    tests (a) the x^2 loadings are used for the slope, (b) the correlated-RE
#    message fires.
# =============================================================================
cat("\n--- 5. brms correlated (x|g) ---\n")
run("correlated (x|g)", {
  G  <- 10
  ng <- rep(c(6L, 30L), 5)
  g  <- factor(rep(seq_len(G), ng))
  N  <- sum(ng)
  x  <- rnorm(N)
  b0 <- rnorm(G, 0, 1.0); b1 <- rnorm(G, 0, 0.7)
  y  <- rnorm(N, b0[as.integer(g)] + b1[as.integer(g)] * x, 1.0)
  dat <- data.frame(y = y, g = g, x = x)

  cap <- with_messages(fit <- fit_brm(y ~ x + (x | g), dat, gaussian()))
  res <- with_messages(pf <- prior_fraction(fit, ndraws = nrow(posterior::as_draws_df(fit))))
  pf  <- res$value
  check("correlated (x|g): correlated-RE message emitted",
        any(grepl("correlated", res$messages, ignore.case = TRUE)))
  check("correlated (x|g): both Intercept and x coefficients present",
        all(c("Intercept", "x") %in% pf$coef))

  dd  <- posterior::as_draws_df(fit)
  pfx <- pf[pf$coef == "x", ]
  sx2 <- tapply(dat$x^2, dat$g, sum)           # slope loading is x, so Z^2 = x^2
  sd_x <- mean(dd$sd_g__x); sy <- mean(dd$sigma)
  pi_exp <- (1 / sd_x^2) / (1 / sd_x^2 + sx2 / sy^2)
  pi_got <- pfx$pi[match(names(pi_exp), pfx$level)]
  check("correlated (x|g): slope pi uses sum(x^2) loadings [brms-internal: Z]",
        all(!is.na(pi_got)) && max(abs(pi_got - pi_exp)) < 1e-5)
})

# =============================================================================
# 6. print / plot smoke tests (reuse a quick gaussian fit)
# =============================================================================
cat("\n--- 6. print / plot ---\n")
run("print/plot", {
  g  <- factor(rep(1:6, each = 10))
  y  <- rnorm(60, rnorm(6)[as.integer(g)], 1)
  fit <- fit_brm(y ~ 1 + (1 | g), data.frame(y = y, g = g), gaussian())
  pf  <- prior_fraction(fit)
  check("print() runs", isTRUE({ capture.output(print(pf)); TRUE }))
  p <- plot(pf)
  check("plot() returns a ggplot", inherits(p, "ggplot"))
})

} # end if(ok_brms)

# =============================================================================
# summary
# =============================================================================
cat("\n================= SUMMARY =================\n")
np <- sum(vapply(.results, `[[`, logical(1), "pass"))
nt <- length(.results)
for (r in .results)
  if (!r$pass) cat(sprintf("  FAIL: %s%s\n", r$name,
                           if (nzchar(r$info)) paste0("   <-- ", r$info) else ""))
cat(sprintf("\n%d / %d checks passed.%s\n", np, nt,
            if (np == nt) "  All good." else "  See FAILs above."))

