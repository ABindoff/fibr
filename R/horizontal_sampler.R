#' Connection-corrected Gibbs sampler for the centred GLMM
#'
#' @description
#' A block Metropolis-within-Gibbs sampler for the centred two-level logistic
#' GLMM that incorporates the **horizontal transport correction**: whenever the
#' base parameters \eqn{(\mu, \sigma)} are proposed to move, the fiber
#' parameters \eqn{\alpha} are pre-displaced along the Fisher metric connection
#' before the accept/reject step.
#'
#' Without the correction the base and fiber proposals are misaligned: moving
#' \eqn{(\mu, \sigma)} while holding \eqn{\alpha} fixed moves the joint state
#' "off the horizontal slice," leading to low acceptance and slow mixing in
#' \eqn{\alpha}.  The correction pre-aligns \eqn{\alpha} with the new base
#' position, dramatically increasing the effective step size in the base
#' direction.
#'
#' Update schedule (per iteration):
#' 1. **Base block**: propose \eqn{(\mu, \log\sigma)} jointly via random walk;
#'    if `use_correction = TRUE`, displace \eqn{\alpha} by
#'    \eqn{A(\theta, \alpha)\,\Delta\theta}; accept/reject jointly.
#' 2. **Fiber block**: propose each \eqn{\alpha_j} independently via random
#'    walk; accept/reject individually.
#' 3. **Fixed-effects block**: propose \eqn{\beta} jointly; accept/reject.
#'
#' @param stan_data   Named list (X, y, group) as passed to Stan.
#' @param n_iter      Number of post-warmup iterations (default 2000).
#' @param n_warmup    Number of warmup iterations (default 1000).
#' @param n_chains    Number of independent chains (default 4).
#' @param init        Optional named list of starting values
#'   (`mu`, `sigma`, `alpha`, `beta`).  Defaults to prior draws.
#' @param use_correction Logical; apply horizontal transport correction when
#'   proposing base moves (default `TRUE`).  Set `FALSE` for the uncorrected
#'   baseline.
#' @param step_base   Initial random-walk SD for the base block
#'   \eqn{(\mu, \log\sigma)} (default 0.1).
#' @param step_alpha  Initial random-walk SD for each \eqn{\alpha_j}
#'   (default 0.3).
#' @param step_beta   Initial random-walk SD for the \eqn{\beta} block
#'   (default 0.15).
#' @param target_rate Target Metropolis acceptance rate for adaptive tuning
#'   (default 0.30).
#' @param seed        Random seed for reproducibility.
#' @param verbose     Print per-chain progress (default `TRUE`).
#'
#' @return A `draws_array` (from the `posterior` package) with the same
#'   variable naming as the centred Stan model: `mu`, `sigma`, `alpha[1]`,
#'   ..., `alpha[J]`, `beta[1]`, `beta[2]`.
#'
#' @export
horizontal_mcmc <- function(stan_data,
                             n_iter       = 2000L,
                             n_warmup     = 1000L,
                             n_chains     = 4L,
                             init         = NULL,
                             use_correction = TRUE,
                             step_base    = 0.10,
                             step_alpha   = 0.30,
                             step_beta    = 0.15,
                             target_rate  = 0.30,
                             seed         = NULL,
                             verbose      = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  J    <- max(stan_data$group)
  N    <- length(stan_data$y)
  n_total <- n_warmup + n_iter
  P    <- 2L + J + 2L   # mu, log_sigma, alpha[1:J], beta[1:2]

  # Parameter names matching Stan output
  par_names <- c("mu", "sigma",
                 paste0("alpha[", seq_len(J), "]"),
                 paste0("beta[",  1:2,         "]"))

  # Storage: n_iter rows x P cols per chain
  all_chains <- array(NA_real_,
                      dim      = c(n_iter, n_chains, P),
                      dimnames = list(NULL, NULL, par_names))

  for (chain_id in seq_len(n_chains)) {
    if (verbose) cat(sprintf("Chain %d/%d ...\n", chain_id, n_chains))

    # ── Initialise ───────────────────────────────────────────────────────────
    if (!is.null(init)) {
      mu_c    <- init$mu
      sigma_c <- init$sigma
      alpha_c <- init$alpha
      beta_c  <- init$beta
    } else {
      mu_c    <- rnorm(1, 0, 1)
      sigma_c <- rexp(1, 2)
      alpha_c <- rnorm(J, mu_c, sigma_c)
      beta_c  <- rnorm(2, 0, 0.5)
    }
    log_sig_c <- log(sigma_c)

    lp_c <- .glmm_log_post(mu_c, log_sig_c, alpha_c, beta_c, stan_data)

    # Per-chain adaptive step sizes
    s_base  <- step_base
    s_alpha <- step_alpha
    s_beta  <- step_beta

    # Acceptance counters for tuning
    acc_base  <- 0L;  n_base  <- 0L
    acc_alpha <- 0L;  n_alpha <- 0L
    acc_beta  <- 0L;  n_beta  <- 0L

    # ── Main loop ─────────────────────────────────────────────────────────────
    for (iter in seq_len(n_total)) {

      # 1. Base block: (mu, log_sigma) + optional horizontal correction to alpha
      prop_mu  <- mu_c    + rnorm(1, 0, s_base)
      prop_lsig <- log_sig_c + rnorm(1, 0, s_base)
      prop_sig  <- exp(prop_lsig)

      # Horizontal correction: pre-displace alpha along the connection.
      # The correction is clamped to ±max_corr per component to prevent
      # runaway displacement when the base step is large relative to sigma.
      if (use_correction) {
        delta_mu  <- prop_mu  - mu_c
        delta_sig <- prop_sig - sigma_c
        G_FF_c <- .glmm_G_FF(sigma_c, alpha_c, stan_data$X, stan_data$group, beta_c)
        G_BF_c <- .glmm_G_BF(sigma_c, mu_c, alpha_c)
        A_c    <- .glmm_connection(G_FF_c, G_BF_c)      # J x 2
        delta_alpha <- as.vector(A_c %*% c(delta_mu, delta_sig))
        max_corr    <- 3 * sigma_c                       # clamp at 3 prior SDs
        delta_alpha <- pmax(pmin(delta_alpha, max_corr), -max_corr)
        alpha_prop  <- alpha_c + delta_alpha
      } else {
        alpha_prop <- alpha_c
      }

      lp_prop <- tryCatch(
        .glmm_log_post(prop_mu, prop_lsig, alpha_prop, beta_c, stan_data),
        error = function(e) -Inf
      )
      if (!is.finite(lp_prop)) lp_prop <- -Inf

      log_r_base <- lp_prop - lp_c
      if (is.finite(log_r_base) && log(runif(1L)) < log_r_base) {
        mu_c      <- prop_mu
        log_sig_c <- prop_lsig
        sigma_c   <- prop_sig
        alpha_c   <- alpha_prop
        lp_c      <- lp_prop
        acc_base  <- acc_base + 1L
      }
      n_base <- n_base + 1L

      # 2. Fiber block: update each alpha_j independently
      for (j in seq_len(J)) {
        a_prop_j <- alpha_c[j] + rnorm(1L, 0, s_alpha)
        alpha_try <- alpha_c
        alpha_try[j] <- a_prop_j
        lp_try <- tryCatch(
          .glmm_log_post(mu_c, log_sig_c, alpha_try, beta_c, stan_data),
          error = function(e) -Inf
        )
        if (!is.finite(lp_try)) lp_try <- -Inf
        log_r_alpha <- lp_try - lp_c
        if (is.finite(log_r_alpha) && log(runif(1L)) < log_r_alpha) {
          alpha_c[j] <- a_prop_j
          lp_c       <- lp_try
          acc_alpha  <- acc_alpha + 1L
        }
        n_alpha <- n_alpha + 1L
      }

      # 3. Fixed-effects block: (beta[1], beta[2]) jointly
      beta_prop <- beta_c + rnorm(2L, 0, s_beta)
      lp_beta <- tryCatch(
        .glmm_log_post(mu_c, log_sig_c, alpha_c, beta_prop, stan_data),
        error = function(e) -Inf
      )
      if (!is.finite(lp_beta)) lp_beta <- -Inf
      log_r_beta <- lp_beta - lp_c
      if (is.finite(log_r_beta) && log(runif(1L)) < log_r_beta) {
        beta_c   <- beta_prop
        lp_c     <- lp_beta
        acc_beta <- acc_beta + 1L
      }
      n_beta <- n_beta + 1L

      # Adaptive step size tuning during warmup
      if (iter <= n_warmup && iter %% 100L == 0L) {
        rate_base  <- acc_base  / n_base
        rate_alpha <- acc_alpha / n_alpha
        rate_beta  <- acc_beta  / n_beta

        s_base  <- .tune_step(s_base,  rate_base,  target_rate)
        s_alpha <- .tune_step(s_alpha, rate_alpha, target_rate)
        s_beta  <- .tune_step(s_beta,  rate_beta,  target_rate)

        acc_base  <- 0L;  n_base  <- 0L
        acc_alpha <- 0L;  n_alpha <- 0L
        acc_beta  <- 0L;  n_beta  <- 0L
      }

      # Store post-warmup draws
      if (iter > n_warmup) {
        t <- iter - n_warmup
        all_chains[t, chain_id, ]  <- c(mu_c, exp(log_sig_c), alpha_c, beta_c)
      }
    }

    if (verbose) {
      cat(sprintf(
        "  step_base=%.3f  step_alpha=%.3f  step_beta=%.3f\n",
        s_base, s_alpha, s_beta
      ))
      cat(sprintf(
        "  acceptance: base=%.2f  alpha=%.2f  beta=%.2f\n",
        acc_base / max(n_base, 1L),
        acc_alpha / max(n_alpha, 1L),
        acc_beta / max(n_beta, 1L)
      ))
    }
  }

  # Return as posterior draws_array
  posterior::as_draws_array(all_chains)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Multiplicative step-size adjustment toward target acceptance rate.
.tune_step <- function(step, rate, target, factor = 1.3, max_step = 0.5) {
  new_step <- if (rate > target) step * factor else step / factor
  min(new_step, max_step)
}
