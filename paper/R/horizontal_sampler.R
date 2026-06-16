#' Block Metropolis-within-Gibbs sampler for the centred GLMM
#'
#' @description
#' A block MwG sampler for the centred two-level logistic GLMM.  The
#' **base block** proposes a random walk on \eqn{(\mu, \log\sigma)} and
#' simultaneously transports the fiber \eqn{\alpha} to the new base point;
#' the transport strategy is controlled by `transport`.
#'
#' **Conditional transport** (`transport = "conditional"`, default): transport
#' uses the per-group Laplace affine map
#' \eqn{T(\alpha)_j = m_j' + (s_j'/s_j)(\alpha_j - m_j)}, which is exactly
#' invertible, so the acceptance ratio is
#' \eqn{\log r = \ell(\theta', \alpha') - \ell(\theta, \alpha) +
#' \sum_j \log(s_j'/s_j)}.  No clamp is needed; the map is well-behaved as
#' \eqn{\sigma \to 0}.
#'
#' **Fisher legacy** (`transport = "fisher_legacy"`): the old Fisher-metric
#' Euler lift with a 3-sigma clamp and no Jacobian correction.  Reproduces
#' the biased behaviour of the original sampler bit-for-bit.  Keep for SBC
#' failure figures.
#'
#' **None** (`transport = "none"`): plain MwG base block; \eqn{\alpha} is not
#' pre-displaced.
#'
#' Update schedule: (1) base block, (2) fiber block (per-\eqn{j} RW),
#' (3) fixed-effects block (\eqn{\beta} RW).  The beta block changes
#' \eqn{m_j, s_j} but transport is only used inside the base block.
#'
#' @param stan_data   Named list (X, y, group) as passed to Stan.
#' @param n_iter      Number of post-warmup iterations (default 2000).
#' @param n_warmup    Number of warmup iterations (default 1000).
#' @param n_chains    Number of independent chains (default 4).
#' @param init        Optional named list of starting values
#'   (`mu`, `sigma`, `alpha`, `beta`).  Defaults to prior draws.
#' @param transport   Character; transport strategy for the base block.
#'   One of `"conditional"` (default, exact), `"fisher_legacy"` (biased,
#'   for comparison), or `"none"` (plain MwG).
#' @param use_correction Deprecated.  If non-`NULL`, overrides `transport`:
#'   `TRUE` maps to `"fisher_legacy"`, `FALSE` to `"none"`.
#' @param step_base   Initial random-walk SD for \eqn{(\mu, \log\sigma)}
#'   (default 0.1).
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
                             n_iter         = 2000L,
                             n_warmup       = 1000L,
                             n_chains       = 4L,
                             init           = NULL,
                             transport      = c("conditional", "fisher_legacy", "none"),
                             use_correction = NULL,
                             step_base      = 0.10,
                             step_alpha     = 0.30,
                             step_beta      = 0.15,
                             target_rate    = 0.30,
                             seed           = NULL,
                             verbose        = TRUE) {

  transport <- match.arg(transport)
  if (!is.null(use_correction)) {
    transport <- if (isTRUE(use_correction)) "fisher_legacy" else "none"
  }

  if (!is.null(seed)) set.seed(seed)

  J       <- max(stan_data$group)
  n_total <- n_warmup + n_iter
  P       <- 2L + J + 2L   # mu, log_sigma, alpha[1:J], beta[1:2]

  par_names <- c("mu", "sigma",
                 paste0("alpha[", seq_len(J), "]"),
                 paste0("beta[",  1:2,         "]"))

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
      mu_c    <- rnorm(1L, 0, 1)
      sigma_c <- rexp(1L, 2L)
      alpha_c <- rnorm(J, mu_c, sigma_c)
      beta_c  <- rnorm(2L, 0, 0.5)
    }
    log_sig_c <- log(sigma_c)
    lp_c      <- .glmm_log_post(mu_c, log_sig_c, alpha_c, beta_c, stan_data)

    s_base  <- step_base
    s_alpha <- step_alpha
    s_beta  <- step_beta

    acc_base   <- 0L;  n_base   <- 0L
    acc_alpha  <- 0L;  n_alpha  <- 0L
    acc_beta   <- 0L;  n_beta   <- 0L
    n_lap_fail <- 0L             # Newton failures in conditional branch

    # ── Main loop ─────────────────────────────────────────────────────────────
    for (iter in seq_len(n_total)) {

      # 1. Base block ──────────────────────────────────────────────────────────
      prop_mu   <- mu_c      + rnorm(1L, 0, s_base)
      prop_lsig <- log_sig_c + rnorm(1L, 0, s_base)
      prop_sig  <- exp(prop_lsig)

      if (transport == "conditional") {

        # Laplace at current point (warm-start Newton from current alpha)
        lap_c <- tryCatch(
          .glmm_laplace(c(mu_c, log_sig_c), beta_c, stan_data, m_init = alpha_c),
          error = function(e) NULL
        )
        # Laplace at proposed point (warm-start from current mode)
        lap_p <- if (!is.null(lap_c))
          tryCatch(
            .glmm_laplace(c(prop_mu, prop_lsig), beta_c, stan_data,
                          m_init = lap_c$m),
            error = function(e) NULL
          )
        else NULL

        if (!is.null(lap_c) && !is.null(lap_p)) {
          tr         <- .glmm_cond_transport(alpha_c, lap_c, lap_p)
          alpha_prop <- tr$alpha_new
          lp_prop    <- tryCatch(
            .glmm_log_post(prop_mu, prop_lsig, alpha_prop, beta_c, stan_data),
            error = function(e) -Inf
          )
          if (!is.finite(lp_prop)) lp_prop <- -Inf
          log_r_base <- lp_prop - lp_c + tr$log_jac
        } else {
          n_lap_fail <- n_lap_fail + 1L
          alpha_prop <- alpha_c
          log_r_base <- -Inf
        }

      } else if (transport == "fisher_legacy") {

        delta_mu  <- prop_mu  - mu_c
        delta_sig <- prop_sig - sigma_c
        G_FF_c    <- .glmm_G_FF(sigma_c, alpha_c, stan_data$X, stan_data$group, beta_c)
        G_BF_c    <- .glmm_G_BF(sigma_c, mu_c, alpha_c)
        A_c       <- .glmm_connection(G_FF_c, G_BF_c)
        d_alpha   <- as.vector(A_c %*% c(delta_mu, delta_sig))
        max_corr  <- 3 * sigma_c
        d_alpha   <- pmax(pmin(d_alpha, max_corr), -max_corr)
        alpha_prop <- alpha_c + d_alpha
        lp_prop    <- tryCatch(
          .glmm_log_post(prop_mu, prop_lsig, alpha_prop, beta_c, stan_data),
          error = function(e) -Inf
        )
        if (!is.finite(lp_prop)) lp_prop <- -Inf
        log_r_base <- lp_prop - lp_c          # no Jacobian: legacy biased behaviour

      } else {   # "none"

        alpha_prop <- alpha_c
        lp_prop    <- tryCatch(
          .glmm_log_post(prop_mu, prop_lsig, alpha_prop, beta_c, stan_data),
          error = function(e) -Inf
        )
        if (!is.finite(lp_prop)) lp_prop <- -Inf
        log_r_base <- lp_prop - lp_c

      }

      if (is.finite(log_r_base) && log(runif(1L)) < log_r_base) {
        mu_c      <- prop_mu
        log_sig_c <- prop_lsig
        sigma_c   <- prop_sig
        alpha_c   <- alpha_prop
        lp_c      <- lp_prop
        acc_base  <- acc_base + 1L
      }
      n_base <- n_base + 1L

      # 2. Fiber block ─────────────────────────────────────────────────────────
      for (j in seq_len(J)) {
        a_prop_j  <- alpha_c[j] + rnorm(1L, 0, s_alpha)
        alpha_try <- alpha_c
        alpha_try[j] <- a_prop_j
        lp_try <- tryCatch(
          .glmm_log_post(mu_c, log_sig_c, alpha_try, beta_c, stan_data),
          error = function(e) -Inf
        )
        if (!is.finite(lp_try)) lp_try <- -Inf
        if (is.finite(lp_try - lp_c) && log(runif(1L)) < lp_try - lp_c) {
          alpha_c[j] <- a_prop_j
          lp_c       <- lp_try
          acc_alpha  <- acc_alpha + 1L
        }
        n_alpha <- n_alpha + 1L
      }

      # 3. Fixed-effects block ─────────────────────────────────────────────────
      beta_prop <- beta_c + rnorm(2L, 0, s_beta)
      lp_beta   <- tryCatch(
        .glmm_log_post(mu_c, log_sig_c, alpha_c, beta_prop, stan_data),
        error = function(e) -Inf
      )
      if (!is.finite(lp_beta)) lp_beta <- -Inf
      if (is.finite(lp_beta - lp_c) && log(runif(1L)) < lp_beta - lp_c) {
        beta_c   <- beta_prop
        lp_c     <- lp_beta
        acc_beta <- acc_beta + 1L
      }
      n_beta <- n_beta + 1L

      # ── Adaptive warmup tuning (every 100 iters) ────────────────────────────
      if (iter <= n_warmup && iter %% 100L == 0L) {
        s_base  <- .tune_step(s_base,  acc_base  / n_base,  target_rate)
        s_alpha <- .tune_step(s_alpha, acc_alpha / n_alpha, target_rate)
        s_beta  <- .tune_step(s_beta,  acc_beta  / n_beta,  target_rate)
        acc_base  <- 0L;  n_base  <- 0L
        acc_alpha <- 0L;  n_alpha <- 0L
        acc_beta  <- 0L;  n_beta  <- 0L
      }

      # ── Store post-warmup draws ──────────────────────────────────────────────
      if (iter > n_warmup) {
        t <- iter - n_warmup
        all_chains[t, chain_id, ] <- c(mu_c, exp(log_sig_c), alpha_c, beta_c)
      }
    }

    if (verbose) {
      cat(sprintf(
        "  transport=%s  step_base=%.3f  step_alpha=%.3f  step_beta=%.3f\n",
        transport, s_base, s_alpha, s_beta
      ))
      cat(sprintf(
        "  acceptance: base=%.2f  alpha=%.2f  beta=%.2f\n",
        acc_base / max(n_base, 1L),
        acc_alpha / max(n_alpha, 1L),
        acc_beta / max(n_beta, 1L)
      ))
      if (n_lap_fail > 0L)
        cat(sprintf("  Newton failures (rejected): %d\n", n_lap_fail))
    }
  }

  posterior::as_draws_array(all_chains)
}

#' Horizontal leapfrog HMC for the centred GLMM
#'
#' @description
#' Two transport strategies, selected by `transport`:
#'
#' **Reparameterised HMC** (`transport = "reparam"`, default): exact HMC in
#' standardised z-coordinates \eqn{z_j = (\alpha_j - m_j)/s_j}, where
#' \eqn{(m_j, s_j)} is the per-group Laplace approximation.  The target
#' density is \eqn{\tilde\ell = \ell(\alpha(z)) + \sum_j \log s_j}.
#' Identity mass matrix; standard leapfrog; gradients via chain rule through
#' the Laplace implicit-function derivatives.  Newton solver warm-started
#' across leapfrog steps.
#'
#' **Fisher legacy** (`transport = "fisher_legacy"`): connection-corrected HMC
#' using a fixed diagonal Fisher metric \eqn{G_d} and horizontal correction
#' \eqn{A(\Delta\mu, \Delta\sigma)} applied at every leapfrog step.  Kept for
#' SBC failure-mode comparison.  Recommended \eqn{L = 1} (MALA regime).
#'
#' @param stan_data Named list with X (N×2), y (N), group (N, 1-indexed).
#' @param n_iter    Post-warmup iterations per chain (default 2000).
#' @param n_warmup  Warmup iterations per chain (default 1000).
#' @param n_chains  Number of independent chains (default 4).
#' @param L         Leapfrog steps per trajectory (default 10 for reparam;
#'   use \eqn{L = 1} for fisher_legacy MALA regime).
#' @param epsilon   Initial step size (default 0.10). Adapted during warmup.
#' @param transport Character; one of `"reparam"` (default) or
#'   `"fisher_legacy"`.
#' @param target_rate Target Metropolis acceptance rate (default 0.65).
#' @param init      Optional named list: \code{mu}, \code{sigma},
#'   \code{alpha}, \code{beta}.  Defaults to prior draws.
#' @param seed      Random seed.
#' @param verbose   Print per-chain progress (default \code{TRUE}).
#'
#' @return A \code{\link[posterior]{draws_array}} with the same variable
#'   naming as the centred Stan model: \code{mu}, \code{sigma},
#'   \code{alpha[1]}, \ldots, \code{beta[2]}.
#'
#' @export
horizontal_hmc <- function(stan_data,
                            n_iter      = 2000L,
                            n_warmup    = 1000L,
                            n_chains    = 4L,
                            L           = 10L,
                            epsilon     = 0.10,
                            transport   = c("reparam", "fisher_legacy"),
                            target_rate = 0.65,
                            init        = NULL,
                            seed        = NULL,
                            verbose     = TRUE) {

  transport <- match.arg(transport)
  if (!is.null(seed)) set.seed(seed)

  J       <- max(stan_data$group)
  P       <- 2L + J + 2L
  n_total <- n_warmup + n_iter

  idx_base  <- 1:2
  idx_fiber <- 3L:(2L + J)
  idx_beta  <- (3L + J):(4L + J)

  par_names <- c("mu", "sigma",
                 paste0("alpha[", seq_len(J), "]"),
                 paste0("beta[",  1:2,         "]"))

  all_chains <- array(NA_real_,
                      dim      = c(n_iter, n_chains, P),
                      dimnames = list(NULL, NULL, par_names))

  for (chain_id in seq_len(n_chains)) {
    if (verbose) cat(sprintf("Chain %d/%d ...\n", chain_id, n_chains))

    # ── Initialise ─────────────────────────────────────────────────────────────
    if (!is.null(init)) {
      mu_c <- init$mu;  ls_c <- log(max(init$sigma, 1e-3))
      al_c <- init$alpha;  be_c <- init$beta
      if (chain_id > 1L) {
        mu_c <- mu_c + rnorm(1L, 0, 0.3);  ls_c <- ls_c + rnorm(1L, 0, 0.3)
        al_c <- al_c + rnorm(J,  0, 0.5);  be_c <- be_c + rnorm(2L, 0, 0.2)
      }
    } else {
      mu_c <- rnorm(1L, 0, 1);  ls_c <- log(rexp(1L, 2L))
      al_c <- rnorm(J, mu_c, exp(ls_c));  be_c <- rnorm(2L, 0, 0.5)
    }

    eps_c  <- epsilon
    n_acc  <- 0L
    n_prop <- 0L

    if (transport == "reparam") {

      # ── Reparameterised HMC in (mu, ls, z, beta) ─────────────────────────────

      # Convert alpha to z at starting point
      lap_init <- tryCatch(
        .glmm_laplace(c(mu_c, ls_c), be_c, stan_data),
        error = function(e) NULL
      )
      if (!is.null(lap_init)) {
        z_c <- (al_c - lap_init$m) / lap_init$s
        m_c <- lap_init$m
      } else {
        z_c <- rep(0, J)
        m_c <- rep(mu_c, J)
      }
      alpha_c <- al_c

      for (iter in seq_len(n_total)) {

        # Gradient at current state; also refreshes m_c warm start and alpha_c
        g_c <- tryCatch(
          .glmm_reparam_grad(mu_c, ls_c, z_c, be_c, stan_data, m_init = m_c),
          error = function(e) NULL
        )
        if (is.null(g_c)) { n_prop <- n_prop + 1L; next }
        m_c     <- g_c$m
        alpha_c <- g_c$alpha

        p0  <- rnorm(P, 0, 1)
        H_s <- -g_c$lp + 0.5 * sum(p0^2)
        if (!is.finite(H_s)) { n_prop <- n_prop + 1L; next }

        # Leapfrog: identity mass matrix; pass gradient through to avoid redundant Newton
        q_new    <- c(mu_c, ls_c, z_c, be_c)
        p_new    <- p0
        grad_cur <- g_c$grad
        m_warm   <- g_c$m
        valid    <- TRUE
        lp_end   <- NA_real_
        alpha_end <- alpha_c

        for (l in seq_len(L)) {
          p_half <- p_new + (eps_c / 2) * grad_cur
          q_new  <- q_new + eps_c * p_half

          g_new <- tryCatch(
            .glmm_reparam_grad(q_new[1L], q_new[2L], q_new[idx_fiber],
                               q_new[idx_beta], stan_data, m_init = m_warm),
            error = function(e) NULL
          )
          if (is.null(g_new)) { valid <- FALSE; break }
          m_warm    <- g_new$m
          grad_cur  <- g_new$grad
          lp_end    <- g_new$lp
          alpha_end <- g_new$alpha

          p_new <- p_half + (eps_c / 2) * grad_cur
        }

        n_prop <- n_prop + 1L

        if (valid && is.finite(lp_end)) {
          H_e   <- -lp_end + 0.5 * sum(p_new^2)
          log_r <- H_s - H_e
          if (is.finite(log_r) && log(runif(1L)) < log_r) {
            mu_c    <- q_new[1L]
            ls_c    <- q_new[2L]
            z_c     <- q_new[idx_fiber]
            be_c    <- q_new[idx_beta]
            m_c     <- m_warm
            alpha_c <- alpha_end
            n_acc   <- n_acc + 1L
          }
        }

        if (iter <= n_warmup && iter %% 100L == 0L) {
          rate  <- n_acc / max(n_prop, 1L)
          eps_c <- .tune_step(eps_c, rate, target = target_rate,
                              factor = 1.3, max_step = 1.0)
          n_acc  <- 0L
          n_prop <- 0L
        }

        if (iter > n_warmup) {
          t <- iter - n_warmup
          all_chains[t, chain_id, ] <- c(mu_c, exp(ls_c), alpha_c, be_c)
        }
      }

    } else {

      # ── Fisher-legacy: connection-corrected HMC in (mu, ls, alpha, beta) ──────

      q_c <- c(mu_c, ls_c, al_c, be_c)

      for (iter in seq_len(n_total)) {

        pars <- .unpack(q_c, J)

        # 1. Diagonal metric at trajectory start — fixed for K and momentum draw
        G_d <- tryCatch(
          .glmm_diag_metric(pars$mu, pars$log_sigma,
                            pars$alpha, pars$beta, stan_data),
          error = function(e) NULL
        )
        if (is.null(G_d) || any(G_d <= 0) || !all(is.finite(G_d))) next

        # 2. Sample momentum p ~ N(0, G_d)
        p0 <- rnorm(P, 0, sqrt(G_d))

        # Kinetic energy — uses FIXED G_d throughout for exact detailed balance
        .K <- function(p) 0.5 * sum(p^2 / G_d)

        # 3. Hamiltonian at start
        lp_s <- tryCatch(
          .glmm_log_post(pars$mu, pars$log_sigma,
                         pars$alpha, pars$beta, stan_data),
          error = function(e) -Inf
        )
        if (!is.finite(lp_s)) { n_prop <- n_prop + 1L; next }
        H_s <- -lp_s + .K(p0)

        # 4. Horizontal leapfrog
        q_new <- q_c
        p_new <- p0
        valid <- TRUE

        for (l in seq_len(L)) {
          pars_l <- .unpack(q_new, J)

          grad_l <- tryCatch(
            .glmm_grad_vec(pars_l$mu, pars_l$log_sigma,
                           pars_l$alpha, pars_l$beta, stan_data),
            error = function(e) NULL
          )
          if (is.null(grad_l) || !all(is.finite(grad_l))) { valid <- FALSE; break }

          # Half-step momentum
          p_half <- p_new + (eps_c / 2) * grad_l

          # Velocity using FIXED G_d
          v <- p_half / G_d

          # Base displacements (in unconstrained parameterisation)
          delta_mu   <- eps_c * v[1L]
          delta_lsig <- eps_c * v[2L]

          # Horizontal correction: A maps (dmu, dsigma) → dalpha.
          # Chain rule: delta_sigma = sigma * delta_log_sigma.
          sigma_l     <- exp(pars_l$log_sigma)
          delta_sigma <- sigma_l * delta_lsig

          G_FF_l <- .glmm_G_FF(sigma_l, pars_l$alpha,
                                stan_data$X, stan_data$group, pars_l$beta)
          G_BF_l <- .glmm_G_BF(sigma_l, pars_l$mu, pars_l$alpha)
          A_l    <- .glmm_connection(G_FF_l, G_BF_l)   # J x 2

          horiz_corr <- as.vector(A_l %*% c(delta_mu, delta_sigma))

          # Full position update: standard velocity + horizontal correction on fiber
          q_new[idx_base]  <- q_new[idx_base]  + c(delta_mu, delta_lsig)
          q_new[idx_fiber] <- q_new[idx_fiber] + eps_c * v[idx_fiber] + horiz_corr
          q_new[idx_beta]  <- q_new[idx_beta]  + eps_c * v[idx_beta]

          if (!all(is.finite(q_new))) { valid <- FALSE; break }

          pars_l2 <- .unpack(q_new, J)
          grad_l2 <- tryCatch(
            .glmm_grad_vec(pars_l2$mu, pars_l2$log_sigma,
                           pars_l2$alpha, pars_l2$beta, stan_data),
            error = function(e) NULL
          )
          if (is.null(grad_l2) || !all(is.finite(grad_l2))) { valid <- FALSE; break }

          p_new <- p_half + (eps_c / 2) * grad_l2
        }

        # 5. Metropolis accept/reject — H uses FIXED G_d from start
        if (valid) {
          pars_e <- .unpack(q_new, J)
          lp_e   <- tryCatch(
            .glmm_log_post(pars_e$mu, pars_e$log_sigma,
                           pars_e$alpha, pars_e$beta, stan_data),
            error = function(e) -Inf
          )
          if (!is.finite(lp_e)) lp_e <- -Inf

          log_r <- H_s - (-lp_e + .K(p_new))
          if (is.finite(log_r) && log(runif(1L)) < log_r) {
            q_c   <- q_new
            n_acc <- n_acc + 1L
          }
        }
        n_prop <- n_prop + 1L

        # 6. Adaptive step size (every 100 warmup iterations)
        if (iter <= n_warmup && iter %% 100L == 0L) {
          rate  <- n_acc / max(n_prop, 1L)
          eps_c <- .tune_step(eps_c, rate, target = target_rate,
                              factor = 1.3, max_step = 1.0)
          n_acc  <- 0L
          n_prop <- 0L
        }

        # 7. Store post-warmup draws (sigma = exp(log_sigma) for output)
        if (iter > n_warmup) {
          t      <- iter - n_warmup
          q_s    <- q_c
          q_s[2L] <- exp(q_c[2L])
          all_chains[t, chain_id, ] <- q_s
        }
      }
    }

    if (verbose) {
      final_rate <- n_acc / max(n_prop, 1L)
      cat(sprintf("  epsilon=%.4f  acceptance=%.2f\n", eps_c, final_rate))
    }
  }

  posterior::as_draws_array(all_chains)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Multiplicative step-size adjustment toward target acceptance rate.
.tune_step <- function(step, rate, target, factor = 1.3, max_step = 0.5) {
  new_step <- if (rate > target) step * factor else step / factor
  min(new_step, max_step)
}
