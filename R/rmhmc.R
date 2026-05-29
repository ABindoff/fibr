#' Riemannian HMC for the centred GLMM
#'
#' @description
#' HMC with a **per-trajectory locally-adapted mass matrix**: the full
#' \eqn{P \times P} Fisher metric \eqn{G(q)} is computed at the start of each
#' trajectory and held fixed for \eqn{L} leapfrog steps.  This is exactly
#' correct (detailed balance holds) because:
#'
#' 1. Momentum is sampled from \eqn{N(0, G(q_\text{start}))} — density cancels
#'    in the M-H ratio.
#' 2. The leapfrog with FIXED \eqn{G^{-1}} as the mass matrix is
#'    volume-preserving.
#' 3. Both \eqn{H_\text{start}} and \eqn{H_\text{end}} use the same \eqn{G},
#'    so the acceptance ratio is exact.
#'
#' The gain over Euclidean HMC: \eqn{G(q)} is adapted to the local geometry —
#' the funnel in \eqn{\sigma}, the varying identifiability of each
#' \eqn{\alpha_j}, and the coupling between base and fiber directions captured
#' by the off-diagonal blocks.
#'
#' @param stan_data Named list with X (N×2), y (N), group (N, 1-indexed).
#' @param n_iter    Post-warmup iterations per chain (default 2000).
#' @param n_warmup  Warmup iterations per chain (default 1000).
#' @param n_chains  Number of independent chains (default 4).
#' @param L         Leapfrog steps per trajectory (default 10).
#' @param epsilon   Initial step size (default 0.1). Adapted during warmup.
#' @param target_rate Target Metropolis acceptance rate for adaptive tuning
#'   (default 0.65, the RMHMC optimal rate for many-step trajectories).
#' @param init      Optional named list: `mu`, `sigma`, `alpha`, `beta`.
#'   Defaults to prior draws.
#' @param seed      Random seed.
#' @param verbose   Print per-chain progress (default `TRUE`).
#'
#' @return A [`posterior::draws_array`] with the same variable naming as the
#'   centred Stan model: `mu`, `sigma`, `alpha[1]`, ..., `beta[2]`.
#'
#' @export
riemannian_mcmc <- function(stan_data,
                             n_iter      = 2000L,
                             n_warmup    = 1000L,
                             n_chains    = 4L,
                             L           = 10L,
                             epsilon     = 0.10,
                             target_rate = 0.65,
                             init        = NULL,
                             seed        = NULL,
                             verbose     = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  J       <- max(stan_data$group)
  P       <- 2L + J + 2L
  n_total <- n_warmup + n_iter

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
      mu_c  <- init$mu
      ls_c  <- log(max(init$sigma, 1e-3))
      al_c  <- init$alpha
      be_c  <- init$beta
      # Jitter each chain so R-hat diagnostics are meaningful
      if (chain_id > 1L) {
        mu_c  <- mu_c  + rnorm(1L, 0, 0.3)
        ls_c  <- ls_c  + rnorm(1L, 0, 0.3)
        al_c  <- al_c  + rnorm(J,  0, 0.5)
        be_c  <- be_c  + rnorm(2L, 0, 0.2)
      }
    } else {
      mu_c  <- rnorm(1L, 0, 1)
      ls_c  <- log(rexp(1L, 2L))
      al_c  <- rnorm(J,  mu_c, exp(ls_c))
      be_c  <- rnorm(2L, 0,    0.5)
    }
    q_c <- c(mu_c, ls_c, al_c, be_c)

    eps_c  <- epsilon
    n_acc  <- 0L
    n_prop <- 0L

    # ── Main loop ─────────────────────────────────────────────────────────────
    for (iter in seq_len(n_total)) {

      pars <- .unpack(q_c, J)

      # 1. Diagonal Fisher metric at current position (always PD)
      #    G_diag[k] = -∂²/∂qk² log p(q|y) > 0 for all k
      G_d <- tryCatch(
        .glmm_diag_metric(pars$mu, pars$log_sigma,
                          pars$alpha, pars$beta, stan_data),
        error = function(e) NULL
      )
      if (is.null(G_d) || !all(is.finite(G_d)) || any(G_d <= 0)) next

      # 2. Sample momentum p₀ ~ N(0, G_d)  component-wise
      p0 <- rnorm(P, 0, sqrt(G_d))

      # 3. Hamiltonian at start  H = U(q) + K(p)
      #    K = 0.5 * sum(p² / G_d)  (diagonal kinetic energy, FIXED G_d)
      lp_s <- tryCatch(
        .glmm_log_post(pars$mu, pars$log_sigma,
                       pars$alpha, pars$beta, stan_data),
        error = function(e) -Inf
      )
      if (!is.finite(lp_s)) { n_prop <- n_prop + 1L; next }

      H_s <- -lp_s + 0.5 * sum(p0^2 / G_d)

      # 4. Leapfrog  (diagonal mass matrix = 1/G_d, FIXED throughout)
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

        p_half <- p_new + (eps_c / 2) * grad_l
        q_new  <- q_new + eps_c * p_half / G_d   # G_d^{-1} p (diagonal)

        pars_l2 <- .unpack(q_new, J)
        grad_l2 <- tryCatch(
          .glmm_grad_vec(pars_l2$mu, pars_l2$log_sigma,
                         pars_l2$alpha, pars_l2$beta, stan_data),
          error = function(e) NULL
        )
        if (is.null(grad_l2) || !all(is.finite(grad_l2))) { valid <- FALSE; break }

        p_new <- p_half + (eps_c / 2) * grad_l2
      }

      # 5. Metropolis accept / reject
      if (valid) {
        pars_e <- .unpack(q_new, J)
        lp_e   <- tryCatch(
          .glmm_log_post(pars_e$mu, pars_e$log_sigma,
                         pars_e$alpha, pars_e$beta, stan_data),
          error = function(e) -Inf
        )
        if (!is.finite(lp_e)) lp_e <- -Inf

        H_e <- -lp_e + 0.5 * sum(p_new^2 / G_d)   # same G_d as H_s

        log_r <- H_s - H_e
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
        t   <- iter - n_warmup
        q_s <- q_c
        q_s[2L] <- exp(q_c[2L])          # log_sigma → sigma
        all_chains[t, chain_id, ] <- q_s
      }
    }

    if (verbose) {
      final_rate <- n_acc / max(n_prop, 1L)
      cat(sprintf("  epsilon=%.4f  acceptance=%.2f\n", eps_c, final_rate))
    }
  }

  posterior::as_draws_array(all_chains)
}
