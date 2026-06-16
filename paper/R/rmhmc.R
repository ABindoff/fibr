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
#' @param method     Metric to use: `"diagonal"` (default, always PD, fast) or
#'   `"softabs"` (SoftAbs regularised full metric; captures off-diagonal coupling).
#' @param softabs_alpha  SoftAbs smoothness parameter alpha (default 1).
#'   Larger values give a sharper approximation to the absolute-value metric.
#' @param target_rate Target Metropolis acceptance rate for adaptive tuning
#'   (default 0.65). For `L = 1` (MALA) use 0.57.
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
                             n_iter        = 2000L,
                             n_warmup      = 1000L,
                             n_chains      = 4L,
                             L             = 10L,
                             epsilon       = 0.10,
                             method        = c("diagonal", "softabs"),
                             softabs_alpha = 1.0,
                             target_rate   = 0.65,
                             init          = NULL,
                             seed          = NULL,
                             verbose       = TRUE) {

  method <- match.arg(method)

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

      # 1. Compute metric at current position
      #    "diagonal": diagonal G, always PD, O(N+J) per call
      #    "softabs" : full SoftAbs G, always PD, captures off-diagonal coupling
      metric <- tryCatch({
        if (method == "diagonal") {
          G_d <- .glmm_diag_metric(pars$mu, pars$log_sigma,
                                   pars$alpha, pars$beta, stan_data)
          if (any(G_d <= 0) || !all(is.finite(G_d))) stop("bad diag")
          list(type = "diagonal", G_d = G_d)
        } else {
          blocks <- .glmm_full_metric(pars$mu, pars$log_sigma,
                                      pars$alpha, pars$beta, stan_data)
          sa     <- .softabs_decomp(blocks, J, softabs_alpha)
          if (!all(is.finite(sa$lambda_sa)) || any(sa$lambda_sa <= 0)) stop("bad sa")
          list(type = "softabs", U = sa$U, lsa = sa$lambda_sa)
        }
      }, error = function(e) NULL)
      if (is.null(metric)) next

      # 2. Sample momentum  p₀ ~ N(0, G)
      #    diagonal: component-wise  softabs: p = U (sqrt(λ_sa) ⊙ z)
      p0 <- if (metric$type == "diagonal") {
        rnorm(P, 0, sqrt(metric$G_d))
      } else {
        as.vector(metric$U %*% (sqrt(metric$lsa) * rnorm(P)))
      }

      # Helper: kinetic energy K = 0.5 p^T G^{-1} p  (same G throughout)
      .K <- function(p) {
        if (metric$type == "diagonal") {
          0.5 * sum(p^2 / metric$G_d)
        } else {
          v <- as.vector(t(metric$U) %*% p)
          0.5 * sum(v^2 / metric$lsa)
        }
      }

      # Helper: velocity v = G^{-1} p
      .vel <- function(p) {
        if (metric$type == "diagonal") {
          p / metric$G_d
        } else {
          v <- as.vector(t(metric$U) %*% p)
          as.vector(metric$U %*% (v / metric$lsa))
        }
      }

      # 3. Hamiltonian at start  H_riem = U(q) + K(p)
      lp_s <- tryCatch(
        .glmm_log_post(pars$mu, pars$log_sigma,
                       pars$alpha, pars$beta, stan_data),
        error = function(e) -Inf
      )
      if (!is.finite(lp_s)) { n_prop <- n_prop + 1L; next }
      H_s <- -lp_s + .K(p0)

      # 4. Leapfrog (FIXED metric throughout the trajectory)
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
        q_new  <- q_new + eps_c * .vel(p_half)

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
