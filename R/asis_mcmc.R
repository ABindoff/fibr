# Ancillarity-Sufficiency Interweaving Strategy (ASIS) MwG for the centred GLMM.
#
# Per iteration, each base-block update is doubled: one proposal with alpha
# held fixed (centred frame) and one proposal with phi = (alpha-mu)/sigma held
# fixed (non-centred frame).  The interleaving ensures the sampler benefits
# from the centering in data-dominated cells and from non-centering in
# prior-dominated cells, without choosing statically.
#
# Reference: Papaspiliopoulos, Roberts & Skold (2007) JRSSB 69(2).
#            Yu & Meng (2011) JCGS 20(3) for the formal ASIS proof.

#' ASIS MwG sampler for the centred GLMM
#'
#' @description
#' Block Metropolis-within-Gibbs with **ancillarity-sufficiency interweaving**:
#' each iteration runs a centred base-block update followed immediately by a
#' non-centred base-block update (phi = (alpha-mu)/sigma held constant, alpha
#' recomputed from the new (mu,sigma)).  The fiber and beta blocks are centred
#' as in \code{\link{horizontal_mcmc}}.
#'
#' @param stan_data   Named list with \code{X} (N×2), \code{y} (N),
#'   \code{group} (N, 1-indexed).
#' @param n_iter      Post-warmup iterations per chain (default 2000).
#' @param n_warmup    Warmup iterations per chain (default 1000).
#' @param n_chains    Number of independent chains (default 4).
#' @param step_base   Initial RW SD for \eqn{(\mu, \log\sigma)} in both blocks
#'   (default 0.10).
#' @param step_alpha  Initial RW SD for each \eqn{\alpha_j} (default 0.30).
#' @param step_beta   Initial RW SD for \eqn{\beta} (default 0.15).
#' @param target_rate Target Metropolis acceptance rate (default 0.30).
#' @param init        Optional named list: \code{mu}, \code{sigma},
#'   \code{alpha}, \code{beta}.
#' @param seed        Random seed.
#' @param verbose     Print per-chain progress (default \code{TRUE}).
#'
#' @return A \code{\link[posterior]{draws_array}} with variable naming matching
#'   the centred Stan model: \code{mu}, \code{sigma}, \code{alpha[1]}, ...,
#'   \code{beta[1]}, \code{beta[2]}.
#'
#' @export
asis_mcmc <- function(stan_data,
                       n_iter      = 2000L,
                       n_warmup    = 1000L,
                       n_chains    = 4L,
                       step_base   = 0.10,
                       step_alpha  = 0.30,
                       step_beta   = 0.15,
                       target_rate = 0.30,
                       init        = NULL,
                       seed        = NULL,
                       verbose     = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  J       <- max(stan_data$group)
  K       <- ncol(stan_data$X)
  P       <- 2L + J + K
  n_total <- n_warmup + n_iter

  par_names <- c("mu", "sigma",
                 paste0("alpha[", seq_len(J), "]"),
                 paste0("beta[",  seq_len(K), "]"))

  all_chains <- array(NA_real_,
                      dim      = c(n_iter, n_chains, P),
                      dimnames = list(NULL, NULL, par_names))

  for (chain_id in seq_len(n_chains)) {
    if (verbose) cat(sprintf("Chain %d/%d ...\n", chain_id, n_chains))

    # ── Initialise ─────────────────────────────────────────────────────────────
    if (!is.null(init)) {
      mu_c  <- init$mu;  ls_c <- log(max(init$sigma, 1e-3))
      al_c  <- init$alpha; be_c <- init$beta
      if (chain_id > 1L) {
        mu_c  <- mu_c  + rnorm(1L, 0, 0.3); ls_c <- ls_c + rnorm(1L, 0, 0.3)
        al_c  <- al_c  + rnorm(J, 0, 0.2);  be_c <- be_c + rnorm(K, 0, 0.2)
      }
    } else {
      mu_c  <- rnorm(1L, 0, 1); ls_c <- log(rexp(1L, 2L))
      al_c  <- rnorm(J, mu_c, exp(ls_c)); be_c <- rnorm(K, 0, 0.5)
    }

    lp_c <- .glmm_log_post(mu_c, ls_c, al_c, be_c, stan_data)
    if (!is.finite(lp_c)) {
      mu_c <- 0; ls_c <- 0; al_c <- rnorm(J, 0, 1); be_c <- rep(0, K)
      lp_c <- .glmm_log_post(mu_c, ls_c, al_c, be_c, stan_data)
    }

    s_base  <- step_base; s_alpha <- step_alpha; s_beta <- step_beta

    acc_base  <- 0L; n_base  <- 0L
    acc_nc    <- 0L; n_nc    <- 0L
    acc_alpha <- 0L; n_alpha <- 0L
    acc_beta  <- 0L; n_beta  <- 0L

    for (iter in seq_len(n_total)) {

      # ── 1. Fiber block: per-group alpha update (centred) ────────────────────
      for (j in seq_len(J)) {
        al_prop      <- al_c
        al_prop[j]   <- al_c[j] + rnorm(1L, 0, s_alpha)
        lp_prop      <- .glmm_log_post(mu_c, ls_c, al_prop, be_c, stan_data)
        log_r        <- lp_prop - lp_c
        n_alpha      <- n_alpha + 1L
        if (is.finite(log_r) && log(runif(1L)) < log_r) {
          al_c  <- al_prop; lp_c <- lp_prop
          acc_alpha <- acc_alpha + 1L
        }
      }

      # ── 2. Base block: (mu, log_sigma) update — centred frame ───────────────
      prop_base <- c(mu_c, ls_c) + rnorm(2L, 0, s_base)
      lp_prop   <- .glmm_log_post(prop_base[1L], prop_base[2L], al_c, be_c, stan_data)
      log_r     <- lp_prop - lp_c
      n_base    <- n_base + 1L
      if (is.finite(log_r) && log(runif(1L)) < log_r) {
        mu_c <- prop_base[1L]; ls_c <- prop_base[2L]; lp_c <- lp_prop
        acc_base <- acc_base + 1L
      }

      # ── 3. NC interweaving block: (mu, log_sigma) update — phi held fixed ───
      # phi_j = (alpha_j - mu) / sigma;  proposed alpha' = mu' + sigma' * phi
      phi <- (al_c - mu_c) / exp(ls_c)
      prop_nc   <- c(mu_c, ls_c) + rnorm(2L, 0, s_base)
      al_nc     <- prop_nc[1L] + exp(prop_nc[2L]) * phi
      lp_prop   <- .glmm_log_post(prop_nc[1L], prop_nc[2L], al_nc, be_c, stan_data)
      log_r     <- lp_prop - lp_c
      n_nc      <- n_nc + 1L
      if (is.finite(log_r) && log(runif(1L)) < log_r) {
        mu_c <- prop_nc[1L]; ls_c <- prop_nc[2L]; al_c <- al_nc; lp_c <- lp_prop
        acc_nc <- acc_nc + 1L
      }

      # ── 4. Beta block ────────────────────────────────────────────────────────
      be_prop <- be_c + rnorm(K, 0, s_beta)
      lp_prop <- .glmm_log_post(mu_c, ls_c, al_c, be_prop, stan_data)
      log_r   <- lp_prop - lp_c
      n_beta  <- n_beta + 1L
      if (is.finite(log_r) && log(runif(1L)) < log_r) {
        be_c <- be_prop; lp_c <- lp_prop
        acc_beta <- acc_beta + 1L
      }

      # ── Step-size adaptation ─────────────────────────────────────────────────
      if (iter <= n_warmup && iter %% 100L == 0L) {
        s_base  <- .tune_step(s_base,  (acc_base + acc_nc) / (n_base + n_nc), target_rate)
        s_alpha <- .tune_step(s_alpha, acc_alpha / max(n_alpha, 1L), target_rate)
        s_beta  <- .tune_step(s_beta,  acc_beta  / max(n_beta,  1L), target_rate)
        acc_base <- 0L; acc_nc <- 0L; n_base <- 0L; n_nc <- 0L
        acc_alpha <- 0L; n_alpha <- 0L
        acc_beta  <- 0L; n_beta  <- 0L
      }

      if (iter > n_warmup) {
        t <- iter - n_warmup
        all_chains[t, chain_id, ] <- c(mu_c, exp(ls_c), al_c, be_c)
      }
    }

    if (verbose) {
      cat(sprintf("  step_base=%.4f  step_alpha=%.4f  step_beta=%.4f\n",
                  s_base, s_alpha, s_beta))
      cat(sprintf("  acceptance: base=%.2f  nc=%.2f  alpha=%.2f  beta=%.2f\n",
                  (acc_base + acc_nc) / max(n_base + n_nc, 1L),
                  acc_nc / max(n_nc, 1L),
                  acc_alpha / max(n_alpha, 1L),
                  acc_beta  / max(n_beta,  1L)))
    }
  }

  posterior::as_draws_array(all_chains)
}
