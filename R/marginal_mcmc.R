# Quadrature-marginalised sampler for the centred GLMM.
#
# Integrates out the random effects alpha_j analytically (per-group adaptive
# Gauss-Hermite quadrature, 15 nodes, centred at the Laplace mode m_j,
# scaled by s_j sqrt(2)) and runs an adaptive random-walk Metropolis sampler
# over (mu, log_sigma, beta). For each stored draw the alpha_j are sampled
# exactly from their conditional posterior via inverse-CDF on a 400-point grid.
#
# Gold standard reference: treat as reference, not competitor, in figures.

# ── Gauss-Hermite nodes and weights ──────────────────────────────────────────

#' Gauss-Hermite nodes and weights (physicist's convention)
#'
#' Computes nodes and weights for \eqn{\int f(x)\exp(-x^2)\,dx \approx
#' \sum_k w_k f(x_k)} via the Golub-Welsch algorithm (tridiagonal eigenproblem
#' for the Hermite recurrence).
#'
#' @param n Integer; number of nodes (default 15).
#' @return Named list: \code{x} (n-vector of nodes, sorted ascending),
#'   \code{w} (n-vector of positive weights; \code{sum(w) = sqrt(pi)}).
#' @keywords internal
.gh_nodes_weights <- function(n = 15L) {
  off  <- sqrt(seq_len(n - 1L) / 2)
  Jmat <- matrix(0, n, n)
  idx  <- cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)
  Jmat[idx]             <- off
  Jmat[idx[, c(2L, 1L)]] <- off
  eig <- eigen(Jmat, symmetric = TRUE)
  x   <- eig$values
  w   <- sqrt(pi) * eig$vectors[1L, ]^2
  ord <- order(x)
  list(x = x[ord], w = w[ord])
}

# ── Marginal log-posterior ────────────────────────────────────────────────────

#' Marginal log-posterior for base parameters (alpha integrated out)
#'
#' Computes \eqn{\log p(\mu, \log\sigma, \beta \mid y)} by numerically
#' integrating out all \eqn{\alpha_j} using adaptive Gauss-Hermite quadrature
#' (nodes centred at the per-group Laplace mode \eqn{m_j}, scaled by
#' \eqn{s_j\sqrt{2}}).
#'
#' Transformation: \eqn{\alpha_j = m_j + s_j\sqrt{2}\,t}, so that
#' \eqn{d\alpha_j = s_j\sqrt{2}\,dt} and the integrand becomes
#' \eqn{\exp(l_j(m_j + s_j\sqrt{2}\,t) + t^2) \cdot \exp(-t^2)} — exactly the
#' form expected by the physicist's GH rule.
#'
#' @param mu    Scalar.
#' @param ls    Scalar; \eqn{\log\sigma}.
#' @param beta  K-vector.
#' @param stan_data Named list with \code{X}, \code{y}, \code{group}.
#' @param gh    Output of \code{\link{.gh_nodes_weights}}.
#' @param m_init Optional J-vector Newton warm start.
#'
#' @return Named list:
#'   \describe{
#'     \item{lp}{Scalar; marginal log-posterior, or \code{-Inf} on failure.}
#'     \item{m}{J-vector Laplace mode (for warm-starting the next call), or
#'       \code{NULL} on failure.}
#'   }
#' @keywords internal
.marginal_log_post <- function(mu, ls, beta, stan_data, gh, m_init = NULL) {
  sigma <- exp(ls)
  theta <- c(mu, ls)

  lap <- tryCatch(
    .glmm_laplace(theta, beta, stan_data, m_init = m_init),
    error = function(e) NULL
  )
  if (is.null(lap)) return(list(lp = -Inf, m = NULL))

  J        <- max(stan_data$group)
  log_ml   <- 0
  eta_base <- as.vector(stan_data$X %*% beta)   # N-vector: x_i' beta

  for (j in seq_len(J)) {
    jj    <- stan_data$group == j
    y_j   <- stan_data$y[jj]
    eb_j  <- eta_base[jj]
    m_j   <- lap$m[j]; s_j <- lap$s[j]
    scale <- s_j * sqrt(2)

    a_nodes <- m_j + scale * gh$x   # length(gh$x) alpha values

    # l_j(a) + x_k^2: log N(a; mu, sigma^2) + lik_j(a) + x_k^2
    lv <- vapply(seq_along(gh$x), function(k) {
      a   <- a_nodes[k]
      eta <- a + eb_j
      dnorm(a, mu, sigma, log = TRUE) +
        sum(y_j * eta - log1p(exp(eta))) +
        gh$x[k]^2
    }, numeric(1L))

    max_lv <- max(lv)
    log_ml <- log_ml + log(scale) + max_lv +
              log(sum(gh$w * exp(lv - max_lv)))
  }

  lp_base <- dnorm(mu, 0, 5, log = TRUE) +
             dexp(sigma, rate = 1, log = TRUE) + ls +
             sum(dnorm(beta, 0, 2, log = TRUE))

  list(lp = lp_base + log_ml, m = lap$m)
}

# ── Conditional alpha sampler ─────────────────────────────────────────────────

#' Sample alpha from the conditional posterior via inverse-CDF on a fine grid
#'
#' For each group \eqn{j}, evaluates the per-group conditional log-density
#' \eqn{l_j(\alpha_j)} on a 400-point equally-spaced grid
#' \eqn{[m_j - 6s_j,\; m_j + 6s_j]}, normalises to a PMF, and draws one
#' sample by inverse-CDF with linear interpolation within each grid cell.
#' Linear interpolation removes the \eqn{O(\delta)} left-endpoint bias that
#' arises from treating the discrete PMF as a step-function CDF.
#'
#' @param theta     2-vector \eqn{(\mu, \log\sigma)}.
#' @param beta      K-vector.
#' @param stan_data Named list with \code{X}, \code{y}, \code{group}.
#' @param lap       Laplace object at \code{(theta, beta)} (output of
#'   \code{\link{.glmm_laplace}}).
#'
#' @return J-vector of \eqn{\alpha_j} draws.
#' @keywords internal
.sample_alpha_cond <- function(theta, beta, stan_data, lap) {
  mu    <- theta[1L]; sigma <- exp(theta[2L])
  J     <- max(stan_data$group)
  eta_base <- as.vector(stan_data$X %*% beta)
  alpha_out <- numeric(J)

  for (j in seq_len(J)) {
    jj    <- stan_data$group == j
    y_j   <- stan_data$y[jj]
    eb_j  <- eta_base[jj]
    m_j   <- lap$m[j]; s_j <- lap$s[j]

    a_grid <- seq(m_j - 6 * s_j, m_j + 6 * s_j, length.out = 400L)

    lv <- vapply(a_grid, function(a) {
      eta <- a + eb_j
      dnorm(a, mu, sigma, log = TRUE) + sum(y_j * eta - log1p(exp(eta)))
    }, numeric(1L))

    lv   <- lv - max(lv)
    prob <- exp(lv); prob <- prob / sum(prob)
    cdf  <- cumsum(prob)
    u    <- runif(1L)
    k    <- which(cdf >= u)

    if (length(k) == 0L) {
      alpha_out[j] <- a_grid[400L]
    } else {
      k_star <- k[1L]
      if (k_star <= 1L) {
        alpha_out[j] <- a_grid[1L]
      } else {
        # linear interpolation within cell [a_grid[k-1], a_grid[k]]
        p_lo  <- cdf[k_star - 1L]
        denom <- cdf[k_star] - p_lo
        if (denom < 1e-15) {
          alpha_out[j] <- a_grid[k_star - 1L]
        } else {
          alpha_out[j] <- a_grid[k_star - 1L] +
            (a_grid[k_star] - a_grid[k_star - 1L]) * (u - p_lo) / denom
        }
      }
    }
  }

  alpha_out
}

# ── Main sampler ──────────────────────────────────────────────────────────────

#' Quadrature-marginalised sampler for the centred GLMM
#'
#' @description
#' Marginalises the random effects \eqn{\alpha_j} analytically via adaptive
#' 15-point Gauss-Hermite quadrature (nodes centred at the per-group Laplace
#' mode \eqn{m_j}, scaled by \eqn{s_j\sqrt{2}}) and runs a scalar-step
#' adaptive random-walk Metropolis sampler over the base parameter vector
#' \eqn{(\mu, \log\sigma, \beta)}.  For each stored post-warmup draw, the
#' \eqn{\alpha_j} are drawn exactly from their conditional posterior via
#' inverse-CDF on a 400-point grid over \eqn{[m_j \pm 6s_j]}.
#'
#' This is the gold-standard reference for the model class.  Treat it as a
#' reference, not a competitor, in benchmark figures.
#'
#' @param stan_data   Named list with \code{X} (N×K), \code{y} (N),
#'   \code{group} (N, 1-indexed).
#' @param n_iter      Post-warmup iterations per chain (default 2000).
#' @param n_warmup    Warmup iterations per chain (default 1000).
#' @param n_chains    Number of independent chains (default 4).
#' @param step        Initial RW step size for all base parameters (default
#'   0.10).  Adapted every 100 warmup iterations toward \code{target_rate}.
#' @param target_rate Target Metropolis acceptance rate (default 0.23;
#'   near-optimal for a \eqn{d \geq 4} Gaussian target).
#' @param init        Optional named list: \code{mu}, \code{sigma},
#'   \code{beta}.  Defaults to prior draws.
#' @param seed        Random seed.
#' @param verbose     Print per-chain progress (default \code{TRUE}).
#'
#' @return A \code{\link[posterior]{draws_array}} with the same variable
#'   naming as the centred Stan model: \code{mu}, \code{sigma},
#'   \code{alpha[1]}, \ldots, \code{beta[1]}, \ldots
#'
#' @export
marginal_mcmc <- function(stan_data,
                           n_iter      = 2000L,
                           n_warmup    = 1000L,
                           n_chains    = 4L,
                           step        = 0.10,
                           target_rate = 0.23,
                           init        = NULL,
                           seed        = NULL,
                           verbose     = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  gh <- .gh_nodes_weights(15L)

  J       <- max(stan_data$group)
  K       <- ncol(stan_data$X)
  P_base  <- 2L + K
  P_out   <- 2L + J + K
  n_total <- n_warmup + n_iter

  par_names <- c("mu", "sigma",
                 paste0("alpha[", seq_len(J), "]"),
                 paste0("beta[",  seq_len(K), "]"))

  all_chains <- array(NA_real_,
                      dim      = c(n_iter, n_chains, P_out),
                      dimnames = list(NULL, NULL, par_names))

  for (chain_id in seq_len(n_chains)) {
    if (verbose) cat(sprintf("Chain %d/%d ...\n", chain_id, n_chains))

    # ── Initialise ─────────────────────────────────────────────────────────────
    if (!is.null(init)) {
      mu_c <- init$mu; ls_c <- log(max(init$sigma, 1e-3)); be_c <- init$beta
      if (chain_id > 1L) {
        mu_c <- mu_c + rnorm(1L, 0, 0.3); ls_c <- ls_c + rnorm(1L, 0, 0.3)
        be_c <- be_c + rnorm(K, 0, 0.2)
      }
    } else {
      mu_c <- rnorm(1L, 0, 1); ls_c <- log(rexp(1L, 2L))
      be_c <- rnorm(K, 0, 0.5)
    }

    r_c  <- .marginal_log_post(mu_c, ls_c, be_c, stan_data, gh)
    lp_c <- r_c$lp; m_c <- r_c$m
    # Fall back to origin if initial point has -Inf log-posterior
    if (!is.finite(lp_c)) {
      mu_c <- 0; ls_c <- 0; be_c <- rep(0, K)
      r_c  <- .marginal_log_post(mu_c, ls_c, be_c, stan_data, gh)
      lp_c <- r_c$lp; m_c <- r_c$m
    }

    eps_c  <- step
    n_acc  <- 0L
    n_prop <- 0L

    for (iter in seq_len(n_total)) {

      q_prop  <- c(mu_c, ls_c, be_c) + rnorm(P_base, 0, eps_c)
      r_prop  <- .marginal_log_post(q_prop[1L], q_prop[2L],
                                     q_prop[3L:(2L + K)],
                                     stan_data, gh, m_init = m_c)
      lp_prop <- r_prop$lp
      m_prop  <- r_prop$m

      log_r <- lp_prop - lp_c
      if (is.finite(log_r) && log(runif(1L)) < log_r) {
        mu_c <- q_prop[1L]; ls_c <- q_prop[2L]
        be_c <- q_prop[3L:(2L + K)]
        lp_c <- lp_prop; m_c <- m_prop
        n_acc <- n_acc + 1L
      }
      n_prop <- n_prop + 1L

      if (iter <= n_warmup && iter %% 100L == 0L) {
        rate  <- n_acc / max(n_prop, 1L)
        eps_c <- .tune_step(eps_c, rate, target = target_rate,
                            factor = 1.3, max_step = 1.0)
        n_acc  <- 0L; n_prop <- 0L
      }

      if (iter > n_warmup) {
        t     <- iter - n_warmup
        lap_c <- tryCatch(
          .glmm_laplace(c(mu_c, ls_c), be_c, stan_data, m_init = m_c),
          error = function(e) NULL
        )
        alpha_c <- if (!is.null(lap_c)) {
          .sample_alpha_cond(c(mu_c, ls_c), be_c, stan_data, lap_c)
        } else {
          rnorm(J, mu_c, exp(ls_c))
        }
        all_chains[t, chain_id, ] <- c(mu_c, exp(ls_c), alpha_c, be_c)
      }
    }

    if (verbose) {
      final_rate <- n_acc / max(n_prop, 1L)
      cat(sprintf("  step=%.4f  acceptance=%.2f\n", eps_c, final_rate))
    }
  }

  posterior::as_draws_array(all_chains)
}
