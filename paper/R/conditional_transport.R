# Laplace conditional transport for the centred GLMM.
#
# Model: y_ij ~ Bernoulli(logit^{-1}(alpha_j + x_i' beta))
#        alpha_j ~ N(mu, sigma^2)
#
# Per group j, the Laplace approximation to alpha_j | theta, beta, data is
# N(m_j, s_j^2) with m_j the log-concave conditional mode and s_j = G_j^{-1/2}.
#
# Conditional (affine) transport between base points:
#   T_{theta -> theta'}(alpha)_j = m_j(theta') + (s_j(theta')/s_j(theta)) * (alpha_j - m_j(theta))
#
# Log-Jacobian: sum_j log(s_j(theta')/s_j(theta)).
# The map is exactly invertible: T_{theta'->theta} = T^{-1}_{theta->theta'}.
#
# All functions operate on a SINGLE parameter vector (one chain draw).
# theta = c(mu, log_sigma) throughout (sampler coordinates).

# ── Laplace approximation ─────────────────────────────────────────────────────

#' Newton mode and scale for the per-group Laplace conditional
#'
#' Finds the mode \eqn{m_j} of \eqn{l_j(\alpha_j)} — the per-group conditional
#' log-density — by Newton's method, and returns the Laplace scale
#' \eqn{s_j = G_j^{-1/2}} where \eqn{G_j = 1/\sigma^2 + s_j^{\rm lik}(m_j)}.
#'
#' Log-concavity of \eqn{l_j} guarantees a unique root and global convergence
#' of Newton with step-halving.
#'
#' @param theta   2-vector \eqn{(\mu, \log\sigma)}.
#' @param beta    K-vector of fixed-effect coefficients.
#' @param stan_data Named list with \code{X} (N×K), \code{y} (N), \code{group}
#'   (N, integer 1..J).
#' @param m_init  Optional J-vector warm start for Newton (default: all \eqn{\mu}).
#'
#' @return Named list:
#'   \describe{
#'     \item{m}{J-vector; per-group conditional modes.}
#'     \item{s}{J-vector; per-group Laplace scales (\eqn{G_j^{-1/2}}).}
#'     \item{G}{J-vector; per-group conditional precisions.}
#'   }
#' @keywords internal
.glmm_laplace <- function(theta, beta, stan_data, m_init = NULL) {
  mu    <- theta[1L]
  sigma <- exp(theta[2L])
  group <- stan_data$group
  X     <- stan_data$X
  y     <- stan_data$y
  J     <- max(group)

  m <- if (!is.null(m_init)) m_init else rep(mu, J)

  for (iter in seq_len(50L)) {
    eta  <- m[group] + as.vector(X %*% beta)
    p    <- plogis(eta)
    r_j  <- as.vector(tapply(y - p,       group, sum))
    s_jm <- as.vector(tapply(p * (1 - p), group, sum))
    G_j  <- 1 / sigma^2 + s_jm
    step <- ((mu - m) / sigma^2 + r_j) / G_j
    al   <- 1.0
    for (k in seq_len(10L)) {
      if (all(is.finite(m + al * step))) break
      al <- al / 2
    }
    m <- m + al * step
    if (max(abs(al * step)) < 1e-12) break
  }

  # Final-pass quantities at the converged mode
  eta  <- m[group] + as.vector(X %*% beta)
  p    <- plogis(eta)
  s_jm <- as.vector(tapply(p * (1 - p), group, sum))
  G_j  <- 1 / sigma^2 + s_jm

  list(m = m, s = 1 / sqrt(G_j), G = G_j)
}

# ── Affine transport ──────────────────────────────────────────────────────────

#' Conditional (Laplace affine) transport between two base points
#'
#' Applies the map
#' \deqn{T(\alpha)_j = m_j^{\rm to} + \frac{s_j^{\rm to}}{s_j^{\rm from}}
#'   (\alpha_j - m_j^{\rm from})}
#' and accumulates the log-Jacobian \eqn{\sum_j \log(s_j^{\rm to}/s_j^{\rm from})}.
#'
#' The map reduces to non-centering when there is no data (\eqn{s_j \to \sigma}),
#' to the identity as \eqn{\sigma \to \infty}, and is well-behaved as
#' \eqn{\sigma \to 0} (unlike the Fisher-metric lift).
#'
#' @param alpha    J-vector; current fiber draw.
#' @param lap_from Laplace object at the FROM base point (output of
#'   \code{\link{.glmm_laplace}}).
#' @param lap_to   Laplace object at the TO base point.
#'
#' @return Named list:
#'   \describe{
#'     \item{alpha_new}{J-vector; transported fiber draw.}
#'     \item{log_jac}{Scalar; \eqn{\sum_j \log(s^{\rm to}/s^{\rm from})}. Add
#'       to the log acceptance ratio.}
#'   }
#' @keywords internal
.glmm_cond_transport <- function(alpha, lap_from, lap_to) {
  ratio     <- lap_to$s / lap_from$s
  alpha_new <- lap_to$m + ratio * (alpha - lap_from$m)
  log_jac   <- sum(log(ratio))
  list(alpha_new = alpha_new, log_jac = log_jac)
}

# ── Implicit-function derivatives ─────────────────────────────────────────────

#' Implicit-function derivatives of the Laplace mode and log-scale
#'
#' At \eqn{a = m_j} (the Laplace mode), the implicit function theorem gives
#' exact partial derivatives of \eqn{m_j(\theta, \beta)} and
#' \eqn{\log s_j(\theta, \beta)} with respect to all base parameters.
#'
#' Derivatives with respect to \eqn{\log\sigma} follow by the chain rule:
#' \eqn{d/d(\log\sigma) = \sigma \cdot d/d\sigma}.
#'
#' @param theta     2-vector \eqn{(\mu, \log\sigma)}.
#' @param beta      K-vector.
#' @param stan_data Named list with \code{X}, \code{y}, \code{group}.
#' @param lap       Output of \code{\link{.glmm_laplace}} at the same
#'   \code{(theta, beta)}.
#'
#' @return Named list of J-vectors (and J×K matrices for beta components):
#'   \code{dm_dmu}, \code{dm_dsig}, \code{dm_db} (J×K),
#'   \code{dls_dmu}, \code{dls_dsig}, \code{dls_db} (J×K).
#'   Multiply \code{*_dsig} by \eqn{\sigma} for the \eqn{\log\sigma} direction.
#' @keywords internal
.glmm_lap_derivs <- function(theta, beta, stan_data, lap) {
  sigma <- exp(theta[2L])
  group <- stan_data$group
  X     <- stan_data$X
  m     <- lap$m
  G_j   <- lap$G

  eta <- m[group] + as.vector(X %*% beta)
  p   <- plogis(eta)
  t_j <- as.vector(tapply(p * (1 - p) * (1 - 2 * p), group, sum))

  # Derivatives of m_j
  dm_dmu  <- (1 / sigma^2) / G_j
  dm_dsig <- 2 * (m - theta[1L]) / (sigma^3 * G_j)

  K     <- length(beta)
  dm_db <- matrix(0, length(m), K)
  for (k in seq_len(K)) {
    wt <- p * (1 - p) * X[, k]
    dm_db[, k] <- -as.vector(tapply(wt, group, sum)) / G_j
  }

  # Derivatives of G_j (needed for dlog s = -(1/2) dG/G)
  dG_dmu  <- t_j * dm_dmu
  dG_dsig <- -2 / sigma^3 + t_j * dm_dsig

  dG_db <- matrix(0, length(m), K)
  for (k in seq_len(K)) {
    wt <- p * (1 - p) * (1 - 2 * p) * X[, k]
    dG_db[, k] <- as.vector(tapply(wt, group, sum)) + t_j * dm_db[, k]
  }

  list(
    dm_dmu   = dm_dmu,
    dm_dsig  = dm_dsig,
    dm_db    = dm_db,
    dls_dmu  = -0.5 * dG_dmu  / G_j,
    dls_dsig = -0.5 * dG_dsig / G_j,
    dls_db   = -0.5 * dG_db   / G_j
  )
}

# ── Reparameterised HMC gradient ─────────────────────────────────────────────

#' Gradient of the reparameterised log-posterior for HMC in z-coordinates
#'
#' Computes the gradient of
#' \deqn{\tilde{\ell}(\mu, \lambda, z, \beta) =
#'   \ell(\mu, \lambda, \alpha(z), \beta) + \sum_j \log s_j(\theta, \beta)}
#' where \eqn{\alpha_j = m_j + s_j z_j}, \eqn{\lambda = \log\sigma}, and
#' \eqn{(m_j, s_j)} are the Laplace mode and scale.
#'
#' Chain rule:
#' \itemize{
#'   \item \eqn{\partial\tilde\ell/\partial z_j = g_j s_j}
#'   \item \eqn{\partial\tilde\ell/\partial\mu  = g_\mu
#'     + \sum_j\bigl[g_j\,\partial m_j/\partial\mu
#'     + (g_j z_j s_j + 1)\,\partial\log s_j/\partial\mu\bigr]}
#'   \item \eqn{\partial\tilde\ell/\partial\lambda = g_\lambda
#'     + \sigma\sum_j\bigl[g_j\,\partial m_j/\partial\sigma
#'     + (g_j z_j s_j + 1)\,\partial\log s_j/\partial\sigma\bigr]}
#'   \item \eqn{\partial\tilde\ell/\partial\beta_k = g_{\beta_k}
#'     + \sum_j\bigl[g_j\,\partial m_j/\partial\beta_k
#'     + (g_j z_j s_j + 1)\,\partial\log s_j/\partial\beta_k\bigr]}
#' }
#' where \eqn{g_j = \partial\ell/\partial\alpha_j} (from the alpha-space gradient)
#' and the partial derivatives of \eqn{m_j}, \eqn{\log s_j} come from
#' \code{\link{.glmm_lap_derivs}}.
#'
#' @param mu    Scalar.
#' @param ls    Scalar; \eqn{\log\sigma}.
#' @param z     J-vector; standardised fiber coordinates.
#' @param beta  K-vector.
#' @param stan_data Named list with \code{X}, \code{y}, \code{group}.
#' @param m_init Optional J-vector warm start for Newton (default: \code{NULL}).
#'
#' @return Named list:
#'   \describe{
#'     \item{grad}{(2+J+K)-vector gradient of \eqn{\tilde\ell} w.r.t.
#'       \eqn{(\mu, \lambda, z, \beta)}.}
#'     \item{lp}{\eqn{\tilde\ell} evaluated at the current point.}
#'     \item{m}{J-vector Laplace mode (for warm-starting the next call).}
#'     \item{alpha}{J-vector \eqn{\alpha = m + s \odot z}.}
#'   }
#'   Returns \code{NULL} if Newton or any downstream computation fails.
#' @keywords internal
.glmm_reparam_grad <- function(mu, ls, z, beta, stan_data, m_init = NULL) {
  sigma  <- exp(ls)
  theta  <- c(mu, ls)

  lap <- tryCatch(
    .glmm_laplace(theta, beta, stan_data, m_init = m_init),
    error = function(e) NULL
  )
  if (is.null(lap)) return(NULL)

  alpha <- lap$m + lap$s * z

  g <- tryCatch(
    .glmm_grad_log_post(mu, ls, alpha, beta, stan_data),
    error = function(e) NULL
  )
  if (is.null(g)) return(NULL)

  d  <- .glmm_lap_derivs(theta, beta, stan_data, lap)
  w  <- g$grad_alpha * z * lap$s + 1   # common weight: g_j z_j s_j + 1

  grad <- c(
    g$grad_mu + sum(g$grad_alpha * d$dm_dmu + w * d$dls_dmu),
    g$grad_log_sigma + sigma * sum(g$grad_alpha * d$dm_dsig + w * d$dls_dsig),
    g$grad_alpha * lap$s,
    g$grad_beta + colSums(g$grad_alpha * d$dm_db + w * d$dls_db)
  )
  if (!all(is.finite(grad))) return(NULL)

  lp_t <- tryCatch(
    .glmm_log_post(mu, ls, alpha, beta, stan_data) + sum(log(lap$s)),
    error = function(e) -Inf
  )
  if (!is.finite(lp_t)) return(NULL)

  list(grad = grad, lp = lp_t, m = lap$m, alpha = alpha)
}

# ── Non-Gaussianity triage statistic ──────────────────────────────────────────

#' Per-group conditional non-Gaussianity statistic
#'
#' \deqn{\kappa_j = -t_j(m_j)\, s_j^3}
#'
#' where \eqn{t_j = \sum_{i \in j} p_{ij}(1-p_{ij})(1-2p_{ij})} is the
#' third-derivative contribution evaluated at the Laplace mode.  The prior
#' contributes zero to the third derivative.
#'
#' \eqn{|\kappa_j|} measures how far the per-group conditional departs from
#' Gaussian: small values indicate the Laplace transport is nearly exact;
#' large values indicate structural non-Gaussianity that the affine map cannot
#' correct.  Use as a triage predictor for residual transport error.
#'
#' @param lap       Output of \code{\link{.glmm_laplace}}.
#' @param beta      K-vector.
#' @param stan_data Named list with \code{X}, \code{y}, \code{group}.
#'
#' @return J-vector of \eqn{\kappa_j} values.
#' @keywords internal
.glmm_kappa <- function(lap, beta, stan_data) {
  group <- stan_data$group
  X     <- stan_data$X
  m     <- lap$m

  eta <- m[group] + as.vector(X %*% beta)
  p   <- plogis(eta)
  t_j <- as.vector(tapply(p * (1 - p) * (1 - 2 * p), group, sum))

  -t_j * lap$s^3
}
