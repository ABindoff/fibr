# Analytic Fisher metric blocks and curvature for the centred GLMM.
#
# Model:  y_ij ~ Bernoulli(logit^{-1}(alpha_j + X_i beta))
#         alpha_j ~ N(mu, sigma^2)
#
# All functions operate on a SINGLE parameter vector (one chain draw).
# They are vectorised over j (groups) but not over draws — the caller
# loops over draws or applies via apply().

# ── Fisher metric blocks ──────────────────────────────────────────────────────

#' Diagonal of the fiber metric block G_FF at one parameter point
#'
#' G_FF\[j,j\] = 1/sigma^2  +  sum_{i in group_j} p_ij * (1 - p_ij)
#'
#' @param sigma  Scalar; group-level SD.
#' @param alpha  J-vector; group intercepts.
#' @param X      N x 2 matrix; fixed-effect design matrix.
#' @param group  N-vector (integer 1..J); group membership.
#' @param beta   2-vector; fixed-effect coefficients.
#' @return J-vector (diagonal of G_FF).
#' @keywords internal
.glmm_G_FF <- function(sigma, alpha, X, group, beta) {
  J   <- length(alpha)
  eta <- as.vector(alpha[group] + X %*% beta)
  p   <- plogis(eta)

  # Likelihood contribution: sum p*(1-p) per group
  lik <- as.vector(tapply(p * (1 - p), group, sum))

  rep(1 / sigma^2, J) + lik
}

#' Off-diagonal block G_BF at one parameter point
#'
#' G_BF\[j, 1\] = -1/sigma^2              (base direction: mu)
#' G_BF\[j, 2\] = -2*(alpha_j - mu)/sigma^3   (base direction: sigma)
#'
#' @param sigma  Scalar.
#' @param mu     Scalar.
#' @param alpha  J-vector.
#' @return J x 2 matrix.
#' @keywords internal
.glmm_G_BF <- function(sigma, mu, alpha) {
  J <- length(alpha)
  cbind(
    rep(-1 / sigma^2, J),            # d/dmu
    -2 * (alpha - mu) / sigma^3      # d/dsigma
  )
}

#' Connection form A at one parameter point
#'
#' A = -G_FF^{-1} G_BF^T   (J x 2)
#'
#' The horizontal lift of a base motion dtheta = (dmu, dsigma) is
#' dalpha = A %*% dtheta.
#'
#' A\[j, 1\] =  1 / (sigma^2 * G_FF\[j\])
#' A\[j, 2\] =  2*(alpha_j - mu) / (sigma^3 * G_FF\[j\])
#'
#' @param G_FF_diag J-vector returned by [.glmm_G_FF()].
#' @param G_BF      J x 2 matrix returned by [.glmm_G_BF()].
#' @return J x 2 matrix.
#' @keywords internal
.glmm_connection <- function(G_FF_diag, G_BF) {
  -G_BF / G_FF_diag   # broadcasts column-wise: each row of G_BF / scalar G_FF[j]
}

# ── Linearised curvature ──────────────────────────────────────────────────────

#' Linearised curvature of the GLMM connection
#'
#' Computes the base-derivative part of the Ehresmann curvature 2-form, with
#' the fiber coordinate \eqn{\alpha} (equivalently \eqn{G_{FF}}) held fixed:
#'
#'   \deqn{F_j^{\text{lin}} = \partial_\mu A_{j,\sigma} - \partial_\sigma A_{j,\mu}
#'         = -2 / (\sigma^5 G_{FF,j}^2)}
#'
#' **This is NOT the full Ehresmann curvature.**  The full curvature of
#' \eqn{A = -G_{FF}^{-1}G_{BF}} is identically zero for this model class:
#' the fiber-derivative (vertical) terms \eqn{A_\mu \partial_\alpha A_\sigma -
#' A_\sigma \partial_\alpha A_\mu = +2/(\sigma^5 G_{FF}^2)} cancel the base
#' terms exactly, so the connection is flat (see Proposition `prop:flat` in the
#' companion paper and `data-raw/verify_flat_connection.R`).
#'
#' This function returns the linearised quantity because that is what
#' [synthetic_holonomy_loop()] integrates (the fiber is frozen at the loop
#' centre), and it equals \eqn{2\pi_j^2/\sigma} in absolute value, making
#' it useful as a proxy for where the prior fraction \eqn{\pi_j} is large.
#'
#' @param G_FF_diag J-vector returned by [.glmm_G_FF()].
#' @param sigma     Scalar; group-level SD.
#' @return J-vector.  All values are negative (\eqn{F_j < 0}).
#' @keywords internal
.glmm_curvature_linearised <- function(G_FF_diag, sigma) {
  -2 / (sigma^5 * G_FF_diag^2)
}

# Back-compat alias — prefer .glmm_curvature_linearised in new code.
.glmm_curvature <- .glmm_curvature_linearised

# ── Prior-vs-likelihood information decomposition ─────────────────────────────

#' Fraction of G_FF\[j\] attributable to the prior (vs likelihood)
#'
#' Returns a J-vector in \[0, 1\].  Values near 1 → prior-dominated (sparse
#' data or large sigma); values near 0 → data-dominated (well-identified).
#' @keywords internal
.glmm_prior_fraction <- function(G_FF_diag, sigma) {
  (1 / sigma^2) / G_FF_diag
}
