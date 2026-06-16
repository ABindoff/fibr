# Fisher metric for the centred GLMM.
#
# Parameter ordering (P = J + 4):
#   q = c(mu, log_sigma, alpha[1:J], beta[1:2])
#
# The OBSERVED Fisher information G = -∇² log p(q|y) is not guaranteed to be
# positive definite away from the posterior mode (the centred GLMM has a
# funnel geometry where G_BB has negative eigenvalues in the tails).
#
# We use the DIAGONAL metric G_diag which IS always positive definite:
#   - Captures the key per-parameter step-size adaptation
#   - Handles the sigma funnel (G_diag[ls] grows with sigma²)
#   - Handles heterogeneous alpha identifiability (G_diag[alpha_j] ∝ n_j)
# The off-diagonal blocks (prior coupling, likelihood coupling) are stored
# separately in .glmm_full_metric_blocks() for research use.

#' Diagonal Fisher metric at one parameter point (always positive definite)
#'
#' @param mu,log_sigma,alpha,beta  Scalar / J-vector / 2-vector of parameters.
#' @param stan_data  Named list with X, y, group.
#' @return P-vector of diagonal metric entries (P = J + 4), in order
#'   (G_mu, G_ls, G_alpha\[1:J\], G_beta\[1:2\]).
#' @keywords internal
.glmm_diag_metric <- function(mu, log_sigma, alpha, beta, stan_data) {

  sigma <- exp(log_sigma)
  J     <- length(alpha)
  X     <- stan_data$X
  group <- stan_data$group

  eta  <- as.vector(alpha[group] + X %*% beta)
  p_i  <- plogis(eta)
  w_i  <- p_i * (1 - p_i)
  w_j  <- as.vector(tapply(w_i, group, sum))
  dev_j <- alpha - mu

  # Diagonal entries of G = -∇² log p
  G_mu    <- 1/25 + J / sigma^2
  G_ls    <- sigma + 2 * sum(dev_j^2) / sigma^2   # always > 0
  G_alpha <- 1 / sigma^2 + w_j                    # always > 0
  G_beta  <- 1/4 + colSums(X^2 * w_i)             # always > 0

  c(G_mu, G_ls, G_alpha, G_beta)
}

#' Full Fisher metric blocks (for research and SoftAbs)
#'
#' Returns the raw metric blocks. NOT guaranteed PD outside the mode; the
#' SoftAbs metric function below regularises the full assembled matrix.
#' @keywords internal
.glmm_full_metric <- function(mu, log_sigma, alpha, beta, stan_data) {

  sigma <- exp(log_sigma)
  J     <- length(alpha)
  X     <- stan_data$X
  group <- stan_data$group

  eta   <- as.vector(alpha[group] + X %*% beta)
  p_i   <- plogis(eta)
  w_i   <- p_i * (1 - p_i)
  w_j   <- as.vector(tapply(w_i, group, sum))
  dev_j <- alpha - mu

  G_FF_diag <- 1 / sigma^2 + w_j

  G_BB <- matrix(c(1/25 + J/sigma^2,
                   2*sum(dev_j)/sigma^2,
                   2*sum(dev_j)/sigma^2,
                   sigma + 2*sum(dev_j^2)/sigma^2), 2L, 2L)

  G_BF <- cbind(rep(-1/sigma^2, J), -2*dev_j/sigma^2)

  G_Beta <- crossprod(X * sqrt(w_i)) + diag(1/4, 2L)

  G_ab <- matrix(0, J, 2L)
  for (k in 1:2) G_ab[, k] <- as.vector(tapply(X[,k]*w_i, group, sum))

  list(G_FF_diag = G_FF_diag, G_BB = G_BB, G_BF = G_BF,
       G_Beta = G_Beta, G_ab = G_ab)
}

#' Assemble the full P×P G matrix from metric blocks
#' @keywords internal
.assemble_G <- function(blocks, J) {
  P     <- 2L + J + 2L
  ib    <- 1:2
  if_   <- 3L:(2L + J)
  ibeta <- (3L + J):(4L + J)

  G <- matrix(0, P, P)
  G[ib,    ib]      <- blocks$G_BB
  G[if_,   if_]     <- diag(blocks$G_FF_diag)
  G[ibeta, ibeta]   <- blocks$G_Beta
  G[if_,   ib]      <- blocks$G_BF
  G[ib,    if_]     <- t(blocks$G_BF)
  G[if_,   ibeta]   <- blocks$G_ab
  G[ibeta, if_]     <- t(blocks$G_ab)
  G
}

# ── SoftAbs metric (Betancourt 2013) ─────────────────────────────────────────

#' SoftAbs regulariser: smooth strictly-positive approximation to |lambda|
#'
#' f(lambda, alpha) = |lambda| / tanh(alpha * |lambda|)
#'
#' Properties:
#'   f(0, alpha)    = 1/alpha    (strictly positive floor)
#'   f(lambda, inf) = |lambda|   (exact absolute value for large alpha)
#'
#' @param lambda  Numeric vector of eigenvalues.
#' @param alpha   Smoothness parameter (default 1; larger = sharper).
#' @return Numeric vector, strictly positive.
#' @keywords internal
.softabs_eval <- function(lambda, alpha = 1.0) {
  al <- alpha * abs(lambda)
  ifelse(al < 1e-8, 1 / alpha, abs(lambda) / tanh(al))
}

#' SoftAbs metric decomposition
#'
#' Computes the SoftAbs Riemannian metric from the full Hessian of the
#' log-posterior.  The Hessian H = ∇² log p(q|y) is indefinite away from the
#' posterior mode; SoftAbs maps each eigenvalue lambda_i of H to the strictly
#' positive value f(lambda_i, alpha), giving a globally PD metric that
#' recovers the absolute-value metric for |lambda| >> 1/alpha.
#'
#' @param blocks  Output of [.glmm_full_metric()].
#' @param J       Number of fiber (group) parameters.
#' @param alpha   SoftAbs smoothness (default 1).
#' @return A list with `U` (eigenvector matrix), `lambda_sa` (SoftAbs
#'   eigenvalues), and `log_det` (log determinant of the metric).
#' @keywords internal
.softabs_decomp <- function(blocks, J, alpha = 1.0) {
  G <- .assemble_G(blocks, J)
  H <- -G                                          # Hessian (not neg-Hessian)
  eig <- eigen(H, symmetric = TRUE)
  lambda_sa <- .softabs_eval(eig$values, alpha)    # always > 0
  list(U         = eig$vectors,
       lambda_sa = lambda_sa,
       log_det   = sum(log(lambda_sa)))
}

#' Gradient of the log-posterior as a flat P-vector
#'
#' @inheritParams .glmm_full_metric
#' @return Numeric vector of length P, in order
#'   (grad_mu, grad_log_sigma, grad_alpha\[1:J\], grad_beta\[1:2\]).
#' @keywords internal
.glmm_grad_vec <- function(mu, log_sigma, alpha, beta, stan_data) {
  g <- .glmm_grad_log_post(mu, log_sigma, alpha, beta, stan_data)
  c(g$grad_mu, g$grad_log_sigma, g$grad_alpha, g$grad_beta)
}

#' Extract parameter components from a flat state vector
#' @keywords internal
.unpack <- function(q, J) {
  list(
    mu        = q[1L],
    log_sigma = q[2L],
    alpha     = q[3L:(2L + J)],
    beta      = q[(3L + J):(4L + J)]
  )
}
