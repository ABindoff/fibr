# Analytic log-posterior and gradient for the centred two-level GLMM.
#
# Parameterisation: (mu, log_sigma, alpha[1:J], beta[1:2])
# Using log_sigma internally removes the sigma > 0 constraint.
# All functions operate on plain R vectors; no Stan or posterior objects.

#' Log-posterior for the centred GLMM (internal parameterisation)
#'
#' @param mu        Scalar hyperparameter mean.
#' @param log_sigma Log of the hyperparameter SD (unconstrained).
#' @param alpha     J-vector of group intercepts.
#' @param beta      2-vector of fixed-effect coefficients.
#' @param stan_data List with elements X (N x 2), y (N), group (N, 1-indexed).
#' @return Scalar log-posterior value.
#' @keywords internal
.glmm_log_post <- function(mu, log_sigma, alpha, beta, stan_data) {
  sigma <- exp(log_sigma)
  eta   <- as.vector(alpha[stan_data$group] + stan_data$X %*% beta)

  # Priors
  lp <- dnorm(mu, 0, 5, log = TRUE)
  lp <- lp + dexp(sigma, rate = 1, log = TRUE) + log_sigma   # Jacobian
  lp <- lp + sum(dnorm(alpha, mu, sigma, log = TRUE))
  lp <- lp + sum(dnorm(beta, 0, 2, log = TRUE))

  # Likelihood (numerically stable log-sum via log1p)
  lp <- lp + sum(stan_data$y * eta - log1p(exp(eta)))

  lp
}

#' Gradient of the log-posterior w.r.t. all parameters
#'
#' @inheritParams .glmm_log_post
#' @return Named list: `grad_mu`, `grad_log_sigma`, `grad_alpha` (J-vector),
#'   `grad_beta` (2-vector).
#' @keywords internal
.glmm_grad_log_post <- function(mu, log_sigma, alpha, beta, stan_data) {
  sigma <- exp(log_sigma)
  J     <- length(alpha)
  N     <- length(stan_data$y)
  eta   <- as.vector(alpha[stan_data$group] + stan_data$X %*% beta)
  p     <- plogis(eta)      # N-vector of event probabilities
  resid <- stan_data$y - p  # N-vector of residuals

  # Gradient wrt mu
  g_mu <- -mu / 25 + sum(alpha - mu) / sigma^2

  # Gradient wrt log_sigma (chain rule: ∂/∂log_sigma = sigma * ∂/∂sigma)
  # ∂ lp / ∂ sigma = -1 (exp prior) + sum[-(1/sigma) + (alpha_j-mu)^2/sigma^3]
  # ∂ lp / ∂ log_sigma = sigma * ∂ lp / ∂ sigma + 1 (Jacobian: d/d(ls) of log_sigma term)
  g_log_sigma <- -sigma + 1 + sum((alpha - mu)^2 / sigma^2 - 1)

  # Gradient wrt alpha_j
  lik_contrib_alpha <- as.vector(tapply(resid, stan_data$group, sum))
  g_alpha <- (mu - alpha) / sigma^2 + lik_contrib_alpha

  # Gradient wrt beta
  g_beta <- as.vector(t(stan_data$X) %*% resid - beta / 4)

  list(
    grad_mu        = g_mu,
    grad_log_sigma = g_log_sigma,
    grad_alpha     = g_alpha,
    grad_beta      = g_beta
  )
}

#' Wrap a named parameter list into a state vector and back
#' @keywords internal
.state_to_list <- function(state, J) {
  list(
    mu        = state[1L],
    log_sigma = state[2L],
    alpha     = state[3L:(2L + J)],
    beta      = state[(3L + J):(4L + J)]
  )
}

#' @keywords internal
.list_to_state <- function(lst) {
  c(lst$mu, lst$log_sigma, lst$alpha, lst$beta)
}
