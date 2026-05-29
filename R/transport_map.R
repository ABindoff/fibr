#' Estimate the holonomy transport map from loop pairs
#'
#' Given loop pairs (start/end indices into the chain) and the fiber draws,
#' fits the linear map \eqn{\hat{H}} such that
#' \deqn{\alpha_{\mathrm{end}} \approx \hat{H}\,\alpha_{\mathrm{start}}}
#' across all detected loops via weighted least squares, then bootstraps the
#' eigenspectrum of \eqn{\hat{H}}.
#'
#' The WLS solution is
#' \deqn{\hat{H} = \bigl(E W S^\top\bigr)\bigl(S W S^\top + \delta I\bigr)^{-1}}
#' where \eqn{S} (\eqn{J \times K}) stacks the start fiber vectors, \eqn{E}
#' (\eqn{J \times K}) stacks the end fiber vectors, and \eqn{W =
#' \mathrm{diag}(w)} is the loop-weight matrix.
#'
#' @param fiber_draws Numeric matrix \[n_iter x J\]. Columns are fiber
#'   (group-level) parameters in the same row order as the base draws that
#'   were passed to [detect_loops()].
#' @param loops Data frame returned by [detect_loops()], with columns
#'   `start`, `end`, `distance`.
#' @param n_bootstrap Number of bootstrap resamples for eigenvalue
#'   uncertainty (default 200).
#' @param ridge Ridge penalty added to the diagonal of the Gram matrix for
#'   numerical stability (default `1e-6`).
#' @param weights `"distance"` (default) weights loops by
#'   \eqn{\exp(-d / \bar{d})} so tighter loops contribute more; `"uniform"`
#'   gives equal weight.
#'
#' @return A list with:
#' \describe{
#'   \item{H}{Estimated \eqn{J \times J} transport matrix.}
#'   \item{eigenvalues}{Complex eigenvalues of `H`, sorted by decreasing modulus.}
#'   \item{boot_eigenvalues}{Complex matrix \[n_bootstrap x J\] of bootstrapped eigenvalues.}
#'   \item{frobenius_dev}{\eqn{\|H - I\|_F}: scalar summary of holonomy magnitude.}
#'   \item{n_loops}{Number of loop pairs used.}
#' }
#'
#' @keywords internal
estimate_transport_map <- function(fiber_draws,
                                   loops,
                                   n_bootstrap = 200L,
                                   ridge       = 1e-6,
                                   weights     = c("distance", "uniform")) {

  weights     <- match.arg(weights)
  fiber_draws <- as.matrix(fiber_draws)
  J <- ncol(fiber_draws)
  K <- nrow(loops)

  if (K < J) {
    warning(sprintf(
      "Only %d loops for J=%d fiber dimensions; H estimate may be unreliable.",
      K, J
    ))
  }

  S <- t(fiber_draws[loops$start, , drop = FALSE])   # J x K
  E <- t(fiber_draws[loops$end,   , drop = FALSE])   # J x K

  w <- .loop_weights(loops$distance, weights)

  H <- .wls_transport(S, E, w, ridge, J)

  evals <- .sorted_eigenvalues(H)

  # Bootstrap by resampling loop pairs
  boot_evals <- matrix(NA_complex_, nrow = n_bootstrap, ncol = J)
  for (b in seq_len(n_bootstrap)) {
    idx_b <- sample.int(K, K, replace = TRUE)
    tryCatch({
      Hb <- .wls_transport(S[, idx_b, drop = FALSE],
                           E[, idx_b, drop = FALSE],
                           w[idx_b], ridge, J)
      boot_evals[b, ] <- .sorted_eigenvalues(Hb)
    }, error = function(e) NULL)
  }

  list(
    H                = H,
    eigenvalues      = evals,
    boot_eigenvalues = boot_evals,
    frobenius_dev    = norm(H - diag(J), type = "F"),
    n_loops          = K
  )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

.loop_weights <- function(distances, type) {
  w <- if (type == "distance") {
    d_mean <- mean(distances)
    if (d_mean == 0) rep(1, length(distances))
    else exp(-distances / d_mean)
  } else {
    rep(1, length(distances))
  }
  # Normalise so sum(w) = K (keeps ridge scale invariant to K)
  w * length(w) / sum(w)
}

# Weighted least squares: H = (E W S')(S W S' + ridge*I)^{-1}
# Uses sqrt-weighting to avoid building a K x K diagonal matrix.
.wls_transport <- function(S, E, w, ridge, J) {
  sw   <- sweep(S, 2L, sqrt(w), "*")   # J x K, columns scaled by sqrt(w_k)
  ew   <- sweep(E, 2L, sqrt(w), "*")

  SWSt <- tcrossprod(sw)               # J x J
  EWSt <- tcrossprod(ew, sw)           # J x J

  diag(SWSt) <- diag(SWSt) + ridge

  EWSt %*% solve(SWSt)
}

.sorted_eigenvalues <- function(H) {
  ev <- eigen(H, only.values = TRUE)$values
  ev[order(Mod(ev), decreasing = TRUE)]
}
