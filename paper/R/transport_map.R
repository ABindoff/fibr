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
#' @param structure `"diagonal"` (default) estimates a per-group scalar
#'   transport \eqn{\alpha_{\mathrm{end},j} \approx h_j\,\alpha_{\mathrm{start},j}},
#'   i.e. \eqn{H = \mathrm{diag}(h_1, \ldots, h_J)}.  This is the correct
#'   model class whenever the fiber metric block \eqn{G_{FF}} is diagonal and
#'   the connection decouples across groups (true for the two-level GLMM,
#'   where parallel transport is per-group contraction and the structure
#'   group is abelian; genuine rotation is impossible).  `"full"` estimates
#'   an unrestricted \eqn{J \times J} matrix; with few loop pairs this
#'   manufactures spurious off-diagonal structure and complex eigenvalues,
#'   so it should be reserved for models with genuinely coupled fibers.
#'
#' @return A list with:
#' \describe{
#'   \item{H}{Estimated \eqn{J \times J} transport matrix (diagonal when
#'     `structure = "diagonal"`).}
#'   \item{eigenvalues}{Eigenvalues of `H`, sorted by decreasing modulus.
#'     Real for `structure = "diagonal"` (returned as complex for a stable
#'     interface); these are the per-group contraction factors \eqn{h_j}.}
#'   \item{boot_eigenvalues}{Complex matrix \[n_bootstrap x J\] of bootstrapped eigenvalues.}
#'   \item{frobenius_dev}{\eqn{\|H - I\|_F}: scalar summary of holonomy magnitude.}
#'   \item{n_loops}{Number of loop pairs used.}
#'   \item{structure}{The structure used (`"diagonal"` or `"full"`).}
#' }
#'
#' @keywords internal
estimate_transport_map <- function(fiber_draws,
                                   loops,
                                   n_bootstrap = 200L,
                                   ridge       = 1e-6,
                                   weights     = c("distance", "uniform"),
                                   structure   = c("diagonal", "full")) {

  weights     <- match.arg(weights)
  structure   <- match.arg(structure)
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

  fit_H <- if (structure == "diagonal") {
    function(S, E, w) .wls_transport_diag(S, E, w, ridge)
  } else {
    function(S, E, w) .wls_transport(S, E, w, ridge, J)
  }

  H <- fit_H(S, E, w)

  # For the diagonal estimator, eigenvalues ARE the per-group contraction
  # factors h_j; keep them in group order (do not sort by modulus) so that
  # bootstrap draws stay aligned to groups.
  extract_evals <- if (structure == "diagonal") {
    function(H) as.complex(diag(H))
  } else {
    .sorted_eigenvalues
  }

  evals <- extract_evals(H)

  # Bootstrap by resampling loop pairs
  boot_evals <- matrix(NA_complex_, nrow = n_bootstrap, ncol = J)
  for (b in seq_len(n_bootstrap)) {
    idx_b <- sample.int(K, K, replace = TRUE)
    tryCatch({
      Hb <- fit_H(S[, idx_b, drop = FALSE],
                  E[, idx_b, drop = FALSE],
                  w[idx_b])
      boot_evals[b, ] <- extract_evals(Hb)
    }, error = function(e) NULL)
  }

  list(
    H                = H,
    eigenvalues      = evals,
    boot_eigenvalues = boot_evals,
    frobenius_dev    = norm(H - diag(J), type = "F"),
    n_loops          = K,
    structure        = structure
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

# Per-group weighted least squares: h_j = sum_k w_k E_jk S_jk / (sum_k w_k S_jk^2 + ridge)
# The correct estimator when the connection decouples across groups (G_FF
# diagonal): each fiber coordinate is transported independently, so H is
# diagonal with real entries.  J separate scalar regressions; far fewer
# parameters than the full J x J fit, hence stable with few loop pairs.
.wls_transport_diag <- function(S, E, w, ridge) {
  num <- as.vector((E * S) %*% w)        # J-vector: sum_k w_k E_jk S_jk
  den <- as.vector((S * S) %*% w) + ridge
  diag(num / den, nrow = nrow(S))
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
