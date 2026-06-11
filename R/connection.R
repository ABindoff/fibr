#' Compute the connection form of a hierarchical MCMC chain
#'
#' @description
#' Estimates the Ehresmann connection on the fiber bundle \eqn{E \to B}
#' induced by the posterior Fisher metric.  At each subsampled chain draw
#' \eqn{t}, the J x K connection matrix
#' \deqn{A(t) = -G_{FF}(t)^{-1} G_{BF}(t)}
#' gives the horizontal lift: if the base moves by \eqn{d\theta}, the fiber
#' must move by \eqn{A(t)\,d\theta} to stay horizontal.
#'
#' Currently supports `method = "analytic_glmm"` for the centred two-level
#' logistic GLMM specified in `inst/stan/glmm_centred.stan`.
#'
#' @section Coordinate convention:
#' All analytic formulas (connection, curvature, prior fraction) are derived
#' with base coordinates \eqn{(\mu, \sigma)} --- **not** \eqn{(\mu, \log\sigma)}.
#' The `sigma` column of `chain` must therefore contain \eqn{\sigma} on its
#' natural (positive) scale, as Stan reports it.  If your chain stores
#' \eqn{\log\sigma} (as the internal samplers in this package do), transform
#' it back with `exp()` before calling this function.  Feeding
#' \eqn{\log\sigma} draws silently produces wrong connection and curvature
#' values: the two parameterisations differ by a chain-rule factor of
#' \eqn{\sigma} in the \eqn{\sigma}-column of \eqn{A} and \eqn{G_{BF}}.
#'
#' @param chain A [`posterior::draws_array`] or named matrix (rows =
#'   iterations, columns = parameters).
#' @param base_vars Character vector of base-space (hyperparameter) column
#'   names; must be `c("mu", "sigma")` for `method = "analytic_glmm"`, with
#'   `sigma` on the natural scale (see *Coordinate convention*).
#' @param fiber_vars Character vector of fiber column names, e.g.
#'   `paste0("alpha[", 1:8, "]")`.
#' @param method `"analytic_glmm"` (default and currently only option).
#' @param stan_data Named list; the original Stan data passed to the model.
#'   Required for `method = "analytic_glmm"`.  Must contain `X`, `y`,
#'   `group` (1-indexed integer vector), and optionally `beta_vars` if the
#'   beta column names differ from `"beta[1]"`, `"beta[2]"`.
#' @param beta_vars Character vector of fixed-effect column names in `chain`
#'   (default `c("beta[1]", "beta[2]")`).
#' @param n_subsample Number of chain draws at which to evaluate A (default
#'   500).  Draws are sampled uniformly without replacement.
#'
#' @return An S3 object of class `fibr_connection`. See
#'   [print.fibr_connection()] and [plot.fibr_connection()].
#'
#' @export
compute_connection <- function(chain,
                               base_vars,
                               fiber_vars,
                               method      = "analytic_glmm",
                               stan_data   = NULL,
                               beta_vars   = c("beta[1]", "beta[2]"),
                               n_subsample = 500L) {

  method <- match.arg(method, "analytic_glmm")

  if (method == "analytic_glmm" && is.null(stan_data))
    stop('stan_data is required for method = "analytic_glmm".')

  # ── Coerce chain to a single pooled matrix ──────────────────────────────────
  chains   <- .split_chains(chain)
  full_mat <- do.call(rbind, chains)

  .check_vars(full_mat, base_vars,  "base_vars")
  .check_vars(full_mat, fiber_vars, "fiber_vars")
  .check_vars(full_mat, beta_vars,  "beta_vars")

  # Guard against log-sigma chains: the analytic formulas require sigma on
  # the natural scale. Negative draws indicate a log-scale column.
  if (method == "analytic_glmm" &&
      any(full_mat[, base_vars[2L]] <= 0)) {
    stop(sprintf(
      paste0("Column '%s' contains non-positive values; the analytic GLMM ",
             "connection requires sigma on the natural scale. If your chain ",
             "stores log(sigma), exp() it first."),
      base_vars[2L]
    ))
  }

  n_total <- nrow(full_mat)
  n_sub   <- min(n_subsample, n_total)
  idx     <- sort(sample.int(n_total, n_sub))
  sub_mat <- full_mat[idx, , drop = FALSE]

  J <- length(fiber_vars)
  K <- length(base_vars)   # = 2

  # Output arrays
  A_arr   <- array(NA_real_, dim = c(n_sub, J, K),
                   dimnames = list(NULL, fiber_vars, base_vars))
  G_FF_mat <- matrix(NA_real_, n_sub, J,
                     dimnames = list(NULL, fiber_vars))
  G_BF_arr <- array(NA_real_, dim = c(n_sub, J, K),
                    dimnames = list(NULL, fiber_vars, base_vars))
  curv_mat <- matrix(NA_real_, n_sub, J,
                     dimnames = list(NULL, fiber_vars))
  pf_mat   <- matrix(NA_real_, n_sub, J,   # prior fraction
                     dimnames = list(NULL, fiber_vars))

  # ── Evaluate at each subsampled draw ────────────────────────────────────────
  X     <- stan_data$X
  y_obs <- stan_data$y
  group <- stan_data$group

  for (t in seq_len(n_sub)) {
    row   <- sub_mat[t, ]
    mu_t  <- row[["mu"]]
    sig_t <- row[["sigma"]]
    alp_t <- as.vector(row[fiber_vars])
    bet_t <- as.vector(row[beta_vars])

    G_FF_t <- .glmm_G_FF(sig_t, alp_t, X, group, bet_t)
    G_BF_t <- .glmm_G_BF(sig_t, mu_t,  alp_t)
    A_t    <- .glmm_connection(G_FF_t, G_BF_t)

    G_FF_mat[t, ]   <- G_FF_t
    G_BF_arr[t, , ] <- G_BF_t
    A_arr[t, , ]    <- A_t
    curv_mat[t, ]   <- .glmm_curvature(G_FF_t, sig_t)
    pf_mat[t, ]     <- .glmm_prior_fraction(G_FF_t, sig_t)
  }

  structure(
    list(
      A           = A_arr,
      G_FF        = G_FF_mat,
      G_BF        = G_BF_arr,
      curvature   = curv_mat,
      prior_frac  = pf_mat,
      base_pts    = sub_mat[, base_vars,  drop = FALSE],
      fiber_pts   = sub_mat[, fiber_vars, drop = FALSE],
      full_mat    = full_mat,         # stored for integrate_transport()
      stan_data   = stan_data,
      beta_vars   = beta_vars,
      base_vars   = base_vars,
      fiber_vars  = fiber_vars,
      method      = method,
      n_subsample = n_sub
    ),
    class = "fibr_connection"
  )
}

#' Print a fibr_connection object
#' @param x A `fibr_connection` object.
#' @param ... Ignored.
#' @export
print.fibr_connection <- function(x, ...) {
  J  <- length(x$fiber_vars)
  cat("fibr connection form\n")
  cat("====================\n")
  cat(sprintf("  Method       : %s\n", x$method))
  cat(sprintf("  Base space   : %s\n", paste(x$base_vars,  collapse = ", ")))
  cat(sprintf("  Fiber dim    : J = %d\n", J))
  cat(sprintf("  Draws eval'd : %d\n\n", x$n_subsample))

  # Connection strength: mean A[j, mu] and A[j, sigma] across draws
  A_mu_mean   <- colMeans(x$A[, , 1L, drop = TRUE])
  A_sig_mean  <- colMeans(x$A[, , 2L, drop = TRUE])
  curv_mean   <- colMeans(x$curvature)
  pf_mean     <- colMeans(x$prior_frac)
  sig_mean    <- mean(x$base_pts[, "sigma"])

  df <- data.frame(
    group        = seq_len(J),
    `A[mu]`      = round(A_mu_mean,  4),
    `A[sigma]`   = round(A_sig_mean, 4),
    `curvature`  = round(curv_mean,  5),
    `prior_frac` = round(pf_mean,    3),
    check.names  = FALSE
  )
  cat("Mean connection coefficients and curvature per group:\n")
  print(df, row.names = FALSE)

  cat(sprintf(
    "\nPosterior mean sigma = %.3f\n",
    sig_mean
  ))
  cat(sprintf(
    "Mean |F[j]| = %.5f  (holonomy per unit base area)\n",
    mean(abs(curv_mean))
  ))

  invisible(x)
}

#' Plot a fibr_connection object
#'
#' @description
#' Two panels:
#' - **Left**: Connection strength `A[j, mu]` vs `sigma` across all draws and
#'   groups, with the theoretical curve `1/(sigma^2 * G_FF[j])` overlaid.
#' - **Right**: Curvature `F[j]` vs `sigma`, showing how the holonomy
#'   magnitude varies across the posterior.
#'
#' @param x A `fibr_connection` object.
#' @param ... Ignored.
#' @return A `ggplot` object (printed invisibly).
#' @export
plot.fibr_connection <- function(x, ...) {
  p <- .plot_connection(x)
  print(p)
  invisible(p)
}

# ── Internal plot helper ──────────────────────────────────────────────────────

.plot_connection <- function(x) {
  n  <- x$n_subsample
  J  <- length(x$fiber_vars)

  sigma_vec <- rep(x$base_pts[, "sigma"], times = J)
  group_vec <- rep(seq_len(J), each = n)

  A_mu_vec  <- as.vector(x$A[, , 1L])
  curv_vec  <- as.vector(x$curvature)
  pf_vec    <- as.vector(x$prior_frac)

  df <- data.frame(
    sigma      = sigma_vec,
    A_mu       = A_mu_vec,
    curvature  = curv_vec,
    prior_frac = pf_vec,
    group      = factor(group_vec)
  )

  # Theoretical reference curve for A[j, mu]: 1/(sigma^2 * G_FF[j])
  # = prior_frac / sigma^2 (since prior_frac = (1/sigma^2)/G_FF)
  df$A_mu_theory <- df$prior_frac / df$sigma^2

  sigma_grid <- seq(min(df$sigma) * 0.9, max(df$sigma) * 1.1, length.out = 100)

  p_left <- ggplot2::ggplot(df, ggplot2::aes(x = sigma, y = A_mu,
                                              colour = group)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = A_mu_theory, group = group),
                       linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
    ggplot2::labs(
      title  = "Connection strength",
      x      = expression(sigma),
      y      = expression(A[j*","~mu]),
      colour = "Group j"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")

  p_right <- ggplot2::ggplot(df, ggplot2::aes(x = sigma, y = curvature,
                                               colour = group)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    ggplot2::labs(
      title  = "Curvature",
      x      = expression(sigma),
      y      = expression(F[j] == -2 / (sigma^5 ~ G[FF][j]^2)),
      colour = "Group j"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right")

  # Combine with patchwork if available, else print separately
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_left + p_right +
      patchwork::plot_annotation(
        title    = "GLMM connection form",
        subtitle = sprintf("J = %d  |  %d subsampled draws  |  dashed = theory",
                           J, x$n_subsample)
      )
  } else {
    # Fallback: facet version
    df_long <- rbind(
      transform(df, panel = "A[j, mu]",  y = A_mu),
      transform(df, panel = "F[j]",      y = curvature)
    )
    ggplot2::ggplot(df_long, ggplot2::aes(x = sigma, y = y, colour = group)) +
      ggplot2::geom_point(alpha = 0.2, size = 0.7) +
      ggplot2::facet_wrap(~ panel, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title  = "GLMM connection form",
        x      = expression(sigma),
        colour = "Group"
      )
  }
}
