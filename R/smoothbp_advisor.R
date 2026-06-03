#' Parameterisation advisor for smoothbp changepoint random effects
#'
#' @description
#' For each random effect on a changepoint location (`omega_k_g` = per-group
#' deviation from the population changepoint at breakpoint k), computes the
#' **Fisher information decomposition** at a subsample of posterior draws.
#'
#' The key quantity is `prior_frac`:
#'
#' \deqn{\text{prior\_frac}_{k,g} =
#'   \frac{G_\text{prior}}{G_\text{prior} + G_\text{lik}}}
#'
#' where
#' \deqn{G_\text{prior} = \sigma_{\text{re},k}^{-2}, \quad
#'   G_\text{lik} = \sigma^{-2}
#'     \sum_{i:\,\text{group}=g}
#'       \!\!\!\left(\frac{\partial\mu_i}{\partial\omega_{k,i}}\right)^{\!2}}
#'
#' \itemize{
#'   \item `prior_frac` \eqn{\to 1}: prior dominates — group changepoints are
#'     poorly identified from data relative to the shrinkage prior.  The sampler
#'     is in the funnel regime and non-centred reparameterisation would help.
#'   \item `prior_frac` \eqn{\to 0}: likelihood dominates — centred
#'     parameterisation is efficient and mixing should be adequate.
#'   \item Mixed: flag individual groups for attention.
#' }
#'
#' **What to do with the results:**
#'
#' When `prior_frac` is high, re-fit with `reparameterise = "omega"`:
#'
#' ```r
#' fit_nc <- smoothbp(..., reparameterise = "omega")
#' fit_nc_ss <- smoothbp_ss(..., reparameterise = "omega")
#' ```
#'
#' This activates the non-centred HMC parameterisation in the Rust sampler:
#' `z[j] = beta_om[j] / sigma_re_om[k]` is sampled with an `N(0,1)` prior,
#' and `beta_om[j] = z[j] * sigma_re_om[k]` is reconstructed automatically.
#' The `sigma_re_om` Gibbs step and the stored draws are unchanged — output
#' is in the original (centred) parameterisation for easy interpretation.
#'
#' Additional options if `reparameterise = "omega"` is insufficient:
#' \enumerate{
#'   \item Increase `warmup` and `iter`.
#'   \item Check `fit$n_divergent` — many divergences confirm a remaining funnel.
#'   \item Fix the changepoint for poorly-identified groups using
#'     `omega = list(fixed(value))`.
#' }
#'
#' The `prior_frac` values quantify the severity: values above 0.8 indicate a
#' serious funnel; 0.6–0.8 suggests moderate difficulty worth addressing.
#'
#' The gradient is computed analytically from the sigmoid smooth-transition
#' likelihood:
#'
#' \deqn{\frac{\partial\mu_i}{\partial\omega_{k,i}} =
#'   -\Bigl[\delta_{k,i}\,\sigma_{ki}\bigl(1 + d_{ki}\,\rho_{ki}(1-\sigma_{ki})\bigr)
#'   + b_{1,i}\,\mathbf{1}_{k=1}\Bigr]}
#'
#' where \eqn{d_{ki} = \tau_i - \omega_{k,i}} and
#' \eqn{\sigma_{ki} = \text{logistic}(d_{ki}\,\rho_{ki})}.
#'
#' @param fit        A `smoothbp_fit` from `smoothbp::smoothbp()` or
#'   `smoothbp::smoothbp_ss()`, with at least one omega random effect
#'   (`omega = list(~ 1 + (1 | group))`).
#' @param n_draws    Number of posterior draws to evaluate the metric at
#'   (default 200; subsampled uniformly).
#' @param threshold_nc Prior fraction above which non-centred is recommended
#'   (default 0.60).
#' @param threshold_c  Prior fraction below which centred is safe
#'   (default 0.40).
#'
#' @return An S3 object of class `fibr_smoothbp_advice`.  Contains one list
#'   element per breakpoint that has omega random effects, each with
#'   `prior_frac_mean`, `prior_frac_q05`, `prior_frac_q95`,
#'   `recommendation`, and `delta` (one entry per group).  `delta[j]` is the
#'   recommended per-group NC mixing fraction
#'   \eqn{\delta_j = \sqrt{\text{prior\_frac}_j}}, which can be passed
#'   directly to `smoothbp()` via `re_fraction = as_smoothbp_re_fraction(advice)`.
#'   Print and plot methods included.
#'
#' @seealso [as_smoothbp_re_fraction()]
#' @export
smoothbp_advisor <- function(fit,
                              n_draws      = 200L,
                              threshold_nc = 0.60,
                              threshold_c  = 0.40) {

  if (!inherits(fit, "smoothbp_fit"))
    stop("`fit` must be a smoothbp_fit object from smoothbp::smoothbp().")

  all_vars   <- posterior::variables(fit$draws)
  re_om_vars <- grep("^omega[0-9]+_re_", all_vars, value = TRUE)

  if (length(re_om_vars) == 0L) {
    message("fibr: no random effects on omega found in this fit.\n",
            "  Refit with e.g. omega = list(~ 1 + (1 | group)).")
    return(invisible(NULL))
  }

  k_vals <- sort(unique(as.integer(
    sub("^omega([0-9]+)_re_.*", "\\1", re_om_vars)
  )))

  # ── Subsample draws ───────────────────────────────────────────────────────
  draws_mat <- as.matrix(posterior::as_draws_matrix(fit$draws))
  n_sub     <- min(as.integer(n_draws), nrow(draws_mat))
  sub_draws <- draws_mat[sort(sample.int(nrow(draws_mat), n_sub)), , drop = FALSE]

  tau <- as.double(fit$data[[fit$time]])
  N   <- length(tau)
  dm  <- fit$dm

  # ── Compute prior_frac per breakpoint ────────────────────────────────────
  bp_results <- Filter(Negate(is.null), lapply(k_vals, function(k) {

    X_om_k   <- dm$X_om[[k]]
    re_mask  <- attr(X_om_k, "re_mask")
    if (is.null(re_mask)) re_mask <- rep(0L, ncol(X_om_k))

    re_cols <- which(re_mask == 1L)
    if (length(re_cols) == 0L) return(NULL)

    col_om  <- dm$col_names_om[[k]]
    col_rho <- dm$col_names_rho[[k]]
    col_del <- dm$col_names_deltas[[k]]
    col_b1  <- dm$col_names_b1

    re_col_names <- col_om[re_cols]
    re_var_names <- paste0("omega", k, "_", re_col_names)
    sigma_re_var <- paste0("sigma_re_omega", k)

    pf_mat <- matrix(NA_real_, n_sub, length(re_cols),
                     dimnames = list(NULL, re_var_names))

    # ── Fully vectorised: extract all draws at once ────────────────────────
    .cv <- function(prefix, cols) paste0(prefix, cols)

    sigma_vec    <- sub_draws[, "sigma",      drop = TRUE]
    sigma_re_vec <- sub_draws[, sigma_re_var, drop = TRUE]

    # Beta draw matrices: n_sub × p_k — use matrix multiplication on all draws
    .bdraw <- function(prefix, cols, X, gamma_prefix = NULL) {
      B <- sub_draws[, .cv(prefix, cols), drop = FALSE]  # n_sub × p
      if (!is.null(gamma_prefix)) {
        G <- .sbp_gamma_mat(sub_draws, gamma_prefix, cols)
        if (!is.null(G)) B <- B * (G > 0.5)
      }
      B %*% t(X)   # n_sub × N
    }

    omega_all <- .bdraw(paste0("omega",  k, "_"), col_om,  X_om_k)
    rho_all   <- .bdraw(paste0("rho",    k, "_"), col_rho, dm$X_rho[[k]])
    delta_all <- .bdraw(paste0("delta",  k, "_"), col_del, dm$X_deltas[[k]],
                        gamma_prefix = paste0("gamma_delta", k, "_"))
    b1_all    <- .bdraw("b1_",                    col_b1,  dm$X_b1,
                        gamma_prefix = "gamma_b1_")

    # Vectorised dmu/domega: n_sub × N
    tau_mat <- matrix(tau, n_sub, N, byrow = TRUE)
    d_all   <- tau_mat - omega_all
    s_all   <- plogis(d_all * rho_all)
    dmu_all <- -(delta_all * s_all * (1 + d_all * rho_all * (1 - s_all)))
    if (k == 1L) dmu_all <- dmu_all - b1_all  # b1 term for first breakpoint

    # G_prior and G_lik for each RE group — all draws at once
    G_prior_vec <- 1 / sigma_re_vec^2     # n_sub

    valid <- is.finite(sigma_vec) & sigma_vec > 0 &
             is.finite(sigma_re_vec) & sigma_re_vec > 0

    for (ji in seq_along(re_cols)) {
      grp_mask <- X_om_k[, re_cols[ji]] == 1L           # N-vector
      # rowSums over group observations: n_sub × n_grp -> n_sub
      G_lik_vec <- rowSums(dmu_all[, grp_mask, drop = FALSE]^2,
                           na.rm = TRUE) / sigma_vec^2  # n_sub
      pf <- G_prior_vec / (G_prior_vec + G_lik_vec)
      pf[!valid] <- NA_real_
      pf_mat[, ji] <- pf
    }

    pf_mean <- colMeans(pf_mat, na.rm = TRUE)
    pf_q05  <- apply(pf_mat, 2L, quantile, 0.05, na.rm = TRUE)
    pf_q95  <- apply(pf_mat, 2L, quantile, 0.95, na.rm = TRUE)
    rec     <- ifelse(pf_mean > threshold_nc, "non-centred",
               ifelse(pf_mean < threshold_c,  "centred (OK)", "borderline"))

    # delta_j = sqrt(prior_frac_j): optimal NC mixing fraction from Gaussian theory.
    # delta = 1 -> fully NC (funnel regime); delta = 0 -> fully centred (data-rich).
    delta <- sqrt(pmax(0, pmin(1, pf_mean)))

    list(breakpoint      = k,
         re_vars         = re_var_names,
         prior_frac_mean = pf_mean,
         prior_frac_q05  = pf_q05,
         prior_frac_q95  = pf_q95,
         recommendation  = rec,
         delta           = delta,
         pf_mat          = pf_mat)
  }))

  structure(list(breakpoints  = bp_results,
                 threshold_nc = threshold_nc,
                 threshold_c  = threshold_c,
                 n_draws      = n_sub),
            class = "fibr_smoothbp_advice")
}

# ── S3 methods ────────────────────────────────────────────────────────────────

#' @export
print.fibr_smoothbp_advice <- function(x, digits = 3L, ...) {

  cat("fibr smoothbp parameterisation advisor\n")
  cat("======================================\n")
  cat(sprintf("Draws evaluated : %d\n", x$n_draws))
  cat(sprintf("Non-centred if prior_frac > %.2f; centred OK if < %.2f\n\n",
              x$threshold_nc, x$threshold_c))

  for (bp in x$breakpoints) {
    k <- bp$breakpoint
    cat(sprintf("Breakpoint %d  (omega%d random effects)\n", k, k))
    cat(rep("-", 50L), "\n", sep = "")

    df <- data.frame(
      group        = bp$re_vars,
      `prior_frac` = round(bp$prior_frac_mean, digits),
      `q05`        = round(bp$prior_frac_q05,  digits),
      `q95`        = round(bp$prior_frac_q95,  digits),
      delta        = round(bp$delta,            digits),
      recommendation = bp$recommendation,
      check.names  = FALSE
    )
    # Shorten group names for display
    df$group <- sub(paste0("omega", k, "_re_"), "", df$group)

    print(df, row.names = FALSE)

    n_nc  <- sum(bp$recommendation == "non-centred")
    n_ok  <- sum(bp$recommendation == "centred (OK)")
    n_tot <- length(bp$recommendation)

    if (n_ok == n_tot) {
      cat("\nSampling geometry OK — centred parameterisation is efficient.\n\n")
    } else {
      severity <- if (n_nc == n_tot) "all groups" else
                  sprintf("%d of %d groups", n_nc, n_tot)
      cat(sprintf(
        "\nFunnel geometry detected in %s at breakpoint %d.\n", severity, k))
      cat("Options (in order of impact):\n")
      cat("  1. Re-fit with reparameterise = \"omega\" (non-centred HMC for omega RE):\n")
      cat("       smoothbp(..., reparameterise = \"omega\")\n")
      cat("  2. Increase warmup/iter if reparameterising does not fully resolve it.\n")
      cat("  3. Check fit$n_divergent to confirm severity.\n")
      cat("  4. Fix changepoints for badly-identified groups: omega = list(fixed(value)).\n\n")
    }
  }

  if (any(sapply(x$breakpoints, function(bp) any(bp$delta > 0 & bp$delta < 1)))) {
    cat("Per-group partial-NC fractions (delta) are available.\n")
    cat("Pass to smoothbp via:\n")
    cat("  fit_nc <- smoothbp(..., re_fraction = as_smoothbp_re_fraction(advice))\n\n")
  }

  invisible(x)
}

#' @export
plot.fibr_smoothbp_advice <- function(x, ...) {
  p <- .plot_sbp_advice(x)
  print(p)
  invisible(p)
}

#' Convert smoothbp advisor output to re_fraction format
#'
#' Extracts the per-group, per-breakpoint NC mixing fractions (`delta`) from a
#' `fibr_smoothbp_advice` object and returns them in the list format expected
#' by the `re_fraction` argument of [smoothbp::smoothbp()] and
#' [smoothbp::smoothbp_ss()].
#'
#' The returned list has one element per breakpoint that has omega random
#' effects.  Each element is a named numeric vector of `delta` values in
#' \eqn{[0, 1]}, one per group.  Groups with `delta > threshold` will use the
#' non-centred parameterisation; others use centred.  The threshold is applied
#' inside the smoothbp sampler (default 0.5).
#'
#' @param advice A `fibr_smoothbp_advice` object returned by
#'   [smoothbp_advisor()].
#' @param ... Unused.
#'
#' @return A named list of numeric vectors, compatible with
#'   `smoothbp(..., re_fraction = ...)`.
#'
#' @examples
#' \dontrun{
#' advice  <- smoothbp_advisor(fit)
#' re_frac <- as_smoothbp_re_fraction(advice)
#' fit_nc  <- smoothbp(..., re_fraction = re_frac)
#' }
#'
#' @export
as_smoothbp_re_fraction <- function(advice, ...) UseMethod("as_smoothbp_re_fraction")

#' @rdname as_smoothbp_re_fraction
#' @export
as_smoothbp_re_fraction.fibr_smoothbp_advice <- function(advice, ...) {
  out <- lapply(advice$breakpoints, function(bp) {
    d <- bp$delta
    names(d) <- sub(paste0("omega", bp$breakpoint, "_re_"), "", bp$re_vars)
    d
  })
  names(out) <- paste0("breakpoint_", sapply(advice$breakpoints, `[[`, "breakpoint"))
  out
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Reconstruct per-observation linear predictor from a single draw row
# (named numeric vector). Used directly in tests.
.sbp_lp <- function(row, prefix, X, col_names, gamma_prefix = NULL) {
  param_names <- paste0(prefix, col_names)
  beta        <- as.numeric(row[param_names])
  beta[is.na(beta)] <- 0

  if (!is.null(gamma_prefix)) {
    g_names <- paste0(gamma_prefix, col_names)
    if (all(g_names %in% names(row))) {
      gamma <- as.numeric(row[g_names])
      gamma[is.na(gamma)] <- 1
      beta <- beta * as.numeric(gamma > 0.5)
    }
  }

  as.vector(X %*% beta)
}

# Extract gamma (spike-and-slab) draw matrix if the gamma variables exist.
# Returns NULL if not a spike-and-slab fit.
.sbp_gamma_mat <- function(draws_mat, gamma_prefix, col_names) {
  g_names <- paste0(gamma_prefix, col_names)
  if (!all(g_names %in% colnames(draws_mat))) return(NULL)
  draws_mat[, g_names, drop = FALSE]
}

# Gradient dmu_i / domega_k_i from the smoothbp sigmoid likelihood.
# From sampler_re.rs: dmu_dom = -(b*s + d*r*s*(1-s)*b) - b1 if is_om1
.sbp_dmu_domega <- function(tau, omega_ki, rho_ki, delta_ki, b1_i,
                             is_first = TRUE) {
  d   <- tau - omega_ki
  s   <- plogis(d * rho_ki)
  dmu <- -(delta_ki * s * (1 + d * rho_ki * (1 - s)))
  if (is_first) dmu <- dmu - b1_i
  dmu
}

# ggplot2 bar chart of prior_frac with coloured recommendation
.plot_sbp_advice <- function(x) {
  rows <- lapply(x$breakpoints, function(bp) {
    k <- bp$breakpoint
    data.frame(
      group       = sub(paste0("omega", k, "_re_"), "", bp$re_vars),
      breakpoint  = paste0("omega", k),
      mean        = bp$prior_frac_mean,
      q05         = bp$prior_frac_q05,
      q95         = bp$prior_frac_q95,
      rec         = bp$recommendation
    )
  })
  df <- do.call(rbind, rows)
  df$label <- paste0(df$breakpoint, "\n", df$group)

  ggplot2::ggplot(df, ggplot2::aes(x = label, y = mean, fill = rec)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = q05, ymax = q95),
                           width = 0.25, colour = "grey30") +
    ggplot2::geom_hline(yintercept = x$threshold_nc, linetype = "dashed",
                        colour = "firebrick",  linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = x$threshold_c, linetype = "dashed",
                        colour = "steelblue",  linewidth = 0.6) +
    ggplot2::scale_fill_manual(
      values = c("non-centred"  = "#B22222",
                 "centred (OK)" = "#2E8B57",
                 "borderline"   = "#FF8C00"),
      drop   = FALSE
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title    = "smoothbp parameterisation advisor",
      subtitle = paste0("Red line = non-centred threshold (", x$threshold_nc,
                        ");  blue line = centred-OK threshold (", x$threshold_c, ")"),
      x        = "Breakpoint × group",
      y        = "prior_frac  (prior / total Fisher information)",
      fill     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(size = 9))
}
