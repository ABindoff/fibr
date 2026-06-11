#' Print a fibr_holonomy object
#'
#' @param x A `fibr_holonomy` object.
#' @param ... Ignored.
#' @export
print.fibr_holonomy <- function(x, ...) {
  cat("fibr holonomy diagnostic\n")
  cat("========================\n")
  cat(sprintf("  Base space   : %s\n", paste(x$base_vars, collapse = ", ")))
  cat(sprintf("  Fiber dim    : J = %d\n", length(x$fiber_vars)))
  cat(sprintf("  Residualised : %s\n", if (isTRUE(x$residualized)) "yes" else "no"))
  cat(sprintf("  Structure    : %s\n", x$structure %||% "full"))
  cat(sprintf("  Loop pairs   : %d\n", x$n_loops))
  cat(sprintf("  ||H - I||_F  : %.4f\n\n", x$frobenius_dev))

  evals <- x$eigenvalues
  if (identical(x$structure, "diagonal")) {
    # Group-aligned per-group contraction factors (real)
    df <- data.frame(
      group       = seq_along(evals),
      h_j         = round(Re(evals), 4),
      check.names = FALSE
    )
    if (!is.null(x$boot_eigenvalues)) {
      ci <- apply(Re(x$boot_eigenvalues), 2L, quantile,
                  probs = c(0.05, 0.95), na.rm = TRUE)
      df$`5%`  <- round(ci[1L, ], 4)
      df$`95%` <- round(ci[2L, ], 4)
    }
    cat("Per-group transport (contraction) factors h_j:\n")
    print(df, row.names = FALSE)
  } else {
    df <- data.frame(
      ` ` = seq_along(evals),
      `|lambda|` = round(Mod(evals),  4),
      `Arg (deg)` = round(Arg(evals) * 180 / pi, 2),
      Re           = round(Re(evals), 4),
      Im           = round(Im(evals), 4),
      check.names  = FALSE
    )
    cat("Eigenvalues of H (sorted by |lambda|):\n")
    print(df, row.names = FALSE)
  }

  if (!is.null(x$boot_eigenvalues)) {
    # Bootstrap 90% CI on ||H_boot - I||_F approximated from eigenvalues
    boot_dev <- apply(x$boot_eigenvalues, 1L, function(ev) {
      sqrt(sum((Mod(ev) - 1)^2))
    })
    ci <- quantile(boot_dev, c(0.05, 0.95), na.rm = TRUE)
    cat(sprintf(
      "\n||H - I||_F  bootstrap 90%% CI: [%.4f, %.4f]\n",
      ci[[1L]], ci[[2L]]
    ))
  }

  invisible(x)
}

#' Plot a fibr_holonomy object
#'
#' @description
#' Three plot types are available:
#' - `"contraction"` (default for diagonal-structure fits): per-group
#'   transport factors \eqn{h_j} with bootstrap intervals, against group
#'   index.  Group identity is preserved, so \eqn{h_j} can be compared
#'   directly with per-group analytic quantities such as the prior fraction
#'   \eqn{\pi_j} or the integrated transport prediction.
#' - `"eigenspectrum"` (default for full-structure fits): eigenvalues of
#'   \eqn{\hat{H}} in the complex plane with the unit circle as reference and
#'   bootstrap clouds.
#' - `"base_loops"`: scatter of detected loop pairs (start vs end iteration),
#'   coloured by base-space distance.
#'
#' @param x A `fibr_holonomy` object.
#' @param type `"contraction"`, `"eigenspectrum"`, or `"base_loops"`.
#'   `NULL` (default) chooses based on the fit structure.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (printed invisibly; assign to capture).
#' @export
plot.fibr_holonomy <- function(x, type = NULL, ...) {
  if (is.null(type)) {
    type <- if (identical(x$structure, "diagonal")) "contraction"
            else "eigenspectrum"
  }
  type <- match.arg(type, c("contraction", "eigenspectrum", "base_loops"))
  p <- switch(type,
    contraction   = .plot_contraction(x),
    eigenspectrum = .plot_eigenspectrum(x),
    base_loops    = .plot_base_loops(x)
  )
  print(p)
  invisible(p)
}

# ── Internal plot helpers ─────────────────────────────────────────────────────

.plot_contraction <- function(x) {
  if (!identical(x$structure, "diagonal")) {
    stop("type = 'contraction' requires a diagonal-structure fit; ",
         "re-run holonomy_diagnostic() with structure = 'diagonal'.")
  }
  J  <- length(x$eigenvalues)
  df <- data.frame(group = factor(seq_len(J)), h = Re(x$eigenvalues))

  p <- ggplot(df, aes(x = group, y = h))

  if (!is.null(x$boot_eigenvalues)) {
    ci <- apply(Re(x$boot_eigenvalues), 2L, quantile,
                probs = c(0.05, 0.95), na.rm = TRUE)
    df$lo <- ci[1L, ]
    df$hi <- ci[2L, ]
    p <- ggplot(df, aes(x = group, y = h)) +
      geom_errorbar(aes(ymin = lo, ymax = hi),
                    width = 0.2, colour = "steelblue")
  }

  p +
    geom_hline(yintercept = c(0, 1), linetype = "dashed", colour = "grey60") +
    geom_point(colour = "firebrick", size = 3) +
    labs(
      title    = "Per-group transport factors",
      subtitle = sprintf("h_j with bootstrap 90%% intervals  |  %d loops",
                         x$n_loops),
      x        = "Group j",
      y        = expression(h[j])
    ) +
    theme_minimal(base_size = 12)
}

.plot_eigenspectrum <- function(x) {
  theta  <- seq(0, 2 * pi, length.out = 300)
  circle <- data.frame(x = cos(theta), y = sin(theta))

  evals_df <- data.frame(
    x     = Re(x$eigenvalues),
    y     = Im(x$eigenvalues),
    label = seq_along(x$eigenvalues)
  )

  p <- ggplot() +
    geom_path(data = circle, aes(x = x, y = y),
              colour = "grey60", linetype = "dashed", linewidth = 0.6)

  if (!is.null(x$boot_eigenvalues)) {
    boot_df <- data.frame(
      x = as.vector(Re(x$boot_eigenvalues)),
      y = as.vector(Im(x$boot_eigenvalues))
    )
    boot_df <- boot_df[complete.cases(boot_df), ]
    p <- p + geom_point(data = boot_df, aes(x = x, y = y),
                        colour = "steelblue", alpha = 0.07, size = 0.7)
  }

  p +
    geom_point(data = evals_df, aes(x = x, y = y),
               colour = "firebrick", size = 3.5) +
    geom_text(data = evals_df, aes(x = x, y = y, label = label),
              nudge_y = 0.05, size = 3, colour = "firebrick") +
    # Mark the identity point
    geom_point(data = data.frame(x = 1, y = 0), aes(x = x, y = y),
               shape = 3, size = 5, colour = "black", stroke = 1) +
    coord_equal() +
    labs(
      title    = "Holonomy eigenspectrum",
      subtitle = sprintf("||H − I||_F = %.4f  |  %d loops  |  J = %d",
                         x$frobenius_dev, x$n_loops, length(x$fiber_vars)),
      x        = expression(Re(lambda)),
      y        = expression(Im(lambda))
    ) +
    theme_minimal(base_size = 12)
}

.plot_base_loops <- function(x) {
  loops <- x$loops
  eps   <- attr(loops, "epsilon")
  sub   <- if (!is.null(eps)) sprintf("%d pairs  |  epsilon = %.4f", nrow(loops), eps)
            else               sprintf("%d pairs", nrow(loops))

  ggplot(loops, aes(x = start, y = end, colour = distance)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_colour_viridis_c(option = "magma", direction = -1) +
    labs(
      title    = "Detected loop pairs",
      subtitle = sub,
      x        = "Start iteration",
      y        = "End iteration",
      colour   = "Distance"
    ) +
    theme_minimal(base_size = 12)
}
