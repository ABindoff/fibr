#' Integrate the connection along detected loops
#'
#' @description
#' For each loop pair (start → end), numerically integrates the connection
#' form A along the chain trajectory between those two iterations:
#'
#' \deqn{\alpha_{\text{transported}} = \alpha_{\text{start}} +
#'   \sum_{t=\text{start}}^{\text{end}-1} A(\theta_t,\alpha_t)\,
#'   \Delta\theta_t}
#'
#' and compares the result against the actual chain value
#' \eqn{\alpha_{\text{end}}}.  This separates the geometric prediction
#' (analytical parallel transport) from the statistical noise of MCMC.
#'
#' @param connection A `fibr_connection` object from [compute_connection()].
#' @param loops Data frame of loop pairs (columns `start`, `end`, `distance`)
#'   as returned by [detect_loops()] or stored in a `fibr_holonomy` object.
#' @param chain Named matrix (rows = iterations, pooled across chains) or
#'   `draws_array`.  Must include `base_vars`, `fiber_vars`, and `beta_vars`.
#' @param n_loops Maximum number of loops to integrate (default 300).
#'   Loops are subsampled by distance (tightest first) before integration.
#'
#' @return A data frame with one row per loop and columns:
#' \describe{
#'   \item{start, end}{Loop boundary indices (global chain iteration).}
#'   \item{length}{Number of steps in the loop (end - start).}
#'   \item{distance}{Base-space closure distance.}
#'   \item{alpha_start_j, alpha_end_j, alpha_transport_j}{Per-group
#'     starting, ending, and transported fiber values.}
#'   \item{transport_error_j}{alpha_transport_j - alpha_start_j.}
#'   \item{empirical_disp_j}{alpha_end_j - alpha_start_j.}
#'   \item{area}{Signed area enclosed by the base trajectory (shoelace).}
#'   \item{holonomy_stokes_j}{Stokes holonomy exp(F_mean_j * area).}
#' }
#'
#' @export
integrate_transport <- function(connection,
                                loops,
                                chain,
                                n_loops = 300L) {

  # ── Setup ──────────────────────────────────────────────────────────────────
  full_mat  <- connection$full_mat
  stan_data <- connection$stan_data
  base_vars <- connection$base_vars
  fib_vars  <- connection$fiber_vars
  bet_vars  <- connection$beta_vars
  J  <- length(fib_vars)
  X  <- stan_data$X
  group <- stan_data$group

  # Select tightest loops up to n_loops
  loops <- loops[order(loops$distance), , drop = FALSE]
  loops <- loops[seq_len(min(n_loops, nrow(loops))), , drop = FALSE]
  K_loops <- nrow(loops)

  # Column indices in full_mat
  base_cols <- base_vars
  fib_cols  <- fib_vars
  bet_cols  <- bet_vars

  # Pre-allocate output
  alpha_start_mat    <- matrix(NA_real_, K_loops, J)
  alpha_end_mat      <- matrix(NA_real_, K_loops, J)
  alpha_trans_mat    <- matrix(NA_real_, K_loops, J)
  area_vec           <- numeric(K_loops)
  stokes_mat         <- matrix(NA_real_, K_loops, J)
  len_vec            <- integer(K_loops)

  for (k in seq_len(K_loops)) {
    i <- loops$start[k]
    j <- loops$end[k]

    # Trajectory: rows i through j (inclusive)
    path <- full_mat[i:j, , drop = FALSE]
    L    <- nrow(path)   # j - i + 1 steps

    alp_i <- as.vector(path[1L, fib_cols])
    alp_j <- as.vector(path[L,  fib_cols])

    # Integrate: sum A(theta_t, alpha_t) * delta_theta over t = 1..L-1
    delta_alpha <- matrix(0, L - 1L, J)
    F_along     <- matrix(NA_real_, L - 1L, J)

    for (t in seq_len(L - 1L)) {
      row_t   <- path[t, ]
      mu_t    <- row_t[[base_vars[1L]]]   # "mu"
      sig_t   <- row_t[[base_vars[2L]]]   # "sigma"
      alp_t   <- as.vector(row_t[fib_cols])
      bet_t   <- as.vector(row_t[bet_cols])
      dtheta  <- as.vector(path[t + 1L, base_cols]) -
                 as.vector(row_t[base_cols])

      G_FF_t  <- .glmm_G_FF(sig_t, alp_t, X, group, bet_t)
      G_BF_t  <- .glmm_G_BF(sig_t, mu_t, alp_t)
      A_t     <- .glmm_connection(G_FF_t, G_BF_t)   # J x 2

      delta_alpha[t, ] <- A_t %*% dtheta
      F_along[t, ]     <- .glmm_curvature_linearised(G_FF_t, sig_t)
    }

    alpha_trans_mat[k, ] <- alp_i + colSums(delta_alpha)
    alpha_start_mat[k, ] <- alp_i
    alpha_end_mat[k, ]   <- alp_j
    len_vec[k]           <- L - 1L

    # Signed area of the base trajectory (shoelace formula)
    base_path <- path[, base_vars, drop = FALSE]
    area_vec[k] <- .shoelace(base_path[, 1L], base_path[, 2L])

    # Stokes holonomy: exp(F_mean * area) per group
    F_mean <- colMeans(F_along, na.rm = TRUE)
    stokes_mat[k, ] <- exp(F_mean * area_vec[k])
  }

  # ── Assemble output data frame ─────────────────────────────────────────────
  result <- data.frame(
    start    = loops$start,
    end      = loops$end,
    length   = len_vec,
    distance = loops$distance,
    area     = area_vec
  )

  for (j in seq_len(J)) {
    jlabel <- gsub("[^[:alnum:]]", "_", fib_vars[j])   # safe column name
    result[[paste0("alpha_start_",     jlabel)]] <- alpha_start_mat[, j]
    result[[paste0("alpha_end_",       jlabel)]] <- alpha_end_mat[,   j]
    result[[paste0("alpha_transport_", jlabel)]] <- alpha_trans_mat[, j]
    result[[paste0("transport_error_", jlabel)]] <- alpha_trans_mat[, j] -
                                                      alpha_start_mat[, j]
    result[[paste0("empirical_disp_",  jlabel)]] <- alpha_end_mat[, j] -
                                                      alpha_start_mat[, j]
    result[[paste0("holonomy_stokes_", jlabel)]] <- stokes_mat[, j]
  }

  attr(result, "fiber_vars") <- fib_vars
  attr(result, "J")          <- J
  attr(result, "connection") <- connection
  class(result)              <- c("fibr_transport", "data.frame")
  result
}

#' Print a fibr_transport object
#' @param x A `fibr_transport` object.
#' @param ... Ignored.
#' @export
print.fibr_transport <- function(x, ...) {
  J      <- attr(x, "J")
  fvars  <- attr(x, "fiber_vars")
  if (is.null(J) || is.null(fvars)) {
    class(x) <- "data.frame"
    return(print(x, ...))
  }
  K      <- nrow(x)

  cat("fibr parallel transport\n")
  cat("=======================\n")
  cat(sprintf("  Loops integrated : %d\n", K))
  cat(sprintf("  Mean loop length : %.1f steps\n", mean(x$length)))
  cat(sprintf("  Mean |area|      : %.4f\n\n", mean(abs(x$area))))

  # Per-group summary: RMSE of transport error vs empirical displacement
  te_rmse <- numeric(J)
  ed_sd   <- numeric(J)
  hs_mean <- numeric(J)

  for (j in seq_len(J)) {
    jlabel <- gsub("[^[:alnum:]]", "_", fvars[j])
    te_rmse[j] <- sqrt(mean(x[[paste0("transport_error_", jlabel)]]^2))
    ed_sd[j]   <- sd(x[[paste0("empirical_disp_",  jlabel)]])
    hs_mean[j] <- mean(x[[paste0("holonomy_stokes_", jlabel)]])
  }

  df <- data.frame(
    group              = seq_len(J),
    `transport RMSE`   = round(te_rmse, 4),
    `empirical SD`     = round(ed_sd,   4),
    `Stokes H (mean)`  = round(hs_mean, 4),
    check.names = FALSE
  )
  cat("Per-group summary:\n")
  print(df, row.names = FALSE)
  cat("\n(transport RMSE: ||alpha_transported - alpha_start||; ")
  cat("empirical SD: sd(alpha_end - alpha_start))\n")
  invisible(x)
}

#' Plot analytical vs empirical fiber displacement or detail loop trajectories
#'
#' @description
#' Two plot types are available:
#' - `"scatter"` (default): compares the connection-predicted displacement
#'   \eqn{\alpha_{\text{transported}} - \alpha_{\text{start}}} against the
#'   observed displacement \eqn{\alpha_{\text{end}} - \alpha_{\text{start}}}
#'   across all loops and fiber dimensions.
#' - `"trajectory"`: plots the detailed step-by-step parallel transport trajectory of
#'   a single loop, contrasting a group with high holonomy against a group with low holonomy.
#'
#' @param x A `fibr_transport` object.
#' @param type `"scatter"` or `"trajectory"`.
#' @param loop_idx Integer; which loop to plot for `"trajectory"` (default 1, the tightest loop).
#' @param groups Integer vector of length 2; group indices to plot for contrast.
#'   If `NULL` (default), automatically selects the groups with the maximum and minimum drift.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (printed invisibly).
#' @export
plot.fibr_transport <- function(x, type = c("scatter", "trajectory"), loop_idx = 1L, groups = NULL, ...) {
  type <- match.arg(type)
  p <- if (type == "scatter") {
    .plot_transport_scatter(x)
  } else {
    .plot_transport_trajectory(x, loop_idx = loop_idx, groups = groups)
  }
  print(p)
  invisible(p)
}

.plot_transport_scatter <- function(x) {
  J     <- attr(x, "J")
  fvars <- attr(x, "fiber_vars")

  rows <- lapply(seq_len(J), function(j) {
    jlabel <- gsub("[^[:alnum:]]", "_", fvars[j])
    data.frame(
      transport = x[[paste0("transport_error_", jlabel)]],
      empirical = x[[paste0("empirical_disp_",  jlabel)]],
      group     = sprintf("j=%d", j)
    )
  })
  df <- do.call(rbind, rows)

  lim <- range(c(df$transport, df$empirical), na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = transport, y = empirical)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = "firebrick", linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, colour = "steelblue") +
    ggplot2::facet_wrap(~ group, nrow = 2L) +
    ggplot2::coord_equal(xlim = lim, ylim = lim) +
    ggplot2::labs(
      title    = "Parallel transport: analytical vs empirical",
      subtitle = "Diagonal = perfect agreement (chain follows horizontal lift)",
      x        = expression(Delta*alpha[transport]~"  (connection integral)"),
      y        = expression(Delta*alpha[empirical]~"  (chain endpoint)  ")
    ) +
    ggplot2::theme_minimal(base_size = 11)

  p
}

.plot_transport_trajectory <- function(x, loop_idx = 1L, groups = NULL) {
  conn <- attr(x, "connection")
  if (is.null(conn)) {
    stop("x must contain a 'connection' attribute to plot trajectories. Please re-run integrate_transport().")
  }

  if (loop_idx < 1L || loop_idx > nrow(x)) {
    stop(sprintf("loop_idx must be between 1 and %d.", nrow(x)))
  }

  loop_row <- x[loop_idx, ]
  start_i  <- loop_row$start
  end_i    <- loop_row$end

  full_mat  <- conn$full_mat
  base_vars <- conn$base_vars
  fib_vars  <- conn$fiber_vars
  bet_vars  <- conn$beta_vars
  stan_data <- conn$stan_data
  J         <- length(fib_vars)
  X         <- stan_data$X
  group     <- stan_data$group
  beta0     <- colMeans(full_mat[, bet_vars, drop = FALSE])

  # Extract path
  path     <- full_mat[start_i:end_i, , drop = FALSE]
  L        <- nrow(path)
  progress <- (seq_len(L) - 1L) / (L - 1L)

  # Parallel transport integration along this loop
  alpha_trans <- matrix(NA_real_, L, J)
  alpha_trans[1L, ] <- as.vector(path[1L, fib_vars])

  for (t in seq_len(L - 1L)) {
    row_t  <- path[t, ]
    mu_t   <- row_t[[base_vars[1L]]]
    sig_t  <- row_t[[base_vars[2L]]]
    alp_t  <- alpha_trans[t, ]
    bet_t  <- as.vector(row_t[bet_vars])
    dtheta <- as.vector(path[t + 1L, base_vars]) - as.vector(row_t[base_vars])

    G_FF_t <- .glmm_G_FF(sig_t, alp_t, X, group, bet_t)
    G_BF_t <- .glmm_G_BF(sig_t, mu_t, alp_t)
    A_t    <- .glmm_connection(G_FF_t, G_BF_t)

    alpha_trans[t + 1L, ] <- alp_t + A_t %*% dtheta
  }

  # Calculate absolute drifts for group selection
  drifts <- alpha_trans[L, ] - alpha_trans[1L, ]

  # Automatic group selection
  if (is.null(groups)) {
    g_high <- which.max(abs(drifts))
    g_low  <- which.min(abs(drifts))
  } else {
    if (length(groups) != 2L) stop("groups must be an integer vector of length 2.")
    g_high <- groups[1L]
    g_low  <- groups[2L]
  }

  # Prepare base space data
  df_base <- data.frame(
    mu       = path[, base_vars[1L]],
    sigma    = path[, base_vars[2L]],
    progress = progress
  )

  # Function to build individual group trajectory plot
  .make_traj_plot <- function(g_idx, label) {
    df_traj <- data.frame(
      progress  = progress,
      transport = alpha_trans[, g_idx],
      empirical = path[, fib_vars[g_idx]]
    )

    start_val <- df_traj$transport[1L]
    end_trans <- df_traj$transport[L]
    end_emp   <- df_traj$empirical[L]

    drift_val  <- end_trans - start_val
    drift_text <- sprintf("Delta*alpha[trans] == %.3f", drift_val)
    y_mid      <- (end_trans + start_val) / 2

    # Y limits for text margin
    y_range <- range(c(df_traj$transport, df_traj$empirical))
    
    col_palette <- ggplot2::scale_colour_viridis_c(option = "plasma", name = "Progress")

    p <- ggplot2::ggplot(df_traj, ggplot2::aes(x = progress)) +
      # Start value reference line
      ggplot2::geom_hline(yintercept = start_val, linetype = "dashed", color = "grey60") +
      # Empirical MCMC path
      ggplot2::geom_path(ggplot2::aes(y = empirical), color = "grey55", linewidth = 1.0, linetype = "dotted") +
      # Parallel transported path
      ggplot2::geom_path(ggplot2::aes(y = transport, color = progress), linewidth = 1.5) +
      col_palette +
      # Start/End points
      ggplot2::geom_point(data = df_traj[1L, ], ggplot2::aes(y = transport), color = "black", size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
      ggplot2::geom_point(data = df_traj[L, ], ggplot2::aes(y = transport), color = "red", size = 2.5, shape = 19) +
      # Drift bracket
      ggplot2::annotate("segment", x = 1.01, xend = 1.01, y = start_val, yend = end_trans,
                        arrow = ggplot2::arrow(ends = "both", angle = 90, length = ggplot2::unit(0.05, "inches")),
                        color = "red", linewidth = 0.8) +
      ggplot2::annotate("text", x = 0.92, y = y_mid, label = drift_text,
                        color = "red", size = 3.2, parse = TRUE, hjust = 1, fontface = "bold") +
      ggplot2::scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                  labels = c("0%", "25%", "50%", "75%", "100%")) +
      ggplot2::labs(
        title    = label,
        x        = "Loop Progress",
        y        = expression(alpha[j])
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(face = "bold", size = 11),
        legend.position    = "none",
        panel.grid.minor   = ggplot2::element_blank()
      )
    p
  }

  p_base <- ggplot2::ggplot(df_base, ggplot2::aes(x = mu, y = sigma, color = progress)) +
    ggplot2::geom_path(linewidth = 2) +
    ggplot2::geom_point(data = df_base[1L, ], color = "black", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    ggplot2::geom_text(data = df_base[1L, ], ggplot2::aes(label = "Start/End"),
                       hjust = -0.2, vjust = -0.5, color = "black", fontface = "bold", size = 3.2) +
    ggplot2::scale_colour_viridis_c(option = "plasma", name = "Progress") +
    ggplot2::labs(
      title    = "Base Space Loop Trajectory",
      subtitle = sprintf("Loop %d (length = %d steps)", loop_idx, L - 1L),
      x        = expression(mu),
      y        = expression(sigma)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle      = ggplot2::element_text(color = "grey40", size = 9),
      legend.position    = "left",
      panel.grid.minor   = ggplot2::element_blank()
    )

  p_high <- .make_traj_plot(g_high, sprintf("Fibre with Holonomy: Group %d (Max Drift)", g_high))
  p_low  <- .make_traj_plot(g_low,  sprintf("Fibre without Holonomy: Group %d (Min Drift)", g_low))

  title_theme <- ggplot2::theme(
    plot.title    = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(color = "grey30", size = 10, hjust = 0.5, margin = ggplot2::margin(b = 10))
  )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_combined <- (p_base | (p_high / p_low)) +
      patchwork::plot_layout(widths = c(1, 1.8)) +
      patchwork::plot_annotation(
        title    = "Fibre Holonomy Trajectory Diagnostic",
        subtitle = "Solid path = parallel transport (analytical); Dotted path = MCMC chain (empirical).",
        theme    = title_theme
      )
    p_combined
  } else {
    print(p_base)
    print(p_high)
    print(p_low)
    p_high
  }
}

#' Linearised holonomy along a synthetic circular loop in base space
#'
#' @description
#' Constructs a circular path in \eqn{(\mu, \sigma)} space and integrates the
#' **linearised** connection form along it (fiber \eqn{\bm{\alpha}} held fixed
#' at `alpha0`, equivalently \eqn{G_{FF}} frozen at the loop centre).  The
#' result is compared with the first-order Stokes approximation
#' \eqn{H_j = 1 + F_j^{\text{lin}} \cdot \pi r^2 / \alpha_{0j}}.
#'
#' **This function does not integrate the true (full) Ehresmann connection.**
#' The true connection is flat (full curvature identically zero;
#' see `data-raw/verify_flat_connection.R`), so its holonomy over any loop is
#' trivial.  What this function computes is the holonomy of the
#' *fiber-frozen* linearisation, which is the object validated in Figure 3
#' of the companion paper and whose curvature is
#' \eqn{F_j^{\text{lin}} = -2/(\sigma^5 G_{FF,j}^2)}.
#'
#' @param connection A `fibr_connection` object.
#' @param mu0,sigma0 Centre of the circle in base space.  Defaults to the
#'   posterior mean of `mu` and `sigma` from the connection's subsampled draws.
#' @param alpha0 Fiber point at which to evaluate the linearised connection
#'   (J-vector).  Defaults to the posterior mean of alpha.
#' @param beta0 Fixed-effect coefficients (2-vector).  Defaults to posterior
#'   mean from the subsampled draws.
#' @param radius Radius of the circle in raw \eqn{(\mu, \sigma)} units
#'   (default 0.3, roughly 0.3 posterior SDs for typical hierarchical models).
#' @param n_steps Number of discretisation steps around the circle (default 200).
#'
#' @return A data frame with one row per group j and columns:
#' \describe{
#'   \item{j}{Group index.}
#'   \item{H_numerical}{Holonomy scalar from integrating the connection:
#'     \eqn{1 + \oint A_j\,d\theta}.}
#'   \item{H_stokes}{Stokes approximation: \eqn{\exp(F_j \cdot \pi r^2)}.}
#'   \item{F_j}{Curvature at the centre point.}
#'   \item{A_mu, A_sigma}{Connection coefficients at the centre.}
#' }
#'
#' @export
synthetic_holonomy_loop <- function(connection,
                                    mu0     = NULL,
                                    sigma0  = NULL,
                                    alpha0  = NULL,
                                    beta0   = NULL,
                                    radius  = 0.3,
                                    n_steps = 200L) {

  sd_obj  <- connection$stan_data
  fvars   <- connection$fiber_vars
  bvars   <- connection$beta_vars
  J       <- length(fvars)

  # Default centre: posterior mean from subsampled draws
  if (is.null(mu0))    mu0    <- mean(connection$base_pts[, "mu"])
  if (is.null(sigma0)) sigma0 <- mean(connection$base_pts[, "sigma"])
  if (is.null(alpha0)) alpha0 <- colMeans(connection$fiber_pts)
  if (is.null(beta0))  beta0  <- colMeans(
    connection$full_mat[, bvars, drop = FALSE]
  )

  # ── Circular path ───────────────────────────────────────────────────────────
  t_grid  <- seq(0, 2 * pi, length.out = n_steps + 1L)
  mu_path    <- mu0    + radius * cos(t_grid)
  sigma_path <- sigma0 + radius * sin(t_grid)

  # Clamp sigma to stay positive
  sigma_path <- pmax(sigma_path, 1e-4)

  # ── Integrate connection along the path (alpha fixed at alpha0) ─────────────
  # The linearised transport: alpha moves by A(theta, alpha0) * dtheta at each step.
  # After a full loop, the net displacement is sum(A * dtheta).
  delta_alpha_total <- numeric(J)

  for (s in seq_len(n_steps)) {
    mu_s    <- mu_path[s]
    sig_s   <- sigma_path[s]
    dmu     <- mu_path[s + 1L]    - mu_s
    dsig    <- sigma_path[s + 1L] - sig_s

    G_FF_s <- .glmm_G_FF(sig_s, alpha0, sd_obj$X, sd_obj$group, beta0)
    G_BF_s <- .glmm_G_BF(sig_s, mu_s,   alpha0)
    A_s    <- .glmm_connection(G_FF_s, G_BF_s)   # J x 2

    delta_alpha_total <- delta_alpha_total + A_s %*% c(dmu, dsig)
  }

  # Shared denominator for both H_numerical and H_stokes: alpha0 floored at
  # 0.05 (in absolute value) to avoid divide-by-zero near the origin.
  # Both scalars must use the same denom so the comparison is consistent.
  denom <- ifelse(abs(alpha0) > 0.05, alpha0, sign(alpha0) * 0.05)

  H_numerical <- 1 + delta_alpha_total / denom

  # ── Stokes approximation ────────────────────────────────────────────────────
  # For a vector bundle (additive parallel transport), the first-order Stokes
  # formula is:
  #   delta_alpha_j = F_j * area      (Stokes; F has units alpha/base^2)
  #   H_j = 1 + delta_alpha_j / alpha0_j = 1 + F_j * area / alpha0_j
  #
  # The exponential form H = exp(F*area) applies to MULTIPLICATIVE holonomy
  # (principal bundles / Lie groups) and is NOT correct for this vector bundle.
  G_FF_c  <- .glmm_G_FF(sigma0, alpha0, sd_obj$X, sd_obj$group, beta0)
  G_BF_c  <- .glmm_G_BF(sigma0, mu0,   alpha0)
  A_c     <- .glmm_connection(G_FF_c, G_BF_c)
  F_c     <- .glmm_curvature_linearised(G_FF_c, sigma0)

  area_circle  <- pi * radius^2
  H_stokes     <- 1 + F_c * area_circle / denom

  data.frame(
    j           = seq_len(J),
    H_numerical = round(as.vector(H_numerical), 5),
    H_stokes    = round(H_stokes,               5),
    F_j         = round(F_c,                    5),
    delta_alpha = round(F_c * area_circle,       5),
    A_mu        = round(A_c[, 1L],              5),
    A_sigma     = round(A_c[, 2L],              5)
  )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Signed area enclosed by a polygon with vertices (x, y) via shoelace formula.
# Returns positive area for counter-clockwise orientation.
.shoelace <- function(x, y) {
  n <- length(x)
  if (n < 3L) return(0)
  i  <- seq_len(n)
  ip <- c(seq(2L, n), 1L)
  0.5 * sum(x[i] * y[ip] - x[ip] * y[i])
}
