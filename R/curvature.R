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
      F_along[t, ]     <- .glmm_curvature(G_FF_t, sig_t)
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

#' Plot analytical vs empirical fiber displacement
#'
#' Scatter plot comparing the connection-predicted displacement
#' \eqn{\alpha_{\text{transported}} - \alpha_{\text{start}}} against the
#' observed displacement \eqn{\alpha_{\text{end}} - \alpha_{\text{start}}}
#' across all loops and fiber dimensions.  Agreement along the diagonal
#' indicates the MCMC chain is tracking the horizontal lift.
#'
#' @param x A `fibr_transport` object.
#' @param ... Ignored.
#' @return A `ggplot` object (printed invisibly).
#' @export
plot.fibr_transport <- function(x, ...) {
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

  print(p)
  invisible(p)
}

#' Holonomy along a synthetic circular loop in base space
#'
#' @description
#' Constructs a circular path in \eqn{(\mu, \sigma)} space, integrates the
#' analytic connection form along it, and compares the result with the Stokes
#' approximation \eqn{H_j \approx \exp(F_j \cdot \pi r^2)}.  This is the
#' clean theoretical validation that bypasses the MCMC loop-area problem.
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

  # H_numerical[j] = 1 + delta_alpha[j] / alpha0[j]
  # (holonomy scalar: how much the fiber is displaced relative to its start)
  # Use absolute displacement when alpha0 is near zero to avoid divide-by-zero
  H_numerical <- 1 + delta_alpha_total / pmax(abs(alpha0), 0.1)

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
  F_c     <- .glmm_curvature(G_FF_c, sigma0)

  area_circle  <- pi * radius^2
  denom        <- ifelse(abs(alpha0) > 0.05, alpha0, sign(alpha0) * 0.05)
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
