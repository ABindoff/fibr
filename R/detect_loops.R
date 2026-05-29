#' Find approximate loops in base (hyperparameter) space
#'
#' Scans a chain for pairs of iterations where the base-space parameters
#' return within `epsilon` of a previous position, separated by at least
#' `min_gap` steps.  These loop pairs are the raw material for estimating
#' the holonomy transport map.
#'
#' @param base_draws Numeric matrix \[n_iter x n_base\]. Each row is one
#'   iteration; columns are the base-space (hyperparameter) variables.
#' @param epsilon Scalar tolerance for loop closure (in standardised units
#'   when `scale = TRUE`). `NULL` auto-selects as the 5th percentile of a
#'   sub-sampled pairwise distance distribution.
#' @param min_gap Minimum number of iterations between loop start and end
#'   (default 50). Screens out autocorrelated neighbours.
#' @param k Number of nearest neighbours to query per point (default 100).
#'   Should exceed the expected number of points within `epsilon`.
#' @param max_loops Maximum loop pairs returned (default 5000). When more
#'   are found, the tightest (smallest distance) are kept.
#' @param scale Logical; standardise base variables before distance
#'   computation so that `epsilon` is in units of posterior SDs (default TRUE).
#'
#' @return A data frame with columns `start`, `end`, `distance` (one row
#'   per loop pair). Attribute `"epsilon"` records the tolerance used.
#'
#' @keywords internal
detect_loops <- function(base_draws,
                         epsilon   = NULL,
                         min_gap   = 50L,
                         k         = 100L,
                         max_loops = 5000L,
                         scale     = TRUE) {

  base_draws <- as.matrix(base_draws)
  n <- nrow(base_draws)

  if (scale) {
    mu_b  <- colMeans(base_draws)
    sd_b  <- apply(base_draws, 2L, sd)
    sd_b[sd_b == 0] <- 1
    base_draws <- sweep(sweep(base_draws, 2L, mu_b), 2L, sd_b, "/")
  }

  if (is.null(epsilon)) {
    n_samp   <- min(800L, n)
    idx      <- sample.int(n, n_samp)
    sub_dist <- as.vector(dist(base_draws[idx, , drop = FALSE]))
    epsilon  <- unname(quantile(sub_dist, 0.05))
    message(sprintf(
      "fibr: auto epsilon = %.4f  (5th pct of %d sampled pairwise distances)",
      epsilon, length(sub_dist)
    ))
  }

  k_use <- min(k, n - 1L)
  # get.knnx returns the k_use+1 nearest neighbours; first is self (dist=0)
  nn <- FNN::get.knnx(base_draws, base_draws, k = k_use + 1L)
  nn_idx  <- nn$nn.index[, -1L, drop = FALSE]   # n x k_use
  nn_dist <- nn$nn.dist[,  -1L, drop = FALSE]

  loop_list <- vector("list", n)
  for (i in seq_len(n)) {
    j    <- nn_idx[i, ]
    d    <- nn_dist[i, ]
    keep <- (j > i + min_gap) & (d < epsilon)
    if (any(keep)) {
      loop_list[[i]] <- data.frame(
        start    = i,
        end      = j[keep],
        distance = d[keep]
      )
    }
  }

  non_null <- !vapply(loop_list, is.null, logical(1L))
  if (!any(non_null)) {
    stop(sprintf(
      "No loops found (epsilon=%.4f, min_gap=%d). Try larger epsilon or smaller min_gap.",
      epsilon, min_gap
    ))
  }

  result <- do.call(rbind, loop_list[non_null])
  rownames(result) <- NULL

  if (nrow(result) > max_loops) {
    result <- result[order(result$distance)[seq_len(max_loops)], ]
  }

  attr(result, "epsilon") <- epsilon
  result
}
