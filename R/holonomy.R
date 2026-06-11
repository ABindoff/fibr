#' Holonomy diagnostic for hierarchical MCMC chains
#'
#' @description
#' Estimates the holonomy transport map \eqn{\hat{H}} from an MCMC chain on
#' a hierarchical model.  Non-trivial holonomy --- eigenvalues of \eqn{\hat{H}}
#' away from 1 --- indicates that the fiber bundle structure of the parameter
#' space is creating a geometric mixing obstruction invisible to standard
#' diagnostics (R-hat, ESS, divergences).
#'
#' Loop detection and transport estimation run **per chain** so that
#' independent chains are never spliced into a single trajectory.  Loop pairs
#' from all chains are pooled for the final regression.
#'
#' @param chain A [`posterior::draws_array`] (recommended), `draws_matrix`,
#'   or a plain numeric matrix with named columns.  Rows are iterations,
#'   columns are parameters.
#' @param base_vars Character vector naming the base-space (hyperparameter)
#'   columns, e.g. `c("mu", "sigma")`.
#' @param fiber_vars Character vector naming the fiber (group-level) parameter
#'   columns, e.g. `paste0("alpha[", 1:8, "]")`.
#' @param epsilon Loop-detection tolerance in standardised base-space units.
#'   `NULL` (default) auto-selects per chain.
#' @param n_bootstrap Bootstrap resamples for eigenvalue uncertainty (default 200).
#' @param min_gap Minimum iterations between loop start and end (default 50).
#' @param k Nearest neighbours queried per point in loop detection (default 100).
#' @param max_loops Maximum loop pairs per chain (default 5000).
#' @param weights Loop weighting for the transport regression; `"distance"`
#'   (default) or `"uniform"`. See [estimate_transport_map()].
#' @param structure Transport map structure; `"diagonal"` (default) fits a
#'   per-group scalar contraction \eqn{h_j} (the correct model class when the
#'   fiber metric is diagonal, as in the two-level GLMM), `"full"` fits an
#'   unrestricted \eqn{J \times J} matrix.  See [estimate_transport_map()].
#' @param ridge Ridge penalty for the transport map Gram matrix (default `1e-6`).
#' @param residualize_fiber If `TRUE` (default), fiber draws are residualised
#'   against the base draws via OLS before transport estimation.  This removes
#'   the linear effect of the base on the fiber (e.g. the mean-shift
#'   `alpha ~ mu` in a centered GLMM), isolating the true vertical component
#'   of the fiber bundle where geometric holonomy lives.  Set to `FALSE` to
#'   operate on raw fiber coordinates.
#'
#' @return An S3 object of class `fibr_holonomy`. See [print.fibr_holonomy()]
#'   and [plot.fibr_holonomy()] for display methods.  The list contains:
#' \describe{
#'   \item{H}{Estimated \eqn{J \times J} transport matrix.}
#'   \item{eigenvalues}{Complex eigenvalues of `H`, sorted by decreasing modulus.}
#'   \item{boot_eigenvalues}{Complex matrix of bootstrapped eigenvalues.}
#'   \item{frobenius_dev}{\eqn{\|H - I\|_F}.}
#'   \item{n_loops}{Total loop pairs used.}
#'   \item{loops}{Data frame of all detected loop pairs.}
#'   \item{residualized}{Logical; whether fiber was residualised against base.}
#'   \item{call, base_vars, fiber_vars}{Metadata.}
#' }
#'
#' @export
holonomy_diagnostic <- function(chain,
                                base_vars,
                                fiber_vars,
                                epsilon            = NULL,
                                n_bootstrap        = 200L,
                                min_gap            = 50L,
                                k                  = 100L,
                                max_loops          = 5000L,
                                weights            = c("distance", "uniform"),
                                structure          = c("diagonal", "full"),
                                ridge              = 1e-6,
                                residualize_fiber  = TRUE) {

  cl        <- match.call()
  weights   <- match.arg(weights)
  structure <- match.arg(structure)

  # ── Split into per-chain matrices ────────────────────────────────────────────
  chains <- .split_chains(chain)   # list of matrices, one per chain

  .check_vars(chains[[1L]], base_vars,  "base_vars")
  .check_vars(chains[[1L]], fiber_vars, "fiber_vars")

  # ── Detect loops per chain, pool results ─────────────────────────────────────
  message("fibr: detecting loops in base space...")
  offsets   <- c(0L, cumsum(vapply(chains, nrow, integer(1L))))
  all_loops <- vector("list", length(chains))

  for (c_idx in seq_along(chains)) {
    base_c <- chains[[c_idx]][, base_vars, drop = FALSE]
    all_loops[[c_idx]] <- tryCatch({
      lp <- detect_loops(base_c, epsilon = epsilon, min_gap = min_gap,
                         k = k, max_loops = max_loops)
      # Offset iteration indices into the global flattened sequence
      lp$start <- lp$start + offsets[c_idx]
      lp$end   <- lp$end   + offsets[c_idx]
      lp
    }, error = function(e) {
      message(sprintf("fibr:   chain %d — %s", c_idx, conditionMessage(e)))
      NULL
    })
  }

  non_null <- !vapply(all_loops, is.null, logical(1L))
  if (!any(non_null)) stop("No loops found in any chain.")

  # Preserve the per-chain epsilon values before rbind drops the attr
  epsilons <- vapply(all_loops[non_null],
                     function(lp) attr(lp, "epsilon") %||% NA_real_,
                     numeric(1L))

  loops <- do.call(rbind, all_loops[non_null])
  rownames(loops) <- NULL
  attr(loops, "epsilon") <- mean(epsilons, na.rm = TRUE)

  message(sprintf("fibr: found %d loop pairs across %d chain(s)",
                  nrow(loops), sum(non_null)))

  # ── Pool all chains into one matrix ──────────────────────────────────────────
  full_mat   <- do.call(rbind, chains)
  full_base  <- full_mat[, base_vars,  drop = FALSE]
  full_fiber <- full_mat[, fiber_vars, drop = FALSE]

  # ── Residualise fiber against base ───────────────────────────────────────────
  # Removes the linear base→fiber effect (e.g. alpha ~ mu in a centred GLMM),
  # isolating the vertical fiber component where holonomy lives.
  if (residualize_fiber) {
    message("fibr: residualising fiber against base...")
    full_fiber <- .residualize(full_fiber, full_base)
  }

  # ── Estimate transport map ───────────────────────────────────────────────────
  message(sprintf("fibr: estimating transport map (n_bootstrap = %d)...", n_bootstrap))
  tm <- estimate_transport_map(full_fiber, loops,
                               n_bootstrap = n_bootstrap,
                               ridge       = ridge,
                               weights     = weights,
                               structure   = structure)

  structure(
    list(
      H                = tm$H,
      eigenvalues      = tm$eigenvalues,
      boot_eigenvalues = tm$boot_eigenvalues,
      frobenius_dev    = tm$frobenius_dev,
      n_loops          = tm$n_loops,
      loops            = loops,
      structure        = tm$structure,
      residualized     = residualize_fiber,
      call             = cl,
      base_vars        = base_vars,
      fiber_vars       = fiber_vars
    ),
    class = "fibr_holonomy"
  )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Returns a list of plain matrices, one per chain.
.split_chains <- function(chain) {
  if (inherits(chain, c("draws_array", "draws_matrix", "draws_df",
                        "draws_list", "draws_rvars"))) {
    arr <- posterior::as_draws_array(chain)
    n_chains <- dim(arr)[2L]
    lapply(seq_len(n_chains), function(c_idx) {
      m <- posterior::as_draws_matrix(
        posterior::subset_draws(arr, chain = c_idx)
      )
      mat <- as.matrix(m)
      colnames(mat) <- posterior::variables(m)
      mat
    })
  } else if (is.matrix(chain) || is.data.frame(chain)) {
    list(as.matrix(chain))
  } else {
    stop("chain must be a draws_array, draws_matrix, matrix, or data frame.")
  }
}

.check_vars <- function(mat, vars, arg) {
  miss <- setdiff(vars, colnames(mat))
  if (length(miss) > 0L)
    stop(sprintf("%s not found in chain: %s", arg, paste(miss, collapse = ", ")))
}

# OLS residuals of fiber against base (with intercept).
# Removes linear base→fiber dependence; result has ~zero column means.
.residualize <- function(fiber_mat, base_mat) {
  X    <- cbind(1, base_mat)
  beta <- solve(crossprod(X), crossprod(X, fiber_mat))
  fiber_mat - X %*% beta
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
