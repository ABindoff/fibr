#' Per-coordinate prior fraction (shrinkage / pooling factor)
#'
#' @description
#' For each group-level (random-effect) coordinate \eqn{\alpha_j} of a fitted
#' hierarchical model, the \emph{prior fraction}
#' \deqn{\pi_j = \frac{\text{prior precision}}{\text{prior precision} +
#'   \text{likelihood information}} = \frac{1/\sigma^2}{G_{FF,j}}}
#' is the share of that coordinate's posterior precision contributed by the
#' prior rather than by its own data. It is the classical shrinkage / pooling
#' factor (Gelman and Pardoe, 2006) and the per-group information ratio of
#' Betancourt and Girolami (2015).
#'
#' \strong{Interpretation.} \eqn{\pi_j \approx 1} means the coordinate is
#' \emph{prior-dominated}: its posterior is essentially the prior pushed through
#' shrinkage, so the estimate is mostly regularisation toward the population and
#' should not be over-interpreted unless the prior is one you would defend.
#' \eqn{\pi_j \approx 0} means the data speak. This is a \emph{prior-influence}
#' report, not a convergence diagnostic, and it is read-only: nothing is
#' reparameterised or refit.
#'
#' \strong{Scope and limits.} The estimate is exact for the common GLM families
#' (gaussian, bernoulli, binomial, poisson, negbinomial) with the standard
#' \code{(... | g)} random-effect structure. For \emph{correlated} random
#' effects it reports the per-marginal fraction (using each coefficient's own
#' \code{sd}); the full story there is the eigenvalues of a matrix pooling
#' factor, and a message is emitted. Coordinates with no data
#' (\code{n_obs == 0}) are flagged with \eqn{\pi = 1}. Smooths and GP terms have
#' correlated coordinates and should be read with that caveat. The diagnostic
#' says nothing about multimodality, aliasing, or likelihood mis-specification.
#'
#' @param x A fitted model. Methods are provided for \code{brmsfit} objects and
#'   a \code{default} method taking the prior precision directly (see Details).
#' @param ... Passed to methods.
#'
#' @return A data frame of class \code{fibr_prior_fraction} with one row per
#'   coordinate and columns \code{group}, \code{coef}, \code{level},
#'   \code{n_obs}, \code{prior_sd}, \code{lik_info}, \code{pi}. Has
#'   \code{print} and \code{plot} methods.
#'
#' @references
#' Gelman and Pardoe (2006), \emph{Technometrics} 48(2):241--251.
#' Betancourt and Girolami (2015), in \emph{Current Trends in Bayesian
#' Methodology with Applications}.
#'
#' @examples
#' ## Manual path (no model fit needed): supply the per-coordinate prior
#' ## precision (1/sigma^2) and likelihood information (sum of per-observation
#' ## Fisher information). This is the closed-form GLMM prior fraction.
#' sigma <- 1.5
#' lik   <- c(0.2, 1.0, 5.0)          # e.g. sum p(1-p) for three groups
#' prior_fraction(1 / sigma^2, lik_information = lik)
#'
#' \dontrun{
#' ## brms path: which group-level estimates are prior-dominated?
#' library(brms)
#' fit <- brm(count ~ 1 + (1 | site), data = my_data, family = poisson())
#' pf  <- prior_fraction(fit)
#' pf                                  # summary: how many coordinates have pi > 0.8
#' plot(pf)                            # pi vs. number of observations
#' }
#'
#' @export
utils::globalVariables(c("n_plot", "group"))

prior_fraction <- function(x, ...) UseMethod("prior_fraction")

#' @describeIn prior_fraction Manual path for any model. Supply the per-coordinate
#'   prior precision \code{x} (\eqn{1/\sigma^2}) and the per-coordinate likelihood
#'   information \code{lik_information} (\eqn{G_{FF,j} - 1/\sigma^2}); optionally a
#'   \code{labels} data frame to carry through. Use this to validate against the
#'   closed-form GLMM or to handle Stan fits this package does not parse.
#' @param lik_information Numeric vector of per-coordinate likelihood information.
#' @param labels Optional data frame of label columns (recycled / bound to output).
#' @export
prior_fraction.default <- function(x, lik_information, labels = NULL, ...) {
  prior_precision <- as.numeric(x)
  if (length(prior_precision) == 1L)
    prior_precision <- rep(prior_precision, length(lik_information))
  if (length(prior_precision) != length(lik_information))
    stop("prior_fraction(): 'x' (prior precision) and 'lik_information' must ",
         "have the same length.")
  pi <- prior_precision / (prior_precision + lik_information)
  pi[!is.finite(pi)] <- 1                      # zero total precision -> all prior
  out <- data.frame(
    prior_sd = 1 / sqrt(prior_precision),
    lik_info = lik_information,
    pi       = pi,
    stringsAsFactors = FALSE
  )
  if (!is.null(labels)) out <- cbind(labels, out)
  class(out) <- c("fibr_prior_fraction", "data.frame")
  out
}

# Per-observation likelihood (Fisher) information w.r.t. the linear predictor:
# the GLM IRLS working weight, (d mu / d eta)^2 / Var(y | mu). Length-N vector.
.glm_information <- function(family, eta, dispersion = 1, trials = 1) {
  switch(family,
    gaussian    = rep(1 / dispersion, length(eta)),          # 1 / sigma^2
    bernoulli   = { p <- stats::plogis(eta); p * (1 - p) },
    binomial    = { p <- stats::plogis(eta); trials * p * (1 - p) },
    poisson     = exp(eta),                                   # mu
    negbinomial = { m <- exp(eta); m / (1 + m / dispersion) },# NB2, dispersion = shape
    stop(sprintf(
      "prior_fraction(): family '%s' is not built in. Compute the likelihood ",
      family),
      "information yourself and use prior_fraction.default().")
  )
}

#' @describeIn prior_fraction Adapter for \pkg{brms} fits. Extracts the
#'   random-effect structure, per-coordinate prior SDs, and the family
#'   information at the posterior mean, and returns the per-coordinate prior
#'   fraction. Requires \pkg{brms}.
#' @param ndraws Number of posterior draws to subsample when forming the
#'   posterior-mean linear predictor (for speed). Default 200.
#' @export
prior_fraction.brmsfit <- function(x, ndraws = 200L, ...) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("prior_fraction(): the 'brms' package is required for brmsfit objects.")
  fit <- x

  fam <- fit$family$family
  if (!fam %in% c("gaussian", "bernoulli", "binomial", "poisson", "negbinomial"))
    stop(sprintf("prior_fraction(): family '%s' not yet supported by the brms ", fam),
         "adapter; use prior_fraction.default() with your own likelihood information.")

  sdata <- brms::standata(fit)
  draws <- posterior::as_draws_df(fit)
  ranef <- fit$ranef                              # term structure (brms internal)
  if (is.null(ranef) || nrow(ranef) == 0L)
    stop("prior_fraction(): the model has no group-level (random) effects.")

  # Linear predictor at the posterior mean (subsampled for speed).
  eta <- colMeans(brms::posterior_linpred(fit, ndraws = ndraws))

  # Dispersion for the family information.
  disp <- 1
  if (fam == "gaussian"    && !is.null(draws$sigma)) disp <- mean(draws$sigma)^2
  if (fam == "negbinomial" && !is.null(draws$shape)) disp <- mean(draws$shape)
  trials <- if (!is.null(sdata$trials)) sdata$trials else 1
  info <- .glm_information(fam, eta, dispersion = disp, trials = trials)

  lvl <- attr(ranef, "levels")                    # named list of level labels
  any_correlated <- FALSE
  rows_out <- list()

  for (id in sort(unique(ranef$id))) {
    term <- ranef[ranef$id == id, , drop = FALSE]
    group <- term$group[1]
    Jname <- paste0("J_", id)
    if (is.null(sdata[[Jname]]))
      stop("prior_fraction(): could not find '", Jname, "' in standata(); ",
           "brms internals may differ in this version. Inspect names(standata(fit)).")
    J     <- sdata[[Jname]]
    nlev  <- as.integer(tapply(rep(1L, length(J)), factor(J, levels = seq_len(max(J))), sum))
    nlev[is.na(nlev)] <- 0L
    labs  <- if (!is.null(lvl) && !is.null(lvl[[group]])) lvl[[group]] else seq_len(max(J))
    if (nrow(term) > 1L) any_correlated <- TRUE   # >1 coef in a term => Cholesky

    for (k in seq_len(nrow(term))) {
      cn       <- term$cn[k]
      coefname <- term$coef[k]
      Z        <- sdata[[paste0("Z_", id, "_", cn)]]
      if (is.null(Z)) Z <- rep(1, length(J))      # intercept-only fallback
      sd_col   <- paste0("sd_", group, "__", coefname)
      if (is.null(draws[[sd_col]]))
        stop("prior_fraction(): could not find draws column '", sd_col, "'.")
      sigma       <- mean(draws[[sd_col]])
      prior_prec  <- 1 / sigma^2

      contrib <- Z^2 * info
      Ilev <- as.numeric(tapply(contrib, factor(J, levels = seq_len(max(J))), sum))
      Ilev[is.na(Ilev)] <- 0
      pi   <- prior_prec / (prior_prec + Ilev)
      pi[!is.finite(pi)] <- 1

      rows_out[[length(rows_out) + 1L]] <- data.frame(
        group = group, coef = coefname,
        level = as.character(labs)[seq_along(Ilev)],
        n_obs = nlev, prior_sd = sigma,
        lik_info = Ilev, pi = pi,
        stringsAsFactors = FALSE
      )
    }
  }

  if (any_correlated)
    message("prior_fraction(): correlated random effects detected; reporting the ",
            "per-marginal prior fraction (each coefficient's own sd). The full ",
            "story is the eigenvalues of the matrix pooling factor.")

  res <- do.call(rbind, rows_out)
  rownames(res) <- NULL
  class(res) <- c("fibr_prior_fraction", "data.frame")
  res
}

#' @export
print.fibr_prior_fraction <- function(x, threshold = 0.8, ...) {
  n <- nrow(x)
  hi <- sum(x$pi > threshold, na.rm = TRUE)
  zero <- sum(x$n_obs == 0L, na.rm = TRUE)
  cat(sprintf("<fibr prior fraction>  %d coordinate(s)\n", n))
  cat(sprintf("  prior-dominated (pi > %.2f): %d (%.0f%%)\n",
              threshold, hi, 100 * hi / max(n, 1L)))
  if (zero > 0)
    cat(sprintf("  unidentified (no data, pi = 1): %d\n", zero))
  cat(sprintf("  pi range: [%.3f, %.3f], median %.3f\n",
              min(x$pi, na.rm = TRUE), max(x$pi, na.rm = TRUE),
              stats::median(x$pi, na.rm = TRUE)))
  cat("  (pi near 1 = estimate is mostly prior/shrinkage; near 0 = data-driven)\n")
  print(utils::head(as.data.frame(x[order(-x$pi), ]), 10), row.names = FALSE)
  if (n > 10) cat(sprintf("  ... %d more rows\n", n - 10))
  invisible(x)
}

#' @export
plot.fibr_prior_fraction <- function(x, threshold = 0.8, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot(): the 'ggplot2' package is required.")
  df <- as.data.frame(x)
  df$n_plot <- pmax(df$n_obs, 0.5)               # keep n_obs = 0 visible on log axis
  ggplot2::ggplot(df, ggplot2::aes(x = n_plot, y = pi, colour = group)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2, colour = "grey50") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "observations loading on the coordinate (log scale)",
      y = expression(pi[j]~"(prior fraction)"),
      colour = "group",
      title = "Prior fraction by coordinate",
      subtitle = "high = prior-dominated (mostly shrinkage); low = data-driven"
    ) +
    ggplot2::theme_minimal()
}
