# fibr <img src="man/figures/logo.svg" align="right" height="140" />

Prior-fraction diagnostics and geometry-guided reparameterisation for hierarchical MCMC.

## Overview

Standard MCMC diagnostics (R̂, ESS, divergence counts) report whether a chain has mixed, not which parameters are causing slow mixing or why. For hierarchical models the centring/non-centring trade-off is the dominant obstruction, and it varies across groups within a single fit.

**fibr** computes the *prior fraction*

```
π_j = (1/σ²) / G_{FF,j}
```

for each group-level coordinate: the share of that coordinate's posterior precision contributed by the prior rather than by the data. This is the classical shrinkage/pooling factor (Gelman & Pardoe 2006) and the per-group information ratio of Betancourt & Girolami (2015). The companion paper shows that once the Ehresmann connection on the base–fiber bundle is known to be flat, π_j is not a heuristic but the unique quantity governing the centring obstruction.

Coordinates with π_j near 1 are prior-dominated: their posteriors reflect regularisation toward the population mean rather than group-specific data. Coordinates with π_j near 0 are data-driven. The same quantity gives the per-group non-centring weight for partial reparameterisation.

The package also provides a chain-level diagnostic (`holonomy_diagnostic`) that estimates base–fiber conditional dependence directly from MCMC output without assuming GLMM structure, and analytic connection/curvature functions for the centred logistic GLMM (`compute_connection`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("ABindoff/fibr")
```

Requires R ≥ 4.1.0. The brms adapter for `prior_fraction()` requires brms and a Stan backend (cmdstanr or rstan).

## Usage

### Prior fraction from a brms fit

```r
library(brms)
library(fibr)

fit <- brm(y ~ 1 + (1 | site), data = dat, family = poisson(),
           chains = 4, cores = 4)

pf <- prior_fraction(fit)
print(pf)   # summary: how many coordinates exceed the 0.8 threshold
plot(pf)    # π_j vs. log(n_obs) per coordinate
```

Supported families: gaussian, bernoulli, binomial, poisson, negbinomial. Correlated random effects are handled per-marginal with a diagnostic message.

### Manual path (any model)

```r
# Supply prior precision and per-group likelihood information directly.
# For a Gaussian outcome: likelihood information = n_j / sigma_y^2.
# For Bernoulli: sum of p_i(1-p_i) within each group.
sigma_re <- 1.5
n_obs    <- c(2, 5, 20, 80)
sigma_y  <- 1.0

prior_fraction(1 / sigma_re^2, lik_information = n_obs / sigma_y^2)
```

This path works for Stan, JAGS, or any other fitting engine. The helper `fibr:::.glm_information()` computes per-observation Fisher weights from the linear predictor for the supported GLM families.

### Analytic connection (centred logistic GLMM)

```r
# chain is a posterior::draws_array from the centred Stan model
conn <- compute_connection(
  chain      = chain,
  base_vars  = c("mu", "sigma"),
  fiber_vars = paste0("alpha[", 1:J, "]"),
  method     = "analytic_glmm",
  stan_data  = stan_data
)
print(conn)   # per-group A[j,mu], A[j,sigma], linearised curvature, π_j
plot(conn)
```

### Chain-level holonomy diagnostic

```r
hd <- holonomy_diagnostic(
  chain      = chain,
  base_vars  = c("mu", "sigma"),
  fiber_vars = paste0("alpha[", 1:J, "]"),
  min_gap    = 50
)
print(hd)
plot(hd)   # eigenvalue spectrum of the estimated transport map
```

`holonomy_diagnostic` detects loops in base-space (hyperparameter) trajectory and estimates the conditional dependence of the fiber at matched loop endpoints. Eigenvalues of the transport map near 1 indicate coupling; the magnitude tracks π_j in prior-dominated groups.

## Background

The joint parameter space of a two-level hierarchical model — hyperparameters (μ, σ) as the base, group parameters (α₁, …, α_J) as the fibers — carries a fiber bundle structure. The Fisher information metric induces an Ehresmann connection A = −G_FF⁻¹ G_BF on this bundle. A natural conjecture is that the centring obstruction arises from the curvature of this connection, manifesting as holonomy when the hyperparameters traverse a closed loop.

The companion paper proves this false for any smooth hierarchical posterior: the connection is flat, its horizontal leaves being the level sets of the fiber score ∂_α log p, so no geometric obstruction exists above the metric level. The non-zero curvature seen when integrating the *linearised* (fiber-frozen) connection is an artefact of that approximation; numerical verification is in `data-raw/verify_flat_connection.R` and `tests/testthat/test-connection.R`.

What remains is statistical: the conditional dependence of fiber on base, governed per group by π_j. Genuine curvature — including rotational holonomy — does appear for connections built from a sampler's working metric (a fixed mass matrix) rather than the true Hessian, making holonomy an algorithmic rather than geometric phenomenon in that setting.

## Reference

Bindoff, A.D. (2026). *A Flat Connection: The Prior Fraction and the Geometry of Centring in Hierarchical MCMC.* Preprint. <https://doi.org/10.5281/zenodo.20687656>

```bibtex
@misc{bindoff2026fibr,
  author = {Bindoff, Aidan D.},
  title  = {A Flat Connection: The Prior Fraction and the Geometry of
             Centring in Hierarchical {MCMC}},
  year   = {2026},
  doi    = {10.5281/zenodo.20687656},
  url    = {https://github.com/ABindoff/fibr}
}
```
