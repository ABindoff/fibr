# fibr <img src="man/figures/logo.svg" align="right" height="140" />

Prior-fraction diagnostics for hierarchical models.

## Overview

Standard MCMC diagnostics (R̂, ESS, divergence counts) report whether a chain has mixed, not which parameters are causing slow mixing or why. For hierarchical models the centring/non-centring trade-off is the dominant obstruction, and it varies across groups within a single fit.

**fibr** computes the *prior fraction*

```
π_j = (1/σ²) / G_{FF,j}
```

for each group-level coordinate: the share of that coordinate's posterior precision contributed by the prior rather than by the data. This is the pooling factor of Gelman & Pardoe (2006) (its complement 1 − π_j is the shrinkage factor); the prior/likelihood balance it captures is the one Betancourt & Girolami (2015) tied to the optimal centred/non-centred parameterisation. The companion paper shows that once the Ehresmann connection on the base–fiber bundle is known to be flat, π_j is not a heuristic but the unique quantity governing the centring obstruction.

Coordinates with π_j near 1 are prior-dominated: their posteriors reflect regularisation toward the population mean rather than group-specific data. Coordinates with π_j near 0 are data-driven. The same quantity gives the per-group non-centring weight for partial reparameterisation: the closed-form optimal weight of Papaspiliopoulos, Roberts & Sköld (2003), given for generalised linear mixed models by Tan & Nott (2013).

The installed package is deliberately small: `prior_fraction()` (with a `brms` method) and the `smoothbp_advisor()` companion. The geometry that motivates it — the analytic connection, curvature, the chain-level coupling diagnostic, and the experimental connection-corrected samplers — is reproduction code for the paper and lives in the source repository's [`paper/`](https://github.com/ABindoff/fibr/tree/master/paper) directory rather than in the installed package, because the paper proves this connection is flat (so those quantities are statistical, not a user-facing geometric diagnostic).

## Installation

```r
# from CRAN (once available)
install.packages("fibr")

# or the development version
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

### Reproducing the paper

The analytic connection, curvature, chain-level coupling diagnostic, and the
connection-corrected samplers are not part of the installed package. They live in
[`paper/R/`](https://github.com/ABindoff/fibr/tree/master/paper/R) and are loaded with:

```r
source("paper/setup.R")   # run from the repository root
```

after which `compute_connection()`, `holonomy_diagnostic()`, and the samplers are
available, and the scripts in `data-raw/` regenerate the paper's figures. See
[`paper/README.md`](https://github.com/ABindoff/fibr/blob/master/paper/README.md). This route needs the heavier dependencies
(Matrix, FNN, deSolve, patchwork, a Stan backend) that the CRAN package does not.

## Background

The joint parameter space of a two-level hierarchical model — hyperparameters (μ, σ) as the base, group parameters (α₁, …, α_J) as the fibers — carries a fiber bundle structure. The Fisher information metric induces an Ehresmann connection A = −G_FF⁻¹ G_BF on this bundle. A natural conjecture is that the centring obstruction arises from the curvature of this connection, manifesting as holonomy when the hyperparameters traverse a closed loop.

The companion paper proves this false for any smooth hierarchical posterior: the connection is flat, its horizontal leaves being the level sets of the fiber score ∂_α log p, so no geometric obstruction exists above the metric level. In the log-concave case this foliation is the classical dual-flat orthogonal foliation of Amari (2001), the information-orthogonal reparameterisation of Cox & Reid (1987), and the gradient mapping of Hessian geometry (Shima 2013); the paper's contribution is to cast it as an Ehresmann connection valid for any smooth posterior (it is the connection that is flat, not the generally-curved Hessian metric). The non-zero curvature seen when integrating the *linearised* (fiber-frozen) connection is an artefact of that approximation; numerical verification is in `data-raw/verify_flat_connection.R` and `paper/tests/test-connection.R`.

What remains is statistical: the conditional dependence of fiber on base, governed per group by π_j. Genuine curvature — including rotational holonomy — does appear for connections built from a sampler's working metric (a fixed mass matrix) rather than the true Hessian, making holonomy an algorithmic rather than geometric phenomenon in that setting.

## Reference

Bindoff, A.D. (2026). *A Flat Connection: The Pooling Factor and the Geometry of Centring in Hierarchical MCMC.* Preprint. <https://doi.org/10.5281/zenodo.20724550>

```bibtex
@misc{bindoff2026fibr,
  author = {Bindoff, Aidan D.},
  title  = {A Flat Connection: The Pooling Factor and the Geometry of
             Centring in Hierarchical {MCMC}},
  year   = {2026},
  doi    = {10.5281/zenodo.20724550},
  url    = {https://github.com/ABindoff/fibr}
}
```
