# Paper reproduction code

This directory holds the geometry and sampler code that reproduces the figures in

> Bindoff, A.D. (2026). *A Flat Connection: The Prior Fraction and the Geometry of
> Centring in Hierarchical MCMC.*

It is **not** part of the installed `fibr` package. The CRAN package ships only the
prior-fraction diagnostic (`prior_fraction()`, with a `brms` method) and the
`smoothbp_advisor()` companion. Everything here is the apparatus behind the paper's
central result, which is that for hierarchical regression the Fisher-metric
connection is *flat*: the analytic connection, its (zero) curvature, the
linearisation artifact, the chain-level base-fiber coupling diagnostic, and the
experimental connection-corrected samplers. Because the geometric holonomy is
identically zero for this model class, these are demonstration tools, not
user-facing diagnostics, so they live with the paper rather than in the package.

## Contents
- `R/` -- the functions: `compute_connection()`, `holonomy_diagnostic()`,
  `synthetic_holonomy_loop()`, curvature/transport, the GLMM Fisher metric and
  log-posterior helpers, and the samplers (`horizontal_hmc()`, `riemannian_mcmc()`,
  `asis_mcmc()`, `marginal_mcmc()`, `integrate_transport()`).
- `tests/` -- the unit tests for that code (the package keeps only the
  prior-fraction and smoothbp-advisor tests).
- `vignettes/` -- the original GLMM worked example.
- `setup.R` -- sources everything in `R/` into the global environment.

## How to reproduce
From the repository root:

```r
source("paper/setup.R")   # defines compute_connection(), holonomy_diagnostic(), ...
library(fibr)             # provides prior_fraction()
# then run the scripts in ../data-raw/, e.g.:
source("data-raw/verify_flat_connection.R")   # curvature ~ machine zero (the result)
source("data-raw/run_simulation_study.R")
```

The `data-raw/` reproduction scripts were written when these functions were
exported from the package; prepend `source("paper/setup.R")` so the now-internal
functions are available.

## Dependencies
Heavier than the CRAN package: `Matrix`, `FNN`, `deSolve`, `patchwork`, and a Stan
backend (`cmdstanr` or `rstan`) for the model fits.
