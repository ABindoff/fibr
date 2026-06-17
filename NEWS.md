# fibr 0.1.1

- Update Zenodo DOI to `10.5281/zenodo.20724550` (version 0.1.1 archive).

# fibr 0.1.0

Initial CRAN release.

## New features

- `prior_fraction()` computes the per-group prior fraction (pooling or
  shrinkage factor) for hierarchical models. Accepts a `brmsfit` object
  directly, or a manual path supplying prior precision and likelihood
  information for any other fitting engine (Stan, JAGS, etc.).
  Supports gaussian, bernoulli, binomial, poisson, and negbinomial families.
  Has `print` and `plot` methods.

- `smoothbp_advisor()` reports the same Fisher information decomposition for
  changepoint random effects in `smoothbp` fits, and recommends centred vs.
  non-centred parameterisation per group. Has `print` and `plot` methods.

## Notes

The underlying geometry (analytic connection, curvature, chain-level coupling
diagnostic, and experimental samplers) is reproduction code for the companion
paper and lives in the source repository's `paper/` directory rather than in
the installed package.
