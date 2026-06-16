# CRAN submission comments — fibr 0.1.0

## Initial submission

## Test environments

- Windows 11 Enterprise (local), R 4.4.x
- win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes on examples

Examples that fit `brms` or `smoothbp` models are wrapped in `\dontrun{}`
because fitting times exceed the CRAN limit and both packages require
external samplers (Stan via cmdstanr or rstan). The manual-path examples
for `prior_fraction()` run unconditionally and cover the core arithmetic.

## Notes on suggested packages

- `brms` and `smoothbp` are in Suggests; the package degrades gracefully
  when they are absent (informative error messages via `requireNamespace()`).
- `cmdstanr` is in Suggests and is not on CRAN; it is only needed by the
  `brms`/`smoothbp` code paths.
- `testthat` is in Suggests; the test suite covers only the manual path
  (no external sampler required).
