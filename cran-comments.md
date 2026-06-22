# CRAN submission comments — fibr 0.1.1

## Resubmission (addressing Konstanze Lauseker's comments)

- Added reference in `Description:` field: `Bindoff (2026) <doi:10.5281/zenodo.20724550>`.
- Added `Additional_repositories: https://stan-dev.r-universe.dev` for `cmdstanr`.
- Replaced `\dontrun{}` with `\donttest{}` in the `prior_fraction()` brms example;
  the block now uses a small synthetic dataset and is fully executable.

## Test environments

- Windows 11 Enterprise (local), R 4.6.0

## R CMD check results

0 errors | 0 warnings | 2 notes

The two notes are expected and acceptable:

1. **CRAN incoming feasibility**: `cmdstanr` is not on CRAN but is available via
   `Additional_repositories: https://stan-dev.r-universe.dev` (confirmed in check
   output: "Availability using Additional_repositories specification: cmdstanr yes").

2. **Top-level files**: `README.md` cannot be checked without pandoc (pandoc is
   available on CRAN's servers); the other flagged files (`fibr_arxiv_submission.zip`,
   `manuscript_bayesian_analysis`) are excluded from the built tarball via
   `.Rbuildignore` and will not appear in CRAN's check.

## Notes on examples

The `brms` example in `prior_fraction()` is wrapped in `\donttest{}` because
`brm()` compiles Stan code (~30 s). It uses a small synthetic dataset and runs
to completion. The manual-path example runs unconditionally and covers the core
arithmetic.

## Notes on suggested packages

- `brms` and `smoothbp` are in Suggests; the package degrades gracefully when
  they are absent (informative error via `requireNamespace()`).
- `cmdstanr` is in Suggests and available via `stan-dev.r-universe.dev`; it is
  only needed by the `brms`/`smoothbp` code paths.
- `testthat` is in Suggests; tests cover only the manual path.
