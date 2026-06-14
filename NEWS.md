# fibr 0.1.0

First release, accompanying the arXiv preprint
"The Footprint of the Connection: Fiber Bundle Geometry and Conditional
Autocorrelation in Hierarchical MCMC" (Bindoff, 2026).

## New features

- `holonomy_diagnostic()` gains `structure = c("diagonal", "full")` argument.
  The default `"diagonal"` returns per-group real contraction factors `h_j`
  (group-aligned, interpretable directly against `pi_j`), replacing the
  full J×J complex transport matrix.  The old `"full"` behaviour is retained.
- `plot.fibr_holonomy()` gains `type = "contraction"` (default for diagonal
  fits), showing per-group `h_j` with bootstrap intervals.
- `compute_connection()` now errors loudly if the sigma column contains
  non-positive values (catches accidental log-sigma inputs).
- `as_smoothbp_re_fraction()` removed (orphaned by smoothbp 0.2.4 decoupling).

## Bug fixes

- `estimate_transport_map()` with `structure = "diagonal"` now preserves
  group identity so that `h_j` can be compared with per-group analytic
  quantities.

## Internal

- SBC validation battery added for the conditional-transport, reparam-HMC,
  and marginal samplers (M3).
- M4a relative-efficiency benchmark (1280-cell, 8-method) added under
  `data-raw/run_m4a.R`.
