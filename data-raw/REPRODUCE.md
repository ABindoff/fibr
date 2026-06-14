# Reproducibility capsule — fibr 0.1.0

Scripts and instructions to regenerate every figure and table in:

> Bindoff, A. D. (2026). *The Footprint of the Connection: Fiber Bundle
> Geometry and Conditional Autocorrelation in Hierarchical MCMC*.
> arXiv preprint.

All scripts are under `data-raw/` and are run from the **package root**
(`fibr/`). Seeds are set inside each script; run order matters where noted.

---

## Software versions used

| Software | Version |
|---|---|
| R | 4.6.0 (2026-04-24 ucrt) |
| posterior | 1.7.0 |
| cmdstanr | 0.9.0.9000 |
| CmdStan | 2.38.0 |
| ggplot2 | see `sessionInfo()` |
| fibr | 0.1.0 |

---

## Run order

### Step 0 — Install the package and compile Stan models

```r
# Install from source (or from GitHub release / Zenodo archive):
# install.packages("fibr_0.1.0.tar.gz", repos = NULL, type = "source")
# or:
# remotes::install_github("ABindoff/fibr@v0.1.0")

library(cmdstanr)
# Compile the four Stan models (one-time; ~30 s each):
cmdstan_model("inst/stan/glmm_centred.stan")
cmdstan_model("inst/stan/glmm_noncentred.stan")
cmdstan_model("inst/stan/glmm_hconnected.stan")
cmdstan_model("inst/stan/glmm_partial_nc.stan")
```

### Step 1 — Sparse GLMM benchmark (base draws; ~2 min)

Generates `data-raw/glmm_sparse_draws.rds` (centred) and
`data-raw/glmm_sparse_nc_draws.rds` (non-centred).
Required by Steps 2, 4, 5.

```bash
Rscript data-raw/simulate_glmm_sparse.R
Rscript data-raw/run_diagnostic_sparse.R   # centred fit
```

### Step 2 — Stokes validation figure (fig:stokes; ~2 min)

```bash
Rscript data-raw/run_connection.R
```
Output: `data-raw/holonomy_vs_radius.png`

### Step 3 — Simulation study (~4 h, checkpointed)

Generates per-cell draws in `data-raw/simstud_draws/` and all three
simulation study figures (fig:heatmap, fig:theory_vs_empirical,
fig:gap_profiles).

```bash
Rscript data-raw/run_simulation_study.R
```
Outputs: `data-raw/simstud_heatmap.png`, `data-raw/simstud_theory_vs_empirical.png`,
`data-raw/simstud_gap_profiles.png`

If interrupted, re-running resumes from the checkpoint
(`data-raw/simstud_results.rds`).

### Step 4 — Centred vs non-centred comparison (fig:comparison; ~5 min)

Requires centred and non-centred sparse draws from Step 1.

```bash
Rscript data-raw/run_comparison.R
```
Output: `data-raw/holonomy_comparison.png`

### Step 5 — Attribution analysis (Section 4.7 numbers; ~10 min)

```bash
Rscript data-raw/run_attribution.R
```
Outputs: `data-raw/attribution_results.rds`,
`data-raw/attribution_scatter.png`, `data-raw/attribution_area.png`,
`data-raw/attribution_scatter_full.png`

### Step 6 — Negative control (Section 4.7, control-pair numbers; ~15 min)

```bash
Rscript data-raw/run_control_pairs.R
```
Outputs: `data-raw/control_pairs.rds`, `data-raw/control_pairs.png`

### Step 7 — ESS table (tab:ess_comparison; ~5 min)

```bash
Rscript data-raw/run_ess_comparison.R
```
Output: `data-raw/ess_comparison.rds` (numbers in manuscript Table 3)

### Step 8 — Relative efficiency benchmark M4a (~8 h, checkpointed)

Optional for the paper (not referenced in any figure or table, but
underlies Section 7 discussion of sampler comparisons).

```bash
Rscript data-raw/run_m4a.R >> data-raw/m4a_run.log 2>&1
```
Outputs: `data-raw/m4a_results.rds`, `data-raw/m4a_summary.csv`,
`data-raw/m4a_rel_eff_pi.png`, `data-raw/m4a_rel_eff_kappa.png`

---

## Pre-computed outputs

The following files are already committed and do not require re-running
the scripts to compile the manuscript PDF:

- All five figure PNGs referenced in `manuscript/fibr_paper.tex`
- `data-raw/glmm_sparse_data.rds`, `data-raw/glmm_sparse_draws.rds`,
  `data-raw/glmm_sparse_nc_draws.rds`
- `data-raw/attribution_results.rds`, `data-raw/control_pairs.rds`,
  `data-raw/holonomy_comparison.rds`

---

## Compiling the manuscript PDF

```bash
cd manuscript
pdflatex fibr_paper.tex
bibtex fibr_paper
pdflatex fibr_paper.tex
pdflatex fibr_paper.tex
```

Expected output: `manuscript/fibr_paper.pdf` (25 pages, zero undefined
references).
