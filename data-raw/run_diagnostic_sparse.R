## Holonomy diagnostic on the sparse-data chain.
## Run from package root: Rscript data-raw/run_diagnostic_sparse.R

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)

invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

library(posterior)
library(FNN)
library(ggplot2)

draws <- readRDS(file.path(pkg_root, "data-raw", "glmm_sparse_draws.rds"))
cat(sprintf("Sparse draws dim: %s\n\n", paste(dim(draws), collapse = " x ")))

base_vars  <- c("mu", "sigma")
fiber_vars <- paste0("alpha[", 1:8, "]")

# ── Run at two min_gap values ─────────────────────────────────────────────────
# min_gap=50: well past the ~3-iter autocorrelation time → fiber ~independent
# min_gap=3:  inside autocorrelation time → fiber retains memory of start

for (gap in c(50L, 3L)) {
  cat(sprintf(
    "\n══════════════════════════════════════════════\n"
  ))
  cat(sprintf("min_gap = %d\n", gap))
  cat(sprintf(
    "══════════════════════════════════════════════\n"
  ))

  hd <- holonomy_diagnostic(
    chain              = draws,
    base_vars          = base_vars,
    fiber_vars         = fiber_vars,
    epsilon            = NULL,
    n_bootstrap        = 200L,
    min_gap            = gap,
    k                  = 100L,
    max_loops          = 5000L,
    residualize_fiber  = TRUE
  )

  print(hd)

  tag <- sprintf("sparse_gap%02d", gap)

  p_eigen <- plot(hd, type = "eigenspectrum")
  ggplot2::ggsave(
    file.path(pkg_root, "data-raw", sprintf("holonomy_%s_eigenspectrum.png", tag)),
    plot = p_eigen, width = 6, height = 6, dpi = 150
  )
  cat(sprintf("Saved: holonomy_%s_eigenspectrum.png\n", tag))
}
