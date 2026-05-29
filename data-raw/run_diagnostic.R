## End-to-end run of the fibr holonomy diagnostic against the GLMM chains.
## Run from the package root: Rscript data-raw/run_diagnostic.R
## Or interactively: source("data-raw/run_diagnostic.R")

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT")
  else "."
)

# Source all R files directly (works before devtools::load_all() is available)
invisible(lapply(
  list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE),
  source
))

library(posterior)
library(FNN)
library(ggplot2)

cat("\n── 1. Load draws ────────────────────────────────────────────────────────\n")
draws <- readRDS(file.path(pkg_root, "data-raw", "glmm_draws.rds"))
cat(sprintf("draws_array dim: %s\n", paste(dim(draws), collapse = " x ")))
cat(sprintf("  [iterations x chains x variables]\n"))

all_vars <- variables(draws)
cat(sprintf("Total variables: %d\n", length(all_vars)))

base_vars  <- c("mu", "sigma")
fiber_vars <- paste0("alpha[", 1:8, "]")

cat(sprintf("\nBase vars  : %s\n", paste(base_vars,  collapse = ", ")))
cat(sprintf("Fiber vars : %s\n", paste(fiber_vars, collapse = ", ")))

# Sanity-check names exist
stopifnot(all(base_vars  %in% all_vars))
stopifnot(all(fiber_vars %in% all_vars))

cat("\n── 2. Run holonomy diagnostic ───────────────────────────────────────────\n")
hd <- holonomy_diagnostic(
  chain       = draws,
  base_vars   = base_vars,
  fiber_vars  = fiber_vars,
  epsilon     = NULL,        # auto
  n_bootstrap = 200L,
  min_gap     = 50L,
  k           = 100L,
  max_loops   = 5000L
)

cat("\n── 3. Results ───────────────────────────────────────────────────────────\n")
print(hd)

cat("\n── 4. Save plots ────────────────────────────────────────────────────────\n")
out_dir <- file.path(pkg_root, "data-raw")

p_eigen <- plot(hd, type = "eigenspectrum")
ggsave(file.path(out_dir, "holonomy_eigenspectrum.png"),
       plot = p_eigen, width = 6, height = 6, dpi = 150)
cat("Saved: data-raw/holonomy_eigenspectrum.png\n")

p_loops <- plot(hd, type = "base_loops")
ggsave(file.path(out_dir, "holonomy_base_loops.png"),
       plot = p_loops, width = 7, height = 5, dpi = 150)
cat("Saved: data-raw/holonomy_base_loops.png\n")

cat("\nDone.\n")
