# paper/setup.R
# Load the geometry/sampler code that reproduces the paper's figures. These
# functions were part of the fibr package R/ during development; they were moved
# here because the paper proves the connection is flat, so they are reproduction
# apparatus rather than user-facing diagnostics, and the CRAN package ships only
# prior_fraction() and smoothbp_advisor().
#
# Run from the repository root:  source("paper/setup.R")
# Then library(fibr) provides prior_fraction(); the functions sourced here provide
# compute_connection(), holonomy_diagnostic(), synthetic_holonomy_loop(), the
# samplers, etc. Needs the heavier deps: Matrix, FNN, deSolve, patchwork, and a
# Stan backend (cmdstanr or rstan) for the model fits in data-raw/.

.fibr_paper_files <- list.files(
  file.path("paper", "R"), pattern = "[.][Rr]$", full.names = TRUE
)
if (length(.fibr_paper_files) == 0L)
  stop("No files found in paper/R/. Run this from the repository root.")
invisible(lapply(.fibr_paper_files, source))
message(sprintf("fibr paper apparatus loaded: %d files from paper/R/",
                length(.fibr_paper_files)))
