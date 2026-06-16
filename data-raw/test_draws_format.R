## Diagnostic: why does as_draws_matrix(fit$draws()) silently fail?
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

pkg_root <- normalizePath(".")
mod_nc <- cmdstan_model(file.path(pkg_root, "inst/stan/glmm_noncentred.stan"))

set.seed(42)
J <- 8L; n_j <- 5L; K <- 2L; N <- J * n_j
alpha_true <- rnorm(J, 0, 1)
group_id   <- rep(seq_len(J), each = n_j)
X          <- matrix(rnorm(N * K), ncol = K)
eta        <- alpha_true[group_id] + X %*% c(0.8, -0.5)
y          <- rbinom(N, 1L, plogis(eta))
stan_data  <- list(N = N, J = J, group = group_id, X = X, y = y)

cat("Running nc_nuts fit...\n")
fit <- mod_nc$sample(
  data = stan_data, chains = 2L, parallel_chains = 2L,
  iter_warmup = 200L, iter_sampling = 500L, seed = 42L, refresh = 0L
)

drws <- fit$draws()
cat("\nClass of fit$draws():", class(drws), "\n")
cat("Dim:", paste(dim(drws), collapse = " x "), "\n")
cat("Dimnames[[3]] (variables):\n")
print(dimnames(drws)[[3]])

cat("\nas_draws_matrix() attempt:\n")
mat <- tryCatch(
  as_draws_matrix(drws),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(mat)) {
  cat("Success — class:", class(mat), "  dim:", paste(dim(mat), collapse = " x "), "\n")
  cat("Colnames (first 10):", paste(head(colnames(mat), 10), collapse = ", "), "\n")
} else {
  cat("FAILED — trying fit$draws('matrix') directly:\n")
  mat2 <- tryCatch(
    fit$draws(format = "matrix"),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(mat2)) {
    cat("fit$draws('matrix') worked — class:", class(mat2),
        "  dim:", paste(dim(mat2), collapse = " x "), "\n")
    cat("Colnames (first 10):", paste(head(colnames(mat2), 10), collapse = ", "), "\n")
  }
}

## Also check the saved centred draws for comparison
draw_file <- file.path(pkg_root, "data-raw/simstud_draws/nj3_s0.5_r01.rds")
if (file.exists(draw_file)) {
  centred_draws <- readRDS(draw_file)
  cat("\nClass of centred_draws (from .rds):", class(centred_draws), "\n")
  mat_c <- tryCatch(as_draws_matrix(centred_draws),
                    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(mat_c))
    cat("as_draws_matrix(centred_draws) worked — dim:",
        paste(dim(mat_c), collapse = " x "), "\n")
}
