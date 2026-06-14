## Replace non-ASCII Unicode characters in R source files with ASCII equivalents.
## Run from the package root: Rscript data-raw/fix_ascii.R
##
## Characters replaced (all in comments or display strings):
##   U+2500 box drawing horizontal  -> -
##   U+2192 rightwards arrow        -> ->
##   U+2014 em dash                 -> --
##   U+2013 en dash                 -> -
##   U+2212 minus sign              -> -
##   U+00D7 multiplication sign     -> x

fix_ascii <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  lines <- gsub("─", "-",  lines)   # box light horizontal
  lines <- gsub("→", "->", lines)   # rightwards arrow
  lines <- gsub("—", "--", lines)   # em dash
  lines <- gsub("–", "-",  lines)   # en dash
  lines <- gsub("−", "-",  lines)   # minus sign
  lines <- gsub("×", "x",  lines)   # multiplication sign
  con <- file(path, open = "wb")
  writeLines(lines, con = con, sep = "\n")
  close(con)
  cat("Fixed:", path, "\n")
}

pkg_root <- normalizePath(
  if (nzchar(Sys.getenv("FIBR_ROOT"))) Sys.getenv("FIBR_ROOT") else "."
)
fix_ascii(file.path(pkg_root, "R", "holonomy.R"))
fix_ascii(file.path(pkg_root, "R", "plot_holonomy.R"))
fix_ascii(file.path(pkg_root, "R", "smoothbp_advisor.R"))
cat("Done.\n")
