setwd("C:/Users/bindoffa/antigravity_projects/fibr")
# Make fibr functions available without installing the package
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
library(posterior)
library(FNN)
library(ggplot2)
library(cmdstanr)

# Make pandoc available (RStudio bundle)
pandoc <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe"
if (file.exists(pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc))
  rmarkdown::find_pandoc(cache = FALSE)
}

cat("Rendering vignette...\n")
rmarkdown::render(
  "vignettes/glmm_example.Rmd",
  output_file = "glmm_example.html",
  output_dir  = "vignettes",
  knit_root_dir = getwd()
)
cat("Done. Output: vignettes/glmm_example.html\n")
