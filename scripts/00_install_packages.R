# 00_install_packages.R
# Installs the R packages used by scripts 01-03. Run once from the project root.
# For a fully reproducible setup, follow with renv::init() and commit renv.lock.

pkgs <- c(
  # data handling and plotting
  "tidyverse", "here", "ragg", "systemfonts", "patchwork", "ggrepel", "viridis",
  # statistical tests (01)
  "car", "rstatix", "openxlsx", "flextable", "officer",
  # multivariate analyses (02)
  "vegan",
  # spatial data (03)
  "sf", "terra", "blockCV", "gstat", "matrixStats",
  # models (03)
  "ranger", "VSURF", "mgcv", "lme4", "lmerTest",
  # reporting
  "knitr", "rmarkdown"
)

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

message("Quarto (>= 1.4) is also required to render 03_*.qmd: https://quarto.org/docs/get-started/")
