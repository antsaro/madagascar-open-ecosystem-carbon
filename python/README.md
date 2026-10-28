# Python / Kaggle edition

`mada_carbon_nested_spatial_RF_GAM_LMM_v4_kaggle.ipynb` is a Python re-implementation of `scripts/03_nested_spatial_RF_GAM_LMM_manuscript_v4.qmd`, built to run on a Kaggle CPU notebook (scikit-learn, statsmodels, Optuna, rasterio).

It is a companion to the R workflow, not a replacement: the **R workflow is the reference used in the paper**. The GAM is a penalised B-spline model and the random forest uses scikit-learn with Optuna tuning, so numbers will not match the R results exactly. To reuse the spatial folds of an R run, set `MADA_RFOLDS` to a `folds_ALL_sites_combined.csv` file of fold assignments (see section 1 of the notebook).

## Paths

The data paths at the top of the notebook (section 1) point to Kaggle dataset locations. Edit them, or set the environment variables `MADA_CSV`, `MADA_COV` and `MADA_OUT`, before running. Use `MADA_FAST=1` for a quick smoke test.
