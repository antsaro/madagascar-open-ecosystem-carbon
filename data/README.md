# Data

The plot-level data and the covariate rasters are **not included** in this repository.
The full dataset is available from the authors on reasonable request (see the *Data availability* statement of the paper). Place your files here as described below and the scripts will find them (paths are resolved with `here::here()`).

```
data/
├── Datamapf_coordo.csv          # plot-level carbon pools, site, land cover, coordinates
└── covariates_30m/              # covariate rasters, one file per covariate and site
    ├── MNT_AF.tif   MNT_AT.tif   MNT_IV.tif
    ├── Slope_AF.tif ...
    └── ...
```

To keep the data elsewhere, set the environment variable `MADA_BASE` to the folder that contains `Datamapf_coordo.csv` and `covariates_30m/`.

## `Datamapf_coordo.csv`

Semicolon-separated (the modelling script also auto-detects `,` and tab). One row per plot (307 plots: 107 AF, 111 AT, 89 IV).

| Column | Description | Unit |
|---|---|---|
| `Name` | Unique plot identifier | |
| `Site` | `AF` (Ankarafantsika), `AT` (Ambohitantely), `IV` (Ivohiboro) | |
| `LU` | Land cover: `Forest`, `Shrubland`, `Grassland`, `Reforestation` | |
| `X`, `Y` | Longitude, latitude (WGS84, EPSG:4326) | degrees |
| `ST30` | Soil organic carbon, 0–30 cm (SOC 0–30) | Mg C ha⁻¹ |
| `ST100` | Soil organic carbon, 0–100 cm (SOC 0–100) | Mg C ha⁻¹ |
| `ST30.100` | Soil organic carbon, 30–100 cm | Mg C ha⁻¹ |
| `AGCGIEC` | Aboveground biomass carbon (AGC) | Mg C ha⁻¹ |
| `BGCGIEC` | Belowground biomass carbon (BGC) | Mg C ha⁻¹ |
| `DWGIEC` | Deadwood carbon | Mg C ha⁻¹ |
| `LITERGIEC` | Litter carbon | Mg C ha⁻¹ |
| `TOTC` | Total carbon (SOC 0–100 + AGC + BGC + deadwood + litter) | Mg C ha⁻¹ |
| `MNT`, `Slope`, `Edge`, `Roughness`, `TPI`, `TRI`, `MAP`, `MAT`, `Cosaspect`, `Sinaspect`, `FF` | Covariate values at the plot (used by scripts 01–02 only) | see below |

Scripts `01` and `02` read the covariates from the CSV. Script `03` ignores the covariate columns and **extracts the covariates from the rasters** at the plot locations.

## Covariate rasters (`covariates_30m/`)

13 single-band rasters per site, EPSG:4326, about 30 m (0.000278°), one shared grid per site. `S` is the site code (`AF`, `AT`, `IV`). Files are searched recursively, so they may sit in one folder or in `AF/`, `AT/`, `IV/` sub-folders.

| File | Model covariate | Type | Description |
|---|---|---|---|
| `MNT_<S>.tif` | Elevation | continuous | Elevation (SRTM) |
| `Slope_<S>.tif` | Slope | continuous | Terrain slope |
| `Edge_<S>.tif` | Edge | continuous | Distance to forest edge |
| `Roughness_<S>.tif` | Roughness | continuous | Terrain roughness |
| `TPI_<S>.tif` | TPI | continuous | Topographic position index |
| `TRI_<S>.tif` | TRI | continuous | Terrain ruggedness index |
| `Cosaspect_<S>.tif` | Cos_aspect | continuous | Cosine of aspect |
| `Sinaspect_<S>.tif` | Sin_aspect | continuous | Sine of aspect |
| `MAP_<S>.tif` | MAP | continuous | Mean annual precipitation |
| `MAT_<S>.tif` | MAT | continuous | Mean annual temperature |
| `ESA_<S>.tif` | LC | categorical | ESA WorldCover land-cover codes |
| `FF_<S>.tif` | Fire | categorical | Fire-frequency classes |
| `Soil_<S>.tif` | Soil | categorical | Soil-map codes |

## Do not commit data

`data/*` and `outputs/*` are git-ignored (except this file). Raw data and rasters are large and may be subject to access restrictions.
