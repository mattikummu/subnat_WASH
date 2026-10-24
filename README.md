# Subnational Drinking Water & Water Quality

Pipeline for estimating subnational (ADM1-level) drinking-water and historical
water-quality indicators for all countries from 1990 to 2023, using XGBoost
models trained on JMP household survey data with population-weighted global predictors.

## Pipeline overview

| Step range | Pipeline | Description |
|---|---|---|
| 1-6 | **Common setup** | Download JMP data, extract regions, link to GADM/DHS/MICS boundaries, select best divisions, extract global predictors, build shared ADM1 template |
| 7-27 | **Drinking water** | Train XGBoost on water indicators, fill ADM0/ADM1 gaps, rasterize, Theil-Sen slopes, population estimates, validation & KPI plots |
| 28-30 | **Water quality** | Classify DynQual WQ rasters (pathogen, salinity, organic), extract to ADM1, Theil-Sen slope maps |
| 31-33 | **Bivariate (WASH x WQ)** | Bivariate maps from scripts 31, 31b, and 32; population trends by combined class |

## Repository structure

```
.
├── 00_setup.R                         # Setup: install packages, pipeline runner
├── 01_download_JMP_data.R
├── 02_extract_JMP_regions.R
├── 03_link_regions_gadm.R
├── 04_select_best_divisions.R
├── 05_extract_global_predictors.R
├── 05b_build_adm1_shared_assets.R
├── 06_prepare_xgboost_data.R
├── 07_train_xgboost_wash.R
├── 08_analyse_xgboost_importance.R
├── 09_prepare_adm0_data.R
├── 09b_prepare_adm0_historical.R
├── 10_fill_adm0_wash_with_xgboost.R
├── 10b_merge_adm0_historical_modern.R
├── 11_fill_adm1_wash_with_xgboost.R
├── 12_rasterize_adm1_wash.R
├── 13_wash_slopes_maps.R
├── 13b_wash_slopes_maps_adm0.R
├── 14_wash_population_estimates.R
├── 15_validate_adm0_wash.R
├── 15b_validate_adm1_second.R
├── 15c_validate_adm0_our_vs_eawag.R
├── 15d_validate_adm1_second_eawag.R
├── 16_validation_adm0_plots_kpi.R
├── 16b_validation_adm1_second_plots_kpi.R
├── 16c_validation_adm0_our_vs_eawag_plots_kpi.R
├── 16d_validation_adm1_second_eawag_plots_kpi.R
├── 17_input_data_maps_wash.R
├── 20_classify_pathogen.R
├── 21_extract_water_quality.R
├── 22_wq_slopes_maps.R
├── 31_wash_wq_bivariate_matrix.R
├── 31b_wash_pathogen_bivariate_matrix.R
├── 32_wq_wash_population_trends.R
├── functions/
│   ├── f_simplify_polyg.R
│   └── f_small_country_training.R
├── data_in/                           # Input data (see below)
├── results/                           # Generated outputs (not in repo)
├── figures_*/                         # Generated figures (not in repo)
├── methods_drinking_water.md
├── LICENSE
└── README.md
```

## Input data

### Included in this repository (`data_in/`)

| File/Directory | Description | Source |
|---|---|---|
| `JMP_2025_WLD.xlsx` | JMP 2025 World Leader Data | [washdata.org](https://washdata.org) |
| `JMP_download_links.xlsx` | JMP subnational data download URLs | [washdata.org](https://washdata.org) |
| `estimates_1980-2015_v1.xlsx` | Historical WASH estimates 1980-2015 | JMP |
| `DHS_bndr/` | DHS survey boundary shapefiles (53 zipped) | [dhsprogram.com](https://dhsprogram.com) |
| `MICS_bndr/` | MICS survey boundary shapefiles | [mics.unicef.org](https://mics.unicef.org) |
| `wq_annual_values/` | Water quality station annual values (4 CSVs) | GEMStat / national monitoring |
| `wq_annual_classes/` | Water quality station annual classes (4 CSVs) | GEMStat / national monitoring |
| `ihme_wat_access/` | IHME water access estimates (ADM0 & ADM1) | [IHME](https://ghdx.healthdata.org/) |
| `EawagPHH-mapping-.../` | Eawag safe drinking water mapping estimates (135 LMICs) | [GitHub](https://github.com/eawag-phh/mapping-safe-drinking-water-use-LMICs) |
| `r_pop_GHS_1985_2025_5arcmin.tif` | GHS-POP population raster (1985-2025, 5-arcmin) | [GHSL](https://human-settlement.emergency.copernicus.eu/) |

### Must be downloaded separately

These files are **not** included in the repository due to size (>100 MB) or
licensing restrictions. Place them in `data_in/` before running the pipeline.

#### 1. GADM administrative boundaries (v4.1)

- **File:** `gadm_410-levels.gpkg` (place in `data_in/` or update path in `05_extract_global_predictors.R`)
- **Download:** <https://gadm.org/download_world.html>
- **Size:** ~1 GB
- **Note:** Required for ADM0/ADM1 polygon geometry and region linking. Version 4.1 is required.

#### 2. JMP subnational Excel files

- **Files:** 47 country Excel files (downloaded automatically by script 01)
- **Downloaded by:** `01_download_JMP_data.R` (requires `JMP_download_links.xlsx`)
- **Destination:** `data_in/JMP_subnational/`

#### 3. DynQual water quality rasters

- **Files:**
  - `data_in/dynqual/pathogen_annualAvg_1980_2019.nc`
  - `data_in/dynqual/salinity_annualAvg_1980_2019.nc`
  - `data_in/dynqual/organic_annualAvg_1980_2019.nc`
- **Download:** <https://doi.org/10.34894/3Q0BVQ> (Utrecht University Dataverse)
- **Size:** ~8.2 GB total
- **Note:** Annual average rasters 1980-2019. Required for the water-quality and bivariate pipelines.

#### 4. Predictor rasters

All predictor rasters are population-weighted and extracted in script 05. They
come from companion projects or public data sources. Paths are hard-coded in
`05_extract_global_predictors.R` and must be updated to match your local setup.

| Raster | Variable | Source | Reference |
|---|---|---|---|
| `r_pop_GHS_1985_2025_5arcmin.tif` | Population | GHSL (included) | [GHSL](https://human-settlement.emergency.copernicus.eu/) |
| `travel_time_to_cities_11.tif` | Travel time to cities | MAP | Weiss et al. (2018), *Nature* |
| `aridity_index.tif` | Aridity index | CGIAR-CSI | Zomer et al. (2022), *Scientific Data* |
| `globgm-wtd-ss.tif` | Groundwater depth | GLOBGM | De Graaf et al. (2024), *Zenodo* |
| `rast_poverty_combined_income_weighted_1990_2022.tif` | Poverty | Companion project | Subnational poverty |
| `rast_msch_1990_2023.tif` | Mean schooling | Companion project | Subnational HDI |
| `rast_gnic_1990_2023.tif` | GNI per capita | Companion project | Subnational GNI |
| `rast_gini_disp_1990_2023.tif` | Gini coefficient | Companion project | Subnational inequality |
| `urbanisationCntryWise_GHS2023a.tif` | Urbanisation | GHSL-derived | GHS Urban Centre Database |
| `rast_gdi_1990_2023.tif` | Gender Development Index | Companion project | Subnational HDI |
| `rast_vdem_adm1_v2x_corr_1990_2024.tif` | Corruption | V-Dem-derived | Coppedge et al. (2024) |

## How to run

### 1. Install required packages

```r
source("00_setup.R")
```

### 2. Run the full pipeline

```r
source("00_setup.R")
run_pipeline()                  # all 33 steps
```

### 3. Run sub-pipelines

```r
source("00_setup.R")
run_common()                    # steps 1-6 (shared setup)
run_water()                     # steps 7-27 (drinking water)
run_wq()                        # steps 28-30 (water quality)
run_bivariate()                 # steps 31-33 (WASH x WQ)
```

### 4. Skip the JMP download (if already present)

```r
run_pipeline(skip_download = TRUE)
```

### 5. Customise predictor paths

Before running, edit the hard-coded file paths in `05_extract_global_predictors.R`
to point to your local copies of the predictor rasters (see table above).

## Output

All generated outputs are written to the following directories (not included
in the repository):

- `results/` - CSVs, GeoPackages, RDS files, model objects, extracted data
- `figures_drinking/` - Drinking water figures and maps
- `figures_bivariate/` - Bivariate WASH x WQ figures
- `figures_wq/` - Water quality figures and maps
- `figures_validation_wash/` - Drinking water validation figures

## License and citation

The code in this repository is licensed under the MIT License (see `LICENSE`).

Data sources are cited above. If you use the pipeline outputs, please cite both
this repository and the relevant underlying data sources (JMP, GHSL, DynQual,
GADM, IHME, Eawag, etc.).
