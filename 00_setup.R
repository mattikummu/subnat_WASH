# 00_setup.R
# ─────────────────────────────────────────────────────────────────────────────
# Master setup script for the subnational drinking-water and water-quality pipeline.
#
# 1. Installs all required packages (if not already installed)
# 2. Optionally runs the full pipeline in sequence
#
# Pipelines:
#   Steps  1-6   Common setup
#   Steps  7-27  Drinking-water (scripts 06-17 + validation/KPI)
#   Steps 28-30  Water quality  (scripts 20-22)
#   Steps 31-33  Bivariate WASH × WQ / pathogen (scripts 31-32)
#
# Usage:
#   source("00_setup.R")           # install packages only
#   source("00_setup.R"); run_pipeline()   # install + run all
# ─────────────────────────────────────────────────────────────────────────────

# ── Required packages ─────────────────────────────────────────────────────────

packages <- c(
  # Data wrangling
  "tidyverse",      # dplyr, tidyr, purrr, stringr, readr, ggplot2
  "readxl",         # reading JMP Excel files
  "httr",           # downloading JMP data files

  # String matching
  "stringdist",     # fuzzy string matching for region name linking
  "stringi",        # Unicode/transliteration helpers

  # Spatial
  "terra",          # raster + vector operations
  "tidyterra",      # tidy interface to terra objects
  "sf",             # simple features for vector data
  "exactextractr",  # fast exact raster extraction (supports terra)

  # Country codes & helpers
  "countrycode",    # ISO3 ↔ numeric country codes
  "zoo",            # interpolation for ADM0 filling
  "plyr",           # round_any() helper

  # Modelling
  "tidymodels",     # modelling framework (includes rsample, recipes, tune, etc.)
  "xgboost",        # gradient boosting engine
  "vip",            # variable importance plots
  "shapviz",        # SHAP value computation and visualisation
  "isotone",        # isotonic regression for bias correction
  "future",         # parallel backend for multisession workers
  "furrr",          # purrr-style parallel mapping

  # Slopes & trends
  "mblm",           # Theil-Sen median-based linear models
  "broom",          # tidy model output

  # Mapping & visualisation
  "tmap",           # thematic map plotting
  "cols4all",       # colour palettes for maps
  "pals",           # additional palettes used in validation scripts
  "rnaturalearth",  # country border polygons
  "ggrepel",        # non-overlapping text labels
  "patchwork",      # combining ggplot panels

  # Geometry simplification helpers
  "qgisprocess",    # QGIS processing used by f_simplify_polyg()
  "geojsonio"       # TopoJSON/GeoJSON helpers used by f_simplify_polyg()
)

new_packages <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(new_packages) > 0) {
  message("Installing missing packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages)
} else {
  message("All packages already installed.")
}

# ── Pipeline runner ───────────────────────────────────────────────────────────

# All scripts in logical execution order.
# The step number (left) is used by run_pipeline(steps = ...) below.

pipeline_scripts <- c(
  # ── Common setup pipeline ────────────────────────────────────────────────
  "01_download_JMP_data.R",                #  1 – Download JMP subnational Excel files
  "02_extract_JMP_regions.R",              #  2 – Extract region-level WASH indicators
  "03_link_regions_gadm.R",                #  3 – Link JMP regions to GADM/DHS/MICS boundaries
  "04_select_best_divisions.R",            #  4 – Rank & select best admin division per country
  "05_extract_global_predictors.R",        #  5 – Global ADM0/ADM1 predictor extraction lookup
  "05b_build_adm1_shared_assets.R",        #  6 – Build shared ADM1 template + simplified geometry

  # ── Drinking-water pipeline ──────────────────────────────────────────────
  "06_prepare_xgboost_data.R",             #  7 – Prepare XGBoost training data (lookup-based)
  "07_train_xgboost_wash.R",               #  8 – Train & evaluate XGBoost ensembles
  "08_analyse_xgboost_importance.R",       #  9 – Variable importance analysis
  "09_prepare_adm0_data.R",                # 10 – Prepare ADM0 JMP water dataset (modern, 2000+)
  "09b_prepare_adm0_historical.R",         # 11 – Prepare ADM0 historical estimates (1990-2015)
  "10_fill_adm0_wash_with_xgboost.R",      # 12 – Fill ADM0 gaps with XGBoost (1990-2023)
  "10b_merge_adm0_historical_modern.R",    # 13 – Merge historical + modern ADM0 (delta-change)
  "11_fill_adm1_wash_with_xgboost.R",      # 14 – Fill ADM1 gaps with XGBoost (1990-2023)
  "12_rasterize_adm1_wash.R",              # 15 – Ratio-adjust to ADM0, normalise, rasterize
  "13_wash_slopes_maps.R",                 # 16 – Theil-Sen slopes + tmap maps
  "13b_wash_slopes_maps_adm0.R",           # 17 – ADM0 Theil-Sen slopes + tmap maps
  "14_wash_population_estimates.R",        # 18 – Population per WASH service level + trends
  "15_validate_adm0_wash.R",               # 19 – Validate ADM0 WASH outputs
  "15b_validate_adm1_second.R",            # 20 – Validate ADM1 WASH (second-stage)
  "15c_validate_adm0_our_vs_eawag.R",      # 21 – Validate ADM0 against Eawag benchmark
  "15d_validate_adm1_second_eawag.R",      # 22 – Validate ADM1 second-stage vs Eawag
  "16_validation_adm0_plots_kpi.R",        # 23 – ADM0 validation KPI plots
  "16b_validation_adm1_second_plots_kpi.R",# 24 – ADM1 second-stage validation KPI plots
  "16c_validation_adm0_our_vs_eawag_plots_kpi.R", # 25 – ADM0 vs Eawag KPI plots
  "16d_validation_adm1_second_eawag_plots_kpi.R", # 26 – ADM1 second-stage vs Eawag KPI plots
  "17_input_data_maps_wash.R",             # 27 – Input-data coverage maps (ADM0/ADM1)

  # ── Water-quality pipeline ───────────────────────────────────────────────
  "20_classify_pathogen.R",                # 28 – Classify DynQual water-quality rasters
  "21_extract_water_quality.R",            # 29 – Extract water-quality class rasters to ADM1
  "22_wq_slopes_maps.R",                   # 30 – WQ Theil-Sen slopes + tmap maps

  # ── Bivariate (WASH × WQ / pathogen) ─────────────────────────────────────
  "31_wash_wq_bivariate_matrix.R",         # 31 – Bivariate maps: WASH × water quality
  "31b_wash_pathogen_bivariate_matrix.R",  # 32 – Bivariate maps: WASH × pathogen quality
  "32_wq_wash_population_trends.R"         # 33 – Population trends by WQ × WASH class

)

run_pipeline <- function(
  steps = seq_along(pipeline_scripts),
  skip_download = FALSE
) {
  if (skip_download) steps <- steps[steps != 1]

  for (i in steps) {
    script <- pipeline_scripts[i]
    message(sprintf("\n══ Step %02d: %s ══════════════════════════════", i, script))
    source(script, echo = FALSE)
    message(sprintf("   Step %d complete.", i))
  }

  message("\n══ Pipeline complete. ══════════════════════════════════════════")
}

# ── Convenience wrappers ──────────────────────────────────────────────────────

run_common    <- function() run_pipeline(steps =  1:6)
run_water     <- function() run_pipeline(steps =  7:27)
run_wq        <- function() run_pipeline(steps = 28:30)
run_bivariate <- function() run_pipeline(steps = 31:33)

message("Run run_pipeline() to execute all steps, or source individual scripts.")
message("")
message("Convenience wrappers:")
message("  run_common()     – steps  1-6   (common setup)")
message("  run_water()      – steps  7-27  (drinking-water pipeline + validation/KPI)")
message("  run_wq()         – steps 28-30  (water-quality classification + extraction + maps)")
message("  run_bivariate()  – steps 31-33  (WASH × WQ bivariate analyses)")
message("")
message("Note: step 1 (download) can be skipped with run_pipeline(skip_download = TRUE)")
message("      if JMP Excel files are already present in data_in/JMP_subnational/")
message("Tip:  run_pipeline(steps = 8:27)  to re-run water from XGBoost training onward.")
message("      run_pipeline(steps = 1:6)   to run only shared setup.")
