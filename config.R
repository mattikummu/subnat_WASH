# Project paths. Set SUBNAT_WASH_DATA_DIR to use data stored outside the repository.

subnat_wash_paths <- function() {
  data_dir <- normalizePath(
    Sys.getenv("SUBNAT_WASH_DATA_DIR", unset = "data_in"),
    mustWork = FALSE
  )
  data_path <- function(...) file.path(data_dir, ...)

  list(
    data_dir = data_dir,
    gadm = data_path("gadm_410-levels.gpkg"),
    population = data_path("r_pop_GHS_1985_2025_5arcmin.tif"),
    dhs = data_path("DHS_bndr"),
    mics = data_path("MICS_bndr"),
    predictors = list(
      poverty = data_path("predictors", "rast_poverty_combined_income_weighted_1990_2022.tif"),
      msch = data_path("predictors", "rast_msch_1990_2023.tif"),
      gnic = data_path("predictors", "rast_gnic_1990_2023.tif"),
      gini = data_path("predictors", "rast_gini_disp_1990_2023.tif"),
      urban = data_path("predictors", "urbanisationCntryWise_GHS2023a.tif"),
      travel_time = data_path("predictors", "travel_time_to_cities_11.tif"),
      aridity = data_path("predictors", "aridity_index.tif"),
      gdi = data_path("predictors", "rast_gdi_1990_2023.tif"),
      corruption = data_path("predictors", "rast_vdem_adm1_v2x_corr_1990_2024.tif"),
      gw_depth = data_path("predictors", "globgm-wtd-ss.tif")
    )
  )
}
