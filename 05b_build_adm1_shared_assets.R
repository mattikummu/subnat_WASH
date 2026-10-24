# 05b_build_adm1_shared_assets.R
# ---------------------------------------------------------------------------
# Build ADM1 static assets for drinking water:
#   1) 5-arcmin template raster (via 1-arcmin rasterize -> 5-arcmin modal)
#   2) simplified geometry-only GPKG
#
# Always (re)creates files on disk regardless of prior existence.
#
# Outputs:
#   results_gis/wash_adm1_template_5arcmin.tif
#   results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg
#   results_gis/polyg_wash_adm0_geomonly_simpl.gpkg
# ---------------------------------------------------------------------------

library(terra)
library(sf)

source("functions/f_simplify_polyg.R")
source("config.R")
paths <- subnat_wash_paths()

dir.create("results_gis", showWarnings = FALSE)

template_wash_path <- "results_gis/wash_adm1_template_5arcmin.tif"

geomonly_wash_simpl_path <- "results_gis/polyg_wash_adm1_1990_2023_geomonly_simpl.gpkg"
geomonly_adm0_wash_simpl_path <- "results_gis/polyg_wash_adm0_geomonly_simpl.gpkg"

gadm_path <- paths$gadm

# ── Load shared ADM1 geometry ──────────────────────────────────────────────
message("== Loading shared ADM1 geometry ==")
adm1_poly <- vect("results/wash_adm1_global.gpkg")
adm1_poly$adm_id <- as.numeric(adm1_poly$adm_id)

if (anyDuplicated(adm1_poly$adm_id) > 0) {
  stop("Duplicate adm_id found in results/wash_adm1_global.gpkg")
}

# ── Template raster ────────────────────────────────────────────────────────
message("== Template raster ==")
message("  Building via 1-arcmin rasterization -> 5-arcmin aggregation ...")

pop_path <- paths$population
r_ref <- rast(pop_path, lyrs = 1)
r_ref_1arcmin <- rast(
  nrows = nrow(r_ref) * 5,
  ncols = ncol(r_ref) * 5,
  ext   = ext(r_ref),
  crs   = crs(r_ref)
)

r_template_1arcmin <- rasterize(adm1_poly, r_ref_1arcmin, field = "adm_id")
r_template <- aggregate(r_template_1arcmin, fact = 5, fun = "modal", na.rm = TRUE)
names(r_template) <- "adm_id"

writeRaster(r_template, template_wash_path,  datatype = "INT4S", gdal = "COMPRESS=LZW", overwrite = TRUE)
message(sprintf("    Saved: %s", template_wash_path))
message(sprintf("  Template dims: %d x %d", nrow(r_template), ncol(r_template)))

# ── Simplified geometry-only GPKGs ─────────────────────────────────────────
message("== Simplified geometry-only GPKGs ==")

geom_only_sv <- adm1_poly[, "adm_id"]

gpkg_geomonly_tmp <- "results_gis/polyg_wash_adm1_1990_2023_geomonly.gpkg"
writeVector(geom_only_sv, gpkg_geomonly_tmp, overwrite = TRUE)

f_simplify_polyg(
  input_gpkg_file  = gpkg_geomonly_tmp,
  output_gpkg_file = geomonly_wash_simpl_path,
  simplify_percentage = 5
)
file.remove(gpkg_geomonly_tmp)

message(sprintf("    Saved: %s", geomonly_wash_simpl_path))

# ── Simplified ADM0 geometry-only GPKGs ────────────────────────────────────
message("== Simplified ADM0 geometry-only GPKGs ==")

lay_info <- st_layers(gadm_path)
adm0_layer <- lay_info$name[grepl("^ADM_0$|adm0|level0|level_0", lay_info$name, ignore.case = TRUE)][1]
if (is.na(adm0_layer)) stop("Could not identify ADM0 layer in: ", gadm_path)

adm0_sf <- st_read(gadm_path, layer = adm0_layer, quiet = TRUE)
adm0_sf$iso3 <- toupper(as.character(adm0_sf$GID_0))
adm0_sf$country_name <- as.character(adm0_sf$COUNTRY)
adm0_sf$adm0_id <- seq_len(nrow(adm0_sf))

adm0_sv <- vect(adm0_sf[, c("adm0_id", "iso3", "country_name")])

gpkg_adm0_geomonly_tmp <- "results_gis/polyg_wash_adm0_geomonly.gpkg"
writeVector(adm0_sv, gpkg_adm0_geomonly_tmp, overwrite = TRUE)

f_simplify_polyg(
  input_gpkg_file  = gpkg_adm0_geomonly_tmp,
  output_gpkg_file = geomonly_adm0_wash_simpl_path,
  simplify_percentage = 5
)
file.remove(gpkg_adm0_geomonly_tmp)

message(sprintf("    Saved: %s", geomonly_adm0_wash_simpl_path))

message("\nDone.")
