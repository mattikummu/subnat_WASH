# 20_classify_pathogen.R
# ---------------------------------------------------------------------------
# Read DynQual annual-average rasters (1980–2019) for pathogen, salinity,
# and organic, apply threshold-based reclassification to a 1–5 class scale,
# and save each result as a multi-layer GeoTIFF (one layer per year).
#
# Classification thresholds:
#   pathogen (cfu/100ml): 0 | 1 | 10 | 100 | 1000 | Inf
#   salinity (mg/L)     : -Inf | 300 | 600 | 1900 | 1200 | Inf
#   organic  (mg O2/L)  : -Inf | 2 | 4 | 8 | 15 | Inf
#
# Inputs:
#   data_in/dynqual/pathogen_annualAvg_1980_2019.nc
#   data_in/dynqual/salinity_annualAvg_1980_2019.nc
#   data_in/dynqual/organic_annualAvg_1980_2019.nc
#
# Outputs:
#   results_gis/pathogen_class_annual_1980_2019.tif
#   results_gis/salinity_class_annual_1980_2019.tif
#   results_gis/organic_class_annual_1980_2019.tif
# ---------------------------------------------------------------------------

library(terra)

### 0. Indicator definitions ----
# Each entry: list(nc = <input path>, out = <output path>, breaks = <class boundaries>)

indicators <- list(
  salinity = list(
    nc = "data_in/dynqual/salinity_annualAvg_1980_2019.nc",
    out = "results_gis/salinity_class_annual_1980_2019.tif",
    breaks = c(-Inf, 300, 600, 900, 1200, Inf)
  ),
  organic = list(
    nc = "data_in/dynqual/organic_annualAvg_1980_2019.nc",
    out = "results_gis/organic_class_annual_1980_2019.tif",
    breaks = c(-Inf, 2, 4, 8, 15, Inf)
  ),
  pathogen = list(
    nc = "data_in/dynqual/pathogen_annualAvg_1980_2019.nc",
    out = "results_gis/pathogen_class_annual_1980_2019.tif",
    breaks = c(-Inf, 1, 10, 100, 1000, Inf)
  )
)

dir.create("results_gis", showWarnings = FALSE, recursive = TRUE)

### 1. Helper: build reclass table for terra::classify ----
#
# `breaks` defines the integer-class anchor points, e.g.:
#   pathogen : c(-Inf, 1, 10, 100, 1000, Inf)  → 0→1, 1→2, 10→3, 100→4, 1000→5
#   salinity : c(-Inf, 300, 600, 900, 1200, Inf) → 0→1, 600→2, …, 2100→5
#   organic  : c(-Inf, 2, 4, 8, 15, Inf)           → 0→1, 2→2, …, 10→5
#
# The lower -Inf break is replaced by `lower` (default 0) so that values
# interpolate linearly from `lower` (= class 1) up to the first finite break
# (= class 2), giving e.g. salinity 300 ≈ class 1.5.  A single catch-all row
# (-Inf, lower) keeps truly negative / missing values at class 1.
# The upper Inf interval maps to a single fixed class 5 row.
#
# With `right = FALSE` in classify(), bins are [from, to).

make_rcl <- function(breaks, class_values = seq_len(length(breaks) - 1),
                     steps = 100, lower = 0) {
  if (length(breaks) != length(class_values) + 1) {
    stop("`breaks` must have exactly one more element than `class_values`.")
  }

  rows <- list()

  # Replace leading -Inf with `lower` for interpolation; add a catch-all row
  eff_breaks <- breaks
  if (is.infinite(breaks[1]) && breaks[1] < 0) {
    eff_breaks[1] <- lower
    rows <- c(rows, list(matrix(c(-Inf, lower, class_values[1]), ncol = 3)))
  }

  n <- length(class_values)
  for (i in seq_len(n)) {
    a   <- eff_breaks[i]
    b   <- eff_breaks[i + 1]
    cv1 <- class_values[i]
    cv2 <- if (i < n) class_values[i + 1] else class_values[i]

    if (is.infinite(b)) {
      # Upper open-ended interval → single fixed row
      rows <- c(rows, list(matrix(c(a, b, cv1), ncol = 3)))
    } else {
      # Finite interval → split into `steps` sub-bins with linear interpolation
      sub_breaks <- seq(a, b, length.out = steps + 1)
      sub_vals   <- seq(cv1, cv2, length.out = steps + 1)[-(steps + 1)]
      rows <- c(rows, list(cbind(
        from    = sub_breaks[-length(sub_breaks)],
        to      = sub_breaks[-1],
        becomes = sub_vals
      )))
    }
  }

  rcl <- do.call(rbind, rows)
  colnames(rcl) <- c("from", "to", "becomes")
  rcl
}

### 2. Classify all indicators ----
ind_name = 'pathogen'
for (ind_name in names(indicators)) {
  ind <- indicators[[ind_name]]
  message(sprintf("\n== %s ==", ind_name))

  # Load
  message("  Loading raster...")
  r_raw <- rast(ind$nc)
  years <- as.integer(format(time(r_raw), "%Y"))
  message(sprintf("  Layers : %d  (%d–%d)", nlyr(r_raw), min(years), max(years)))
  message(sprintf("  CRS    : %s", crs(r_raw, describe = TRUE)$name))

  # Classify
  message("  Classifying...")
  rcl     <- make_rcl(ind$breaks)
  r_class <- classify(r_raw, rcl = rcl, right = FALSE, include.lowest = TRUE)

  # Restore metadata
  names(r_class) <- paste0(ind_name, "_class_", years)
  time(r_class)  <- time(r_raw)

  message(sprintf("  Value range: %.2f – %.2f  (step 0.1 within each class)",
                  global(r_class, "min", na.rm = TRUE)[[1]],
                  global(r_class, "max", na.rm = TRUE)[[1]]))

  # Save
  message("  Saving...")
  writeRaster(
    r_class,
    ind$out,
    gdal      = c("COMPRESS=LZW", "PREDICTOR=3"),
    overwrite = TRUE
  )
  message(sprintf("  Saved: %s", ind$out))
}

message("\nDone.")
