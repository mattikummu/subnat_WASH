### Rasterize ADM1-level sanitation predictions
# Creates rasters and geopackages of the four sanitation indicators at
# ADM1 resolution for 1990-2023.
#
# Key steps:
#   1. Scale ADM1 values to match JMP ADM0 national totals
#   2. Normalise so s_bas + s_lim + s_unimp + s_od = 100
#   3. Read shared 5-arcmin ADM1 template (built in script 05b)
#   4. For each indicator x year, classify() template -> value raster
#
# Inputs:
#   results_sanit/sanit_adm1_global.gpkg                    - script 50
#   results_sanit/JMP_adm1_sanitation_interpPredExtrap.csv   - script 50
#   results_sanit/JMP_adm0_sanitation_merged_1990_2023.csv   - script 49b
#   results/global_adm1_predictors_1990_2023.csv             - shared ADM1 pop
#   GHS_POP raster (for 5-arcmin grid template)
#   results_gis_sanit/sanit_adm1_template_5arcmin.tif        - script 05b
#
# Outputs:
#   results_sanit/sanit_adm1_ratio_adjusted.csv
#   results_gis_sanit/sanit_adm1_values_1990_2023_long.csv
#   results_gis_sanit/sanit_adm1_values_1990_2023_wide.csv
#   results_gis_sanit/rast_sanit_{indicator}_1990_2023.tif

library(terra)
library(dplyr)
library(tidyr)
library(readr)

library(purrr)

### 0. Setup ----
years_in   <- 1990:2023
indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")
dir.create("results_gis_sanit", showWarnings = FALSE)

### 1. Load ADM1 data and population ----
message("== Loading inputs ==")

adm1_poly <- vect("results_sanit/sanit_adm1_global.gpkg")
adm1_poly$adm_id <- as.numeric(adm1_poly$adm_id)
message(sprintf("  ADM1 polygons: %d", nrow(adm1_poly)))

if (anyDuplicated(adm1_poly$adm_id) > 0) {
  stop("Duplicate adm_id found in sanit_adm1_global.gpkg; rasterization requires unique adm_id per polygon")
}

sanit_adm1 <- read_csv("results_sanit/JMP_adm1_sanitation_interpPredExtrap.csv",
                        show_col_types = FALSE) |>
  filter(year %in% years_in)
message(sprintf("  ADM1 data: %d rows", nrow(sanit_adm1)))

# Load ADM0 filled data (from script 49b - merged historical + modern)
sanit_adm0 <- read_csv("results_sanit/JMP_adm0_sanitation_merged_1990_2023.csv",
                        show_col_types = FALSE) |>
  filter(year %in% years_in)
message(sprintf("  ADM0 data: %d rows", nrow(sanit_adm0)))

# Load population from shared current ADM1 predictor table
pop_long <- read_csv("results/global_adm1_predictors_1990_2023.csv",
                      show_col_types = FALSE) |>
  filter(year %in% years_in) |>
  select(adm_id, year, pop)
message(sprintf("  Population data: %d adm_id x year rows", nrow(pop_long)))

### 2. Join population to ADM1 predictions ----
message("== Computing ratio adjustment ==")

adm1_pop <- sanit_adm1 |>
  left_join(pop_long, by = c("adm_id", "year"))

n_no_pop <- sum(is.na(adm1_pop$pop))
if (n_no_pop > 0) {
  message(sprintf("  WARNING: %d adm_id-year rows missing population", n_no_pop))
}

### 3. Pop-weighted ADM0 from ADM1 ----
adm0_from_adm1 <- adm1_pop |>
  filter(!is.na(pop)) |>
  group_by(iso3, year) |>
  summarise(
    s_bas_adm0_fromADM1   = sum(s_bas_filled   * pop, na.rm = TRUE) /
                            sum(pop[!is.na(s_bas_filled)],   na.rm = TRUE),
    s_lim_adm0_fromADM1   = sum(s_lim_filled   * pop, na.rm = TRUE) /
                            sum(pop[!is.na(s_lim_filled)],   na.rm = TRUE),
    s_unimp_adm0_fromADM1 = sum(s_unimp_filled * pop, na.rm = TRUE) /
                            sum(pop[!is.na(s_unimp_filled)], na.rm = TRUE),
    s_od_adm0_fromADM1    = sum(s_od_filled    * pop, na.rm = TRUE) /
                            sum(pop[!is.na(s_od_filled)],    na.rm = TRUE),
    .groups = "drop"
  )

message(sprintf("  Pop-weighted ADM0: %d country-years", nrow(adm0_from_adm1)))

### 4. Load JMP ADM0 filled values ----
jmp_adm0 <- sanit_adm0 |>
  select(iso3, year,
         s_bas_jmp_adm0   = s_bas_filled,
         s_lim_jmp_adm0   = s_lim_filled,
         s_unimp_jmp_adm0 = s_unimp_filled,
         s_od_jmp_adm0    = s_od_filled)

### 5. Apply scaling: adjusted_ADM1 = ADM1 × (JMP_ADM0 / pop-weighted ADM0) ----
# This preserves relative differences between regions while matching national totals

adm1_adjusted <- adm1_pop |>
  left_join(adm0_from_adm1, by = c("iso3", "year")) |>
  left_join(jmp_adm0, by = c("iso3", "year"))

for (ind in indicators) {
  filled_col <- paste0(ind, "_filled")
  fromADM1_col <- paste0(ind, "_adm0_fromADM1")
  jmp_col <- paste0(ind, "_jmp_adm0")

  scaling_factor <- if_else(
    !is.na(adm1_adjusted[[fromADM1_col]]) & adm1_adjusted[[fromADM1_col]] > 0,
    adm1_adjusted[[jmp_col]] / adm1_adjusted[[fromADM1_col]],
    1.0
  )

  adm1_adjusted[[paste0(ind, "_adjusted")]] <-
    adm1_adjusted[[filled_col]] * scaling_factor
}

### 7. Normalise so 4 indicators sum to 100 ----
message("== Normalising to sum = 100 ==")

adj_cols <- paste0(indicators, "_adjusted")

adm1_adjusted <- adm1_adjusted |>
  mutate(
    total_adj = s_bas_adjusted + s_lim_adjusted +
                s_unimp_adjusted + s_od_adjusted,
    across(all_of(adj_cols),
           ~ if_else(total_adj > 0, . / total_adj * 100, NA_real_),
           .names = "{.col}_norm")
  )

for (ind in indicators) {
  adm1_adjusted[[paste0(ind, "_final")]] <-
    adm1_adjusted[[paste0(ind, "_adjusted_norm")]]
}

# Fallback: if ratio adjustment yielded NA but original filled is available
adm1_adjusted <- adm1_adjusted |>
  mutate(
    has_final = !is.na(s_bas_final) & !is.na(s_lim_final) &
                !is.na(s_unimp_final) & !is.na(s_od_final),
    total_orig = s_bas_filled + s_lim_filled +
                 s_unimp_filled + s_od_filled
  )

for (ind in indicators) {
  final_col  <- paste0(ind, "_final")
  filled_col <- paste0(ind, "_filled")

  adm1_adjusted[[final_col]] <- if_else(
    !adm1_adjusted$has_final & !is.na(adm1_adjusted[[filled_col]]) &
      adm1_adjusted$total_orig > 0,
    adm1_adjusted[[filled_col]] / adm1_adjusted$total_orig * 100,
    adm1_adjusted[[final_col]]
  )
}

# Clamp to [0, 100]
for (ind in indicators) {
  final_col <- paste0(ind, "_final")
  adm1_adjusted[[final_col]] <- pmax(0, pmin(100, adm1_adjusted[[final_col]]))
}

# Verify
total_check <- adm1_adjusted |>
  filter(!is.na(s_bas_final)) |>
  mutate(total_final = s_bas_final + s_lim_final +
                       s_unimp_final + s_od_final)
message(sprintf("  Final sum range: %.2f - %.2f (mean %.2f)",
                min(total_check$total_final, na.rm = TRUE),
                max(total_check$total_final, na.rm = TRUE),
                mean(total_check$total_final, na.rm = TRUE)))

### 8. Save ratio-adjusted CSV ----
out_csv <- adm1_adjusted |>
  select(iso3, adm_id, region_name, year, pop,
         s_bas_original   = s_bas_filled,
         s_lim_original   = s_lim_filled,
         s_unimp_original = s_unimp_filled,
         s_od_original    = s_od_filled,
         s_bas   = s_bas_final,
         s_lim   = s_lim_final,
         s_unimp = s_unimp_final,
         s_od    = s_od_final)

write_csv(out_csv, "results_sanit/sanit_adm1_ratio_adjusted.csv")
message(sprintf("  Saved: results_sanit/sanit_adm1_ratio_adjusted.csv (%d rows)", nrow(out_csv)))

### 9. Load shared raster template (from script 05b) ----
message("== Loading shared 5-arcmin raster template ==")

template_path <- "results_gis_sanit/sanit_adm1_template_5arcmin.tif"
if (!file.exists(template_path)) {
  stop("Missing shared template raster. Run script 05b_build_adm1_shared_assets.R first.")
}

r_template <- rast(template_path)
message(sprintf("  Loaded: %s (%d x %d)", template_path, nrow(r_template), ncol(r_template)))

template_ids <- na.omit(unique(values(r_template)))
data_ids     <- unique(out_csv$adm_id)
id_overlap <- sum(data_ids %in% template_ids)
message(sprintf("  Template admIDs: %d | Data admIDs: %d | Overlap: %d",
                length(template_ids), length(data_ids), id_overlap))

if (id_overlap < 0.5 * length(data_ids)) {
  stop("Low adm_id overlap between data and template (<50%): likely adm_id coding mismatch")
}

### 10. Save attribute CSVs (long + wide) ----
message("== Saving attribute CSVs ==")

csv_long_path <- sprintf("results_gis_sanit/sanit_adm1_values_%d_%d_long.csv",
                         min(years_in), max(years_in))
csv_wide_path <- sprintf("results_gis_sanit/sanit_adm1_values_%d_%d_wide.csv",
                         min(years_in), max(years_in))

write_csv(out_csv, csv_long_path)

wide_list <- lapply(indicators, function(ind) {
  out_csv |>
    select(adm_id, year, value = all_of(ind)) |>
    pivot_wider(names_from = year, values_from = value,
                names_prefix = paste0(ind, "_"))
})

wide_all <- reduce(wide_list, left_join, by = "adm_id")
write_csv(wide_all, csv_wide_path)
message(sprintf("  Saved: %s", csv_long_path))
message(sprintf("  Saved: %s", csv_wide_path))
rm(wide_list, wide_all)
gc(verbose = FALSE)

### 11. Rasterization function ----

sanit_data2raster <- function(indicator, inData, inTemplate, inYears) {

  temp_data <- inData |>
    select(adm_id, year, value = all_of(indicator)) |>
    mutate(
      adm_id = as.numeric(adm_id),
      value  = case_when(
        is.na(value) ~ NA_real_,
        value < 0    ~ 0,
        value > 100  ~ 100,
        TRUE         ~ value
      )
    )

  used_ids    <- unique(temp_data$adm_id)
  all_ids     <- as.numeric(na.omit(unique(values(inTemplate))))
  missing_ids <- all_ids[!all_ids %in% used_ids]
  tmpl <- inTemplate
  if (length(missing_ids) > 0) tmpl[tmpl %in% missing_ids] <- NA

  coll_raster <- rast()
  for (yr in inYears) {
    yr_data   <- temp_data |> filter(year == yr)
    yr_raster <- classify(tmpl, cbind(yr_data$adm_id, yr_data$value))
    terra::add(coll_raster) <- yr_raster
  }

  names(coll_raster) <- paste0(indicator, "_", inYears)
  coll_raster
}

### 12. Rasterize all indicators ----
message("== Rasterizing indicators ==")

for (ind in indicators) {
  message(sprintf("  Rasterizing %s ...", ind))

  r_out <- sanit_data2raster(ind, out_csv, r_template, years_in)

  out_file <- sprintf("results_gis_sanit/rast_sanit_%s_%d_%d.tif",
                       ind, min(years_in), max(years_in))
  writeRaster(r_out, out_file, gdal = "COMPRESS=LZW", overwrite = TRUE)
  message(sprintf("    Saved: %s (%d layers)", out_file, nlyr(r_out)))
  rm(r_out)
  gc(verbose = FALSE)
}

### 13. Summary ----
message("\n== Summary ==")
message(sprintf("Indicators: %s", paste(indicators, collapse = ", ")))
message(sprintf("Years: %d-%d (%d layers each)", min(years_in), max(years_in),
                length(years_in)))
message(sprintf("ADM1 polygons: %d", nrow(adm1_poly)))
message("\nOutputs:")
message("  results_sanit/sanit_adm1_ratio_adjusted.csv")
message(sprintf("  %s", csv_long_path))
message(sprintf("  %s", csv_wide_path))
for (ind in indicators) {
  message(sprintf("  results_gis_sanit/rast_sanit_%s_%d_%d.tif", ind, min(years_in), max(years_in)))
}
message("\nDone.")
