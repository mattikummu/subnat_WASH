# select_best_divisions.R
# ─────────────────────────────────────────────────────────────────────────────
# For each country, rank every available subnational boundary-source by:
#   1. all_covered      – TRUE when every (region, year) has at least one WASH value
#   2. spatial_coverage – fraction of unique polygons vs best boundary source
#   3. n_years          – number of distinct survey years
#   4. wash_coverage    – fraction of WASH indicator cells that are non-NA
#   5. n_regions        – number of distinct spatial units (finer = better tie-break)
#
# Then write geometry and data separately:
#   results/wash_best_division.gpkg         – best boundary geometries (keys only)
#   results/wash_best_division_data.csv     – WASH data for best boundary
#   results/wash_second_division.gpkg       – 2nd-best boundary (countries with ≥2)
#   results/wash_second_division_data.csv   – WASH data for 2nd-best boundary
#
# Join later using: iso3, region_name, year, boundary_source
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(countrycode)
  library(purrr)
  library(sf)
})

dir.create("results", showWarnings = FALSE)

source("config.R")
paths <- subnat_wash_paths()
dhs_dir <- paths$dhs
mics_dir <- paths$mics
gadm_path <- paths$gadm

# ── 1. Load matched WASH data ─────────────────────────────────────────────────
message("Loading results/JMP_regions_matched_data.csv ...")
df <- read_csv("results/JMP_regions_matched_data.csv", show_col_types = FALSE) |>
  mutate(.row = row_number())

key_cols <- c("iso3", "region_name", "year", "data_source", "survey_type",
              "survey_family", "boundary_source", "match_type")

wash_data_cols <- intersect(
  c("Improved", "Not.Improved", "Basic", "Limited", "Unimproved", "Surface.water"),
  names(df)
)
message(sprintf("  %d rows | %d countries | WASH cols: %s",
                nrow(df), n_distinct(df$iso3),
                paste(wash_data_cols, collapse = ", ")))

# ── 1a. Load geometry for polygon-based region counting ──────────────────────
# Different region name spellings (e.g., "Asante" / "Ashanti") can map to the
# same physical polygon.  Counting unique polygon centroids instead of unique
# name strings gives the true number of distinct regions per boundary source.
message("Loading geometry from results/JMP_regions_matched.gpkg ...")
wash_sv   <- vect("results/JMP_regions_matched.gpkg")
geom_df   <- as.data.frame(wash_sv) |> mutate(.row = row_number())
cent_crds <- crds(centroids(wash_sv))
geom_df$poly_id <- paste(round(cent_crds[, 1], 3), round(cent_crds[, 2], 3))

poly_per_source <- geom_df |>
  group_by(iso3, boundary_source) |>
  summarise(n_polys = n_distinct(poly_id), .groups = "drop")

max_polys <- poly_per_source |>
  group_by(iso3) |>
  summarise(max_polys = max(n_polys), .groups = "drop")

message(sprintf("  %d features | polygon counts computed", nrow(wash_sv)))

# ── 1a′. Merge boundary sources with identical polygon sets ──────────────────
# Different boundary files (e.g., DHS 2013 vs DHS 2017, or geo1_sr2010 vs
# geo1_sr2018) can produce the exact same set of admin polygons.  When they
# are kept as separate boundary_sources, each gets credit only for its own
# years — splitting what should be a single time series.  Detect identical
# polygon fingerprints (sorted centroids) and relabel to a canonical source
# so all years are counted together.
source_fingerprint <- geom_df |>
  group_by(iso3, boundary_source) |>
  summarise(poly_fp = paste(sort(unique(poly_id)), collapse = "|"), .groups = "drop")

# Within each (iso3, fingerprint) group, pick the source with the most rows
# in the data as the canonical name (proxy for "most years")
source_rows <- df |> count(iso3, boundary_source, name = "n_rows")

canon_map <- source_fingerprint |>
  left_join(source_rows, by = c("iso3", "boundary_source")) |>
  mutate(n_rows = coalesce(n_rows, 0L)) |>
  group_by(iso3, poly_fp) |>
  mutate(
    canon = boundary_source[which.max(n_rows)],
    n_in_group = n()
  ) |>
  ungroup()

merged_groups <- canon_map |> filter(n_in_group > 1)
if (nrow(merged_groups) > 0) {
  n_merged <- n_distinct(paste(merged_groups$iso3, merged_groups$poly_fp))
  n_sources_merged <- nrow(merged_groups) -
    n_distinct(paste(merged_groups$iso3, merged_groups$canon))
  message(sprintf("  Merging %d equivalent boundary sources into %d canonical groups",
                  n_sources_merged + n_distinct(paste(merged_groups$iso3, merged_groups$canon)),
                  n_merged))
  for (iso in unique(merged_groups$iso3)) {
    grps <- merged_groups |> filter(iso3 == iso)
    for (fp in unique(grps$poly_fp)) {
      g <- grps |> filter(poly_fp == fp)
      message(sprintf("    %s: %s → %s",
                      iso, paste(g$boundary_source, collapse = " + "), unique(g$canon)))
    }
  }
}

# Build remap lookup and apply to data + geometry
remap <- canon_map |> select(iso3, boundary_source, canon)
df <- df |>
  left_join(remap, by = c("iso3", "boundary_source")) |>
  mutate(boundary_source = coalesce(canon, boundary_source)) |>
  select(-canon)
geom_df <- geom_df |>
  left_join(remap, by = c("iso3", "boundary_source")) |>
  mutate(boundary_source = coalesce(canon, boundary_source)) |>
  select(-canon)

# Keep full copies (before spatial/year filtering) so output best/second
# divisions retain all admin areas even when some areas have no WASH data.
df_all   <- df
geom_all <- geom_df

# Recompute polygon counts after merging
poly_per_source <- geom_df |>
  group_by(iso3, boundary_source) |>
  summarise(n_polys = n_distinct(poly_id), .groups = "drop")

max_polys <- poly_per_source |>
  group_by(iso3) |>
  summarise(max_polys = max(n_polys), .groups = "drop")

# ── 1b. Filter out poor spatial matches ──────────────────────────────────────
# Drop boundary sources whose unique polygon count is <80% of the country's
# achievable maximum (e.g., KEN GADM matching only Nairobi, TJK GADM missing
# GBAO and DRS, GHA GADM missing Brong Ahafo).
message("Filtering poor spatial matches ...")

adequate_sources <- poly_per_source |>
  left_join(max_polys, by = "iso3") |>
  filter(max_polys < 3 |                              # tiny countries: keep all
           (n_polys >= 3 & n_polys >= 0.8 * max_polys))

n_sources_before <- n_distinct(paste(df$iso3, df$boundary_source))
df <- df |> semi_join(adequate_sources, by = c("iso3", "boundary_source"))
n_sources_after <- n_distinct(paste(df$iso3, df$boundary_source))
message(sprintf("  Dropped %d boundary sources with <80%% spatial coverage",
                n_sources_before - n_sources_after))

# ── 1c. Filter out incomplete years ──────────────────────────────────────────
# Different JMP surveys may use different admin-level divisions across years
# (e.g., MDG 2004 uses 6 provinces, 2009+ uses 22 regions). When a boundary
# matches one level, it gets only 1-2 stray hits from the other level.
# Drop years where a boundary source matched <80% of its max regions.
message("Filtering incomplete years ...")

year_region <- df |>
  group_by(iso3, boundary_source, year) |>
  summarise(n_regions_yr = n_distinct(region_name), .groups = "drop")

max_regions <- year_region |>
  group_by(iso3, boundary_source) |>
  summarise(max_regions = max(n_regions_yr), .groups = "drop")

complete_years <- year_region |>
  left_join(max_regions, by = c("iso3", "boundary_source")) |>
  mutate(coverage = n_regions_yr / max_regions) |>
  filter(coverage >= 0.8)

n_before <- nrow(df)
df <- df |> semi_join(complete_years, by = c("iso3", "boundary_source", "year"))
n_dropped <- n_before - nrow(df)
message(sprintf("  Dropped %d rows from incomplete years (%d → %d rows)",
                n_dropped, n_before, nrow(df)))

# ── 2. Score each (iso3, boundary_source) ────────────────────────────────────
message("Scoring divisions ...")

wash_data_matrix <- df |> select(all_of(wash_data_cols))

df_scored <- df |>
  mutate(
    row_coverage = rowMeans(!is.na(wash_data_matrix))
  )

division_scores <- df_scored |>
  group_by(iso3, boundary_source) |>
  summarise(
    n_years         = n_distinct(year),
    n_regions       = n_distinct(region_name),
    n_obs           = n(),
    wash_coverage   = mean(row_coverage),
    all_covered     = all(row_coverage > 0),
    survey_families = paste(sort(unique(survey_family)), collapse = "+"),
    years_range     = paste0(min(year), "–", max(year)),
    .groups         = "drop"
  ) |>
  # Add polygon-based spatial coverage
  left_join(poly_per_source, by = c("iso3", "boundary_source")) |>
  left_join(max_polys, by = "iso3") |>
  mutate(spatial_coverage = n_polys / max_polys) |>
  # Rank: completeness → full spatial coverage → most years → data coverage
  arrange(iso3, desc(all_covered), desc(spatial_coverage), desc(n_years),
          desc(wash_coverage), desc(n_regions)) |>
  group_by(iso3) |>
  mutate(division_rank = row_number()) |>
  ungroup()

# ── 3. Print ranking summary ──────────────────────────────────────────────────
message("\n── Division ranking summary ────────────────────────────────────────")
division_scores |>
  filter(division_rank <= 2) |>
  mutate(bs_short = if_else(
    nchar(boundary_source) > 45,
    paste0(substr(boundary_source, 1, 42), "..."),
    boundary_source
  )) |>
  select(iso3, division_rank, all_covered, spatial_coverage, n_years,
         n_regions, wash_coverage, years_range, survey_families, bs_short) |>
  print(n = 80)

message(sprintf(
  "\nCountries with 1  division: %d",
  sum(division_scores$division_rank == 1 &
      !(division_scores$iso3 %in%
          filter(division_scores, division_rank == 2)$iso3))
))
message(sprintf(
  "Countries with 2+ divisions: %d",
  n_distinct(filter(division_scores, division_rank == 2)$iso3)
))

# Save full ranking to CSV for reference
division_scores |>
  select(iso3, boundary_source, n_polys, max_polys, spatial_coverage,
         n_years, n_regions, n_obs,
         wash_coverage, all_covered, survey_families, years_range, division_rank) |>
  write_csv("results/division_ranking.csv")
message("Saved: results/division_ranking.csv")

load_boundary_sv <- function(path, family) {
  if (family != "DHS") {
    shp <- file.path(path, paste0(basename(path), ".shp"))
    return(tryCatch(vect(shp), error = function(e) NULL))
  }

  shp_paths <- if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    zfiles <- tryCatch(unzip(path, list = TRUE)$Name, error = function(e) character(0))
    matches <- zfiles[grepl("shps/sdr_subnational_boundaries[^/]*\\.shp$", zfiles)]
    paste0("/vsizip/", path, "/", matches)
  } else {
    list.files(file.path(path, "shps"),
               pattern = "^sdr_subnational_boundaries.*\\.shp$",
               full.names = TRUE)
  }

  if (length(shp_paths) == 0) {
    return(NULL)
  }

  candidates <- map(shp_paths, ~ tryCatch(vect(.x), error = function(e) NULL)) |>
    compact()

  if (length(candidates) == 0) {
    return(NULL)
  }

  candidates[[which.max(map_int(candidates, nrow))]]
}

load_complete_boundary <- function(iso3, boundary_source) {
  if (identical(boundary_source, "GADM_ADM1")) {
    sv <- tryCatch(vect(gadm_path, layer = "ADM_1"), error = function(e) NULL)
    if (is.null(sv)) return(NULL)
    sv <- sv[sv$GID_0 == iso3, ]
    if (nrow(sv) == 0) return(NULL)
    sv$iso3 <- iso3
    sv$region_name <- as.character(sv$NAME_1)
    sv$boundary_source <- boundary_source
    return(sv)
  }

  dhs_path <- file.path(dhs_dir, boundary_source)
  if (file.exists(dhs_path)) {
    sv <- load_boundary_sv(dhs_path, "DHS")
    if (is.null(sv) || nrow(sv) == 0) return(NULL)
    sv$iso3 <- iso3
    sv$region_name <- as.character(sv$DHSREGEN)
    sv$boundary_source <- boundary_source
    return(sv)
  }

  mics_path <- file.path(mics_dir, boundary_source)
  if (dir.exists(mics_path)) {
    sv <- load_boundary_sv(mics_path, "MICS")
    if (is.null(sv) || nrow(sv) == 0) return(NULL)
    sv$iso3 <- iso3
    sv$region_name <- as.character(sv$ADMIN_NAME)
    sv$boundary_source <- boundary_source
    return(sv)
  }

  NULL
}

align_to_schema <- function(sv, template_cols) {
  out_df <- as.data.frame(matrix(NA, nrow = nrow(sv), ncol = length(template_cols)))
  names(out_df) <- template_cols
  present <- intersect(template_cols, names(sv))
  if (length(present) > 0) {
    out_df[present] <- as.data.frame(sv)[present]
  }
  values(sv) <- out_df
  sv
}

# ── 4. Subset data and geometry ──────────────────────────────────────────────
best_keys   <- division_scores |> filter(division_rank == 1) |>
  select(iso3, boundary_source)
second_keys <- division_scores |> filter(division_rank == 2) |>
  select(iso3, boundary_source)

# Subset WASH data
best_data   <- df_all |> semi_join(best_keys, by = c("iso3", "boundary_source"))
second_data <- df_all |> semi_join(second_keys, by = c("iso3", "boundary_source"))

message("Building best-source geometry directly from full boundary sources ...")
best_schema_cols <- names(wash_sv)
best_geom_list <- vector("list", nrow(best_keys))

for (i in seq_len(nrow(best_keys))) {
  country <- best_keys$iso3[i]
  source_name <- best_keys$boundary_source[i]

  sv_full <- load_complete_boundary(country, source_name)

  if (is.null(sv_full)) {
    matched_fallback <- geom_all |>
      filter(iso3 == country, boundary_source == source_name) |>
      group_by(iso3, boundary_source, poly_id) |>
      slice_head(n = 1) |>
      ungroup()

    if (nrow(matched_fallback) == 0) {
      message(sprintf("  %s (%s): no raw source or matched fallback geometry found",
                      country, source_name))
      next
    }

    best_geom_list[[i]] <- st_as_sf(wash_sv[matched_fallback$.row, ])
    message(sprintf("  %s (%s): raw boundary not found; using matched fallback (%d regions)",
                    country, source_name, nrow(matched_fallback)))
    next
  }

  cent <- crds(centroids(sv_full))
  sv_full$poly_id <- paste(round(cent[, 1], 3), round(cent[, 2], 3))
  sv_full$iso3 <- country
  sv_full$boundary_source <- source_name

  sv_full <- sv_full[!duplicated(sv_full$poly_id), ]
  sv_full <- align_to_schema(sv_full, best_schema_cols)

  best_geom_list[[i]] <- st_as_sf(sv_full)

  message(sprintf("  %s (%s): loaded full geometry (%d regions)",
                  country, source_name, nrow(sv_full)))
}

best_geom_list <- compact(best_geom_list)

if (length(best_geom_list) == 0) {
  stop("No best-source geometries could be built")
}

best_sf <- bind_rows(best_geom_list)
best_sv <- vect(best_sf)
best_df <- as.data.frame(best_sv)
best_cent <- crds(centroids(best_sv))
best_dedup <- best_df |>
  mutate(poly_id = paste(round(best_cent[, 1], 3), round(best_cent[, 2], 3))) |>
  select(iso3, boundary_source, poly_id) |>
  distinct()

second_dedup <- geom_all |>
  semi_join(second_keys, by = c("iso3", "boundary_source")) |>
  group_by(iso3, boundary_source, poly_id) |>
  slice_head(n = 1) |>
  ungroup()
second_sv <- if (nrow(second_dedup) > 0) wash_sv[second_dedup$.row, ] else NULL

# ── 4b. Assign stable adm_id ─────────────────────────────────────────────────
# Encoding:
#   Best division:   iso3n * 1000 + row          (e.g. 4001 .. 4099 for AFG)
#   Second division: iso3n * 1000 + 500 + row  (e.g. 4501 .. 4599 for AFG)
#
# Rule: adm_id % 1000 < 500  =>  best division
#       adm_id % 1000 >= 500  =>  second division
message("Assigning adm_id ...")

make_adm_ids <- function(dedup_df, offset) {
  if (!"iso3" %in% names(dedup_df)) {
    stop("make_adm_ids(): required column 'iso3' is missing from dedup_df")
  }

  cntry_id_vec <- suppressWarnings(
    countrycode::countrycode(as.character(dedup_df$iso3), "iso3c", "iso3n")
  )

  dedup_df |>
    mutate(
      cntry_id = cntry_id_vec,
      cntry_id = as.integer(cntry_id),
      cntry_id = if_else(is.na(cntry_id), 900L + dense_rank(iso3), cntry_id)
    ) |>
    arrange(iso3, poly_id) |>
    group_by(iso3) |>
    mutate(adm_id = as.integer(cntry_id) * 1000L + offset + row_number()) |>
    ungroup() |>
    select(iso3, poly_id, adm_id)
}

best_lookup   <- make_adm_ids(best_dedup,   0L)
second_lookup <- if (nrow(second_dedup) > 0)
  make_adm_ids(second_dedup, 500L) else tibble(iso3 = character(), poly_id = character(), adm_id = integer())

adm_id_lookup <- bind_rows(
  best_lookup   |> mutate(division_rank = 1L),
  second_lookup |> mutate(division_rank = 2L)
)

best_dedup  <- best_dedup  |> left_join(best_lookup,   by = c("iso3", "poly_id"))
best_sv$adm_id <- best_dedup$adm_id

if (!is.null(second_sv)) {
  second_dedup <- second_dedup |> left_join(second_lookup, by = c("iso3", "poly_id"))
  second_sv$adm_id <- second_dedup$adm_id
}

# Add adm_id to tabular outputs via natural key join (robust to name variants)
# The CSV and GPKG have independent row numbering, so .row cannot be used as a
# join key. Instead join on (iso3, boundary_source, region_name, year,
# survey_family) to match survey observations to their polygons, then map
# polygons to adm_id via poly_id.
name_poly_map <- geom_all |>
  select(iso3, boundary_source, region_name, year, survey_family, poly_id) |>
  distinct()

best_id_map <- best_dedup |>
  select(iso3, boundary_source, poly_id, adm_id) |>
  distinct()

best_data <- best_data |>
  left_join(name_poly_map, by = c("iso3", "boundary_source", "region_name", "year", "survey_family")) |>
  left_join(best_id_map, by = c("iso3", "boundary_source", "poly_id")) |>
  select(-poly_id)

if (nrow(second_dedup) > 0) {
  second_id_map <- second_dedup |>
    select(iso3, boundary_source, poly_id, adm_id) |>
    distinct()

  second_data <- second_data |>
    left_join(name_poly_map, by = c("iso3", "boundary_source", "region_name", "year", "survey_family")) |>
    left_join(second_id_map, by = c("iso3", "boundary_source", "poly_id")) |>
    select(-poly_id)
}

write_csv(adm_id_lookup, "results/adm_id_lookup.csv")
message(sprintf("  %d best + %d second polygons assigned adm_id; saved results/adm_id_lookup.csv",
                nrow(best_lookup), nrow(second_lookup)))

# ── 5. Write outputs ──────────────────────────────────────────────────────────
message("\n── Writing output ──────────────────────────────────────────────────")

writeVector(best_sv, "results/wash_best_division.gpkg", overwrite = TRUE)
message(sprintf("Saved: results/wash_best_division.gpkg   (%d features, %d countries)",
                nrow(best_sv), n_distinct(values(best_sv)$iso3)))

if (!is.null(second_sv)) {
  writeVector(second_sv, "results/wash_second_division.gpkg", overwrite = TRUE)
  message(sprintf("Saved: results/wash_second_division.gpkg (%d features, %d countries)",
                  nrow(second_sv), n_distinct(values(second_sv)$iso3)))
}

write_csv(best_data |> select(-.row), "results/wash_best_division_data.csv")
message(sprintf("Saved: results/wash_best_division_data.csv   (%d rows)",
                nrow(best_data)))

write_csv(second_data |> select(-.row), "results/wash_second_division_data.csv")
message(sprintf("Saved: results/wash_second_division_data.csv (%d rows)",
                nrow(second_data)))

# ── 6. Summary stats ─────────────────────────────────────────────────────────
message("\n── Best-division breakdown ─────────────────────────────────────────")
best_data |>
  group_by(iso3) |>
  summarise(
    n_years   = n_distinct(year),
    n_regions = n_distinct(region_name),
    years     = paste0(min(year), "–", max(year)),
    .groups   = "drop"
  ) |>
  arrange(desc(n_years)) |>
  print(n = 30)
