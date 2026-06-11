suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(terra)
})

haqi_large_countries <- function() {
  c("CHN", "RUS", "CAN", "AUS", "ARG", "KAZ",
    "DZA", "COD", "SAU", "LBY", "MNG", "SDN", "BRA")
}

get_small_countries_no_subnat <- function(countries_with_subnat,
                                          gadm_path,
                                          max_adm1 = 15L,
                                          large_countries = haqi_large_countries()) {
  gadm1 <- terra::vect(gadm_path, layer = "ADM_1")
  gadm1_df <- tibble(iso3 = as.character(gadm1$GID_0))

  adm1_counts <- gadm1_df |>
    filter(!is.na(iso3)) |>
    count(iso3, name = "n_adm1")

  adm1_counts |>
    filter(
      n_adm1 <= max_adm1,
      !iso3 %in% countries_with_subnat,
      !iso3 %in% large_countries
    ) |>
    pull(iso3) |>
    sort() |>
    unique()
}

make_country_adm_id_lookup <- function(iso3_vec) {
  iso_tbl <- tibble(iso3 = sort(unique(iso3_vec)))

  iso_tbl |>
    mutate(
      iso3n = suppressWarnings(countrycode::countrycode(iso3, "iso3c", "iso3n")),
      iso3n = as.integer(iso3n),
      iso3n = if_else(is.na(iso3n), 900L + dense_rank(iso3), iso3n),
      adm_id = as.integer(iso3n) * 100000L
    ) |>
    select(iso3, adm_id)
}
