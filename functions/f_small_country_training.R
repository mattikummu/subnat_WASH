suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(terra)
  library(countrycode)
})

# Countries excluded from being represented as a single ADM0 row in training.
# Rationale: in these countries, the within-country variance of the predictors
# (and of WASH outcomes) is so much larger than the between-country variance
# that the population-weighted ADM0 average is not a meaningful point in
# predictor space — it would mislead the model rather than inform it.
#
# The selection combines:
#   - very large land area (CAN, RUS, AUS, KAZ, DZA, SAU, LBY, MNG, SDN, ARG)
#   - very high population (CHN, IND, USA, IDN, PAK, NGA, MEX, BRA)
#   - large + heterogeneous (COD)
#
# When a country in this list HAS subnational training data, it stays in via
# the ADM1 path. When it doesn't, it is simply omitted from training and from
# the ADM0 prediction grid. Adjust the list to taste — e.g. drop CAN if
# Canadian ADM0 data is informative for adjacent ADM0 predictions.
haqi_large_countries <- function() {
  c(
    # Very large by area
    "CHN", "RUS", "CAN", "AUS", "ARG", "KAZ",
    "DZA", "COD", "SAU", "LBY", "MNG", "SDN", "BRA",
    # Very large by population (added)
    "IND", "USA", "IDN", "PAK", "NGA", "MEX"
  )
}

# Original function — kept for backward compatibility, in case other scripts
# call it. New script-05 logic should call get_adm0_eligible_countries()
# instead, which uses a different (data-driven) eligibility rule.
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

# New eligibility rule: any country with JMP ADM0 observations that is NOT
# already represented in the subnational training set, excluding the very
# large countries above.
#
# This is broader than get_small_countries_no_subnat() — it includes:
#   (a) small countries without any GADM ADM1 data (the original case), AND
#   (b) countries whose subnational data was REJECTED by script 04
#       (e.g. spatial coverage <80%, incomplete years, name-matching failures).
#
# Case (b) is the new addition you asked for: countries like high-income
# states whose subnational data didn't pass quality gates can re-enter
# training as ADM0 rows.
get_adm0_eligible_countries <- function(countries_with_subnat,
                                        adm0_obs_iso3,
                                        large_countries = haqi_large_countries()) {

  eligible <- setdiff(unique(adm0_obs_iso3), countries_with_subnat)
  eligible <- setdiff(eligible, large_countries)
  sort(unique(eligible))
}

# adm_id scheme for ADM0 country rows.
#
# Subnational scheme (set in script 04): iso3n * 1000 + offset + row_idx
#   - best division:    iso3n * 1000 + 1     ..  iso3n * 1000 + 499
#   - second division:  iso3n * 1000 + 501 ..  iso3n * 1000 + 999
#
# ADM0 scheme (this function): plain iso3n (3-digit integer).
#   - Range: roughly 4 .. 894 for assigned ISO codes
#   - Custom 9xx codes for non-standard iso3 strings
#   - Cannot collide with subnational scheme (smallest subnat id is
#     iso3n * 1000 + 1, which is at least 4001 for iso3n = 4)
#
# This guarantees each adm_id is globally unique and identifies its level
# unambiguously: id < 1000 ==> ADM0; id >= 1000 ==> subnational.
make_country_adm_id_lookup <- function(iso3_vec) {
  iso_tbl <- tibble(iso3 = sort(unique(iso3_vec)))

  iso_tbl <- iso_tbl |>
    mutate(
      iso3n = suppressWarnings(countrycode::countrycode(iso3, "iso3c", "iso3n")),
      iso3n = as.integer(iso3n)
    )

  na_idx <- which(is.na(iso_tbl$iso3n))
  if (length(na_idx) > 0) {
    if (length(na_idx) > 99L) {
      stop("Too many non-standard ISO3 country codes to map into ADM0 fallback range 901-999")
    }
    iso_tbl$iso3n[na_idx] <- 900L + seq_along(na_idx)
  }

  iso_tbl |>
    transmute(iso3, adm_id = as.integer(iso3n))
}
