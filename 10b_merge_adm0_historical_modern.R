# 10b_merge_adm0_historical_modern.R
# ---------------------------------------------------------------------------
# Merge historical (1990-2015, from 08b) and modern (1990-2023, from 09)
# ADM0 water data using delta-change harmonisation.
#
# The modern JMP_2025_WLD data is the authoritative source. Historical
# estimates use partly different reporting methods, so jumps may exist
# at the junction year (2000).
#
# Delta-change method (for countries WITH historical observations):
#   1. For each country & indicator, compute the mean difference between
#      modern_filled and hist_filled over the overlap window (2000-2002):
#        offset = mean(modern − hist,  over 2000-2002)
#   2. Shift the entire historical series by that offset:
#        hist_adjusted = hist_filled + offset
#   3. Clamp to [0, 100]
#   4. For 1990-1999: use adjusted historical values
#   5. Proportionally rescale all four indicators to sum to 100
#
# For countries WITHOUT historical observations:
#   - Use 09's XGBoost-filled values for 1990-1999 directly
#
# For 2000 onwards: always use 09's filled values (the authoritative source)
#
# Inputs:
#   results/JMP_adm0_water_interpPredExtrap_hist.csv   (from 08b)
#   results/JMP_adm0_water_interpPredExtrap.csv        (from 09, now 1990-2023)
#   results/JMP_adm0_water_historical.csv              (from 08b, raw observed)
#
# Outputs:
#   results/JMP_adm0_water_merged_1990_2023.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

clamp_pct <- function(x) pmax(0, pmin(100, x))

# ---- 1. Read datasets ----
message("== Reading modern and historical ADM0 data ==")
modern <- read_csv("results/JMP_adm0_water_interpPredExtrap.csv",
                   show_col_types = FALSE)
hist   <- read_csv("results/JMP_adm0_water_interpPredExtrap_hist.csv",
                   show_col_types = FALSE)

# Read raw historical observations to identify countries with real data
hist_obs <- read_csv("results/JMP_adm0_water_historical.csv",
                     show_col_types = FALSE)

max_modern_year <- max(modern$year, na.rm = TRUE)
message(sprintf("  Modern (09):     %d-%d (%d rows, %d countries)",
                min(modern$year), max_modern_year,
                nrow(modern), n_distinct(modern$iso3)))
message(sprintf("  Historical (08b): %d-%d (%d rows, %d countries)",
                min(hist$year), max(hist$year),
                nrow(hist), n_distinct(hist$iso3)))

# Countries with at least one historical observation (Basic or Surface.water)
countries_with_hist <- hist_obs |>
  filter(!is.na(Basic) | !is.na(Surface.water)) |>
  distinct(iso3) |>
  pull(iso3)

message(sprintf("  Countries WITH historical observations: %d",
                length(countries_with_hist)))
message(sprintf("  Countries WITHOUT (will use 09 xgboost fill): %d",
                n_distinct(modern$iso3) - length(countries_with_hist)))

wash_indicators <- c("Basic", "Limited", "Unimproved", "Surface.water")
overlap_years   <- 2000:2002

# ---- 2. Build long-format tables with _filled values ----
pivot_filled <- function(df) {
  df |>
    select(iso3, year, country_name,
           ends_with("_filled"), ends_with("_fill_type")) |>
    pivot_longer(
      cols = -c(iso3, year, country_name),
      names_to = c("indicator", ".value"),
      names_pattern = "(.+)_(filled|fill_type)"
    )
}

modern_long <- pivot_filled(modern) |>
  rename(modern_filled = filled, modern_fill_type = fill_type)

hist_long <- pivot_filled(hist) |>
  rename(hist_filled = filled, hist_fill_type = fill_type)

# ---- 3. Calculate additive offset per country × indicator ----
#      (only for countries WITH historical observations)
message("== Computing delta-change offsets (overlap: 2000-2002) ==")

offsets <- modern_long |>
  filter(year %in% overlap_years, iso3 %in% countries_with_hist) |>
  inner_join(
    hist_long |> filter(year %in% overlap_years, iso3 %in% countries_with_hist) |>
      select(iso3, year, indicator, hist_filled),
    by = c("iso3", "year", "indicator")
  ) |>
  filter(!is.na(modern_filled), !is.na(hist_filled)) |>
  summarise(
    offset = mean(modern_filled - hist_filled, na.rm = TRUE),
    n_overlap = n(),
    .by = c(iso3, indicator)
  )

# Summary of offsets
message("  Offset statistics (countries with historical data only):")
offset_stats <- offsets |>
  summarise(
    mean_abs_offset = mean(abs(offset)),
    median_abs_offset = median(abs(offset)),
    max_abs_offset = max(abs(offset)),
    .by = indicator
  )
for (i in seq_len(nrow(offset_stats))) {
  message(sprintf("    %s: mean |offset| = %.2f, median = %.2f, max = %.2f",
                  offset_stats$indicator[i],
                  offset_stats$mean_abs_offset[i],
                  offset_stats$median_abs_offset[i],
                  offset_stats$max_abs_offset[i]))
}

# ---- 4. Adjust historical series for countries WITH observations ----
message("== Adjusting historical series (1990-1999) for countries with observations ==")

hist_pre2000 <- hist_long |>
  filter(year < 2000, iso3 %in% countries_with_hist) |>
  left_join(offsets |> select(iso3, indicator, offset), by = c("iso3", "indicator")) |>
  mutate(
    offset = replace_na(offset, 0),
    adjusted = clamp_pct(hist_filled + offset),
    fill_type_out = case_when(
      is.na(hist_filled) ~ NA_character_,
      offset == 0        ~ hist_fill_type,
      TRUE               ~ paste0(hist_fill_type, "_delta_adj")
    )
  )

# ---- 5. Rescale so indicators sum to 100 ----
message("== Rescaling adjusted historical to maintain sum = 100 ==")

hist_adjusted_wide <- hist_pre2000 |>
  select(iso3, year, country_name, indicator, adjusted, fill_type_out) |>
  pivot_wider(
    names_from  = indicator,
    values_from = c(adjusted, fill_type_out),
    names_glue  = "{indicator}_{.value}"
  ) |>
  rename(
    Basic_adj         = Basic_adjusted,
    Limited_adj       = Limited_adjusted,
    Unimproved_adj    = Unimproved_adjusted,
    Surface.water_adj = Surface.water_adjusted,
    Basic_ft         = Basic_fill_type_out,
    Limited_ft       = Limited_fill_type_out,
    Unimproved_ft    = Unimproved_fill_type_out,
    Surface.water_ft = Surface.water_fill_type_out
  )

# Proportionally scale all four indicators to sum to 100
hist_rescaled <- hist_adjusted_wide |>
  mutate(
    raw_sum = Basic_adj + Limited_adj + Unimproved_adj + Surface.water_adj,

    Basic_final = case_when(
      is.na(Basic_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(Basic_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(Basic_adj)
    ),
    Limited_final = case_when(
      is.na(Limited_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(Limited_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(Limited_adj)
    ),
    Unimproved_final = case_when(
      is.na(Unimproved_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(Unimproved_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(Unimproved_adj)
    ),
    Surface.water_final = case_when(
      is.na(Surface.water_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(Surface.water_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(Surface.water_adj)
    )
  )

# ---- 6. Assemble final output ----
message("== Assembling merged output ==")

# A) Historical portion (1990-1999) for countries WITH historical data
hist_out <- hist_rescaled |>
  transmute(
    iso3, year, country_name,
    Basic_filled         = Basic_final,
    Basic_fill_type      = Basic_ft,
    Limited_filled       = Limited_final,
    Limited_fill_type    = Limited_ft,
    Unimproved_filled    = Unimproved_final,
    Unimproved_fill_type = Unimproved_ft,
    Surface.water_filled    = Surface.water_final,
    Surface.water_fill_type = Surface.water_ft
  )

# B) Pre-2000 portion for countries WITHOUT historical data → use 09's output
modern_pre2000 <- modern |>
  filter(year < 2000, !iso3 %in% countries_with_hist) |>
  select(
    iso3, year, country_name,
    Basic_filled,         Basic_fill_type,
    Limited_filled,       Limited_fill_type,
    Unimproved_filled,    Unimproved_fill_type,
    Surface.water_filled, Surface.water_fill_type
  )

message(sprintf("  Historical (delta-adjusted): %d rows (%d countries)",
                nrow(hist_out), n_distinct(hist_out$iso3)))
message(sprintf("  XGBoost fallback (pre-2000):  %d rows (%d countries)",
                nrow(modern_pre2000), n_distinct(modern_pre2000$iso3)))

# C) Modern portion (2000+) for ALL countries
modern_out <- modern |>
  filter(year >= 2000) |>
  select(
    iso3, year, country_name,
    Basic_filled,         Basic_fill_type,
    Limited_filled,       Limited_fill_type,
    Unimproved_filled,    Unimproved_fill_type,
    Surface.water_filled, Surface.water_fill_type
  )

merged <- bind_rows(hist_out, modern_pre2000, modern_out) |>
  arrange(iso3, year)

# ---- 7. Verify sum constraint ----
merged <- merged |>
  mutate(
    total = Basic_filled + Limited_filled + Unimproved_filled + Surface.water_filled
  )

n_complete <- sum(!is.na(merged$total))
n_bad <- sum(!near(merged$total, 100, tol = 0.5) & !is.na(merged$total))
n_na  <- sum(is.na(merged$total))
message(sprintf("  Complete rows: %d, NA rows: %d", n_complete, n_na))
if (n_bad > 0) {
  message(sprintf("WARNING: %d rows where sum != 100 (tol=0.5)", n_bad))
} else {
  message("  Sum constraint OK: all complete rows sum to ~100")
}

merged <- merged |> select(-total)

# ---- 8. Check for jumps at 2000 boundary ----
message("== Checking continuity at year 2000 ==")

for (ind in wash_indicators) {
  filled_col <- paste0(ind, "_filled")

  vals_99 <- merged |> filter(year == 1999) |> select(iso3, val_99 = all_of(filled_col))
  vals_00 <- merged |> filter(year == 2000) |> select(iso3, val_00 = all_of(filled_col))

  paired <- inner_join(vals_99, vals_00, by = "iso3") |>
    filter(!is.na(val_99), !is.na(val_00)) |>
    mutate(jump = abs(val_00 - val_99))

  if (nrow(paired) > 0) {
    message(sprintf("  %s: mean |jump| = %.2f, median = %.2f, max = %.2f (n=%d)",
                    ind, mean(paired$jump), median(paired$jump),
                    max(paired$jump), nrow(paired)))
  } else {
    message(sprintf("  %s: no paired 1999/2000 data", ind))
  }
}

# ---- 9. Save ----
out_path <- sprintf("results/JMP_adm0_water_merged_1990_%d.csv", max_modern_year)
write_csv(merged, out_path)
message(sprintf("\nSaved: %s (%d rows, %d countries, years %d-%d)",
                out_path, nrow(merged), n_distinct(merged$iso3),
                min(merged$year), max(merged$year)))

message("\nDone.")
