# 50b_merge_adm0_historical_modern_sanit.R
# ---------------------------------------------------------------------------
# Merge historical (1990-2015, from 48b) and modern (1990-2023, from 49)
# ADM0 sanitation data using delta-change harmonisation.
#
# The modern JMP_2025_WLD data is the authoritative source. Historical
# estimates use partly different reporting methods, so jumps may exist
# at the junction year (2000).
#
# Delta-change method (for countries WITH historical observations):
#   1. For each country & indicator, compute the mean difference between
#      modern_filled and hist_filled over the overlap window (2000-2002):
#        offset = mean(modern - hist, over 2000-2002)
#   2. Shift the entire historical series by that offset:
#        hist_adjusted = hist_filled + offset
#   3. Clamp to [0, 100]
#   4. For 1990-1999: use adjusted historical values
#   5. Proportionally rescale all four indicators to sum to 100
#
# For countries WITHOUT historical observations:
#   - Use 49's XGBoost-filled values for 1990-1999 directly
#
# For 2000 onwards: always use 49's filled values (the authoritative source)
#
# Inputs:
#   results_sanit/JMP_adm0_sanitation_interpPredExtrap_hist.csv   (from 48b)
#   results_sanit/JMP_adm0_sanitation_interpPredExtrap.csv        (from 49)
#   results_sanit/JMP_adm0_sanitation_historical.csv              (from 48b, raw observed)
#
# Outputs:
#   results_sanit/JMP_adm0_sanitation_merged_1990_YYYY.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

clamp_pct <- function(x) pmax(0, pmin(100, x))

# ---- 1. Read both datasets ----
message("== Reading modern and historical ADM0 sanitation data ==")
modern <- read_csv("results_sanit/JMP_adm0_sanitation_interpPredExtrap.csv",
                   show_col_types = FALSE)
hist   <- read_csv("results_sanit/JMP_adm0_sanitation_interpPredExtrap_hist.csv",
                   show_col_types = FALSE)

# Raw historical observations are needed to identify which countries
# actually have data — used to gate the delta-change adjustment below.
hist_obs <- read_csv("results_sanit/JMP_adm0_sanitation_historical.csv",
                     show_col_types = FALSE)

max_modern_year <- max(modern$year, na.rm = TRUE)
message(sprintf("  Modern: %d-%d (%d rows, %d countries)",
                min(modern$year), max_modern_year,
                nrow(modern), n_distinct(modern$iso3)))
message(sprintf("  Historical: %d-%d (%d rows, %d countries)",
                min(hist$year), max(hist$year),
                nrow(hist), n_distinct(hist$iso3)))

sanit_indicators <- c("s_bas", "s_lim", "s_unimp", "s_od")
overlap_years    <- 2000:2002

# Sanitation 'with historical observations' = country has at least one
# non-NA value of either anchor indicator (s_bas or s_od).
countries_with_hist <- hist_obs |>
  filter(!is.na(s_bas) | !is.na(s_od)) |>
  distinct(iso3) |>
  pull(iso3)

message(sprintf("  Countries WITH historical observations: %d",
                length(countries_with_hist)))
message(sprintf("  Countries WITHOUT (will use 49 xgboost fill): %d",
                n_distinct(modern$iso3) - length(countries_with_hist)))

# ---- 2. Build long-format tables with _filled values ----
pivot_filled <- function(df, suffix) {
  df |>
    select(iso3, year, country_name,
           ends_with("_filled"), ends_with("_fill_type")) |>
    pivot_longer(
      cols = -c(iso3, year, country_name),
      names_to = c("indicator", ".value"),
      names_pattern = "(.+)_(filled|fill_type)"
    )
}

modern_long <- pivot_filled(modern, "modern") |>
  rename(modern_filled = filled, modern_fill_type = fill_type)

hist_long <- pivot_filled(hist, "hist") |>
  rename(hist_filled = filled, hist_fill_type = fill_type)

# ---- 3. Calculate additive offset per country × indicator ----
message("== Computing delta-change offsets (overlap: 2000-2002) ==")

offsets <- modern_long |>
  filter(year %in% overlap_years, iso3 %in% countries_with_hist) |>
  inner_join(
    hist_long |> filter(year %in% overlap_years,
                        iso3 %in% countries_with_hist) |>
      select(iso3, year, indicator, hist_filled),
    by = c("iso3", "year", "indicator")
  ) |>
  filter(!is.na(modern_filled), !is.na(hist_filled)) |>
  summarise(
    offset = mean(modern_filled - hist_filled, na.rm = TRUE),
    n_overlap = n(),
    .by = c(iso3, indicator)
  )

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

# ---- 4. Adjust historical series ----
message("== Adjusting historical series (1990-1999) ==")

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
  )

hist_adjusted_wide <- hist_adjusted_wide |>
  rename(
    s_bas_adj       = s_bas_adjusted,
    s_lim_adj       = s_lim_adjusted,
    s_unimp_adj     = s_unimp_adjusted,
    s_od_adj        = s_od_adjusted,
    s_bas_ft        = s_bas_fill_type_out,
    s_lim_ft        = s_lim_fill_type_out,
    s_unimp_ft      = s_unimp_fill_type_out,
    s_od_ft         = s_od_fill_type_out
  )

hist_rescaled <- hist_adjusted_wide |>
  mutate(
    raw_sum = s_bas_adj + s_lim_adj + s_unimp_adj + s_od_adj,

    s_bas_final = case_when(
      is.na(s_bas_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(s_bas_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(s_bas_adj)
    ),
    s_lim_final = case_when(
      is.na(s_lim_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(s_lim_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(s_lim_adj)
    ),
    s_unimp_final = case_when(
      is.na(s_unimp_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(s_unimp_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(s_unimp_adj)
    ),
    s_od_final = case_when(
      is.na(s_od_adj) ~ NA_real_,
      !is.na(raw_sum) & raw_sum > 0 ~ clamp_pct(s_od_adj * 100 / raw_sum),
      TRUE ~ clamp_pct(s_od_adj)
    )
  )

# ---- 6. Assemble final output ----
message("== Assembling merged output ==")

hist_out <- hist_rescaled |>
  transmute(
    iso3, year, country_name,
    s_bas_filled      = s_bas_final,
    s_bas_fill_type   = s_bas_ft,
    s_lim_filled      = s_lim_final,
    s_lim_fill_type   = s_lim_ft,
    s_unimp_filled    = s_unimp_final,
    s_unimp_fill_type = s_unimp_ft,
    s_od_filled       = s_od_final,
    s_od_fill_type    = s_od_ft
  )

# B) Pre-2000 portion for countries WITHOUT historical data -> use 49's output
modern_pre2000 <- modern |>
  filter(year < 2000, !iso3 %in% countries_with_hist) |>
  select(
    iso3, year, country_name,
    s_bas_filled,   s_bas_fill_type,
    s_lim_filled,   s_lim_fill_type,
    s_unimp_filled, s_unimp_fill_type,
    s_od_filled,    s_od_fill_type
  )

message(sprintf("  Historical (delta-adjusted): %d rows (%d countries)",
                nrow(hist_out), n_distinct(hist_out$iso3)))
message(sprintf("  XGBoost fallback (pre-2000): %d rows (%d countries)",
                nrow(modern_pre2000), n_distinct(modern_pre2000$iso3)))

modern_out <- modern |>
  filter(year >= 2000) |>
  select(
    iso3, year, country_name,
    s_bas_filled,      s_bas_fill_type,
    s_lim_filled,      s_lim_fill_type,
    s_unimp_filled,    s_unimp_fill_type,
    s_od_filled,       s_od_fill_type
  )

merged <- bind_rows(hist_out, modern_pre2000, modern_out) |>
  arrange(iso3, year)

# ---- 7. Verify sum constraint ----
merged <- merged |>
  mutate(
    total = s_bas_filled + s_lim_filled + s_unimp_filled + s_od_filled
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

for (ind in sanit_indicators) {
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
out_path <- sprintf("results_sanit/JMP_adm0_sanitation_merged_1990_%d.csv", max_modern_year)
write_csv(merged, out_path)
message(sprintf("\nSaved: %s (%d rows, %d countries, years %d-%d)",
                out_path, nrow(merged), n_distinct(merged$iso3),
                min(merged$year), max(merged$year)))

message("\nDone.")
