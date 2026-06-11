# 09_prepare_adm0_data.R
# ---------------------------------------------------------------------------
# Prepare ADM0 (country-level) JMP data for modelling.
#
# Input:
#   data_in/JMP_2025_WLD.xlsx  (sheet "wat")
#
# Output:
#   results/JMP_adm0_water.csv
#
# Notes:
#   - Uses the 'wat' sheet (machine-readable layout)
#   - Reads TOTAL service-level columns:
#       wat_basal_t, wat_lim_t, wat_unimp_t, wat_ns_t
#   - Keeps four indicators: Basic, Limited, Unimproved, Surface.water
#   - Constraint: Basic + Limited + Unimproved + Surface.water = 100
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(countrycode)
})

dir.create("results", showWarnings = FALSE)

xlsx_path <- "data_in/JMP_2025_WLD.xlsx"
if (!file.exists(xlsx_path)) {
  stop("Input file not found: ", xlsx_path)
}

message("== Reading ADM0 JMP workbook (sheet 'wat') ==")

# Row 1 = header names, row 2+ = data
# Columns: 1=name, 2=year, 15=wat_basal_t, 16=wat_lim_t,
#           17=wat_unimp_t, 18=wat_ns_t, 46=iso3
wat_raw <- read_excel(xlsx_path, sheet = "wat", skip = 1, col_names = FALSE)

# First row is the machine-readable header; rest is data
wat_data <- wat_raw[-1, ]

adm0_water <- tibble(
  country_name  = as.character(wat_data[[1]]),
  year          = as.numeric(wat_data[[2]]),
  iso3          = as.character(wat_data[[46]]),
  Basic         = as.numeric(wat_data[[15]]),
  Limited       = as.numeric(wat_data[[16]]),
  Unimproved    = as.numeric(wat_data[[17]]),
  Surface.water = as.numeric(wat_data[[18]])
) |>
  filter(!is.na(country_name), !is.na(year), !is.na(iso3)) |>
  filter(year >= 2000, year <= 2025)

# Verify sum constraint
adm0_water <- adm0_water |>
  mutate(total = Basic + Limited + Unimproved + Surface.water)

sum_ok <- all(near(adm0_water$total, 100, tol = 0.01), na.rm = TRUE)
if (!sum_ok) {
  n_bad <- sum(!near(adm0_water$total, 100, tol = 0.01), na.rm = TRUE)
  message(sprintf("WARNING: %d rows where Basic+Limited+Unimproved+Surface.water != 100", n_bad))
}

adm0_water <- adm0_water |>
  select(-total) |>
  relocate(iso3, country_name, year, Basic, Limited, Unimproved, Surface.water)

write_csv(adm0_water, "results/JMP_adm0_water.csv")

message(sprintf(
  "Saved: results/JMP_adm0_water.csv (%d rows, %d countries)",
  nrow(adm0_water), n_distinct(adm0_water$iso3, na.rm = TRUE)
))

message("Done.")
