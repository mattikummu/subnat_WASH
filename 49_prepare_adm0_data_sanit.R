# 49_prepare_adm0_data_sanit.R
# ---------------------------------------------------------------------------
# Prepare ADM0 (country-level) JMP sanitation data for modelling.
#
# Input:
#   data_in/JMP_2025_WLD.xlsx  (sheet "san")
#
# Output:
#   results_sanit/JMP_adm0_sanitation.csv
#
# Notes:
#   - Uses the 'san' sheet (machine-readable layout)
#   - Reads TOTAL service-level columns:
#       san_basal_t (col 17), san_lim_t (col 18),
#       san_unimp_t (col 19), san_ns_t (col 20)
#       iso3 (col 52)
#   - Keeps four indicators: s_bas, s_lim, s_unimp, s_od
#   - Constraint: s_bas + s_lim + s_unimp + s_od = 100
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(countrycode)
})

dir.create("results_sanit", showWarnings = FALSE)

xlsx_path <- "data_in/JMP_2025_WLD.xlsx"
if (!file.exists(xlsx_path)) {
  stop("Input file not found: ", xlsx_path)
}

message("== Reading ADM0 JMP workbook (sheet 'san') ==")

san_raw <- read_excel(xlsx_path, sheet = "san", skip = 1, col_names = FALSE)

san_data <- san_raw[-1, ]

adm0_sanit <- tibble(
  country_name = as.character(san_data[[1]]),
  year         = as.numeric(san_data[[2]]),
  iso3         = as.character(san_data[[52]]),
  s_bas        = as.numeric(san_data[[17]]),
  s_lim        = as.numeric(san_data[[18]]),
  s_unimp      = as.numeric(san_data[[19]]),
  s_od         = as.numeric(san_data[[20]])
) |>
  filter(!is.na(country_name), !is.na(year), !is.na(iso3)) |>
  filter(year >= 2000, year <= 2025)

# Verify sum constraint
adm0_sanit <- adm0_sanit |>
  mutate(total = s_bas + s_lim + s_unimp + s_od)

sum_ok <- all(near(adm0_sanit$total, 100, tol = 0.01), na.rm = TRUE)
if (!sum_ok) {
  n_bad <- sum(!near(adm0_sanit$total, 100, tol = 0.01), na.rm = TRUE)
  message(sprintf("WARNING: %d rows where s_bas+s_lim+s_unimp+s_od != 100", n_bad))
}

adm0_sanit <- adm0_sanit |>
  select(-total) |>
  relocate(iso3, country_name, year, s_bas, s_lim, s_unimp, s_od)

write_csv(adm0_sanit, "results_sanit/JMP_adm0_sanitation.csv")

message(sprintf(
  "Saved: results_sanit/JMP_adm0_sanitation.csv (%d rows, %d countries)",
  nrow(adm0_sanit), n_distinct(adm0_sanit$iso3, na.rm = TRUE)
))

message("Done.")
