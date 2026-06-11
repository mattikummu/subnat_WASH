library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

in_dir <- "data_in/JMP_subnational"

# Target var1/var2 combinations to extract
target_vars <- tibble(
  var1 = c("Water",        "Water",        "Water ladder", "Water ladder",
           "Water ladder", "Water ladder"),
  var2 = c("Improved",     "Not Improved", "Basic",        "Limited",
           "Unimproved",   "Surface water")
)

read_country <- function(path) {
  ds <- read_excel(path, sheet = "data_summary")

  df <- ds |>
    semi_join(target_vars, by = c("var1", "var2")) |>
    select(iso3, year, data_source, var2,
           matches("^region_\\d+_name$"),
           matches("^region_\\d+$"))

  # Pivot region name columns to long
  names_long <- df |>
    select(iso3, year, data_source, var2, matches("^region_\\d+_name$")) |>
    pivot_longer(
      cols      = matches("^region_\\d+_name$"),
      names_to  = "region_num",
      values_to = "region_name",
      names_pattern = "region_(\\d+)_name"
    )

  # Pivot region value columns to long
  values_long <- df |>
    select(iso3, year, data_source, var2, matches("^region_\\d+$")) |>
    pivot_longer(
      cols      = matches("^region_\\d+$"),
      names_to  = "region_num",
      values_to = "value",
      names_pattern = "region_(\\d+)"
    )

  left_join(names_long, values_long, by = c("iso3", "year", "data_source", "var2", "region_num")) |>
    filter(!is.na(region_name)) |>
    select(-region_num) |>
    pivot_wider(names_from = var2, values_from = value) |>
    select(iso3, region_name, year, data_source,
           any_of(c("Improved", "Not Improved", "Basic",
                    "Limited", "Unimproved", "Surface water")))
}

# Process all files
files <- list.files(in_dir, pattern = "\\.xlsx$", full.names = TRUE)

all_regions <- map(files, \(f) {
  tryCatch(
    read_country(f),
    error = function(e) {
      warning("Failed: ", basename(f), " - ", conditionMessage(e))
      NULL
    }
  )
}, .progress = TRUE) |>
  list_rbind()

# Save results
dir.create("results", showWarnings = FALSE)
write.csv(all_regions, "results/JMP_regions_water.csv")

message("\nExtracted ", nrow(all_regions), " rows from ", length(files), " files")
message("Countries: ", n_distinct(all_regions$iso3))
message("Regions:   ", n_distinct(all_regions$region_name))
message("Years:     ", min(all_regions$year), "–", max(all_regions$year))
