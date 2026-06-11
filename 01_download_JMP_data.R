library(readxl)
library(httr)

# Read download links
df <- read_excel("data_in/JMP_download_links.xlsx")

# Create output directory
out_dir <- "data_in/JMP_subnational"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Download each file
results <- lapply(seq_len(nrow(df)), function(i) {
  iso3 <- df$ISO3[i]
  country <- df$Country[i]
  url <- df$Link[i]

  out_file <- file.path(out_dir, paste0(iso3, ".xlsx"))

  # Skip if already downloaded
  if (file.exists(out_file)) {
    message(sprintf("[%d/%d] Skipping %s (%s) - already exists", i, nrow(df), country, iso3))
    return(list(iso3 = iso3, country = country, status = "skipped"))
  }

  message(sprintf("[%d/%d] Downloading %s (%s)...", i, nrow(df), country, iso3))

  resp <- tryCatch(
    GET(url, write_disk(out_file, overwrite = TRUE), timeout(60),
        add_headers("User-Agent" = "Mozilla/5.0")),
    error = function(e) e
  )

  if (inherits(resp, "error")) {
    message("  ERROR: ", conditionMessage(resp))
    return(list(iso3 = iso3, country = country, status = "error", message = conditionMessage(resp)))
  }

  if (http_error(resp)) {
    message("  HTTP error: ", status_code(resp))
    file.remove(out_file)
    return(list(iso3 = iso3, country = country, status = "http_error", code = status_code(resp)))
  }

  message("  Done (", file.size(out_file), " bytes)")
  list(iso3 = iso3, country = country, status = "ok")
})

# Summary
statuses <- sapply(results, `[[`, "status")
message("\n--- Download summary ---")
message("OK:         ", sum(statuses == "ok"))
message("Skipped:    ", sum(statuses == "skipped"))
message("HTTP error: ", sum(statuses == "http_error"))
message("Error:      ", sum(statuses == "error"))
