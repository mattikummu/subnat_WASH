library(terra)
library(tidyterra)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(stringdist)
library(readr)
library(countrycode)

gadm_path <- "/Users/mkummu/R/GIS_data_common/gadm_410-levels.gpkg"
dhs_dir   <- "data_in/DHS_bndr"
mics_dir  <- "data_in/MICS_bndr"

# ── Helpers ───────────────────────────────────────────────────────────────────

normalise <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    str_to_lower() |>
    str_replace_all("[-\u2013\u2014]", " ") |>   # hyphens/dashes → space
    str_replace_all("[^a-z0-9 ]", "") |>
    str_squish()
}

# French → English translations for DHS region names in francophone countries
# (normalised form: lowercase, ASCII, no punctuation, squished)
fr_to_en <- c(
  "adamaoua"     = "adamawa",
  "est"          = "east",
  "extreme nord" = "far north",
  "nord ouest"   = "northwest",
  "nordouest"    = "northwest",
  "north west"   = "northwest",   # en hyphenated → solid
  "nord est"     = "northeast",
  "nordest"      = "northeast",
  "north east"   = "northeast",
  "ouest"        = "west",
  "sud"          = "south",
  "sud ouest"    = "southwest",
  "sudouest"     = "southwest",
  "south west"   = "southwest",   # en hyphenated → solid
  "sud est"      = "southeast",
  "sudest"       = "southeast",
  "south east"   = "southeast",
  "nord"         = "north",
  # Philippines ARMM / BARMM
  "autonomous region in muslim mindanao"                 = "armm",
  "bangsamoro autonomous region in muslim mindanao"      = "barmm",
  # Iraq: Arabic vs Kurdish governorate name
  "taamem"       = "kirkuk",
  "tameem"       = "kirkuk",
  "tamim"        = "kirkuk",
  # Iraq: transliteration variants
  "muthanna"     = "muthana",
  "nenava"       = "nainawa",
  "nainava"      = "nainawa",
  "ninawa"       = "nainawa",
  "nineva"       = "nainawa",
  "qadisiya"     = "qadisyah",
  "qadisyia"     = "qadisyah",
  "qadisya"      = "qadisyah",
  "missan"       = "misan",
  "suleimaniya"  = "sulaimaniya",
  "sulimania"    = "sulaimaniya",
  "sulaimania"   = "sulaimaniya",
  "sulaymanyia"  = "sulaimaniya",
  # KAZ: Kazakh names (after 'oblysy' suffix stripped by translate_fr_en)
  "almati city"         = "almaty city",   # typo almati → almaty (geo1_kz has 'almaty city')
  "karagandy"           = "karaganda",     # GADM: 'karaganda'
  "kostanai"            = "kostanay",      # geo1_kz: 'kostanay'
  "ontustik kazakhstan" = "south kazakhstan",  # GADM: 'south kazakhstan'
  "shygys kazakhstan"   = "east kazakhstan",   # GADM: 'east kazakhstan'
  # IDN: old/variant province names → canonical GADM Indonesian names
  "nanggroe aceh darussalam" = "aceh",
  "irian jaya"               = "papua",
  "cenrtal sulawesi"         = "sulawesi tengah",  # typo 'cenrtal'
  "central sulawesi"         = "sulawesi tengah",
  "north sumatra"            = "sumatera utara",
  "north sumatera"           = "sumatera utara",
  "south sumatra"            = "sumatera selatan",
  "south sumatera"           = "sumatera selatan",
  "west sumatra"             = "sumatera barat",
  "west sumatera"            = "sumatera barat",
  "west java"                = "jawa barat",
  "central java"             = "jawa tengah",
  "east java"                = "jawa timur",
  "west kalimantan"          = "kalimantan barat",
  "central kalimantan"       = "kalimantan tengah",
  "south kalimantan"         = "kalimantan selatan",
  "east kalimantan"          = "kalimantan timur",
  "north kalimantan"         = "kalimantan utara",
  "north sulawesi"           = "sulawesi utara",
  "south sulawesi"           = "sulawesi selatan",
  "southeast sulawesi"       = "sulawesi tenggara",
  "west sulawesi"            = "sulawesi barat",
  "north maluku"             = "maluku utara",
  "west nusa tenggara"       = "nusa tenggara barat",
  "east nusa tenggara"       = "nusa tenggara timur",
  "west papua"               = "papua barat",
  "dki jakarta"              = "jakarta raya",
  "jakarta"                  = "jakarta raya",
  "bangka belitung"          = "kep bangka belitung",
  "riau islands"             = "kepulauan riau",
  # General: strip '-ern' suffix from compass directions (boundary names like
  # 'eastern', 'northeastern' → 'east', 'northeast' to match JMP bare forms)
  "eastern"      = "east",
  "western"      = "west",
  "northern"     = "north",
  "southern"     = "south",
  "northeastern" = "northeast",
  "northwestern" = "northwest",
  "southeastern" = "southeast",
  "southwestern" = "southwest",
  # MKD: Macedonian transliterated region names → geo1_mk boundary names
  "pelagoniski"  = "pelagonia",
  "pelagonija"   = "pelagonia",
  "poloski"      = "polog",
  "skopski"      = "skopje",
  "vardarski"    = "vardar",
  # DZA: French compound names → geo1_dz boundary English names
  # (boundary has: south, northwest, north center, northeast,
  #  haut plateau east/center/west)
  "nord centre"           = "north center",
  "haut plateau centre"   = "haut plateau center",
  "haut plateau est"      = "haut plateau east",
  "haut plateau ouest"    = "haut plateau west",
  "hauts plateaux centre" = "haut plateau center",
  "hauts plateaux est"    = "haut plateau east",
  "hauts plateaux ouest"  = "haut plateau west",
  # COL: Spanish region names
  "pacifica"              = "pacific",
  "atlantica"             = "atlantic",
  "territorios nacionales"= "orinoquia and amazon",
  "orinoquiaamazonia"     = "orinoquia and amazon",
  # STP/Portuguese compass words
  "norte"="north", "sul"="south", "centro"="central",
  # UKR 2005 MICS: Crimea abbreviation and oblast name mismatch
  "ar krym" = "autonomous republic of crimea",
  # MDG DHS: variant spellings of Amoroni'i Mania region
  # "AMORON I MANIA" / "Amoron I Mania" → normalise → "amoron i mania"
  # "Anamoroni'i Mania" → normalise strips apostrophe → "anamoronii mania"
  "amoron i mania"  = "amoroni mania",
  "anamoronii mania"= "amoroni mania",
  # NPL 2019 MICS: named provinces (province 4 = Gandaki, 6 = Karnali, 7 = Sudurpashchim)
  "gandaki province"      = "province 4",
  "karnali province"      = "province 6",
  "sudoorpaschim province"= "province 7"
)

translate_fr_en <- function(x) {
  # Strip leading Roman numeral region codes (e.g. "i ", "iva ", "xii ")
  x <- str_remove(x, "^m{0,3}(cm|cd|d?c{0,3})(xc|xl|l?x{0,3})(ix|iv|v?i{0,3})[ab]? (?=[a-z])")
  # Strip trailing generic 'region' suffix (e.g. PHL "ilocos region" -> "ilocos")
  # Only strip when other words remain (avoid e.g. bare "region" becomes "")
  x <- if_else(str_detect(x, "^.+ region$"), str_remove(x, " region$"), x)
  # KAZ: strip Kazakh oblast suffix 'oblysy' (e.g. "akmola oblysy" -> "akmola")
  x <- str_remove(x, " oblysy$")
  # IDN: strip leading BPS province codes (e.g. "11 aceh" -> "aceh", "32 west java" -> "west java")
  x <- str_remove(x, "^\\d{2,3} (?=[a-z])")
  # GUY: normalise zero-padded region numbers (e.g. "region 01" -> "region 1")
  x <- str_replace(x, "^region 0(\\d)$", "region \\1")
  # DOM 2013: strip "region " prefix when followed only by a Roman numeral
  # (e.g. "region i" -> "i", "region vii" -> "vii")
  x <- str_remove(x, "^region (?=[ivxlcdm]+$)")
  # NPL 2019: strip 'no' from numbered provinces (e.g. "province no 1" -> "province 1")
  x <- str_replace(x, "^province no (\\d)$", "province \\1")
  # DZA: strip 'ept N ' programme-area prefix (e.g. "ept 4 hauts plateaux centre" -> "hauts plateaux centre")
  x <- str_remove(x, "^ept \\d+ ")
  # STP: strip Portuguese 'regiao' prefix (e.g. "regiao norte" -> "norte" -> "north")
  x <- str_remove(x, "^regiao (autonoma de |do |da |de |autónoma do )?")
  # Word-level compass -ern suffix stripping (e.g. "mid northern" -> "mid north", "south western" -> "south west")
  x <- str_replace_all(x, "\\beastern\\b",      "east")
  x <- str_replace_all(x, "\\bwestern\\b",      "west")
  x <- str_replace_all(x, "\\bnorthern\\b",     "north")
  x <- str_replace_all(x, "\\bsouthern\\b",     "south")
  x <- str_replace_all(x, "\\bnortheastern\\b", "northeast")
  x <- str_replace_all(x, "\\bnorthwestern\\b", "northwest")
  x <- str_replace_all(x, "\\bsoutheastern\\b", "southeast")
  x <- str_replace_all(x, "\\bsouthwestern\\b", "southwest")
  # SRB geo1_rs: strip "part of ..." parenthetical suffix from ADMIN_NAME
  x <- str_remove(x, "\\s*part of .*$")
  ifelse(x %in% names(fr_to_en), fr_to_en[x], x)
}

# Key columns for linking geometry to data (no spaces — safe for terra)
key_cols <- c("iso3", "region_name", "year", "data_source", "survey_type",
              "survey_family", "boundary_source", "match_type")

# WASH indicator columns (dot-separated, R-safe names)
wash_data_cols <- c("Improved", "Not.Improved", "Basic", "Limited",
                    "Unimproved", "Surface.water")

# ── 1. Load & clean JMP data ──────────────────────────────────────────────────

source("02_extract_JMP_regions.R")

# Normalize column names: "Not Improved" → "Not.Improved", "Surface water" → "Surface.water"
# This prevents terra::rbind() from creating duplicate suffixed columns
names(all_regions) <- make.names(names(all_regions))
message("WASH columns in all_regions: ",
        paste(intersect(wash_data_cols, names(all_regions)), collapse = ", "))

all_regions <- all_regions |>
  mutate(
    survey_type = str_extract(data_source, "[^_]+$") |>
      str_remove_all("\\d+$") |>
      str_to_upper(),
    survey_family = case_when(
      survey_type %in% c("DHS", "MIS", "AIS", "OTH") ~ "DHS",
      survey_type == "MICS"                            ~ "MICS",
      TRUE                                             ~ "OTHER"
    )
  ) |>
  filter(
    !str_detect(region_name, "^\\.\\."),
    !(iso3 == "EGY" & str_detect(region_name, regex("\\b(rural|urban)\\b", ignore_case = TRUE))),
    !str_detect(region_name, regex("^(rural|urban)$", ignore_case = TRUE))
  ) |>
  mutate(
    region_name = str_remove(region_name, "^\\[[A-Z0-9]+\\]\\s*"),
    region_name = str_remove(region_name, "\\s*\\([^)]+\\)\\s*$") |> str_squish()
  )

iso3_jmp <- unique(all_regions$iso3)

# ── 2. Build boundary catalogues ──────────────────────────────────────────────

# --- DHS: zip files + unzipped folders ---
# Filename pattern: "{Country Name} {Year} {TYPE} - sdr_subnational_boundaries_*.zip"
dhs_catalogue <- tibble(path = list.files(dhs_dir, full.names = TRUE)) |>
  mutate(
    fname         = basename(path),
    year          = str_match(fname, "\\b(\\d{4})\\b")[, 2] |> as.integer(),
    survey_type   = str_match(fname, "\\d{4} ([A-Z]+)")[, 2],
    country_raw   = str_match(fname, "^(.+?) \\d{4}")[, 2],
    iso3          = countrycode(country_raw, "country.name", "iso3c",
                                custom_match = c(
                                  "CAR"                       = "CAF",
                                  "Congo Democratic Republic" = "COD",
                                  "Cote d'Ivoire"             = "CIV",
                                  "The Gambia"                = "GMB",
                                  "Sao Tome and Principe"     = "STP"
                                )),
    survey_family = "DHS"
  ) |>
  filter(!is.na(iso3), !is.na(year))

# --- MICS: subfolders ---
# Folder pattern: "geo1_{iso2}", "geo1_{iso2}{year}", "geo1_{iso2}{year1}_{year2}"
mics_catalogue <- tibble(
  path   = list.dirs(mics_dir, recursive = FALSE, full.names = TRUE),
  folder = basename(list.dirs(mics_dir, recursive = FALSE, full.names = TRUE))
) |>
  filter(str_starts(folder, "geo1_") | str_starts(folder, "geo2_")) |>
  mutate(
    suffix    = str_remove(folder, "^geo[12]_"),
    iso2      = str_extract(suffix, "^[a-z]+") |> str_to_upper(),
    year_str  = str_remove(suffix, "^[a-z]+"),
    # year_from/year_to: NA means boundary applies to all years
    year_from = if_else(nchar(year_str) >= 4,
                        str_sub(year_str, 1, 4) |> as.integer(),
                        NA_integer_),
    year_to   = case_when(
      nchar(year_str) == 9 ~ str_sub(year_str, 6, 9) |> as.integer(),  # e.g. "2017_2023"
      !is.na(year_from)    ~ year_from,
      TRUE                 ~ NA_integer_
    ),
    iso3          = countrycode(iso2, "iso2c", "iso3c"),
    survey_family = "MICS"
  ) |>
  filter(!is.na(iso3))

message("DHS boundaries:  ", nrow(dhs_catalogue), " files for ",
        n_distinct(dhs_catalogue$iso3), " countries")
message("MICS boundaries: ", nrow(mics_catalogue), " folders for ",
        n_distinct(mics_catalogue$iso3), " countries")

# ── 3. Cache boundary region names ───────────────────────────────────────────

load_dhs_names <- function(path) {
  shp_paths <- if (str_ends(path, "\\.zip")) {
    # list all sdr_subnational_boundaries*.shp inside the zip
    zfiles <- tryCatch(unzip(path, list = TRUE)$Name, error = function(e) character(0))
    matches <- zfiles[str_detect(zfiles, "shps/sdr_subnational_boundaries[^/]*\\.shp$")]
    paste0("/vsizip/", path, "/", matches)
  } else {
    list.files(file.path(path, "shps"),
               pattern = "^sdr_subnational_boundaries.*\\.shp$",
               full.names = TRUE)
  }
  if (length(shp_paths) == 0) return(character(0))
  map(shp_paths, function(shp) {
    v <- tryCatch(vect(shp), error = function(e) NULL)
    if (is.null(v)) return(character(0))
    df <- as.data.frame(v)
    normalise(df$DHSREGEN[!is.na(df$DHSREGEN)]) |> translate_fr_en() |> unique()
  }) |> unlist() |> unique()
}

load_mics_names <- function(path) {
  shp <- file.path(path, paste0(basename(path), ".shp"))
  v <- tryCatch(vect(shp), error = function(e) NULL)
  if (is.null(v)) return(character(0))
  df <- as.data.frame(v)
  normalised <- normalise(df$ADMIN_NAME[!is.na(df$ADMIN_NAME)])
  translated <- translate_fr_en(normalised)
  unique(translated)
}

message("Caching boundary region names...")
dhs_name_cache  <- set_names(map(dhs_catalogue$path,  load_dhs_names),  dhs_catalogue$path)
mics_name_cache <- set_names(map(mics_catalogue$path, load_mics_names), mics_catalogue$path)

# ── 4. Select best boundary per (iso3, survey_family, jmp_year) ───────────────
#
# For each (iso3, survey_family, year) in JMP:
#   - Score all candidate boundary files: fraction of JMP region names matched
#   - Pick highest score; tie-break by year proximity
#   - MICS boundaries without year info: year distance = 0 (applies universally)

best_boundary <- function(iso, family, yr) {
  empty <- tibble(best_path = NA_character_, best_score = 0,
                  n_match = 0L, n_total = 0L)

  score_cands <- function(cands, name_cache, yr_dist_fn) {
    if (nrow(cands) == 0) return(NULL)
    cands |>
      mutate(
        bnames      = map(path, ~ name_cache[[.x]] %||% character(0)),
        n_match     = map_int(bnames, ~ as.integer(sum(jmp_norms %in% .x))),
        n_total     = length(jmp_norms),
        score_exact = n_match / n_total,
        score_fuzzy = map_dbl(bnames, function(bn) {
          if (length(bn) == 0) return(0)
          mean(map_dbl(jmp_norms, function(jn) {
            1 - min(stringdist(jn, bn, method = "jw"))
          }))
        }),
        score   = score_exact + 0.01 * score_fuzzy,
        yr_dist = map_int(seq_len(n()), yr_dist_fn)
      ) |>
      arrange(desc(score), yr_dist) |>
      slice_head(n = 1)
  }

  if (family == "DHS") {
    cands_primary      <- dhs_catalogue  |> filter(iso3 == iso)
    cache_primary      <- dhs_name_cache
    yr_dist_primary    <- function(i) abs(cands_primary$year[i] - yr)
    cands_fallback     <- mics_catalogue |> filter(iso3 == iso)
    cache_fallback     <- mics_name_cache
    yr_dist_fallback   <- function(i) {
      from <- cands_fallback$year_from[i]; to <- cands_fallback$year_to[i]
      if (is.na(from)) return(0L)
      if (yr >= from && yr <= to) return(0L)
      as.integer(min(abs(yr - from), abs(yr - to)))
    }
  } else if (family == "MICS") {
    cands_primary      <- mics_catalogue |> filter(iso3 == iso)
    cache_primary      <- mics_name_cache
    yr_dist_primary    <- function(i) {
      from <- cands_primary$year_from[i]; to <- cands_primary$year_to[i]
      if (is.na(from))            return(0L)
      if (yr >= from && yr <= to) return(0L)
      as.integer(min(abs(yr - from), abs(yr - to)))
    }
    cands_fallback     <- dhs_catalogue  |> filter(iso3 == iso)
    cache_fallback     <- dhs_name_cache
    yr_dist_fallback   <- function(i) abs(cands_fallback$year[i] - yr)
  } else {
    return(empty)
  }

  if (nrow(cands_primary) == 0 && nrow(cands_fallback) == 0) return(empty)

  jmp_norms <- all_regions |>
    filter(iso3 == iso, survey_family == family, year == yr) |>
    pull(region_name) |> normalise() |> translate_fr_en() |> unique()

  if (length(jmp_norms) == 0) return(empty)

  best_primary  <- score_cands(cands_primary,  cache_primary,  yr_dist_primary)
  best_fallback <- score_cands(cands_fallback, cache_fallback, yr_dist_fallback)

  # Use cross-family fallback only when primary has poor coverage (< 50%)
  scored <- if (!is.null(best_primary) && best_primary$score_exact >= 0.5) {
    best_primary
  } else if (!is.null(best_fallback) &&
             (is.null(best_primary) || best_fallback$score_exact > best_primary$score_exact)) {
    best_fallback
  } else {
    best_primary
  }

  if (is.null(scored)) return(empty)

  tibble(best_path  = scored$path,  best_score = scored$score_exact,
         n_match    = as.integer(scored$n_match),
         n_total    = as.integer(scored$n_total))
}

jmp_keys <- all_regions |>
  filter(survey_family %in% c("DHS", "MICS")) |>
  distinct(iso3, survey_family, year)

message("Selecting best boundary per (iso3, survey_family, year) [",
        nrow(jmp_keys), " combinations]...")

boundary_selection <- jmp_keys |>
  mutate(result = pmap(list(iso3, survey_family, year), best_boundary)) |>
  unnest(result)

message("\n── Boundary match quality ──────────────────────")
boundary_selection |>
  mutate(quality = case_when(
    best_score >= 0.75 ~ ">=75% names match",
    best_score >= 0.50 ~ "50-74% match",
    best_score >  0    ~ "<50% match",
    TRUE               ~ "no boundary found"
  )) |>
  count(survey_family, quality) |>
  arrange(survey_family, desc(n)) |>
  print()

message("\n── Boundary selected per (country, year) ───────")
boundary_selection |>
  filter(!is.na(best_path)) |>
  mutate(boundary_file = basename(best_path)) |>
  distinct(iso3, survey_family, year, boundary_file, best_score) |>
  arrange(iso3, survey_family, year) |>
  print(n = 50)

# ── 5. Join JMP rows to matched boundary polygons ─────────────────────────────

load_boundary_sv <- function(path, family, jmp_norms = NULL) {
  if (family != "DHS") {
    shp <- file.path(path, paste0(basename(path), ".shp"))
    return(tryCatch(vect(shp), error = function(e) NULL))
  }
  # DHS: may have multiple sdr_subnational_boundaries*.shp
  shp_paths <- if (str_ends(path, "\\.zip")) {
    zfiles <- tryCatch(unzip(path, list = TRUE)$Name, error = function(e) character(0))
    matches <- zfiles[str_detect(zfiles, "shps/sdr_subnational_boundaries[^/]*\\.shp$")]
    paste0("/vsizip/", path, "/", matches)
  } else {
    list.files(file.path(path, "shps"),
               pattern = "^sdr_subnational_boundaries.*\\.shp$",
               full.names = TRUE)
  }
  if (length(shp_paths) == 0) return(NULL)
  vects <- map(shp_paths, ~ tryCatch(vect(.x), error = function(e) NULL))
  vects <- compact(vects)
  if (length(vects) == 0) return(NULL)
  if (length(vects) == 1 || is.null(jmp_norms)) {
    # Only one candidate or no JMP names supplied — take the largest
    return(vects[[which.max(map_int(vects, nrow))]])
  }
  # Multiple shapefiles: pick the one whose names best overlap with JMP names
  scores <- map_dbl(vects, function(v) {
    bn <- normalise(as.data.frame(v)$DHSREGEN) |> translate_fr_en()
    sum(jmp_norms %in% bn) / length(jmp_norms)
  })
  vects[[which.max(scores)]]
}

jmp_with_boundary <- all_regions |>
  filter(survey_family %in% c("DHS", "MICS")) |>
  left_join(
    boundary_selection |> select(iso3, survey_family, year, best_path, best_score),
    by = c("iso3", "survey_family", "year")
  ) |>
  filter(!is.na(best_path), best_score > 0)

message("\nBuilding survey-boundary spatial layer...")
sv_results <- jmp_with_boundary |>
  group_by(best_path, iso3, survey_family, year) |>
  group_map(function(rows, key) {
    rows <- rows |> mutate(iso3 = key$iso3, survey_family = key$survey_family,
                           year = key$year)
    path_family  <- if (key$best_path %in% dhs_catalogue$path) "DHS" else "MICS"
    jmp_rn_vals  <- normalise(rows$region_name) |> translate_fr_en() |> unique()
    sv           <- load_boundary_sv(key$best_path, path_family, jmp_rn_vals)
    if (is.null(sv)) return(NULL)

    name_col     <- if (path_family == "DHS") "DHSREGEN" else "ADMIN_NAME"
    sv_df        <- as.data.frame(sv)
    sv_df$row_id <- seq_len(nrow(sv_df))
    sv_df$rn     <- normalise(sv_df[[name_col]]) |> translate_fr_en()

    joined <- rows |>
      mutate(rn = normalise(region_name) |> translate_fr_en()) |>
      left_join(sv_df |> select(row_id, rn), by = "rn") |>
      filter(!is.na(row_id))

    if (nrow(joined) == 0) return(NULL)

    # Geometry with key columns only (WASH data never enters SpatVector)
    out_sv <- sv[joined$row_id]
    values(out_sv) <- joined |>
      mutate(boundary_source = basename(key$best_path),
             match_type      = "survey_boundary") |>
      select(all_of(key_cols)) |>
      as.data.frame()

    # WASH data as plain data frame (avoids terra column-name mangling)
    out_data <- joined |>
      mutate(boundary_source = basename(key$best_path),
             match_type      = "survey_boundary") |>
      select(all_of(c(key_cols, wash_data_cols)))

    list(geom = out_sv, data = out_data)
  }) |>
  compact()

survey_sv   <- if (length(sv_results) > 0) do.call(rbind, map(sv_results, "geom")) else NULL
survey_data <- map(sv_results, "data") |> bind_rows()
message("Survey-boundary matched rows: ", if (!is.null(survey_sv)) nrow(survey_sv) else 0)

# ── 6. GADM as alternative boundary source ───────────────────────────────────
# GADM provides a universal ADM1 boundary set for all countries.
# Treat it as just another boundary option (alongside DHS/MICS/OTHER survey boundaries)
# so that select_best_divisions.R can rank GADM against survey boundaries per country.

gadm_rows <- all_regions

message("\nRows for GADM boundary: ", nrow(gadm_rows),
        " (", n_distinct(gadm_rows$iso3), " countries)")

message("Loading GADM ADM_1...")
gadm1 <- vect(gadm_path, layer = "ADM_1") |>
  filter(GID_0 %in% iso3_jmp)

gadm_lookup <- gadm1 |>
  as_tibble() |>
  select(GID_0, GID_1, NAME_1, VARNAME_1) |>
  mutate(alt_names = str_split(VARNAME_1, "\\|")) |>
  unnest_longer(alt_names, keep_empty = TRUE) |>
  pivot_longer(c(NAME_1, alt_names), names_to = NULL, values_to = "gadm_name") |>
  filter(!is.na(gadm_name), gadm_name != "") |>
  mutate(gadm_name_norm = normalise(gadm_name) |> translate_fr_en()) |>
  distinct(GID_0, GID_1, gadm_name_norm)

gadm1_names <- gadm1 |> as_tibble() |> select(GID_1, NAME_1)

gadm_regions <- gadm_rows |>
  distinct(iso3, region_name) |>
  mutate(region_norm = normalise(region_name) |> translate_fr_en())

# Exact match
gadm_exact <- gadm_regions |>
  left_join(gadm_lookup,
            by = c("iso3" = "GID_0", "region_norm" = "gadm_name_norm"),
            relationship = "many-to-many") |>
  slice_head(n = 1, by = c(iso3, region_name)) |>
  mutate(match_type = if_else(!is.na(GID_1), "exact", NA_character_))

gadm_matched_exact <- gadm_exact |> filter(!is.na(GID_1))
gadm_unmatched     <- gadm_exact |> filter(is.na(GID_1)) |> select(-GID_1, -match_type)

# Fuzzy match for remaining
fuzzy_best <- function(iso, name_norm) {
  cands <- gadm_lookup |> filter(GID_0 == iso)
  if (nrow(cands) == 0)
    return(tibble(GID_1 = NA_character_, dist = NA_real_))
  dists <- stringdist(name_norm, cands$gadm_name_norm, method = "jw")
  i     <- which.min(dists)
  tibble(GID_1 = cands$GID_1[i], dist = dists[i])
}

gadm_fuzzy <- gadm_unmatched |>
  rowwise() |>
  mutate(fuzzy_best(iso3, region_norm)) |>
  ungroup() |>
  mutate(match_type = case_when(
    dist <= 0.10 ~ "fuzzy_high",
    dist <= 0.20 ~ "fuzzy_medium",
    TRUE         ~ NA_character_
  )) |>
  filter(!is.na(match_type))

gadm_region_map <- bind_rows(
  gadm_matched_exact |> select(iso3, region_name, GID_1, match_type),
  gadm_fuzzy         |> select(iso3, region_name, GID_1, match_type)
) |>
  slice_head(n = 1, by = c(iso3, region_name))

message("GADM exact:   ", nrow(gadm_matched_exact))
message("GADM fuzzy:   ", nrow(gadm_fuzzy))
message("GADM no match:", nrow(gadm_regions) - nrow(gadm_region_map))

gadm_joined <- gadm_rows |>
  left_join(gadm_region_map, by = c("iso3", "region_name")) |>
  left_join(gadm1_names, by = "GID_1") |>
  filter(!is.na(GID_1))

gadm_sv   <- NULL
gadm_data <- NULL
if (nrow(gadm_joined) > 0) {
  # WASH data as plain data frame (never stored in SpatVector)
  gadm_data <- gadm_joined |>
    mutate(boundary_source = "GADM_ADM1") |>
    select(all_of(c(key_cols, wash_data_cols)), GID_1)

  # Geometry with key columns only
  gadm_sv <- gadm1 |>
    filter(GID_1 %in% gadm_joined$GID_1) |>
    left_join(
      gadm_data |> select(all_of(key_cols), GID_1) |> as.data.frame(),
      by = "GID_1"
    )
  # Strip GADM-native columns, keep only keys
  values(gadm_sv) <- values(gadm_sv)[, key_cols, drop = FALSE]
  gadm_data <- gadm_data |> select(-GID_1)
}
message("GADM spatial rows: ", if (!is.null(gadm_sv)) nrow(gadm_sv) else 0)

# ── 7. Combine and save ───────────────────────────────────────────────────────
# Geometry: key columns only (safe for terra::rbind — no space-named columns)
# Data: plain data frame with key columns + WASH indicators

dir.create("results", showWarnings = FALSE)

geom_parts <- compact(list(survey_sv, gadm_sv))
wash_geom <- if (length(geom_parts) > 1) {
  do.call(rbind, geom_parts)
} else if (length(geom_parts) == 1) {
  geom_parts[[1]]
} else {
  NULL
}

wash_data <- bind_rows(survey_data, gadm_data)

if (!is.null(wash_geom)) {
  message("\n── Final output ────────────────────────────────")
  message("Geometry features: ", nrow(wash_geom))
  message("Data rows:         ", nrow(wash_data))
  message("Countries covered: ", n_distinct(wash_data$iso3))

  for (col in wash_data_cols) {
    n_val <- sum(!is.na(wash_data[[col]]))
    n_na  <- sum(is.na(wash_data[[col]]))
    message(sprintf("  %s: %d values, %d NA", col, n_val, n_na))
  }

  writeVector(wash_geom, "results/JMP_regions_matched.gpkg", overwrite = TRUE)
  message("Saved: results/JMP_regions_matched.gpkg  (geometry + keys)")

  write_csv(wash_data, "results/JMP_regions_matched_data.csv")
  message("Saved: results/JMP_regions_matched_data.csv  (keys + WASH data)")
}

# Boundary selection report for review
boundary_selection |>
  mutate(boundary_file = basename(coalesce(best_path, ""))) |>
  select(iso3, survey_family, year, boundary_file, best_score, n_match, n_total) |>
  write_csv("results/boundary_match_report.csv")
message("Match report: results/boundary_match_report.csv")
