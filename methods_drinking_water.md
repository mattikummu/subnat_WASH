# Methods – Subnational Drinking Water Access (1990–2023)

## Overview

This document describes the full methodology for producing global, subnational
estimates of drinking water service levels for 1990–2023. The pipeline produces
four mutually exclusive and exhaustive service-level indicators that sum to
100 % for every administrative unit and year:

| Indicator | JMP ladder rung |
|-----------|----------------|
| **Basic** | Improved source, ≤30 min round trip |
| **Limited** | Improved source, >30 min |
| **Unimproved** | Unprotected well/spring, vendor, tanker |
| **Surface water** | River, canal, lake, pond |

The spatial output covers ~3,040 first-level administrative divisions (ADM1) in
240 countries/territories at 5-arcminute (~10 km) resolution, 1990–2023.

---

## 1. Observed Subnational Data

### 1.1 Source

Subnational drinking water data come from the WHO/UNICEF Joint Monitoring
Programme (JMP) country-specific Excel workbooks (downloaded individually per
country from the JMP website). Each workbook contains a `data_summary` worksheet
with regional survey observations.

### 1.2 Indicators extracted

Six indicator combinations are extracted from `var1`/`var2` columns:

- Water / Improved
- Water / Not Improved
- Water ladder / Basic
- Water ladder / Limited
- Water ladder / Unimproved
- Water ladder / Surface water

The four ladder-rung columns (Basic, Limited, Unimproved, Surface water) are
the primary modelling targets. "Improved" and "Not Improved" serve as quality
checks.

### 1.3 Coverage

| Dimension | Value |
|-----------|-------|
| Countries with subnational data | **102** |
| Unique subnational regions (country × name) | **2,339** |
| Survey-region-year observation rows | **4,280** |
| Survey year range | **1997–2020** |
| Data sources | DHS (incl. MIS, AIS, OTH), MICS, and other surveys |

### 1.4 Survey family classification

Each observation row is classified according to its source code suffix:

- `survey_family = "DHS"` for survey types DHS, MIS, AIS, OTH
- `survey_family = "MICS"` for MICS surveys
- `survey_family = "OTHER"` for all other sources

### 1.5 Data cleaning steps

1. Rows with region names starting with `..` (country-level summaries) are
   removed.
2. Egypt rural/urban rows and bare "rural"/"urban" region rows are removed.
3. Square-bracketed prefixes (e.g., `[N]`) and parenthetical suffixes are
   stripped from region names.
4. Region names are Unicode-normalised to ASCII, lowercased, and
   punctuation-stripped for subsequent string matching.

---

## 2. Linking Subnational Data to Administrative Boundaries

### 2.1 Boundary sources

Three boundary datasets are used:

| Source | Description |
|--------|-------------|
| **DHS boundaries** | Country-specific shapefiles from the DHS Programme (`sdr_subnational_boundaries*.shp`); one file per survey × year; field `DHSREGEN` |
| **MICS boundaries** | Country-specific shapefiles from UNICEF MICS (`geo1_*` and `geo2_*` folders); field `ADMIN_NAME` |
| **GADM 4.1** | Global Administrative Database layer `ADM_1`; fields `NAME_1`, `VARNAME_1` (pipe-separated aliases); used as universal fallback |

DHS catalogued: 1 file per survey-year. MICS: 1 folder per country (with
optional year ranges in the folder name, e.g., `geo1_sn2017_2023`).

### 2.2 Region-name matching

Region names are harmonised using:

1. **Normalisation**: transliterate to ASCII (`stringi`), lowercase, replace
   hyphens/dashes with spaces, remove non-alphanumeric characters, squish
   whitespace.
2. **Translation dictionary**: ~100 stored mappings convert French/indigenous
   language names to English equivalents as used in the boundary files (e.g.,
   French compass words `nord→north`, `ouest→west`; country-specific overrides
   for Iraq, Indonesia, Kazakhstan, Madagascar, Nepal, etc.).
3. **Suffix stripping**: Roman numeral region codes, trailing "region"/"oblysy"
   suffixes, BPS province codes.
4. **Exact matching**: normalised JMP name matched to normalised boundary name.
5. **Fuzzy matching** (Jaro–Winkler distance, via `stringdist`): for unmatched
   names, accept fuzzy matches with distance ≤ 0.10 (`fuzzy_high`) or ≤ 0.20
   (`fuzzy_medium`).

### 2.3 Boundary selection per (country, survey family, year)

For each unique (iso3, survey family, year), candidate boundary files are scored:

- **Score** = (fraction of JMP region names with exact match) + 0.01 × (mean
  Jaro–Winkler similarity across names).
- Ties broken by year proximity.
- Cross-family fallback (e.g., MICS boundary used for a DHS survey) only if
  primary family score < 50 %.
- GADM is always a parallel candidate for all observations.

### 2.4 Best administrative division selection

After matching all observations to all candidate boundaries, the boundary that
best represents each country's subnational data is selected via a ranking
procedure (script `04_select_best_divisions.R`):

1. **Merge equivalent geometries**: boundary files that produce identical polygon
   centroids (fingerprinting) are collapsed to a canonical source — preventing a
   single admin division from being treated as multiple independent sources.
2. **Spatial coverage filter**: boundary sources whose unique polygon count is
   < 80 % of the country's maximum achievable are dropped (prevents partial
   national coverage).
3. **Year completeness filter**: survey years where a boundary matched < 80 % of
   its maximum region count are dropped.
4. **Ranking criteria** (in order):
   - `all_covered` — all (region × year) observations have ≥ 1 WASH value
   - `spatial_coverage` — polygon count / maximum polygons
   - `n_years` — number of distinct survey years
   - `wash_coverage` — mean fraction of non-NA WASH indicators
   - `n_regions` — number of distinct spatial units (finer is better)
5. Each country's **best division** (rank 1) and optionally a **second division**
   (rank 2) are retained. The second division provides additional training
   observations for the machine learning model.

**Best-division summary:**

| Metric | Value |
|--------|-------|
| Countries with a best division | **99** |
| Countries also with a second division | **29** |
| Best-division ADM1 polygons | **1,293** |

Stable integer identifiers (`adm_id`) are assigned to all polygons:
`adm_id = ISO3_numeric × 1000 + polygon_row` (best division: rows 1–499;
second division: rows 500–999).

---

## 3. Global ADM1 Coverage (Training + Prediction)

Countries without any JMP subnational data receive GADM 4.1 ADM1 polygons as
fallback ("gadm_fallback" source), so that XGBoost predictions can be
computed for all territories.

| Source | ADM1 polygons | Countries |
|--------|--------------|-----------|
| Best division (DHS/MICS/GADM-matched) | 1,293 | 99 |
| GADM fallback | 1,747 | 145 |
| **Total global ADM1** | **3,040** | **240** |

All 3,040 units receive population-weighted predictor extraction and XGBoost
predictions for 1990–2023.

From the final ADM1 interpolation tables (`results/JMP_adm1_water_interpPredExtrap.csv`):
- Countries with any observed ADM1 data used in the final dataset: **95**
- Countries without observed ADM1 data (fully model-predicted at ADM1): **145**

The 145-country set corresponds to ISO3s in `results/adm1_global_lookup.csv`
with `source = "gadm_fallback"`.

---

## 4. National (ADM0) Observed Data

### 4.1 Modern data (2000–2024)

Source: **JMP_2025_WLD.xlsx**, sheet `wat` (machine-readable format).
Columns used: `wat_basal_t` (Basic), `wat_lim_t` (Limited),
`wat_unimp_t` (Unimproved), `wat_ns_t` (Surface water) — all total
(not urban/rural disaggregated).

| Metric | Value |
|--------|-------|
| Country-year rows | **5,660** |
| Countries / territories | **232** |
| Year range | **2000–2024** |

Data constraint checked: Basic + Limited + Unimproved + Surface water = 100
(tolerance 0.01 pp).

Observed ADM0 country coverage check (from interpolated ADM0 tables):
- Countries with any observed modern JMP (2000+): **231**
- Countries with any observed historical JMP: **211**
- Countries with any observed ADM0 data in either source (modern OR historical): **232**

### 4.2 Historical data (1990–2015)

Source: **estimates_1980-2015_v1.xlsx** (older JMP release).
Only two indicators are available historically:
- `Basic` (= `wat_total_improved_t`)
- `Surface water` (= `wat_surface_water_t`)

`Limited` and `Unimproved` are not directly reported for 1990–1999; they are
derived from XGBoost predictions (see §6.2).

---

## 5. Predictor Variables

Ten population-weighted predictors are extracted for all ADM0 and ADM1 units
for 1990–2023 using `exactextractr` (population weights from GHS-POP 5-arcmin,
1985–2025):

| Predictor | Description | Source / Notes |
|-----------|-------------|----------------|
| `travel_time` | Log₁₀ travel time to nearest city (min) | Weiss et al. global friction surface; static |
| `poverty` | Income poverty rate (%) | Subnational poverty raster; 1990–2022; last year padded |
| `msch` | Mean schooling years | Subnational HDI project; 1990–2023 |
| `gnic` | Log GNI per capita (natural log) | Subnational GNI raster; 1990–2023 |
| `gini` | Gini disposable income inequality | Subnational Gini raster; 1990–2023 |
| `urban` | Urbanisation rate (country-normalised) | GHS-SMOD; 1985–2023 (GHS 2023a) |
| `aridity` | Aridity index | Static raster |
| `gdi` | Gender Development Index | Subnational HDI project; 1990–2023 |
| `corruption` | V-Dem control of corruption (`v2x_corr`) | 1990–2024; last year padded |
| `gw_depth` | Groundwater table depth (m) | GLOBGM steady-state |

Population weighting formula:
$$\text{predictor}_{a,t} = \frac{\sum_{\text{pixel}} \text{var}_{t} \cdot \text{pop}_{t}}{\sum_{\text{pixel}} \text{pop}_{t}}$$

Extraction is parallelised across predictor rasters using `future`/`furrr`
(5 workers, up to 4 GB GDAL cache per worker). ADM0 and all ADM1 sets are
extracted in a single combined pass to minimise I/O.

**Predictor handling in XGBoost**: predictor NAs are retained and passed
directly to XGBoost, which handles missing values natively by learning
the optimal default direction per split. Only rows with NA in the *response*
variable are dropped.

---

## 6. XGBoost Ensemble Model

### 6.1 Training data

| Dimension | Value |
|-----------|-------|
| Total training rows | **7,269** |
| Countries | **223** |
| Admin units (adm_id) | **1,623** |
| ADM0 rows (countries without subnational data) | 2,837 |
| Best-division ADM1 rows | 3,384 |
| Second-division ADM1 rows | 1,048 |

Year sampling: to avoid excessive temporal autocorrelation and give each
admin area equal spatial weight, **5 years are sampled per admin area** (with
replacement) from their observed survey years. Admin areas with fewer than 5
observed years receive duplicated rows.

Four separate models are trained — one per WASH indicator (Basic, Limited,
Unimproved, Surface water) — in parallel using `future`/`furrr` (4 workers).

### 6.2 Model specification

Gradient boosting via the `xgboost` engine (R package `tidymodels`/`parsnip`):

- **Trees**: fixed at 1,500; actual number selected by early stopping
- **Validation slice**: 20 % internal validation for early stopping
- **`nthread`**: 1 per model (parallelism at indicator level, not within model)
- **No feature normalisation**: XGBoost is scale-invariant; omitting
  `step_normalize` keeps SHAP plots on original predictor scales

**Tuned hyperparameters** (space-filling grid, 40 combinations per indicator):

| Parameter | Range |
|-----------|-------|
| `tree_depth` | 2–6 |
| `learn_rate` | 10⁻³–10⁻⁰·⁷ (~0.001–0.2) |
| `loss_reduction` | 10⁻²–10¹·⁵ |
| `min_n` | 2–30 |
| `mtry` (`colsample_bytree`) | 0.5–1.0 (proportion) |
| `sample_size` (`subsample`) | 0.5–1.0 (proportion) |
| `stop_iter` | 20–80 rounds |
| `lambda` (L2 regularisation) | 10⁻²–10² (sampled on log scale) |

### 6.3 Cross-validation and train/test splits

**Grouped train/test split** (80/20): split at the admin-area level, stratified
by country:
- Countries with ≥ 2 admin areas: 80 % of admin areas to train, 20 % to test.
- Countries with exactly 1 admin area: always assigned to train (to avoid
  fully losing single-region countries from the test set).

**10-fold group cross-validation** (grouped by `adm_id`) is used within the
training set to select best hyperparameters (metric: RMSE).

**Country-level hold-out** (independent evaluation): 80 % of *countries* in
training, 20 % held out — assesses cross-country transferability directly
relevant to the 145 territories with no subnational training data.

### 6.4 Bootstrap bagging ensemble

Using the tuned hyperparameters, **50 bootstrap models** are trained on
admin-area-grouped bootstraps of the training data. Final predictions are the
**median** across the 50 models; the **standard deviation** quantifies
per-prediction uncertainty.

### 6.5 Bias correction

After bagging, a **binned calibration** correction is applied:

1. Ensemble median predictions on the test set are divided into 10 quantile
   bins.
2. Per-bin mean bias (actual − predicted) is computed.
3. Each prediction is shifted by its bin's mean bias, then clamped to [0, 100].

---

## 7. National (ADM0) Gap Filling (2000–2023)

For each country × indicator × year where the JMP national value is missing,
the value is filled using the following priority hierarchy:

| Priority | Method | Applies when |
|----------|--------|--------------|
| 1 | **Reported** | JMP value present |
| 2 | **Linear interpolation** | ≥ 2 reported years bracket the gap |
| 3 | **XGBoost ratio extrapolation (backward)** | Before first reported year; anchored to first observation × (XGBoost_t / XGBoost_first) |
| 4 | **XGBoost ratio extrapolation (forward)** | After last reported year; anchored to last observation |
| 5 | **Direct XGBoost prediction** | No observations at all |

Interpolation uses `zoo::na.approx` (linear, rule = 1, no extrapolation).
All values are clamped to [0, 100].

---

## 8. Historical ADM0 Extension (1990–1999) and Delta-Change Harmonisation

Because the modern JMP 2025 data begins at 2000, historical estimates
(1990–1999) are built from the legacy JMP estimates file (§4.2) and harmonised
using a **delta-change method** to avoid discontinuities at the year-2000 join:

1. For each country and indicator with historical observations, compute the
   mean additive offset over the 2000–2002 overlap window:
   $$\text{offset}_{i,c} = \overline{(\text{modern}_{t} - \text{historical}_{t})} \quad t \in \{2000, 2001, 2002\}$$
2. Shift the historical series: $\text{hist\_adj}_{t} = \text{hist\_filled}_{t} + \text{offset}_{i,c}$, clamped to [0, 100].
3. Proportionally rescale all four indicators to ensure they sum to 100 in
   every country-year.

For `Limited` and `Unimproved` (not in the historical observation file):
XGBoost predictions are scaled proportionally so that
Limited + Unimproved = 100 − (Basic + Surface water):
$$\text{Limited}_{\text{hist}} = \text{residual} \times \frac{\hat{y}_{\text{Limited}}}{\hat{y}_{\text{Limited}} + \hat{y}_{\text{Unimproved}}}$$

Countries without any historical observations use the modern dataset's XGBoost
back-filled values for 1990–1999 directly.

---

## 9. Subnational (ADM1) Gap Filling (1990–2023)

The same interpolation/extrapolation hierarchy as ADM0 (§7) is applied per
`adm_id`:

1. Reported (observed survey value for that region × year)
2. Linear interpolation (requires ≥ 2 observed years)
3. XGBoost ratio extrapolation (backward from first observation)
4. XGBoost ratio extrapolation (forward from last observation)
5. Direct XGBoost prediction

After filling, each ADM1 series for each indicator covers 1990–2023 as a
continuous annual time series.

---

## 10. ADM0 Ratio Adjustment and Normalisation

To ensure national consistency, ADM1 values are scaled so that their
population-weighted mean matches the ADM0 national total:

$$\text{ADM1}^*_{r,t} = \text{ADM1}_{r,t} \times \frac{\text{JMP\_ADM0}_{c,t}}{\overline{\text{ADM1}}_{c,t}^{\text{pop-wtd}}}$$

where $\overline{\text{ADM1}}_{c,t}^{\text{pop-wtd}}$ is the GHS-POP
population-weighted mean of filled ADM1 values for country $c$ in year $t$.
If the denominator is zero or missing, the scaling factor is set to 1 (no
adjustment).

Following scaling, the four indicator values are normalised within each ADM1 ×
year to sum exactly to 100:
$$\text{ADM1}^{**}_{r,t,i} = \frac{\text{ADM1}^*_{r,t,i}}{\sum_i \text{ADM1}^*_{r,t,i}} \times 100$$

Values are clamped to [0, 100] throughout.

---

## 11. Rasterisation

A **5-arcminute (~10 km at equator) global raster template** is built once
(script `05b`) by rasterising the global ADM1 polygon set using `terra`.
Each pixel encodes the integer `adm_id` of the polygon it falls in.

For each indicator × year, the template is classified from `adm_id` values to
indicator values using `terra::classify()`. This produces:

- 4 multi-layer GeoTIFF files:
  `rast_wash_{Basic|Limited|Unimproved|Surface.water}_1990_2023.tif`
  (34 annual layers each, LZW compressed)
- 4 slope rasters (Theil-Sen annual trends, §12)

---

## 12. Trend Analysis

**Theil-Sen slopes** are computed per ADM1 unit (`adm_id`) for each WASH
indicator over 1990–2023 using the `mblm` package (median-based linear
models). Statistical significance is assessed from the regression p-value.

Slope output rasters encode percentage-point change per year; non-significant
units are masked. Maps are produced with `tmap`.

---

## 13. Population Estimates per Service Level

For each ADM1 unit × year, the number of people on each water service level is:
$$\text{pop}_{r,t,i} = \text{ADM1}^{**}_{r,t,i} / 100 \times \text{pop\_GHS}_{r,t}$$

Population data: GHS-POP at 5-arcminute resolution (1985–2025), extracted via
`exactextractr`. Country and macro-region summaries are aggregated and Theil-Sen
trend slopes are computed for the population without basic water access
(Limited + Unimproved + Surface water).

---

## 14. Software and Reproducibility

All analyses are conducted in R 4.5. Key packages:

| Package | Role |
|---------|------|
| `terra`, `exactextractr` | Raster / vector operations, pop-weighted extraction |
| `sf` | Vector I/O and spatial joins |
| `tidymodels`, `xgboost` | Model specification, tuning, fitting |
| `future`, `furrr` | Parallel processing |
| `zoo` | Linear interpolation of time series |
| `mblm`, `broom` | Theil-Sen slope estimation |
| `stringdist`, `stringi` | Fuzzy string matching for region names |
| `tmap`, `ggplot2`, `patchwork` | Mapping and visualisation |
| `countrycode` | ISO3 ↔ numeric country codes |

The pipeline is fully scripted (scripts `01`–`14`) and reproducible from
raw JMP Excel downloads. Intermediate results are cached to disk.

---

## 15. Key Assumptions and Limitations

1. **Subnational survey representativeness**: JMP subnational data rely on
   DHS and MICS surveys, which may not cover all years or all admin units
   within a country. Regions not surveyed receive XGBoost predictions only.

2. **Single best division per country**: The pipeline selects one primary
   administrative division per country. Countries that changed administrative
   boundaries over the study period may have residual inconsistencies.

3. **Boundary fuzzy matching**: Region names matched with Jaro–Winkler
   distance ≤ 0.20 may introduce occasional mismatches, particularly for
   countries with many similarly named regions.

4. **Predictor availability**: Several predictor rasters end before 2023
   (e.g., poverty ends in 2022); missing trailing years are padded with the
   last available value, assuming no change.

5. **Historical indicators (1990–1999)**: `Limited` and `Unimproved` are not
   directly observed in the historical JMP file. They are derived from XGBoost
   predictions rescaled to the observed residual (100 − Basic − Surface.water),
   which introduces additional uncertainty for these categories before 2000.

6. **Delta-change harmonisation**: The additive offset is estimated over
   2000–2002 and applied uniformly across 1990–1999. Structural breaks in
   reporting methods before 2000 are not modelled explicitly.

7. **Countries without subnational data**: 145 territories use GADM ADM1
   boundaries with fully model-predicted values. Predictions for these areas
   rely entirely on the cross-country transferability of the XGBoost models,
   which is assessed via the country-level hold-out evaluation but is inherently
   less reliable than areas with survey data.

8. **XGBoost does not enforce the sum-to-100 constraint**: The four indicators
   are modelled independently and then normalised. Normalisation may shift
   individual indicators slightly from their raw model predictions.

9. **5-arcminute spatial resolution**: Sub-ADM1 spatial heterogeneity is not
   captured; all pixels within an ADM1 polygon receive identical indicator values.
