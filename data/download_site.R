## download_site.R
## ---------------
## Parameterized site download script for FluxCourseModelCalib.Rmd.
## Run this script ONCE before running prepare_site.R or the exercise Rmd.
## Downloads the AmeriFlux FLUXNET FULLSET data for any valid site ID and
## saves the half-hourly (HH) and daily (DD) CSVs to data/{site_id}/.
##
## PREREQUISITES
## -------------
## 1. Python 3.11–3.13 installed on your system.
## 2. ~/.Renviron contains your AmeriFlux attribution info:
##      AMERIFLUX_USER_NAME=Your Full Name
##      AMERIFLUX_USER_EMAIL=you@institution.edu
##      AMERIFLUX_INTENDED_USE=5    (5 = Education)
## 3. R packages: pak (CRAN), tidyverse, lubridate.
##    Run install.packages(c("pak","tidyverse","lubridate")) if missing.
##
## HOW TO RUN (from the project root directory)
## ----------------------------------------------
##   Rscript data/download_site.R
##   -- or inside RStudio --
##   source("data/download_site.R")
##
## OUTPUT
##   data/{site_id}/{site_id}_HH.csv   half-hourly FLUXMET (30-min or 60-min)
##   data/{site_id}/{site_id}_DD.csv   daily FLUXMET
##   data/{site_id}/{site_id}_MM.csv   monthly FLUXMET
##   data/{site_id}/{site_id}_ERA5_HR.csv  ERA5 half-hourly drivers
##   data/{site_id}/{site_id}_ERA5_DD.csv  ERA5 daily drivers
## The raw ZIP and unzipped files go to data/{site_id}/raw/ (gitignored).

## ── Set your site here ────────────────────────────────────────────────────────
## When sourced from FluxCourseSiteExplorer.Rmd, site_id and year are set in
## the calling environment from params; the guards below preserve those values.
if (!exists("site_id")) site_id <- "US-NR1"  # Change to any valid AmeriFlux or ICOS site ID
if (!exists("year"))    year    <- 2015       # Year you plan to use for the exercise
##              (Mike Dietze's FluxCourseModelCalib.Rmd uses US-NR1, 2015)
## ─────────────────────────────────────────────────────────────────────────────

## ---------------------------------------------------------------------------
## Metadata header — timestamped for reproducibility
## ---------------------------------------------------------------------------
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Commit:", tryCatch(system("git rev-parse --short HEAD", intern=TRUE),
                        error=function(e) "unknown"), "\n\n")
cat("Downloading site:", site_id, "| Exercise year:", year, "\n\n")

## ---------------------------------------------------------------------------
## 0. Directory setup
## ---------------------------------------------------------------------------
## All paths are relative to the project root so the script works whether
## called via Rscript or source() from inside RStudio.

out_dir   <- file.path("data", site_id)      # processed CSVs (git-tracked)
raw_dir   <- file.path(out_dir, "raw")        # downloaded ZIPs (gitignored)
unzip_dir <- file.path(raw_dir, "unzipped")  # extracted files (gitignored)

dir.create(out_dir,   showWarnings=FALSE, recursive=TRUE)
dir.create(raw_dir,   showWarnings=FALSE, recursive=TRUE)
dir.create(unzip_dir, showWarnings=FALSE, recursive=TRUE)

## ---------------------------------------------------------------------------
## 1. Set Python virtualenv path BEFORE loading fluxnet
## ---------------------------------------------------------------------------
## The fluxnet package uses Python under the hood via reticulate.
## CRITICAL: set RETICULATE_PYTHON before any library(fluxnet) or
## library(reticulate) call so reticulate locks in the correct virtualenv.
## This virtualenv is created by flux_install_shuttle() in §2.

Sys.setenv(RETICULATE_PYTHON = "~/.virtualenvs/fluxnet/bin/python")
cat("RETICULATE_PYTHON:", Sys.getenv("RETICULATE_PYTHON"), "\n")

## ---------------------------------------------------------------------------
## 2. Install and load the fluxnet package (pinned to v0.3.2)
## ---------------------------------------------------------------------------
## v0.3.2 is the first release that correctly extracts the _HR_ sub-daily files
## produced by the AmeriFlux FLUXNET v1.3_r1 product.  Earlier versions
## silently skipped them by searching only for the old "_HH_" pattern.

if (!requireNamespace("pak", quietly=TRUE)) install.packages("pak")
message("Ensuring EcosystemEcologyLab/fluxnet-package @ v0.3.2 is installed ...")
pak::pak("EcosystemEcologyLab/fluxnet-package@v0.3.2")

library(fluxnet)
library(tidyverse)
library(lubridate)

## ---------------------------------------------------------------------------
## 3. Install the fluxnet-shuttle Python library
## ---------------------------------------------------------------------------
## flux_install_shuttle() creates the "fluxnet" virtualenv (if absent) and
## installs the fluxnet-shuttle Python package.  This must be called ONCE
## before any flux_download() or flux_extract() call.
## On subsequent runs it checks and reuses the existing virtualenv.

message("Setting up fluxnet-shuttle Python environment ...")
flux_install_shuttle()
message("fluxnet-shuttle ready.\n")

## ---------------------------------------------------------------------------
## 4. AmeriFlux attribution check
## ---------------------------------------------------------------------------
## AmeriFlux sites are shared under CC-BY-4.0 with attribution requirements.
## Check that attribution variables are set; warn but do not stop if absent
## (the download may proceed but you should fill these in for data reuse).

for (v in c("AMERIFLUX_USER_NAME", "AMERIFLUX_USER_EMAIL",
            "AMERIFLUX_INTENDED_USE")) {
  if (!nzchar(Sys.getenv(v)))
    message("NOTE: ", v, " not set in ~/.Renviron. ",
            "Add it for AmeriFlux attribution compliance.")
}

## ---------------------------------------------------------------------------
## 5. Confirm the site exists in the AmeriFlux manifest
## ---------------------------------------------------------------------------
## flux_listall() queries the AmeriFlux API and returns metadata for every
## registered site.  We verify the requested site_id exists and that the
## requested year falls within the published record span.

message("Querying AmeriFlux site catalogue ...")
site_manifest <- flux_listall()
site_meta     <- dplyr::filter(site_manifest, site_id == !!site_id)

if (nrow(site_meta) == 0)
  stop(site_id, " not found in the AmeriFlux manifest.\n",
       "Check the site ID at https://ameriflux.lbl.gov/sites/site-list-and-info/")

cat("\n--- ", site_id, " manifest entry ---\n", sep="")
print(as.data.frame(site_meta))
cat(sprintf("Manifest record span: %s – %s\n",
            site_meta$first_year, site_meta$last_year))

## Warn if the requested year is outside the published record
if (!is.na(site_meta$first_year) && year < site_meta$first_year)
  warning("Year ", year, " is before the site's first year (", site_meta$first_year, ").")
if (!is.na(site_meta$last_year) && year > site_meta$last_year)
  warning("Year ", year, " is after the site's last year (", site_meta$last_year, ").")

## ---------------------------------------------------------------------------
## 6. Download the full FLUXNET record for this site
## ---------------------------------------------------------------------------
## flux_download() wraps the fluxnet-shuttle Python library, which handles
## AmeriFlux authentication.  It creates one ZIP file per site in raw_dir.
##
## Bug in v0.3.2: flux_download() raises an error if the ZIP already exists
## in download_dir.  Workaround: skip the call when the ZIP is already present.

existing_zips <- list.files(raw_dir, pattern="\\.zip$", full.names=TRUE)

if (length(existing_zips) > 0) {
  message("Existing ZIP(s) found — skipping download:")
  message("  ", paste(basename(existing_zips), collapse="\n  "))
} else {
  message("Downloading ", site_id, " FULLSET data (may take several minutes) ...")
  flux_download(site_ids=site_id, download_dir=raw_dir)
  message("Download complete.")
}

## ---------------------------------------------------------------------------
## 7. Extract the ZIP
## ---------------------------------------------------------------------------
## flux_extract() unzips everything in raw_dir to unzip_dir.
## v0.3.2 correctly handles _HR_ filenames (the v0.3.1 bug is fixed).
## If it fails (network issue, corrupted ZIP, etc.) we fall back to raw unzip.

message("Extracting ZIP files ...")
extract_ok <- tryCatch({
  flux_extract(zip_dir=raw_dir, output_dir=unzip_dir)
  TRUE
}, error=function(e) {
  message("flux_extract() failed: ", e$message)
  FALSE
})

if (!extract_ok) {
  ## Fallback: use the zip package to extract manually
  message("Trying manual fallback extraction with the zip package ...")
  if (!requireNamespace("zip", quietly=TRUE))
    stop("Install the 'zip' package first: install.packages('zip')")
  zip_files <- list.files(raw_dir, pattern="\\.zip$", full.names=TRUE)
  if (length(zip_files) == 0) stop("No ZIP files in: ", raw_dir)
  zip_file  <- tail(sort(zip_files), 1)  # use the most recent version
  message("Extracting: ", basename(zip_file))
  zip::unzip(zip_file, exdir=unzip_dir, junkpaths=TRUE)
  message("Fallback extraction complete.")
}

## ---------------------------------------------------------------------------
## 8. Discover extracted files and detect sub-daily label (HR vs HH)
## ---------------------------------------------------------------------------
## AmeriFlux FLUXNET v1.3_r1 uses '_HR_' for the 30-min file.
## FLUXNET-2015 used '_HH_'.  flux_discover_files() returns both; we detect
## whichever is present so flux_read() uses the correct resolution code.

message("Discovering extracted files ...")
file_manifest <- flux_discover_files(data_dir=unzip_dir)

cat("\n--- File manifest ---\n")
print(as.data.frame(file_manifest))

subdaily_label <- if ("HR" %in% file_manifest$time_resolution) "HR" else "HH"
message("Sub-daily manifest label detected: '", subdaily_label, "'")

## Build a map from flux_read() resolution code to manifest label
RES_MAP <- c(h=subdaily_label, d="DD", w="WW", m="MM", y="YY")

## ---------------------------------------------------------------------------
## 9. Read all FLUXMET and ERA5 resolutions
## ---------------------------------------------------------------------------
## flux_read() returns a single tibble per call (not a list) and already
## replaces -9999 with NA.  It renames the timestamp columns:
##   HR/HH → DATETIME_START / DATETIME_END (POSIXct)
##   DD    → DATE (Date)
##   WW    → DATE_START / DATE_END (Date)
##   MM/YY → TIMESTAMP (character) or YEAR depending on product version

read_one <- function(res_code, dataset_name) {
  manifest_label <- RES_MAP[res_code]
  present <- unique(file_manifest$time_resolution[file_manifest$dataset==dataset_name])
  if (!(manifest_label %in% present)) {
    message("  ", dataset_name, " ", manifest_label, ": absent from delivery — skipping")
    return(NULL)
  }
  message("  Reading ", dataset_name, " ", manifest_label, " ...")
  tryCatch(
    flux_read(file_manifest, resolution=res_code,
              datasets=dataset_name, site_ids=site_id),
    error=function(e) { message("  ERROR: ", e$message); NULL }
  )
}

message("\n--- Reading FLUXMET (all resolutions) ---")
fluxmet <- setNames(
  lapply(names(RES_MAP), read_one, "FLUXMET"),
  paste0("fluxmet_", names(RES_MAP))
)

message("\n--- Reading ERA5 (all resolutions) ---")
era5 <- setNames(
  lapply(names(RES_MAP), read_one, "ERA5"),
  paste0("era5_", names(RES_MAP))
)

hr      <- fluxmet$fluxmet_h   # sub-daily FLUXMET
dd      <- fluxmet$fluxmet_d   # daily FLUXMET
era5_hr <- era5$era5_h         # ERA5 sub-daily
era5_dd <- era5$era5_d         # ERA5 daily

## ---------------------------------------------------------------------------
## 10. Save processed CSVs
## ---------------------------------------------------------------------------
## The sub-daily output is always named {site_id}_HH.csv regardless of the
## product's internal _HR_ vs _HH_ convention, following our repo's naming
## standard (matches data/US-NR1/US-NR1_HH.csv etc.).

message("\n--- Saving CSVs ---")

save_csv <- function(df, fname, label) {
  if (is.null(df)) {
    message("  Skipping ", label, " (not available from this product)")
    return(invisible(NULL))
  }
  path <- file.path(out_dir, fname)
  readr::write_csv(df, path)
  message(sprintf("  %-40s  %d rows × %d cols  (%.1f MB)",
                  fname, nrow(df), ncol(df), file.size(path)/1e6))
}

save_csv(hr,              paste0(site_id, "_HH.csv"),      "FLUXMET sub-daily")
save_csv(dd,              paste0(site_id, "_DD.csv"),      "FLUXMET daily")
save_csv(fluxmet$fluxmet_m, paste0(site_id, "_MM.csv"),   "FLUXMET monthly")
save_csv(era5_hr,         paste0(site_id, "_ERA5_HR.csv"), "ERA5 sub-daily")
save_csv(era5_dd,         paste0(site_id, "_ERA5_DD.csv"), "ERA5 daily")

## ---------------------------------------------------------------------------
## 11. Confirmation summary
## ---------------------------------------------------------------------------

cat("\n", strrep("=", 65), "\n", "CONFIRMATION SUMMARY\n", strrep("=", 65), "\n", sep="")
cat(sprintf("Site:       %s\n", site_id))
cat(sprintf("Year:       %d\n", year))
cat(sprintf("Output dir: %s\n\n", out_dir))

if (!is.null(hr)) {
  ## Identify the best NEE column available
  nee_col <- if ("NEE_VUT_REF" %in% names(hr)) "NEE_VUT_REF"
             else if ("NEE_CUT_REF" %in% names(hr)) "NEE_CUT_REF"
             else NA_character_

  ## Check presence and completeness of key SSEM driver and flux variables
  key_vars <- c(nee_col, "TA_ERA", "SW_IN_ERA", "TA_F", "SW_IN_F")
  key_vars <- key_vars[!is.na(key_vars) & key_vars %in% names(hr)]

  cat("Key variable availability in sub-daily file:\n")
  for (v in key_vars) {
    n_total <- nrow(hr); n_nonNA <- sum(!is.na(hr[[v]]))
    cat(sprintf("  %-30s  %d/%d non-NA  (%.1f%%)\n",
                v, n_nonNA, n_total, 100*n_nonNA/n_total))
  }

  ## Report available years and check whether the exercise year is present
  if ("DATETIME_START" %in% names(hr)) {
    years_avail <- sort(unique(lubridate::year(hr$DATETIME_START)))
    cat(sprintf("\nYears available: %d – %d  (%d years)\n",
                min(years_avail), max(years_avail), length(years_avail)))
    n_yr_rows <- sum(lubridate::year(hr$DATETIME_START) == year)
    if (n_yr_rows > 0) {
      cat(sprintf("Year %d confirmed present: %d rows\n", year, n_yr_rows))
    } else {
      cat(sprintf("WARNING: year %d NOT found in downloaded data!\n", year))
    }
  }

  ## Note if NEE_VUT_REF is absent (relevant for US-NR1 users)
  if (!is.na(nee_col) && nee_col == "NEE_CUT_REF")
    cat("\nNote: NEE_VUT_REF absent — prepare_site.R will use NEE_CUT_REF.\n",
        "See Papale et al. (2006) for CUT vs VUT differences.\n")
}

cat(sprintf("\nFiles written to %s/:\n", out_dir))
for (f in list.files(out_dir, pattern="\\.csv$", full.names=TRUE))
  cat(sprintf("  %-40s  %.1f MB\n", basename(f), file.size(f)/1e6))

cat("\nNext step: run data/prepare_site.R to create model-ready objects.\n")
cat("\ndownload_site.R complete.\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
sessionInfo()
