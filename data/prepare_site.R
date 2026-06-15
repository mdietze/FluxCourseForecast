## prepare_site.R
## ---------------
## Converts the downloaded FLUXNET CSV for any site into the exact R objects
## that Mike Dietze's FluxCourseModelCalib.Rmd expects.  Run this script
## (or source it in Mike's setup chunk) BEFORE knitting FluxCourseModelCalib.Rmd.
##
## OBJECTS PRODUCED IN THE GLOBAL ENVIRONMENT
## -------------------------------------------
##   flux        data.frame  FLUXNET data, year-filtered, -9999→NA, +NEE_REF cols
##   date        POSIXct     timestamp vector (same length as nrow(flux))
##   inputs      data.frame  SSEM drivers: date (POSIXct), temp (°C), PAR (umol/m2/s)
##   nep         numeric     Net Ecosystem Production = -NEE (umol/m2/s)
##   nep.qc      numeric     QC flag for NEE (0=measured, 1-3=gap-fill quality)
##   nep.unc     numeric     Joint uncertainty on NEE (umol/m2/s)
##   X           matrix      [ne × 3] initial C pool ensemble (Mg/ha)
##   X.orig      matrix      copy of X before any calibration changes it
##   params      data.frame  [ne × 9] prior parameter ensemble for SSEM
##
## PREREQUISITES
## -------------
##   1. Run data/download_site.R first to obtain {site_id}_HH.csv.
##   2. R packages: tidyverse, lubridate, compiler (base R).
##   3. R/functions.R must be accessible from the working directory.
##
## HOW TO RUN (from the project root)
## ------------------------------------
##   source("data/prepare_site.R")   — from inside R or RStudio
##   Rscript data/prepare_site.R     — from the terminal
##
## TESTING
##   After running, verify the checklist at the end of this script.
##   All checks should report [PASS].

## ── Set your site and year here ───────────────────────────────────────────────
## When sourced from FluxCourseSiteExplorer.Rmd, site_id, year, and data_dir
## are set in the calling environment from params; the guards below preserve them.
if (!exists("site_id"))  site_id  <- "US-NR1"  # AmeriFlux site ID matching download_site.R
if (!exists("year"))     year     <- 2015       # year to filter to (Mike's Rmd uses US-NR1 2015)
if (!exists("data_dir")) data_dir <- "data"     # path to the data folder (relative to project root)
## ─────────────────────────────────────────────────────────────────────────────

## ---------------------------------------------------------------------------
## Metadata header
## ---------------------------------------------------------------------------
cat("prepare_site.R\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Commit:", tryCatch(system("git rev-parse --short HEAD", intern=TRUE),
                        error=function(e) "unknown"), "\n")
cat("Site:", site_id, "| Year:", year, "\n\n")

## ---------------------------------------------------------------------------
## 0. Load model functions and packages
## ---------------------------------------------------------------------------
## rdirichlet.orig and beta.match (used in §11 parameter priors) live in
## R/functions.R.  We source it only if these functions are not already in
## the session, so sourcing prepare_site.R inside Mike's setup chunk (which
## also sources functions.R) does not reload functions unnecessarily.

if (!exists("rdirichlet.orig", mode="function")) {
  fpath <- file.path("R", "functions.R")
  if (!file.exists(fpath))
    stop("Cannot find R/functions.R.\n",
         "Run this script from the project root directory (not from data/).")
  source(fpath)
  message("Loaded R/functions.R")
}

## Also load utils.R if present (provides helper functions for plotting)
if (!exists("plot_params", mode="function") &&
    file.exists(file.path("R", "utils.R"))) {
  source(file.path("R", "utils.R"))
  message("Loaded R/utils.R")
}

library(tidyverse)   # read_csv, dplyr verbs, pipe
library(lubridate)   # year(), as_date(), etc.
library(compiler)    # cmpfun() — used to byte-compile rdirichlet for speed

## ---------------------------------------------------------------------------
## 1. Resolve ensemble size (ne)
## ---------------------------------------------------------------------------
## ne is defined in Mike's setup chunk (ne = 100 or up to 5000).
## If prepare_site.R is run standalone before ne is defined, default to 100.

if (!exists("ne") || !is.numeric(ne) || length(ne) != 1) {
  ne <- 100
  message("'ne' not found in environment — defaulting to ne = 100")
}
cat("Ensemble size (ne):", ne, "\n\n")

## ---------------------------------------------------------------------------
## 2. Locate the sub-daily CSV
## ---------------------------------------------------------------------------
## download_site.R always writes the sub-daily output as {site_id}_HH.csv,
## following our repository convention (see also data/US-NR1/US-NR1_HH.csv).
## We also check for {site_id}_HR.csv produced by the older
## data/US-MMS/extract_hr_workaround.R for backward compatibility.

hh_path <- file.path(data_dir, site_id, paste0(site_id, "_HH.csv"))
hr_path <- file.path(data_dir, site_id, paste0(site_id, "_HR.csv"))

if (file.exists(hh_path)) {
  csv_path    <- hh_path
  product_res <- "HH"   # 30-min half-hourly
} else if (file.exists(hr_path)) {
  csv_path    <- hr_path
  product_res <- "HR"   # 60-min hourly (or 30-min filed under HR naming)
} else {
  stop(
    "No sub-daily CSV found for site '", site_id, "'.\n",
    "Expected:  ", hh_path, "\n",
    "Or:        ", hr_path, "\n",
    "Run data/download_site.R first to download the site data."
  )
}

cat("CSV source:      ", csv_path, "\n")
cat("(Large files may take 30–90 seconds to read.)\n\n")

## ---------------------------------------------------------------------------
## 3. Read the CSV and replace -9999 with NA
## ---------------------------------------------------------------------------
## FLUXNET products use -9999 as the universal missing-value sentinel.
## readr::read_csv() does NOT auto-convert -9999 to NA, so we pass
## na = c("NA", "-9999") to handle both the original -9999 and any
## previously-written "NA" strings.
##
## Note: flux_read() / readr::write_csv() already converts -9999 → NA before
## saving, so files from download_site.R typically arrive clean.  The na=
## argument is a defensive safeguard for files from other pipelines.

flux_raw <- readr::read_csv(
  csv_path,
  na             = c("NA", "-9999"),   # both missing-value encodings → NA
  show_col_types = FALSE               # suppress readr column-type output
)

cat(sprintf("Loaded: %d rows × %d columns\n", nrow(flux_raw), ncol(flux_raw)))

## ---------------------------------------------------------------------------
## 4. Filter to the requested year
## ---------------------------------------------------------------------------
## Two timestamp formats exist depending on which pipeline produced the file:
##
##   DATETIME_START:  ISO 8601 string "YYYY-MM-DDTHH:MM:SSZ" saved by
##                    readr::write_csv() from flux_read() output (download_site.R).
##                    readr auto-parses this as POSIXct on read-back.
##   TIMESTAMP_START: 12-digit integer YYYYMMDDHHMM saved by
##                    data.table::fwrite() (extract_hr_workaround.R).
##                    Stored as numeric; we extract the year from the first 4 chars.

if ("DATETIME_START" %in% names(flux_raw)) {
  ## ISO 8601 path
  yr_vec     <- lubridate::year(flux_raw$DATETIME_START)
  ts_format  <- "DATETIME_START"
  cat("Timestamp column: DATETIME_START (ISO 8601 → POSIXct)\n")
} else if ("TIMESTAMP_START" %in% names(flux_raw)) {
  ## Integer YYYYMMDDHHMM path — prevent scientific notation on large ints
  yr_vec     <- as.integer(
    substr(format(flux_raw$TIMESTAMP_START, scientific=FALSE, digits=12), 1, 4)
  )
  ts_format  <- "TIMESTAMP_START"
  cat("Timestamp column: TIMESTAMP_START (integer YYYYMMDDHHMM)\n")
} else {
  stop(
    "Cannot find DATETIME_START or TIMESTAMP_START in:\n  ", csv_path,
    "\nCheck that download_site.R completed successfully for this site."
  )
}

## Apply the year filter
flux <- flux_raw[yr_vec == year, ]
cat(sprintf("Filtered to year %d: %d rows (%.1f%% of full %d-row record)\n",
            year, nrow(flux), 100*nrow(flux)/nrow(flux_raw), nrow(flux_raw)))

if (nrow(flux) == 0)
  stop("No rows found for year ", year,
       ". Check that the downloaded record covers this year.")

## Free memory — the full-record file can be ~450 MB for US-NR1
rm(flux_raw)

## ---------------------------------------------------------------------------
## 5. Parse timestamps into a POSIXct date vector
## ---------------------------------------------------------------------------
## Mike's Rmd requires `date` as a POSIXct vector of the same length as flux.
##
## DATETIME_START (US-NR1 style): already POSIXct from readr auto-detection,
##   but we call as.POSIXct() defensively in case readr stored it as character.
## TIMESTAMP_START (US-MMS style): parse with strptime() using the
##   YYYYMMDDHHMM format string.  format(..., scientific=FALSE) prevents
##   large integers from being coerced to "1.23e+11" before strptime().

if (ts_format == "DATETIME_START") {
  date <- as.POSIXct(flux$DATETIME_START, tz="UTC")
} else {
  date <- as.POSIXct(
    strptime(
      format(flux$TIMESTAMP_START, scientific=FALSE, digits=12),
      format = "%Y%m%d%H%M",
      tz     = "UTC"
    )
  )
}

n_na_date <- sum(is.na(date))
if (n_na_date > 0)
  warning(n_na_date, " NA values in parsed date vector — inspect timestamp column.")

cat(sprintf("date: %d POSIXct values  [%s  to  %s]\n",
            length(date),
            format(min(date, na.rm=TRUE), "%Y-%m-%d %H:%M"),
            format(max(date, na.rm=TRUE), "%Y-%m-%d %H:%M")))

## ---------------------------------------------------------------------------
## 6. Auto-detect model timestep (seconds per step)
## ---------------------------------------------------------------------------
## Mike's SSEM.orig detects timestep internally from inputs$date:
##   timestep <- as.numeric(median(diff(as.numeric(inputs$date))))
## We replicate this detection here so that the parameter prior calculations
## in §11 use the correct per-timestep scaling factor (timestep/86400/365).
##
## Expected values:
##   HH (30-min):  timestep = 1800 s
##   HR (60-min):  timestep = 3600 s

timestep <- as.integer(
  median(diff(as.numeric(date[!is.na(date)])), na.rm=TRUE)
)
cat(sprintf("Detected timestep: %d s  (%s)\n",
            timestep,
            if      (timestep == 1800) "half-hourly (HH)"
            else if (timestep == 3600) "hourly (HR)"
            else    paste0("non-standard — expected 1800 or 3600")))

## ---------------------------------------------------------------------------
## 7. Build the `inputs` data frame (SSEM meteorological drivers)
## ---------------------------------------------------------------------------
## Mike's SSEM expects a data frame named `inputs` with three columns:
##   date  POSIXct  timestamp (same as the `date` vector above)
##   temp  numeric  air temperature (°C), gap-filled
##   PAR   numeric  photosynthetically active radiation (umol m-2 s-1)
##
## PAR is derived from incoming shortwave radiation (SW_IN, W m-2) using the
## Campbell and Norman (1998, p.151) conversion: PAR ≈ SW_IN / 0.486.
## This converts the full shortwave band to the 400–700 nm photosynthetic band.
##
## ERA5 reanalysis columns (TA_ERA, SW_IN_ERA) are PREFERRED over the
## tower-based gap-filled columns (TA_F, SW_IN_F) because ERA5 provides a
## continuous record even when tower measurements are missing.  If ERA5
## columns are absent, we fall back to TA_F and SW_IN_F with an approximate
## PAR conversion (× 2.1 umol/J ≈ W m-2 → umol m-2 s-1 broadband).

use_era <- "TA_ERA" %in% names(flux) && "SW_IN_ERA" %in% names(flux)

if (use_era) {
  temp_col   <- "TA_ERA"
  swin_col   <- "SW_IN_ERA"
  par_factor <- 1 / 0.486   # Campbell and Norman (1998) p.151
  cat("Driver columns: TA_ERA, SW_IN_ERA (ERA5; preferred)\n")
} else {
  temp_col   <- "TA_F"
  swin_col   <- "SW_IN_F"
  par_factor <- 2.1          # approximate W m-2 → umol m-2 s-1 for broadband SW
  message(
    "Note: ERA5 columns absent. Falling back to TA_F / SW_IN_F.\n",
    "PAR conversion uses approximate factor 2.1 (not Campbell & Norman 0.486).\n",
    "For best accuracy, use a product version that includes TA_ERA and SW_IN_ERA."
  )
}

if (!temp_col %in% names(flux))
  stop("Column '", temp_col, "' not found in ", csv_path,
       ".\nCheck that download_site.R ran successfully.")
if (!swin_col %in% names(flux))
  stop("Column '", swin_col, "' not found in ", csv_path, ".")

inputs <- data.frame(
  date = date,
  temp = flux[[temp_col]],               # air temperature (°C)
  PAR  = flux[[swin_col]] * par_factor   # PAR (umol m-2 s-1)
)

cat(sprintf("inputs: %d rows | temp [%.1f, %.1f] °C | PAR [%.0f, %.0f] umol/m2/s\n",
            nrow(inputs),
            min(inputs$temp, na.rm=TRUE), max(inputs$temp, na.rm=TRUE),
            min(inputs$PAR,  na.rm=TRUE), max(inputs$PAR,  na.rm=TRUE)))

## ---------------------------------------------------------------------------
## 8. Select the best available NEE column
## ---------------------------------------------------------------------------
## Mike's Rmd hardcodes flux$NEE_VUT_REF (variable u*-threshold partitioning).
## Some site products (including US-NR1 v1.3_r1) provide only NEE_CUT_REF
## (constant u*-threshold).  Both estimate NEE but can differ by ±10–20 gC/m2/yr.
##
## We create three new columns so Mike's Rmd needs only ONE change
## (replace NEE_VUT_REF with NEE_REF):
##   flux$NEE_REF         — whichever NEE variant is available
##   flux$NEE_REF_QC      — matching QC flag
##   flux$NEE_REF_JOINTUNC — matching joint uncertainty
##
## u* threshold methods (Papale et al. 2006, Biogeosciences):
##   VUT: threshold varies year-to-year; generally preferred for long records.
##   CUT: single threshold for the full record; simpler, comparably valid.

if ("NEE_VUT_REF" %in% names(flux)) {
  ## VUT method is available — use Mike's default
  flux$NEE_REF <- flux$NEE_VUT_REF
  flux$NEE_REF_QC <- if ("NEE_VUT_REF_QC" %in% names(flux))
                       flux$NEE_VUT_REF_QC else NA_real_
  flux$NEE_REF_JOINTUNC <- if ("NEE_VUT_REF_JOINTUNC" %in% names(flux))
                             flux$NEE_VUT_REF_JOINTUNC else NA_real_
  cat("NEE column: NEE_VUT_REF (VUT method; Mike's default)\n")

} else if ("NEE_CUT_REF" %in% names(flux)) {
  ## VUT absent — fall back to CUT method
  flux$NEE_REF <- flux$NEE_CUT_REF
  flux$NEE_REF_QC <- if ("NEE_CUT_REF_QC" %in% names(flux))
                       flux$NEE_CUT_REF_QC else NA_real_
  flux$NEE_REF_JOINTUNC <- if ("NEE_CUT_REF_JOINTUNC" %in% names(flux))
                             flux$NEE_CUT_REF_JOINTUNC else NA_real_
  message(
    "Note: NEE_VUT_REF not available for ", site_id, ".\n",
    "Using NEE_CUT_REF (constant u* threshold method).\n",
    "CUT and VUT can differ by 10-20 gC m-2 yr-1.\n",
    "See Papale et al. (2006) for details."
  )

} else {
  stop(
    "Neither NEE_VUT_REF nor NEE_CUT_REF found in ", csv_path, ".\n",
    "Check that download_site.R produced a complete FLUXMET file."
  )
}

## If JOINTUNC is entirely absent, fill with NA (column must exist for the checklist)
if (all(is.na(flux$NEE_REF_JOINTUNC)))
  message("Note: NEE_REF_JOINTUNC is all NA (absent in this product version).")

## ---------------------------------------------------------------------------
## 9. Create observation vectors (nep, nep.qc, nep.unc)
## ---------------------------------------------------------------------------
## NEP (Net Ecosystem Production) = -NEE by sign convention.
## Positive NEP = net carbon uptake by the ecosystem (photosynthesis > respiration).
## Negative NEP = net carbon release (respiration > photosynthesis).
##
## Mike's Assessment chunk uses these three variable names.
## The ONE change to Mike's Rmd is replacing NEE_VUT_REF with NEE_REF here:
##   ORIGINAL:  nep     = -flux$NEE_VUT_REF
##              nep.qc  = flux$NEE_VUT_REF_QC
##              nep.unc = flux$NEE_VUT_REF_JOINTUNC
##   MODIFIED:  nep     = -flux$NEE_REF
##              nep.qc  = flux$NEE_REF_QC
##              nep.unc = flux$NEE_REF_JOINTUNC
##
## All other code in Mike's Rmd that references nep, nep.qc, nep.unc works
## unchanged — only the source column name changes.

nep     <- -flux$NEE_REF           # sign-flip: NEE uptake is negative, NEP positive
nep.qc  <- flux$NEE_REF_QC        # 0 = directly measured, 1-3 = gap-fill quality tier
nep.unc <- flux$NEE_REF_JOINTUNC  # joint uncertainty (umol m-2 s-1); may be all NA

cat(sprintf("nep:     %d values | %.1f%% non-NA\n",
            length(nep), 100*mean(!is.na(nep))))
cat(sprintf("nep.qc:  %d values | QC=0 (measured): %d rows (%.1f%%)\n",
            length(nep.qc), sum(nep.qc==0, na.rm=TRUE),
            100*mean(nep.qc==0, na.rm=TRUE)))

## ---------------------------------------------------------------------------
## 10. Initial conditions  X  and  X.orig
## ---------------------------------------------------------------------------
## X is the [ne × 3] ensemble matrix of starting C pool sizes (Mg/ha):
##   Column 1: Leaf biomass
##   Column 2: Structural biomass (wood + coarse roots + fine roots)
##   Column 3: Soil organic matter (litter + CWD + mineral soil C)
##
## For US-NR1: exact values from Mike's Rmd, derived from the AmeriFlux BADM
##   metadata file (AMF_AA-Flx_FLUXNET-BIF_CCBY4_20220606).
## For all other sites: reasonable temperate forest defaults.
##   Replace with site-specific BADM values for best results.

pool.lab <- c("leaf", "wood", "SOM")  # column labels for diagnostic output

if (site_id == "US-NR1") {

  ## ── US-NR1: BADM-derived pool sizes ──────────────────────────────────────
  ## Source: AmeriFlux BADM — https://ameriflux.lbl.gov/sites/siteinfo/US-NR1
  ## Unit conversion: g/m2 × (10000 m2/ha) / (1e6 g/Mg) = g/m2 × 0.01 = Mg/ha
  ## (equivalently: 1e-6 × 1e4 = 1e-2 = 0.01)

  Bwood <- 14500 * 1e-6 * 10000   # aboveground tree biomass: ~145 Mg/ha
  Bleaf <- 2950  * 0.01  * 4.2/14 # leaf biomass with empirical LAI correction factor
  ## SOM = litter (Mg/ha) + CWD (Mg/ha) + soil (kg/m2 converted: × 1e-3 × 10000)
  ## The two-element vectors give two replicate BADM measurements per pool.
  SOM   <- c(1.57, 1.58) + c(0.49, 1.39) + c(2.06, 2.59) * 1e-3 * 10000

  ## BADM-derived trait data used in the parameter priors below (§11)
  ## LMA (g/m2) from BADM; SLA = 1000/LMA because SLA is in m2/kg
  SLA_badm       <- 1e3 / c(193.7, 205.1, 237.7)  # m2/kg for US-NR1 conifers
  litterfall_badm <- c(215.8) * 0.01 * 3            # gC/m2/yr → Mg/ha/yr

  ## Ensemble sampling: leaf and wood have only one BADM point each,
  ## so we assume 10% CV as the prior uncertainty.  SOM uses observed
  ## variability between the two replicate measurements.
  X_mean <- c(mean(Bleaf), mean(Bwood), mean(SOM))
  X <- as.matrix(X_mean)
  if (ne > 1) {
    X <- as.matrix(cbind(
      rnorm(ne, X_mean[1], Bleaf    * 0.1),  # 10% CV on leaf biomass
      rnorm(ne, X_mean[2], Bwood    * 0.1),  # 10% CV on wood biomass
      rnorm(ne, X_mean[3], sd(SOM))          # natural variability across BADM reps
    ))
  }

} else {

  ## ── Other sites: temperate forest defaults ──────────────────────────────
  ## These are reasonable starting points for a temperate deciduous or
  ## mixed forest.  Replace with BADM values from:
  ##   https://ameriflux.lbl.gov/sites/siteinfo/{site_id}#related
  ## or from the BADM download in download_site.R.

  message(
    "Note: Using temperate forest DEFAULT initial conditions for ", site_id, ".\n",
    "  Bwood = 100 Mg/ha (±50%),  Bleaf = 3 Mg/ha (±50%),  SOM = 100 Mg/ha (±50%)\n",
    "  For best results, replace with site-specific BADM values:\n",
    "  https://ameriflux.lbl.gov/sites/siteinfo/", site_id, "#related"
  )

  ## Note: The US-NR1 initial conditions above were derived from BADM biomass
  ## and soil carbon measurements. For other sites, different data may be
  ## available (e.g. forest inventory, soil surveys, remote sensing) that
  ## would support different or more direct prior specifications. The generic
  ## defaults here are intentionally broad to reflect this uncertainty.

  Bwood <- 100   # structural biomass (Mg/ha)
  Bleaf <- 3     # leaf biomass (Mg/ha)
  SOM   <- 100   # total soil organic matter (Mg/ha)

  ## Generic SLA prior: temperate forest mean ~10 m2/kg, SD ~3 m2/kg
  ## Replace with site-specific trait data where available (e.g. TRY database)
  SLA_mean <- 10   # m2/kg
  SLA_sd   <- 3    # m2/kg
  SLA_badm <- rnorm(10, SLA_mean, SLA_sd)  # synthetic sample for mean/sd calls below
  litterfall_badm <- c(150) * 0.01 * 3       # ~150 gC/m2/yr → Mg/ha/yr

  X_mean <- c(Bleaf, Bwood, SOM)
  X <- as.matrix(X_mean)
  if (ne > 1) {
    X <- as.matrix(cbind(
      rnorm(ne, Bleaf, Bleaf * 0.5),
      rnorm(ne, Bwood, Bwood * 0.5),
      rnorm(ne, SOM,   SOM   * 0.5)
    ))
  }

}

## Enforce positivity — biomass pools cannot be negative
X[X < 0] <- 0

## Save a copy before calibration overwrites it (needed by Mike's EKI section)
X.orig <- X

for (i in 1:3) {
  cat(sprintf("  X col %d (%s): mean=%.2f, sd=%.2f Mg/ha\n",
              i, pool.lab[i],
              mean(X[, i]),
              if (ne > 1) sd(X[, i]) else 0))
}

## ---------------------------------------------------------------------------
## 11. Parameter priors
## ---------------------------------------------------------------------------
## This section reproduces Mike's parameter prior distributions exactly,
## drawn from R/functions.R utilities and BADM-anchored SSEM defaults.
## See FluxCourseModelCalib.Rmd §Parameters for ecological interpretations.
##
## Site-varying inputs (set in §10):
##   SLA_badm        — specific leaf area from BADM (m2/kg)
##   litterfall_badm — annual litterfall flux from BADM (Mg/ha/yr)
##   X[,1]           — leaf biomass ensemble (Mg/ha)
##   timestep        — seconds per model step (from §6)

params <- list()

## --- Light use efficiency (alpha) ---
## Lognormal prior: alpha is strictly positive.  Mean 0.02 umol C / umol PAR.
params$alpha <- rlnorm(ne, log(0.02), 0.05)

## --- Soil respiration Q10 ---
## Temperature sensitivity: Rh multiplied by Q10 for every 10°C increase.
## Normal prior centred on 2.1 (widely observed in soils; ranges ~1.5–3.5).
params$Q10 <- rnorm(ne, 2.1, 0.5)

## --- Basal soil respiration (Rbasal) ---
## Rh = Rbasal × SOM × Q10^(T/10) in the model.
## Divided by Q10^2.5 to express Rbasal at 0°C baseline (25°C reference).
params$Rbasal <- rlnorm(ne, log(0.2), 1.5) / (params$Q10^2.5)

## --- Fractional carbon allocation (falloc) ---
## Dirichlet prior over (Ra, wood NPP, leaf NPP) — must sum to 1.
## Ra ≈ 0.5 of GPP (Litton et al. 2007); leaf NPP ≈ 31.5% of total NPP (Quaife et al. 2008).
Ra    <- 0.5
alloc <- matrix(c(Ra, (1 - 0.315)*(1 - Ra), 0.315*(1 - Ra)), 1)
## Neff controls how tightly the ensemble clusters around `alloc`
Neff  <- matrix(rnbinom(ne, 1000, mu=50), ne)
## Byte-compile rdirichlet.orig for speed (the ne×1000 matrix multiply is slow otherwise)
rdirichlet    <- compiler::cmpfun(rdirichlet.orig)
params$falloc <- rdirichlet(ne, Neff %*% alloc)
## Note: when converted to data.frame below, falloc becomes falloc.1, falloc.2, falloc.3

## --- Specific leaf area (SLA) ---
## Links leaf biomass to LAI in the model: LAI = Bleaf × SLA.
## Normal prior centred on site BADM values.
params$SLA <- rnorm(ne, mean(SLA_badm), sd(SLA_badm))

## --- Leaf litterfall turnover rate ---
## Fraction of leaf C lost per timestep via senescence and abscission.
## Estimated from BADM litterfall flux (Mg C/ha/yr) ÷ leaf biomass (Mg C/ha).
## Monte Carlo approach propagates uncertainty in both numerator and denominator.
if (length(litterfall_badm) > 1) {
  sdlitterfall <- sd(litterfall_badm)
} else {
  sdlitterfall <- 0.1 * litterfall_badm  # assume 10% uncertainty if only one value
}

## sd(X[,1]) needs at least 2 elements; guard against ne=1 edge case
sd_leaf <- if (ne > 1) sd(X[, 1]) else X_mean[1] * 0.1

## Ratio: annual turnover rate = litterfall / leaf biomass
lit <- rnorm(10000, mean(litterfall_badm), sdlitterfall/sqrt(2)) /
       rnorm(10000, mean(X[, 1]),          sd_leaf/sqrt(2))

## Scale from yr-1 to fraction per timestep using the detected timestep
lit.mu <- rnorm(ne, mean(lit), sd(lit)) * timestep/86400/365
lit.sd <- 1/sqrt(rgamma(ne, 10, 10*var(lit))) * timestep/86400/365

## Moment-match to a Beta distribution (turnover must be in [0, 1])
litterfall.param <- beta.match(lit.mu, lit.sd^2)
params$litterfall <- rbeta(ne, litterfall.param$a, litterfall.param$b)

## --- Background tree mortality rate ---
## Prior mean: 1/142 yr-1 (Dietze et al. 2011, Global Change Biology).
## This represents the probability that any given tree dies in a given year.
mortality.mu    <- 1/rpois(ne, 142) * timestep/86400/365
mortality.sd    <- rbeta(ne, 4, 4) * mortality.mu * timestep/86400/365
mortality.param <- beta.match(mortality.mu, mortality.sd^2)
params$mortality <- rbeta(ne, mortality.param$a, mortality.param$b)

## Flatten to a data frame; enforce positivity (all params must be ≥ 0)
params <- as.data.frame(params)
params[params < 0] <- 0

cat(sprintf("params: %d rows × %d columns\n", nrow(params), ncol(params)))
cat(sprintf("Column names: %s\n", paste(names(params), collapse=", ")))

## ---------------------------------------------------------------------------
## 12. Final checklist
## ---------------------------------------------------------------------------
## Verify that every object Mike's Rmd expects exists and has the right shape.

cat("\n", strrep("=", 65), "\n", sep="")
cat("PREPARE_SITE.R OBJECT CHECKLIST\n")
cat(strrep("=", 65), "\n\n")
cat(sprintf("Site: %s | Year: %d | ne: %d\n\n", site_id, year, ne))

OK  <- "[PASS]"
BAD <- "[FAIL]"
chk <- function(label, test) {
  pass <- isTRUE(test)
  cat(sprintf("  %s  %s\n", if(pass) OK else BAD, label))
  pass
}

results <- c(
  ## flux
  chk("flux is a data.frame",                       is.data.frame(flux)),
  chk("flux has rows > 0",                          nrow(flux) > 0),
  chk(sprintf("flux has %d rows (year %d)",
              nrow(flux), year),                    nrow(flux) > 0),
  chk("flux$NEE_REF column present",                "NEE_REF" %in% names(flux)),
  chk("flux$NEE_REF_QC column present",             "NEE_REF_QC" %in% names(flux)),
  chk("flux$NEE_REF_JOINTUNC column present",       "NEE_REF_JOINTUNC" %in% names(flux)),
  chk("flux has no -9999 sentinel in NEE_REF",
      !any(flux$NEE_REF == -9999, na.rm=TRUE)),

  ## date
  chk("date is POSIXct",                            inherits(date, "POSIXct")),
  chk("date length == nrow(flux)",                  length(date) == nrow(flux)),
  chk("date has no NAs",                            sum(is.na(date)) == 0),

  ## inputs
  chk("inputs is a data.frame",                     is.data.frame(inputs)),
  chk("inputs$date column exists",                  "date" %in% names(inputs)),
  chk("inputs$temp column exists",                  "temp" %in% names(inputs)),
  chk("inputs$PAR column exists",                   "PAR"  %in% names(inputs)),
  chk("inputs has same rows as flux",               nrow(inputs) == nrow(flux)),

  ## nep, nep.qc, nep.unc
  chk("nep is numeric",                             is.numeric(nep)),
  chk("nep length == nrow(flux)",                   length(nep) == nrow(flux)),
  chk("nep.qc is numeric or integer",               is.numeric(nep.qc) || is.integer(nep.qc)),
  chk("nep.qc length == nrow(flux)",                length(nep.qc) == nrow(flux)),
  chk("nep.unc is numeric",                         is.numeric(nep.unc)),
  chk("nep.unc length == nrow(flux)",               length(nep.unc) == nrow(flux)),

  ## X and X.orig
  chk("X is a matrix",                              is.matrix(X)),
  chk("X has 3 columns (leaf, wood, SOM)",          ncol(X) == 3),
  chk(sprintf("X has ne=%d rows", ne),              nrow(X) == ne || (ne == 1 && length(X) == 3)),
  chk("X has no negative values",                   !any(X < 0, na.rm=TRUE)),
  chk("X.orig has same dimensions as X",            identical(dim(X.orig), dim(X))),

  ## params
  chk("params is a data.frame",                     is.data.frame(params)),
  chk(sprintf("params has ne=%d rows", ne),         nrow(params) == ne),
  chk("params has exactly 9 columns",               ncol(params) == 9),
  chk("params has 'alpha' column",                  "alpha"     %in% names(params)),
  chk("params has 'Q10' column",                    "Q10"       %in% names(params)),
  chk("params has 'Rbasal' column",                 "Rbasal"    %in% names(params)),
  chk("params has 'falloc.1' column",               "falloc.1"  %in% names(params)),
  chk("params has 'falloc.2' column",               "falloc.2"  %in% names(params)),
  chk("params has 'falloc.3' column",               "falloc.3"  %in% names(params)),
  chk("params has 'SLA' column",                    "SLA"       %in% names(params)),
  chk("params has 'litterfall' column",             "litterfall" %in% names(params)),
  chk("params has 'mortality' column",              "mortality"  %in% names(params)),
  chk("params has no NAs",                          !any(is.na(params))),
  chk("params has no negative values",              !any(params < 0, na.rm=TRUE))
)

n_pass <- sum(results)
n_fail <- sum(!results)

cat(strrep("-", 65), "\n")
if (n_fail == 0) {
  cat(sprintf("\nREADY — all %d checks passed.\n\n", n_pass))
  cat("Objects created in global environment:\n")
  cat("  flux, date, inputs, nep, nep.qc, nep.unc, X, X.orig, params\n\n")
  cat("Required change in FluxCourseModelCalib.Rmd Assessment chunk:\n")
  cat("  nep     = -flux$NEE_VUT_REF       -->  nep     = -flux$NEE_REF\n")
  cat("  nep.qc  = flux$NEE_VUT_REF_QC    -->  nep.qc  = flux$NEE_REF_QC\n")
  cat("  nep.unc = flux$NEE_VUT_REF_JOINTUNC --> nep.unc = flux$NEE_REF_JOINTUNC\n")
} else {
  cat(sprintf("\nWARNING — %d check(s) FAILED. Review the [FAIL] lines above.\n", n_fail))
}

cat("\n", strrep("=", 65), "\n\n", sep="")
cat("prepare_site.R complete.\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
