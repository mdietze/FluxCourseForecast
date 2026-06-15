# CLAUDE.md — FluxCourseForecast (multisite prep)

## What this repo is and what these files do

FluxCourseForecast is a teaching repository for the Fluxcourse model-calibration exercise. It provides the R Markdown notebook `FluxCourseModelCalib.Rmd` (authored by Mike Dietze, BU) together with supporting R code in `R/functions.R` and `R/utils.R`. The three files added on this branch extend the exercise to work with any AmeriFlux site, not just the originally hard-coded one. `data/download_site.R` downloads FLUXNET FULLSET data for a user-specified site ID and year via the EcosystemEcologyLab fluxnet package, extracts half-hourly and daily CSVs, and saves ERA5 driver files alongside the tower data. `data/prepare_site.R` reads those CSVs and constructs the R objects (`flux`, `date`, `inputs`, `nep`, `nep.qc`, `nep.unc`, `X`, `X.orig`, `params`) that Mike's Rmd expects in its global environment — including site-specific initial conditions and parameter priors. `dashboard.html` is a pre-built student-facing summary artifact generated from the davidjpmoore/FluxCourseForecast repository; rebuild instructions are in that repo.

## Hard rules

- **Do not modify `FluxCourseModelCalib.Rmd`** except for the three NEE_REF lines described below.
- **Do not modify `R/functions.R` or `R/utils.R`.**
- **Git operations** (add, commit, push, merge, reset) require explicit human approval before running.

## Data sources

- FLUXNET data: EcosystemEcologyLab fluxnet package v0.3.2 (`EcosystemEcologyLab/fluxnet-package@v0.3.2`). Raw ZIPs land in `data/{site_id}/raw/` (gitignored); processed CSVs are written to `data/{site_id}/` (tracked).
- ERA5 drivers (`TA_ERA`, `SW_IN_ERA`) are preferred over tower gap-filled columns (`TA_F`, `SW_IN_F`). `prepare_site.R` selects ERA5 automatically when present and falls back to tower met otherwise.
- Do not download or commit raw FLUXNET ZIPs; use the fluxnet package only.

## Required change to `FluxCourseModelCalib.Rmd`

Mike's Assessment chunk hardcodes `flux$NEE_VUT_REF`. Some site products (including US-NR1 v1.3_r1) provide only `NEE_CUT_REF` (constant u\* threshold) instead of `NEE_VUT_REF` (variable u\* threshold). `prepare_site.R` resolves this by creating a unified column `flux$NEE_REF` (whichever variant is available) along with `NEE_REF_QC` and `NEE_REF_JOINTUNC`. The only required Rmd edit is three lines in the Assessment chunk:

```r
# Before
nep     = -flux$NEE_VUT_REF
nep.qc  = flux$NEE_VUT_REF_QC
nep.unc = flux$NEE_VUT_REF_JOINTUNC

# After
nep     = -flux$NEE_REF
nep.qc  = flux$NEE_REF_QC
nep.unc = flux$NEE_REF_JOINTUNC
```

No other code in the Rmd needs to change; everything downstream already uses `nep`, `nep.qc`, and `nep.unc`.

## `dashboard.html`

This file is a pre-built artifact. Do not regenerate it here. Rebuild instructions and source code are in the `davidjpmoore/FluxCourseForecast` repository.
