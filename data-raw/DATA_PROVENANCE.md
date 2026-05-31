# Dataset provenance and redistribution checklist

Date: 2026-05-31
Package version: 0.5.0

This directory records what can be verified from the current repository before
CRAN submission. It is excluded from the CRAN source package through
`.Rbuildignore`.

The current repository contains processed `.rda` files in `data/`, but it does
not contain the original raw files or complete reconstruction scripts for all
datasets. Therefore, this file is a maintainer-facing provenance checklist, not
a raw-data reconstruction workflow.

## Required maintainer action before CRAN upload

Before submitting to CRAN, confirm that every dataset included in `data/` can be
redistributed with the package under the package license or under compatible
terms. If this cannot be confirmed for a dataset, remove that dataset from the
CRAN submission or replace it with a fully documented redistributable example.

The redistribution rights listed below cannot be independently verified from the
repository alone.

## bison

- Package file: `data/bison.rda`
- Documentation: `R/bison.R`, `man/bison.Rd`
- Described source: reduced package dataset derived from the bison movement
  application in Rivest, Duchesne, Nicosia and Fortin (2016).
- Known transformations from package documentation: reduced to 5,696 rows and
  seven variables for angular-regression examples. Variables include movement
  direction, lagged movement directions, angles toward meadows and canopy gaps,
  and corresponding distance or weight variables.
- Raw-data reconstruction status: raw GPS files and a reconstruction script are
  not present in the repository.
- Redistribution status: maintainer confirmation required.

## multiplebison

- Package file: `data/multiplebison.rda`
- Documentation: `R/multiplebison.R`, `man/multiplebison.Rd`
- Described source: GPS-derived movement metrics for two plains bison
  (`1044-a` and `1045-a`) monitored in Prince Albert National Park,
  Saskatchewan, Canada, during July to October 2013.
- Known transformations from package documentation: individual hourly GPS tracks
  were merged, rows were retained when the time gap between successive fixes was
  exactly one hour, and data were reshaped to a wide format with one row per
  timestamp. Turning angles were computed from consecutive movement directions.
  Direction and turning-angle variables are stored in degrees.
- Raw-data reconstruction status: raw GPS files and a reconstruction script are
  not present in the repository.
- Redistribution status: maintainer confirmation required.

## Sandhopper

- Package file: `data/Sandhopper.rda`
- Documentation: `R/Sandhopper.r`, `man/Sandhopper.Rd`
- Described source: sandhopper escape-orientation example discussed by
  D'Elia (2001) and Rivest and Kato (2019).
- Known transformations from package documentation: processed package data
  contain 72 rows and 22 variables, including individual identifiers, repeated
  escape orientations, environmental variables and morphometric variables.
- Raw-data reconstruction status: a raw-data reconstruction script is not
  present in the repository.
- Redistribution status: maintainer confirmation required.

## noshiro

- Package file: `data/noshiro.rda`
- Documentation: `R/noshiro.R`, `man/noshiro.Rd`
- Described source: package data file containing direction of steepest descent
  and direction of movement for 678 locations associated with the Noshiro
  earthquake.
- Known transformations from package documentation: no transformations are
  documented beyond the processed two-column data frame.
- Raw-data reconstruction status: raw data and a reconstruction script are not
  present in the repository.
- Original bibliographic provenance: not documented in the repository.
- Redistribution status: maintainer confirmation required.

## Suggested future reconstruction scripts

When raw sources and permissions are available, add explicit scripts such as:

- `data-raw/bison.R`
- `data-raw/multiplebison.R`
- `data-raw/Sandhopper.R`
- `data-raw/noshiro.R`

Each script should read only raw files stored under `data-raw/`, produce the
corresponding object in `data/`, document all unit conversions and filtering,
and call `usethis::use_data(..., overwrite = TRUE)` only for the intended
dataset.
