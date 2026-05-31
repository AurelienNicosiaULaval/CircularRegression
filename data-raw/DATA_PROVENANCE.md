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
- Angular units: radians.
- Raw-data reconstruction status: raw GPS files and a reconstruction script are
  not present in the repository.
- Redistribution status: maintainer confirmation required.

## multiplebison

- Package file: `data/multiplebison.rda`
- Documentation: `R/multiplebison.R`, `man/multiplebison.Rd`
- Described source: GPS-derived movement metrics for two plains bison
  (`1044-a` and `1045-a`) monitored in Prince Albert National Park,
  Saskatchewan, Canada, during July to October 2013.
- Scientific provenance status: partially documented. The maintainer identifies
  these data as associated with doctoral work and/or an article on bison
  movement. A related article has been identified: Nicosia, Duchesne, Rivest and
  Fortin (2017), "A multi-state conditional logistic regression model for the
  analysis of animal movement", The Annals of Applied Statistics, 11(3),
  1537-1560, doi:10.1214/17-AOAS1045. This article analyzes GPS-collared bison
  movement in Prince Albert National Park using an individual winter trajectory
  from November 2013 to April 2014. The repository does not verify whether the
  current two-bison July-October 2013 data file is a direct subset or
  preprocessing product of that analysis.
- Known transformations from package documentation: individual hourly GPS tracks
  were merged, rows were retained when the time gap between successive fixes was
  exactly one hour, and data were reshaped to a wide format with one row per
  timestamp. Turning angles were computed from consecutive movement directions.
  Direction and turning-angle variables are stored in degrees.
- Angular units: degrees.
- Raw-data reconstruction status: raw GPS files and a reconstruction script are
  not present in the repository.
- Redistribution status: redistribution permission should be explicitly
  confirmed or archived before CRAN submission.

## Sandhopper

- Package file: `data/Sandhopper.rda`
- Documentation: `R/Sandhopper.r`, `man/Sandhopper.Rd`
- Described source: sandhopper escape-orientation example discussed by
  D'Elia (2001) and Rivest and Kato (2019).
- Known transformations from package documentation: processed package data
  contain 72 rows and 22 variables, including individual identifiers, repeated
  escape orientations, environmental variables and morphometric variables.
- Angular units: degrees for `Azimuth`, `LN1`-`LN5` and `DirW`.
- Raw-data reconstruction status: a raw-data reconstruction script is not
  present in the repository.
- Redistribution status: maintainer confirmation required.

## noshiro

- Package file: `data/noshiro.rda`
- Documentation: `R/noshiro.R`, `man/noshiro.Rd`
- Described source: processed Noshiro earthquake direction data used in Rivest
  (1997), containing direction of steepest descent (`DIRDSC`) and direction of
  movement (`DIRMV`) for 678 locations.
- Scientific provenance status: documented by maintainer notes. Louis-Paul
  Rivest informed the maintainer that the data are treated in Rivest (1997),
  that duplicate observations were removed before inclusion, and that the
  resulting processed dataset has `n = 678`.
- Original analysis: `DIRMV` is the response variable to be predicted from
  `DIRDSC`.
- Bibliographic provenance: Rivest, L.-P. (1997). A decentered predictor for
  circular-circular regression. Biometrika, 84, 717-726.
- Known transformations from package documentation: duplicate observations were
  removed before inclusion in this package. No additional transformations are
  documented in the repository.
- Angular units: not documented in the repository or in the provenance notes
  provided to the maintainer.
- Raw-data reconstruction status: raw data and a reconstruction script are not
  present in the repository.
- Redistribution status: redistribution permission should be explicitly
  confirmed or archived before CRAN submission.

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
