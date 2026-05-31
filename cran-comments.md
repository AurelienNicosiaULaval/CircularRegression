## Test environments

- Local macOS 26.5, R 4.5.0, aarch64-apple-darwin20.
- GitHub Actions R-CMD-check workflow on macOS, Windows and Ubuntu
  (latest observed run on 2026-05-31: success).
- win-builder R-release and R-devel submissions were sent on 2026-05-31;
  replace this line with the email results before CRAN upload.

## R CMD check results

This is a new submission.

`devtools::test()`:

- 0 failures
- 0 warnings
- 0 skipped tests
- 107 passed expectations

`devtools::build_vignettes()`:

- completed successfully

`devtools::check(args = "--as-cran", error_on = "never")`:

- 0 errors
- 0 warnings
- 0 notes

`spelling::spell_check_package(".")`:

- no spelling errors found

`R CMD check --as-cran` on a clean external source copy:

- 0 errors
- 0 warnings
- 2 notes

Notes:

- New submission.
- Local HTML validation skipped because the installed `tidy` executable is not
  recent enough.

The HTML validation note is environment-specific. No WARNING or ERROR was
observed.

Pending before upload:

- Copy the win-builder R-release and R-devel email results into this file.
- Run or record R-hub results if required for the final submission notes.

## Changes in this version

- Added `circular_regression()` as the main fixed-effect interface.
- Added S3 methods for printing, summaries, coefficients, fitted values,
  residuals, predictions, plots, log-likelihoods and information criteria.
- Added `angular_two_step()` as an explicit consensus-then-homogeneous workflow.
- Improved numerical stability for the consensus likelihood through stable
  Bessel-ratio and log-Bessel computations.
- Added validation for finite angles, non-negative modifiers, valid weights,
  control values and initial values.
- Replaced the old draft vignette with reproducible HTML vignettes and added a
  package-data workflow vignette for reviewer support.
- Expanded the test suite for estimation, predictions, NA handling, weights,
  small samples, modulo invariance, package datasets, wrapper outputs and
  random-effect methods.
- Added pkgdown configuration and a GitHub Actions R CMD check workflow.
- Added a maintainer-facing `data-raw/DATA_PROVENANCE.md` checklist for dataset
  sources, transformations and redistribution-permission confirmation.

## Downstream dependencies

No reverse dependencies are known at this stage.

## Dataset provenance

The package includes four documented datasets: `bison`, `multiplebison`,
`Sandhopper` and `noshiro`. Their Rd files include source/provenance notes. The
repository also tracks a maintainer-facing checklist in
`data-raw/DATA_PROVENANCE.md`; this directory is excluded from the CRAN source
package. The checklist records the transformations currently documented in the
package and identifies redistribution permissions that must be confirmed by the
maintainer before upload.
