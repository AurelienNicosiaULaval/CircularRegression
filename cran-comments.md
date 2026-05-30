## Test environments

- Local macOS 26.5, R 4.5.0, aarch64-apple-darwin20.

## R CMD check results

`devtools::test()`:

- 0 failures
- 0 warnings
- 0 skipped tests
- 107 passed expectations

`devtools::build_vignettes()`:

- completed successfully

`devtools::check(args = "--no-manual", error_on = "never")`:

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

## Downstream dependencies

No reverse dependencies are known at this stage.
