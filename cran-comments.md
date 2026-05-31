## Test environments

- Local macOS 26.5, R 4.5.0, aarch64-apple-darwin20.
- GitHub Actions R-CMD-check workflow on macOS, Windows and Ubuntu
  (latest observed run on 2026-05-31: success).

## R CMD check results

This is a new submission.

`devtools::test()`:

- 0 failures
- 0 warnings
- 0 skipped tests
- 118 passed expectations

`devtools::build_vignettes()`:

- completed successfully

`devtools::check(args = "--as-cran", error_on = "never")`:

- 0 errors
- 0 warnings
- 0 notes

`spelling::spell_check_package(".")`:

- no spelling errors found

`urlchecker::url_check()`:

- all URLs are correct

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

## External checks

win-builder R-release: submitted on 2026-05-31 for the current 0.5.0 source.
Result is not recorded in this repository.

win-builder R-devel: submitted on 2026-05-31 for the current 0.5.0 source.
Result is not recorded in this repository.

R-hub: not run in this session.

## Dataset redistribution

On 2026-05-31, the package maintainer confirmed that the processed datasets
included in the package can be redistributed with the package under the package
license or under compatible terms.

## Changes in this version

- Added `circular_regression()` as the main fixed-effect interface.
- Added S3 methods for printing, summaries, coefficients, fitted values,
  residuals, predictions, plots, log-likelihoods and information criteria.
- Added `summary.angular_re()` and `print.summary.angular_re()`.
- Added `angular_two_step()` as an explicit consensus-then-homogeneous workflow.
- Improved numerical stability for the consensus likelihood through stable
  Bessel-ratio and log-Bessel computations.
- Aligned consensus `logLik`, `AIC` and `BIC` with the full von Mises
  likelihood by including the normalizing constant.
- Added validation for finite angles, non-negative modifiers, valid weights,
  control values and initial values.
- Preserved model-frame `na.action` information in `angular()` objects.
- Replaced the old draft vignette with reproducible HTML vignettes and added a
  package-data workflow vignette for reviewer support.
- Expanded the test suite for estimation, predictions, NA handling, weights,
  small samples, modulo invariance, package datasets, wrapper outputs and
  random-effect methods.
- Added pkgdown configuration and a GitHub Actions R CMD check workflow.
- Added a maintainer-facing `data-raw/DATA_PROVENANCE.md` checklist for dataset
  sources, transformations and redistribution-permission confirmation.
- Clarified `noshiro` provenance using maintainer-supplied information from
  Louis-Paul Rivest and documented the remaining unverified provenance items for
  `multiplebison`.
- Added Nicosia et al. (2017) as a related scientific reference for the
  `multiplebison` study context while documenting that the exact mapping to the
  processed package dataset remains unverified from the repository.

## Downstream dependencies

No reverse dependencies are known at this stage.

## Dataset provenance

The package includes four documented datasets: `bison`, `multiplebison`,
`Sandhopper` and `noshiro`. Their Rd files include source/provenance notes. The
repository also tracks a maintainer-facing checklist in
`data-raw/DATA_PROVENANCE.md`; this directory is excluded from the CRAN source
package. The checklist records the transformations currently documented in the
package and records maintainer confirmation of redistribution status for the
processed package datasets.

For `noshiro`, the Rd file now identifies Rivest (1997) as the scientific
source and states that the processed data were provided by Louis-Paul Rivest to
the maintainer after duplicate observations had been removed. The repository
does not include raw Noshiro data, a reconstruction script, documented angular
units, or an archived redistribution-permission record. For `multiplebison`,
Nicosia et al. (2017) is documented as a related scientific reference for
GPS-collared bison movement in Prince Albert National Park. The exact mapping
between that article and the processed two-bison July-October 2013 package file
remains not verifiable from the repository.
