# CRAN readiness checklist

Date: 2026-05-30

## Package metadata

- DESCRIPTION includes a maintainer, license, URL and BugReports field.
- DESCRIPTION includes `Language: en-US`.
- Strong dependencies are limited to packages used by core functionality.
- Plotting dependencies are kept in Suggests and used conditionally.
- .Rbuildignore excludes local project files, audit files, benchmark scripts,
  CRAN comments and temporary vignette logs.

## API and documentation

- Existing exported functions are preserved.
- `circular_regression()` provides the main fixed-effect interface.
- S3 methods are implemented for model objects where appropriate.
- Documentation is regenerated with roxygen2.
- Examples are short and reproducible.
- Four vignettes provide an introduction, an applied bison example, diagnostics
  and a package-data workflow.
- `inst/CITATION` returns version 0.5.0 from a temporary source installation.

## Statistical and numerical checks

- Angles are treated modulo `2*pi` through sine and cosine operations.
- Circular residuals are wrapped with `atan2(sin(x), cos(x))`.
- Consensus likelihood calculations use stable `logI0()` and `A1()` helpers.
- Invalid responses, angle covariates, modifiers, weights, controls and initial
  values are rejected with explicit errors.

## Test and check commands

Commands run successfully on 2026-05-30:

```r
devtools::document()
devtools::test()
devtools::build_vignettes()
pkgdown::build_site(preview = FALSE)
spelling::spell_check_package(".")
```

Commands run successfully on 2026-05-30:

```sh
R CMD build <clean temporary source copy>
R CMD check --as-cran CircularRegression_0.5.0.tar.gz
```

Observed result in a clean external check directory:

- 0 errors
- 0 warnings
- 2 notes

The two notes were the expected CRAN incoming "New submission" note and a local
HTML Tidy version note.

`devtools::check(args = "--no-manual", error_on = "never")` returned 0 errors,
0 warnings and 0 notes.

`pkgdown::build_site(preview = FALSE)` returned successfully with a non-blocking
pkgdown metadata warning that no canonical site URL is configured.
