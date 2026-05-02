# CircularRegression package audit

Date: 2026-05-02
Branch: `amelioration-1`

## Sources consulted

- CRAN Repository Policy, CRAN Repository Maintainers, accessed 2026-05-02:
  https://cran.r-project.org/web/packages/policies.html
- Writing R Extensions, R Core Team, version for R 4.6.0, 2026:
  https://cran.r-project.org/doc/manuals/r-release/R-exts.html
- The R Journal submission guidance, The R Journal, accessed 2026-05-02:
  https://rjournal.github.io/submissions.html
- The R Journal guidelines for papers about R packages, The R Journal,
  accessed 2026-05-02:
  https://rjournal.github.io/R_package_guidelines.html

## Executive summary

The package now has a clearer fixed-effect API centered on
`circular_regression()`, while preserving the previous model-specific functions
`angular()`, `consensus()`, `angular_two_step()`, `angular_re()` and the
special-case wrappers. The main changes address CRAN and R Journal readiness:
input validation, S3 method consistency, numerical stability, reproducible
vignettes, a broader test suite and package metadata cleanup.

The statistical methodology has not been changed except for justified
numerical safeguards. In particular, the consensus log-likelihood still uses
the same model structure, but the Bessel computations are now stabilized for
large concentration-resultant values.

## Critical issues found and addressed

1. Stale branch state and partial checkout risk

   The requested branch was `amelioration-1`, without accent. The branch was
   checked out and synchronized with `origin/amelioration-1` using
   `git pull --ff-only` before package edits.

2. Missing main interface

   The branch did not provide a single high-level fixed-effect modeling
   function comparable to the user experience of `glm()`. The package now
   exports `circular_regression()`, with methods
   `method = "two_step"`, `"homogeneous"` and `"consensus"`.

3. Unstable consensus likelihood for large linear predictors

   The previous `log(besselI(k, 0, expon.scaled = TRUE)) + k` calculation can
   become non-finite for very large `k`. This caused failed fits on realistic
   `bison` subsets. The implementation now uses centralized stable helpers for
   `logI0(k)` and `A1(k) = I1(k) / I0(k)`, with asymptotic fallbacks. The
   step-halving logic now treats non-finite proposals as failed proposals.

4. Weak validation of model inputs

   The package now validates finite numeric responses and angular covariates,
   finite non-negative modifiers in `x:z` terms, non-negative finite weights,
   non-zero total weights, valid controls and initial-value lengths.

5. S3 method gaps

   The package now includes S3 methods for `print()`, `summary()`, `coef()`,
   `fitted()`, `residuals()`, `predict()`, `plot()`, `logLik()`, `AIC()` and
   `BIC()` where appropriate for fixed-effect model objects.

6. Non-CRAN-friendly vignette material

   The old draft overview vignette contained broken citation syntax, heavy
   dependencies and auxiliary TeX and BibTeX files. It was replaced with two
   short reproducible HTML vignettes.

## Important improvements made

1. Package structure

   Internal helpers were centralized in `R/model_core.R`: formula parsing,
   model-frame handling, design-matrix construction, angle wrapping, stable
   Bessel helpers, reference-angle selection and prediction design for
   `newdata`.

2. API design

   The new main interface is:

   ```r
   circular_regression(
     formula,
     data,
     method = c("two_step", "homogeneous", "consensus")
   )
   ```

   The previous functions remain available for users who need direct control
   of a specific model.

3. Documentation

   Roxygen documentation was extended for the new interface, S3 methods,
   `angular_two_step()` and updated validation behavior. The README was
   rewritten as a concise package entry point.

4. Tests

   Tests were expanded to cover simulation-based recovery, predictions,
   `circular_regression()` S3 methods, `bison` fitting, NA and weight handling,
   invalid inputs, modulo `2*pi` invariance and small-sample reference-only
   fits.

5. Dependency hygiene

   Plotting packages are now suggested rather than imported by the namespace.
   Plot methods call them conditionally.

## Minor improvements made

- Removed stale generated files and project artifacts from the package source.
- Added `.Rbuildignore` patterns for local files, audit material, benchmarks
  and CRAN submission notes.
- Added `cran-comments.md` and `cran-readiness.md`.
- Added a lightweight benchmark script under `benchmarks/`, excluded from the
  CRAN build.

## Remaining limitations and recommendations

1. The package is closer to CRAN readiness, but final submission should only be
   made after repeating a successful tarball check with `R CMD check --as-cran`
   on the final submitted source archive.

2. For an R Journal submission, the package must be available on CRAN or
   Bioconductor, include unit tests and a full workflow vignette. The current
   package now addresses these points, but the manuscript should still include
   broader comparison with related circular-regression software and a clear
   explanation of the intended user workflow.

3. The random-effects model `angular_re()` remains separate from
   `circular_regression()` in this version, as requested. A future version
   could add a parallel high-level interface only if the API and statistical
   interpretation are fully justified.

## Verification status

The following checks were run on 2026-05-02:

- `devtools::test()`: 0 failures, 0 warnings, 73 passed expectations.
- `devtools::build_vignettes()`: completed successfully.
- `R CMD check --as-cran` from a clean external source copy: 0 errors,
  0 warnings, 2 notes.

The two notes were the expected CRAN incoming "New submission" note and a local
HTML Tidy version note.
