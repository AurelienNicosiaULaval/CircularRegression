# CircularRegression 0.5.0

* Add `circular_regression()` as the main fixed-effect modeling interface.
* Add `angular_two_step()` as an explicit consensus-then-homogeneous workflow.
* Add S3 methods for coefficients, fitted values, residuals, predictions,
  summaries, plots, log-likelihoods and information criteria.
* Improve consensus numerical stability for large Bessel-function arguments.
* Add validation for finite angles, non-negative modifiers, weights, controls
  and initial values.
* Align `logLik.consensus()`, `AIC.consensus()` and `BIC.consensus()` with the
  full von Mises log-likelihood by including the normalizing constant. This
  changes absolute consensus likelihood and information-criterion values but
  does not change the fitted estimates.
* Add `summary.angular_re()` and `print.summary.angular_re()`.
* Preserve model-frame `na.action` information in `angular()` objects.
* Replace the draft overview vignette with two reproducible HTML vignettes.
* Add pkgdown configuration and a workflow diagnostics article.
* Add a package-data workflow vignette for R Journal reviewer support.
* Add a minimal GitHub Actions R CMD check workflow.
* Clarify random-effects and special-wrapper documentation.
* Expand tests for simulation recovery, predictions, NA handling, weights,
  modulo invariance and small-sample fits.

# CircularRegression 0.4.0

* See `CHANGELOG_0.4.0.md` in the development repository for development
  notes. That file is not included in the CRAN build.

# CircularRegression 0.1.1

* Improve stability of angular and consensus estimators with QR-based updates and better handling of reference-only models.
* Implement observation weights in `consensus()` and remove the deprecated `model` argument across the API.
* Make `angular_re()` usable without workarounds, add Hessian conditioning checks, and provide more robust predictions.
* Refresh documentation and vignette examples, and add an automated test suite covering key model features.
