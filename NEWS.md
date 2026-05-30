# CircularRegression 0.5.0

* Add `circular_regression()` as the main fixed-effect modeling interface.
* Add `angular_two_step()` as an explicit consensus-then-homogeneous workflow.
* Add S3 methods for coefficients, fitted values, residuals, predictions,
  summaries, plots, log-likelihoods and information criteria.
* Improve consensus numerical stability for large Bessel-function arguments.
* Add validation for finite angles, non-negative modifiers, weights, controls
  and initial values.
* Replace the draft overview vignette with two reproducible HTML vignettes.
* Add pkgdown configuration and a workflow diagnostics article.
* Expand tests for simulation recovery, predictions, NA handling, weights,
  modulo invariance and small-sample fits.

# CircularRegression 0.1.1

* Improve stability of angular and consensus estimators with QR-based updates and better handling of reference-only models.
* Implement observation weights in `consensus()` and remove the deprecated `model` argument across the API.
* Make `angular_re()` usable without workarounds, add Hessian conditioning checks, and provide more robust predictions.
* Refresh documentation and vignette examples, and add an automated test suite covering key model features.
