# CircularRegression 0.4.0

- Added natural-parameter outputs to all special-case wrappers:
  - `meanDirectionModel()`
  - `decentredPredictorModel()`
  - `presnellModel()`
  - `jammalamadakaModel()`
- Wrappers now return:
  - `natural_parameters` (estimate, model SE, robust SE)
  - `natural_vcov_model`
  - `natural_vcov_robust`
  - `natural_jacobian`
- Implemented delta-method uncertainty propagation for wrapper back-transformations from homogeneous-model `beta` estimates.
- Wrapper fitting now defaults to `reference = "first"` when no explicit `reference` is supplied, to keep special-case parameterization stable.
- Expanded tests to validate natural-parameter recovery and associated SE computations.

# CircularRegression 0.3.0

- Reworked the core API to align with Rivest et al. (2016):
  - `angular()` now supports explicit reference strategies (`auto`, `first`, `name`).
  - `consensus()` now estimates kappa parameters directly (`initkappa`).
  - `coef.consensus()` now supports `type = "kappa"` and `type = "beta"`.
- Added `select_reference_angle()` and `angular_two_step()` implementing the recommended two-step workflow (consensus then homogeneous).
- Deprecated `pick_reference_angle()` in favor of `select_reference_angle()`.
- Removed the package vignette from the release tarball pending a full R Journal companion article.
- Fixed critical convergence bugs (`maxiter = 0`) in both `angular()` and `consensus()`.
- Added explicit identifiability checks for homogeneous fits.
- Stabilized `logLik.angular()` numerically and ensured `nobs`/`df` attributes for model-comparison generics.
- Made `predict.angular_re()` robust to missing calling-environment objects by storing model structures inside fitted objects.
- Updated wrappers (`meanDirectionModel`, `decentredPredictorModel`, `presnellModel`, `jammalamadakaModel`) to return two-step workflow outputs.
- Expanded tests for non-regression and new API behavior.

# CircularRegression 0.1.1

- Improve stability of angular and consensus estimators with QR-based updates and better handling of reference-only models.
- Implement observation weights in `consensus()` and remove the deprecated `model` argument across the API.
- Make `angular_re()` usable without workarounds, add Hessian conditioning checks, and provide more robust predictions.
- Refresh documentation and vignette examples, and add an automated test suite covering key model features.
