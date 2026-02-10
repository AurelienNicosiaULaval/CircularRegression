# CircularRegression 0.4.0 - Release Notes

This document summarizes user-facing changes introduced in `v0.4.0`.

## Highlights

- Added natural-parameter outputs for all special-case wrappers.
- Added delta-method standard errors for natural parameters, using covariance from homogeneous `angular()` fits.
- Added both model-based and robust uncertainty propagation in wrapper outputs.
- Stabilized special-case parameterization by defaulting wrapper calls to `reference = "first"` when `reference` is not explicitly provided.

## New Wrapper Outputs

The following wrappers now include natural-parameter inference outputs:

- `meanDirectionModel()`
- `decentredPredictorModel()`
- `presnellModel()`
- `jammalamadakaModel()`

Each wrapper output now includes:

- `natural_parameters`: matrix with columns
  - `estimate`
  - `se_model`
  - `se_robust`
- `natural_vcov_model`: delta-method covariance using `homogeneous_fit$varcov0`
- `natural_vcov_robust`: delta-method covariance using `homogeneous_fit$varcov1`
- `natural_jacobian`: Jacobian used for the delta-method transformation

## Natural-Parameter Mappings

### Mean Direction Wrapper

- Natural parameter:
  - `mu`
- Transformation:
  - recovered from the synthetic two-angle representation as `atan2(beta_2, beta_1)` (wrapped to `[-pi, pi)`).

### Decentred Predictor Wrapper

- Natural parameters:
  - `r_hat`
  - `alpha_hat`
  - `beta_hat`
- Transformation:
  - recovered from the synthetic `(w, w + pi/2, 0, pi/2)` parameterization.

### Presnell Wrapper

- Natural parameters:
  - `beta2`
  - `beta3`
  - `beta4`
- Transformation:
  - ratio-based recovery from the synthetic `(0, pi/2, 0:w, pi/2:w)` form.

### Jammalamadaka-Sengupta Wrapper

- Natural parameters:
  - `gamma2`
  - `gamma3`
  - `gamma4`
  - `gamma5`
  - `gamma6`
- Transformation:
  - linear map from the synthetic six-term representation.

## Inference Details

For each wrapper:

1. Fit homogeneous model through the two-step pipeline.
2. Extract `beta_hat` and covariance (`varcov0` and `varcov1`) from `homogeneous_fit`.
3. Apply transformation `eta = g(beta_hat)` to natural parameters.
4. Compute Jacobian `J = d g / d beta`.
5. Propagate uncertainty:
   - model-based: `J %*% varcov0 %*% t(J)`
   - robust: `J %*% varcov1 %*% t(J)`
6. Report `se_model` and `se_robust` from diagonal terms.

## Behavior Change

Wrappers now use `reference = "first"` by default if `reference` is not supplied explicitly.

Rationale:

- keeps special-case synthetic parameterizations stable and interpretable,
- ensures deterministic natural-parameter recovery across runs.

If users pass `reference` explicitly, that value is respected.

## Documentation and Tests

- Wrapper documentation updated to reflect new outputs.
- Test suite expanded with:
  - non-regression checks for natural-parameter outputs,
  - consistency checks versus manual back-transformations.

## Validation Summary

Validated in repository:

- `devtools::test()` -> all tests pass.
- `R CMD build .` -> source tarball built.
- `R CMD check CircularRegression_0.4.0.tar.gz --no-manual` -> `Status: OK`.

## Versioning

- Package version bumped to `0.4.0`.
- Release date updated in `DESCRIPTION`.
- README version badge updated.
- Citation note updated to `R package version 0.4.0`.
