# CircularRegression 0.7.1

* Correct the consensus score for modified directions and its information limit
  at zero resultant length. Convergence now requires a small scaled score;
  step-halving stagnation alone is not convergence. Store convergence diagnostics.
* Correct homogeneous sandwich covariance using the full observed derivative
  in Rivest et al. (2016), equations 13-14. Recompute natural-parameter robust
  covariances using this corrected matrix.
* Implement `vcov.consensus(robust = TRUE)` using the weighted observation
  scores; preserve the legacy model covariance in `varcov1`.
* Retain existing affine and harmonic ratios and add `natural_orientation`.
  Both reconstructed vector components must be multiplied by this sign.
  Support negative continuous covariates in the affine wrapper.
* Solve the homogeneous concentration likelihood equation numerically. Report
  infinite concentration and unavailable Wald inference at perfect alignment.
* Require compatible, converged, nested fits on identical data for T2 tests.
* Include the normalizing constant in grouped log-likelihoods, flag boundary
  inference, and expose cluster scores, observed information and convergence.
  Replace singleton-group Sandhopper examples by repeated trials.
* Validate numeric inputs before coercion, long lags, integer controls, ignored
  arguments, and new prediction groups. Return NA for undefined mean directions.
* Preserve direction/modifier roles when R reorders interaction labels; share
  this parser with grouped fitting and new-data prediction.
* Preserve PIT endpoints and use scaled von Mises quadrature at high concentration.
* Use the public ggcircular 0.1.0 API; remove the development-revision dependency.

# CircularRegression 0.7.0

* Replace the normal residual histogram and normal Q-Q plot with circular
  diagnostics for homogeneous and consensus angular models: a residual rose,
  a fitted von Mises density overlay when appropriate, a probability-integral
  transform (PIT) Q-Q plot, and residual plots against fitted direction and
  observation order.
* Reimplement the fixed and conditional residual roses of `angular_re()` with
  `ggplot2` and `ggcircular`, while retaining the established S3 `plot()`
  interface and the distinction documented by Rivest and Kato (2019).
* Make paired `angular_re()` roses use one pooled kernel bandwidth, identical
  angular bins, and a common radial limit; expose the optional `bandwidth`
  argument for reproducible comparisons.
* Warn when a special-case natural-parameter transformation is numerically
  singular and return the evaluated threshold in `natural_stability`.
* Guard the consensus-to-homogeneous normalization with a scale-relative
  reference-coefficient threshold and expose it in `reference_stability`.
* Preserve the orientation of a negative consensus reference coefficient by
  fitting the homogeneous model in a signed reference chart. This prevents an
  unintended pi shift in two-step fits and special-case wrappers.
* Apply the same signed, scale-relative normalization to
  `coef.consensus(type = "beta")` and record the reference sign and stability
  diagnostics as attributes.

# CircularRegression 0.6.0

* Replace numerical differentiation in the special-case wrappers with exact
  analytic Jacobians for transformed parameters and their delta-method
  covariance matrices.
* Correct the first-harmonic circular-circular wrapper to return named cosine
  and sine components of the fitted mean curve.
* Clarify that `presnellModel()` fits an affine resultant-vector mean curve,
  rather than the projected-normal likelihood of Presnell, Morrison and
  Littell (1998).
* Document the baseline interpretation of `meanDirectionModel()` and the
  non-identification boundary of `decentredPredictorModel()`.

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
* Clarify documented provenance for `noshiro` and the remaining provenance
  limitations for `multiplebison`.
* Add a related scientific reference for the `multiplebison` study context.
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
