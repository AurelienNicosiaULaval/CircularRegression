# CircularRegression 0.3.0 - Migration Changelog

This document summarizes user-facing changes introduced in `v0.3.0` and provides migration guidance from `v0.2.x`.

## Scope

Version `0.3.0` aligns package behavior with Rivest et al. (2016):
- explicit separation between consensus and homogeneous error models,
- reference-angle selection based on empirical mean cosine,
- two-step workflow (`consensus` then `angular`) as recommended strategy.

## Breaking Changes

### 1) `angular()` API changed

Old:
```r
angular(formula, data, initbeta = NULL, control = list())
```

New:
```r
angular(
  formula,
  data,
  reference = c("auto", "first", "name"),
  initbeta = NULL,
  control = list()
)
```

What changed:
- reference handling is now explicit.
- default `reference = "auto"` selects the plain angle term maximizing
  `mean(cos(y - x_j))`.
- model identifiability is checked explicitly; non-identifiable designs now error early.

---

### 2) `consensus()` initialization changed

Old:
```r
consensus(formula, data, weights = NULL, initbeta = NULL, control = list())
```

New:
```r
consensus(formula, data, weights = NULL, initkappa = NULL, control = list())
```

What changed:
- consensus is parameterized directly in terms of `kappa_j` (article form).
- `initbeta` was removed; use `initkappa`.

---

### 3) `coef.consensus()` changed

Old behavior: single implicit coefficient extraction.

New behavior:
```r
coef(object, type = c("kappa", "beta"), reference = c("auto", "first", "name"))
```

- `type = "kappa"` (default): returns estimated `kappa_j`.
- `type = "beta"`: returns scale-free ratios relative to selected reference.

---

### 4) Reference-selection API replaced

- New: `select_reference_angle(formula, data, tie_method = c("first", "random"))`.
- Legacy: `pick_reference_angle()` remains available as deprecated wrapper for one transition version.

---

### 5) Wrappers now return two-step outputs

The following wrappers now return a standardized two-step object payload instead of a single homogeneous fit:
- `meanDirectionModel()`
- `decentredPredictorModel()`
- `presnellModel()`
- `jammalamadakaModel()`

Returned list keys:
- `consensus_fit`
- `homogeneous_fit`
- `reference`
- `data_aug`
- `formula`

---

## New Features

### 1) Two-step workflow

New function:
```r
angular_two_step(
  formula,
  data,
  weights = NULL,
  reference = c("auto", "first", "name"),
  control_consensus = list(),
  control_angular = list()
)
```

Pipeline:
1. fit consensus model,
2. select reference angle,
3. initialize homogeneous coefficients from consensus ratios,
4. fit homogeneous model,
5. return both fits in one object.

### 2) New S3 classes/methods

- `select_reference_angle`
- `angular_two_step`
- `summary.angular_two_step`
- `print.angular_two_step`

---

## Bug Fixes

### 1) `maxiter = 0` no longer crashes

Fixed in both:
- `angular()`
- `consensus()`

### 2) `angular_re` prediction no longer depends on caller environment

`predict.angular_re()` now uses model structures stored in fitted objects; calls remain valid even if original formula/data symbols are removed from the calling environment.

### 3) NA handling alignment for `cluster` in `angular_re`

`cluster` is now aligned with model-frame filtering behavior to avoid silent length mismatches.

### 4) Model-comparison consistency

- stabilized `logLik.angular()` with scaled Bessel evaluation,
- guaranteed `df` and `nobs` attributes for compatibility with `AIC/BIC/logLik` generics,
- `BIC.consensus` and stored `fit$BIC` are consistent.

---

## Migration Guide

### A) `consensus(initbeta = ...)` -> `consensus(initkappa = ...)`

Before:
```r
fit <- consensus(y ~ x1 + x2:z2, data = d, initbeta = c(0.5, 0.1))
```

After:
```r
fit <- consensus(y ~ x1 + x2:z2, data = d, initkappa = c(0.5, 0.1))
```

### B) Explicit reference control in `angular()`

Before:
```r
fit <- angular(y ~ x1 + x2:z2, data = d)
```

After (same call still valid, but now explicit if desired):
```r
fit <- angular(y ~ x1 + x2:z2, data = d, reference = "auto")
# or
fit <- angular(y ~ x1 + x2:z2, data = d, reference = c("name", "x1"))
```

### C) Getting consensus beta-like coefficients

Before:
```r
coef(fit)
```

After:
```r
coef(fit, type = "kappa")
coef(fit, type = "beta", reference = c("name", "x1"))
```

### D) Replace legacy reference helper

Before:
```r
pick_reference_angle(y ~ x1 + x2, data = d)
```

After:
```r
select_reference_angle(y ~ x1 + x2, data = d)
```

### E) Prefer two-step workflow

New recommended usage:
```r
fit2 <- angular_two_step(y ~ x1 + x2:z2, data = d)
fit2$consensus_fit
fit2$homogeneous_fit
```

---

## Validation Summary

Validated in repository:
- local test suite passes,
- package checks pass (with expected environment-related NOTES only when pandoc is unavailable).

---

## Deprecation Timeline

- `pick_reference_angle()` is currently soft-deprecated.
- Planned removal window: after one transition release cycle (next major/minor policy decision).
