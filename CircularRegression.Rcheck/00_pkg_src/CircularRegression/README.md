<img src="figure/hex-logo.png" align="right" width="180" alt="CircularRegression hex logo" />

# CircularRegression

[![Version](https://img.shields.io/badge/version-0.4.0-2A9D8F)](https://github.com/AurelienNicosiaULaval/CircularRegression)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-0B2239.svg)](https://opensource.org/license/gpl-3-0)
[![R >= 3.5](https://img.shields.io/badge/R-%3E%3D%203.5-276DC3.svg)](https://www.r-project.org/)

`CircularRegression` provides regression methods for circular responses (angles in radians), with an API aligned to Rivest, Duchesne, Nicosia and Fortin (2016, JRSS C). It is designed for movement ecology and other directional-data settings where standard linear models are inappropriate.

## Scope

The package implements:

- The homogeneous angular regression model via `angular()`.
- The consensus model via `consensus()`.
- Article-based reference-angle selection via `select_reference_angle()`.
- The recommended two-step workflow via `angular_two_step()`:
  consensus fit, reference selection, then homogeneous fit.
- Specialized wrappers (`meanDirectionModel()`, `decentredPredictorModel()`, `presnellModel()`, `jammalamadakaModel()`).
- Core S3 methods (`print`, `summary`, `coef`, `residuals`, `fitted`, `logLik`, `AIC`, `BIC`, model comparison utilities).

## Installation

Development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/CircularRegression")
```

## Quick Start (Recommended Workflow)

```r
library(CircularRegression)
data(bison)

# Keep this run short in interactive sessions
d <- bison[seq_len(600), ]

form <- y.dir ~ y.prec + y.prec2 + x.meadow + x.meadow:z.meadow + x.gap + x.gap:z.gap

fit <- angular_two_step(
  formula = form,
  data = d,
  control_consensus = list(maxiter = 100),
  control_angular = list(maxiter = 100)
)

fit$reference
coef(fit$consensus_fit, type = "kappa")
coef(fit$homogeneous_fit)
```

## Fit Models Separately

```r
library(CircularRegression)
data(bison)

d <- bison[seq_len(600), ]
form <- y.dir ~ y.prec + y.prec2 + x.meadow + x.meadow:z.meadow + x.gap + x.gap:z.gap

fit_consensus <- consensus(formula = form, data = d)
fit_homogeneous <- angular(formula = form, data = d, reference = "auto")

summary(fit_consensus)
summary(fit_homogeneous)
```

## Included Data

- `bison`: bison movement directions with landscape covariates.
- `multiplebison`: synchronized movement features for two tracked bison.
- `noshiro`: earthquake-related directional variables.
- `Sandhopper`: repeated-orientation escape experiment data.

All modeling functions assume angular quantities are expressed in radians.

## Migration Notes for 0.4.0

- Special-case wrappers now return natural-parameter summaries and delta-method SEs:
  - `natural_parameters` (`estimate`, `se_model`, `se_robust`)
  - `natural_vcov_model`, `natural_vcov_robust`, `natural_jacobian`
- Wrappers default to `reference = "first"` if `reference` is not provided explicitly.

## Migration Notes for 0.3.0

- `consensus(initbeta = ...)` is replaced by `consensus(initkappa = ...)`.
- `angular()` now uses explicit reference strategies: `"auto"`, `"first"`, or `c("name", "<angle>")`.
- `select_reference_angle()` is the primary reference-selection helper.
- `pick_reference_angle()` is retained as a deprecated compatibility wrapper.
- `angular_two_step()` is the documented default analysis path.

## Documentation

Use package help pages directly:

```r
?angular_two_step
?angular
?consensus
?select_reference_angle
```

## Citation

Methodological reference:

Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016). A general angular regression model for the analysis of data on animal movement in ecology. *Journal of the Royal Statistical Society: Series C (Applied Statistics)*, 65(3), 445-463.

Package citation:

```r
citation("CircularRegression")
```

## Issues and Contributions

- Bug reports and feature requests:
  <https://github.com/AurelienNicosiaULaval/CircularRegression/issues>
- Contributions are welcome through pull requests with a minimal reproducible example and tests when relevant.
