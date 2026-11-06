# CircularRegression <img src="man/figures/logo.png" align="right" width="120" alt="CircularRegression hex logo" />

CircularRegression fits regression models for circular response data, such as
movement directions or angles measured in radians. The package implements the
general angular regression framework of Rivest, Duchesne, Nicosia and Fortin
(2016), including homogeneous angular regression, consensus regression, a
two-step workflow, and selected special-case wrappers.

[![CRAN version](https://img.shields.io/cran/v/CircularRegression?label=CRAN)](https://cran.r-project.org/package=CircularRegression)
[![GitHub release](https://img.shields.io/github/v/release/AurelienNicosiaULaval/CircularRegression?label=GitHub%20release)](https://github.com/AurelienNicosiaULaval/CircularRegression/releases/latest)
[![R-CMD-check](https://github.com/AurelienNicosiaULaval/CircularRegression/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/AurelienNicosiaULaval/CircularRegression/actions/workflows/R-CMD-check.yaml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-236660)](https://aureliennicosiaulaval.github.io/CircularRegression/)

## Installation

### CRAN

[CRAN](https://cran.r-project.org/package=CircularRegression) currently provides
version 0.5.1. The CRAN page has the reference manual and vignettes for that
version.

```r
install.packages("CircularRegression")
```

### GitHub release

Use the tagged [0.7.3 release](https://github.com/AurelienNicosiaULaval/CircularRegression/releases/tag/v0.7.3)
to follow the examples below and the guides on the package website. Install
`remotes` once with `install.packages("remotes")` before using either GitHub
command.

```r
remotes::install_github("AurelienNicosiaULaval/CircularRegression@v0.7.3")
```

### Development version

Install from the `main` branch to try the latest changes. The code on this
branch can change between installations.

```r
remotes::install_github("AurelienNicosiaULaval/CircularRegression@main")
```

## Main interface

```r
library(CircularRegression)

data(bison)
d <- bison[seq_len(100), ]

fit <- circular_regression(
  y.dir ~ y.prec + x.meadow:z.meadow,
  data = d
)

summary(fit)
coef(fit)
head(predict(fit))
vcov(fit, robust = TRUE)
fit$fit_status
fit$inference_available
```

The formula syntax uses angular variables as terms. A term of the form `x`
adds a direction directly. A term of the form `x:z` adds a direction `x`
weighted by a finite non-negative modifier `z`.

## Model-specific interfaces

The original interfaces remain available:

- `angular()` fits the homogeneous angular regression model.
- `consensus()` fits the consensus angular regression model.
- `angular_two_step()` fits the consensus model, selects a reference direction,
  and then fits the homogeneous model.
- `angular_re()` fits the random-intercept extension for clustered circular
  outcomes.

The package provides S3 methods for printing, summarising, coefficients,
fitted values, residuals, predictions, plots, information criteria, and
log-likelihoods where appropriate.

## Checking a fit

All response and explanatory angles must be in radians. Modifiers must be
finite and non-negative. The direction in `x:z` is `x`; reversing direction
and modifier changes the model.

Check `convergence`, `fit_status` and `inference_available` before interpreting
standard errors. A small score does not establish a finite maximum at a
concentration boundary. Unbounded likelihoods and unconverged fits do not
provide ordinary Wald standard errors. Fixed-effect sandwiches assume
independent observations; grouped sandwiches assume independent groups.

The [vignettes and quick guides](https://aureliennicosiaulaval.github.io/CircularRegression/articles/)
cover a first fit, model checking and comparisons, bison movement, special mean
curves, and repeated observations.

## Citation

To cite `CircularRegression` from R, use:

```r
citation("CircularRegression")
```

Citation metadata is also provided in `inst/CITATION`.

## References

Rivest, L.-P., Duchesne, T., Nicosia, A. and Fortin, D. (2016). A general
angular regression model for the analysis of data on animal movement in
ecology. Journal of the Royal Statistical Society: Series C (Applied
Statistics), 65(3), 445-463.

Rivest, L.-P. and Kato, S. (2019). A random-effects model for clustered
circular data. Canadian Journal of Statistics, 47(4), 712-728.

Related software: [ggcircular](https://github.com/AurelienNicosiaULaval/ggcircular)
for circular graphics and [gmov](https://github.com/AurelienNicosiaULaval/gmov)
for movement models.
