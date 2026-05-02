# CircularRegression

CircularRegression fits regression models for circular response data, such as
movement directions or angles measured in radians. The package implements the
general angular regression framework of Rivest, Duchesne, Nicosia and Fortin
(2016), including homogeneous angular regression, consensus regression, a
two-step workflow, and selected special-case wrappers.

## Installation

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/CircularRegression")
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

## References

Rivest, L.-P., Duchesne, T., Nicosia, A. and Fortin, D. (2016). A general
angular regression model for the analysis of data on animal movement in
ecology. Journal of the Royal Statistical Society: Series C (Applied
Statistics), 65(3), 445-463.

Rivest, L.-P. and Kato, S. (2019). A random-effects model for clustered
circular data. Canadian Journal of Statistics, 47(4), 712-728.
