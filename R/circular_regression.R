#' Fit Circular Regression Models with a Common Interface
#'
#' \code{circular_regression()} is the recommended high-level interface for
#' fixed-effect circular regression models in this package. It keeps the
#' existing model-specific functions available while providing a common,
#' \code{glm()}-like entry point for the homogeneous angular model, the
#' consensus model, and the two-step workflow.
#'
#' @param formula Model formula with a circular response on the left-hand side
#'   and explanatory direction terms on the right-hand side. Terms must be of
#'   the form \code{x} or \code{x:z}, where \code{x} is an angular covariate and
#'   \code{z} is a finite non-negative modifier.
#' @param data Data frame, list, or environment containing model variables.
#' @param method Model fitting strategy. \code{"two_step"} first fits
#'   \code{\link{consensus}()} and then \code{\link{angular}()} with the selected
#'   reference angle. \code{"homogeneous"} calls \code{\link{angular}()}.
#'   \code{"consensus"} calls \code{\link{consensus}()}.
#' @param weights Optional non-negative observation weights. Used by
#'   \code{"two_step"} and \code{"consensus"}.
#' @param reference Reference-angle rule for homogeneous fits. Use
#'   \code{"auto"}, \code{"first"}, or \code{c("name", "x_ref")}.
#' @param control Optional control list passed to the selected model. For
#'   \code{method = "two_step"}, use \code{control_consensus} and
#'   \code{control_angular} for separate control.
#' @param control_consensus Control list for the consensus step.
#' @param control_angular Control list for the homogeneous step.
#' @param na.action Function used by \code{\link[stats]{model.frame}} to handle
#'   missing values. The default is \code{\link[stats]{na.omit}}.
#' @param ... Additional arguments passed to the underlying model-specific
#'   function.
#'
#' @return An object of class \code{"circular_regression"} containing the fitted
#'   model, the method used, the original call, and the formula.
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @examples
#' data(bison)
#' d <- bison[seq_len(100), ]
#' fit <- circular_regression(
#'   y.dir ~ y.prec + x.meadow:z.meadow,
#'   data = d
#' )
#' coef(fit)
#' head(predict(fit))
#' @export
circular_regression <- function(
  formula,
  data,
  method = c("two_step", "homogeneous", "consensus"),
  weights = NULL,
  reference = c("auto", "first", "name"),
  control = list(),
  control_consensus = list(),
  control_angular = list(),
  na.action = stats::na.omit,
  ...
) {
  call <- match.call()
  method <- match.arg(method)

  fit <- switch(
    method,
    two_step = angular_two_step(
      formula = formula,
      data = data,
      weights = weights,
      reference = reference,
      control_consensus = control_consensus,
      control_angular = control_angular,
      na.action = na.action,
      ...
    ),
    homogeneous = angular(
      formula = formula,
      data = data,
      reference = reference,
      control = control,
      na.action = na.action,
      ...
    ),
    consensus = consensus(
      formula = formula,
      data = data,
      weights = weights,
      control = control,
      na.action = na.action,
      ...
    )
  )

  out <- list(
    call = call,
    formula = formula,
    method = method,
    fit = fit,
    nobs = .fit_nobs(.primary_circular_fit(fit))
  )
  class(out) <- c(
    "circular_regression",
    paste0("circular_regression_", method)
  )
  out
}

.primary_circular_fit <- function(object) {
  if (inherits(object, "circular_regression")) {
    object <- object$fit
  }
  if (inherits(object, "angular_two_step")) {
    return(object$homogeneous_fit)
  }
  object
}

.fit_nobs <- function(object) {
  object$nobs %||% length(object$y)
}

#' Methods for Circular Regression Objects
#'
#' @param x,object A \code{"circular_regression"} object.
#' @param newdata Optional data frame for prediction.
#' @param type Prediction type passed to the underlying model.
#' @param k Numeric penalty used by \code{AIC()}.
#' @param ... Additional arguments passed to the underlying method.
#' @rdname circular_regression-methods
#' @export
print.circular_regression <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\n")
  cat("Number of observations:", x$nobs, "\n\n")
  print(.primary_circular_fit(x))
  invisible(x)
}

#' @rdname circular_regression-methods
#' @export
summary.circular_regression <- function(object, ...) {
  primary <- .primary_circular_fit(object)
  out <- list(
    call = object$call,
    method = object$method,
    nobs = object$nobs,
    model = summary(primary, ...)
  )
  class(out) <- "summary.circular_regression"
  out
}

#' @rdname circular_regression-methods
#' @export
print.summary.circular_regression <- function(x, ...) {
  cat("Circular regression summary\n")
  cat("Method:", x$method, "\n")
  cat("Number of observations:", x$nobs, "\n\n")
  print(x$model, ...)
  invisible(x)
}

#' @rdname circular_regression-methods
#' @export
coef.circular_regression <- function(object, ...) {
  coef(.primary_circular_fit(object), ...)
}

#' @rdname circular_regression-methods
#' @export
fitted.circular_regression <- function(object, ...) {
  fitted(.primary_circular_fit(object), ...)
}

#' @rdname circular_regression-methods
#' @export
residuals.circular_regression <- function(object, ...) {
  residuals(.primary_circular_fit(object), ...)
}

#' @rdname circular_regression-methods
#' @export
predict.circular_regression <- function(
  object,
  newdata = NULL,
  type = c("response", "components"),
  ...
) {
  type <- match.arg(type)
  predict(.primary_circular_fit(object), newdata = newdata, type = type, ...)
}

#' @rdname circular_regression-methods
#' @export
plot.circular_regression <- function(x, ...) {
  plot(.primary_circular_fit(x), ...)
  invisible(x)
}

#' @rdname circular_regression-methods
#' @export
logLik.circular_regression <- function(object, ...) {
  stats::logLik(.primary_circular_fit(object), ...)
}

#' @rdname circular_regression-methods
#' @export
AIC.circular_regression <- function(object, ..., k = 2) {
  stats::AIC(.primary_circular_fit(object), ..., k = k)
}

#' @rdname circular_regression-methods
#' @export
BIC.circular_regression <- function(object, ...) {
  stats::BIC(.primary_circular_fit(object), ...)
}
