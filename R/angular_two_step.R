#' Two-step angular regression workflow (consensus then homogeneous)
#'
#' Implements the strategy recommended by Rivest et al. (2016): fit the
#' consensus model first, choose the reference angle, then fit the homogeneous
#' model with the reference fixed and initialized from consensus ratios.
#'
#' @param formula Model formula with terms of the form \code{x} or \code{x:z}.
#' @param data Data frame containing model variables.
#' @param weights Optional non-negative weights passed only to \code{consensus()}
#'   for initialization. The final homogeneous fit is unweighted.
#' @param reference Reference-angle rule for the homogeneous step:
#'   \code{"auto"}, \code{"first"}, or \code{c("name", "x_ref")}.
#' @param control_consensus Control list for \code{consensus()}.
#' @param control_angular Control list for \code{angular()}.
#' @param na.action Function used by \code{\link[stats]{model.frame}} to handle
#'   missing values. The default is \code{\link[stats]{na.omit}}.
#' @return An object of class \code{"angular_two_step"}. The
#'   \code{reference_stability} component records the selected consensus
#'   coefficient, its scale-relative threshold, and the signed reference chart
#'   used to preserve the orientation of the consensus resultant.
#'   \code{convergence}, \code{fit_status}, \code{inference_available} and
#'   \code{score_norm} describe the final homogeneous fit. An unbounded
#'   consensus initialization stops the two-step procedure.
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @examples
#' data(bison)
#' fit <- angular_two_step(y.dir ~ y.prec + x.meadow:z.meadow, bison[1:100, ])
#' fit$homogeneous_fit$convergence
#' head(predict(fit))
#' @export
angular_two_step <- function(
  formula,
  data,
  weights = NULL,
  reference = c("auto", "first", "name"),
  control_consensus = list(),
  control_angular = list(),
  na.action = stats::na.omit
) {
  call <- match.call()

  fit_cons <- consensus(
    formula = formula,
    data = data,
    weights = weights,
    control = control_consensus,
    na.action = na.action
  )
  if (isTRUE(fit_cons$alignment_boundary)) {
    stop("The consensus step has no finite maximum; use angular() directly to specify a homogeneous fit.", call. = FALSE)
  }

  ref <- .resolve_reference(
    reference = reference,
    y = fit_cons$y,
    angle_data = fit_cons$angle_data,
    candidates = fit_cons$plain_angles,
    tie_method = "first"
  )
  ref_name <- ref$ref_name
  ref_idx <- .reference_term_index(fit_cons$term_info, ref_name)

  reference_stability <- .reference_coefficient_stability(
    kappa = fit_cons$kappa,
    ref_idx = ref_idx
  )
  kappa_ref <- reference_stability$coefficient
  if (!isTRUE(reference_stability$stable)) {
    stop(
      "Cannot initialize homogeneous step from consensus: the absolute selected ",
      "reference coefficient (",
      format(abs(kappa_ref), digits = 4),
      ") does not exceed the scale-relative threshold (",
      format(reference_stability$threshold, digits = 4),
      "). Choose another reference angle.",
      call. = FALSE
    )
  }
  reference_sign <- reference_stability$reference_sign
  initbeta <- as.numeric(fit_cons$kappa[-ref_idx] / abs(kappa_ref))
  names(initbeta) <- fit_cons$paramname[-ref_idx]

  fit_hom <- angular(
    formula = formula,
    data = data,
    reference = c("name", ref_name),
    reference_sign = reference_sign,
    initbeta = initbeta,
    control = control_angular,
    na.action = na.action
  )

  out <- list(
    call = call,
    formula = formula,
    reference = ref_name,
    reference_scores = ref$scores,
    reference_stability = reference_stability,
    consensus_fit = fit_cons,
    homogeneous_fit = fit_hom,
    initbeta_from_consensus = initbeta,
    nobs = fit_hom$nobs,
    convergence = fit_hom$convergence,
    fit_status = fit_hom$fit_status,
    inference_available = fit_hom$inference_available,
    score_norm = fit_hom$score_norm
  )
  class(out) <- "angular_two_step"
  out
}

.reference_coefficient_stability <- function(
  kappa,
  ref_idx,
  tolerance = sqrt(.Machine$double.eps)
) {
  if (!is.numeric(kappa) || length(kappa) < 1L || any(!is.finite(kappa))) {
    stop("'kappa' must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (length(ref_idx) != 1L || !is.numeric(ref_idx) || !is.finite(ref_idx) ||
      ref_idx != as.integer(ref_idx) || ref_idx < 1L || ref_idx > length(kappa)) {
    stop("'ref_idx' must identify one element of 'kappa'.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("'tolerance' must be a single finite positive number.", call. = FALSE)
  }

  coefficient <- unname(kappa[[as.integer(ref_idx)]])
  coefficient_scale <- sqrt(sum(kappa^2))
  threshold <- tolerance * coefficient_scale

  list(
    stable = coefficient_scale > 0 && abs(coefficient) > threshold,
    reference_sign = if (coefficient < 0) -1 else 1,
    coefficient = coefficient,
    coefficient_scale = coefficient_scale,
    threshold = threshold,
    tolerance = tolerance
  )
}

#' Methods for two-step angular regression objects
#'
#' @param x,object An object of class \code{"angular_two_step"}.
#' @param step Which fitted step to use for \code{coef()} or \code{vcov()}.
#' @param robust Logical; request the covariance based on independent observation
#'   scores rather than the model covariance.
#' @param newdata Optional data frame for prediction.
#' @param type Prediction type passed to \code{\link{predict.angular}()}.
#' @param k Numeric penalty used by \code{AIC()}.
#' @param ... Additional arguments passed to the underlying method.
#' @rdname angular_two_step-methods
#' @examples
#' data(bison)
#' fit <- angular_two_step(y.dir ~ y.prec, bison[1:100, ])
#' coef(fit, step = "consensus")
#' coef(fit)
#' summary(fit)
#' @export
print.angular_two_step <- function(x, ...) {
  cat("Two-step angular workflow (Rivest et al., 2016)\n")
  .print_fit_status(x)
  cat("Reference angle:", x$reference, "\n")
  if (!is.null(x$reference_stability)) {
    cat(
      "Reference coefficient:",
      format(x$reference_stability$coefficient, digits = 6),
      "| Relative threshold:",
      format(x$reference_stability$threshold, digits = 3),
      "\n"
    )
  }
  cat(
    "Consensus logLik:",
    format(x$consensus_fit$logLik, digits = 6),
    "| Homogeneous max cosine:",
    format(x$homogeneous_fit$MaxCosine, digits = 6),
    "\n"
  )
  invisible(x)
}

#' @rdname angular_two_step-methods
#' @export
summary.angular_two_step <- function(object, ...) {
  out <- list(
    call = object$call,
    reference = object$reference,
    reference_stability = object$reference_stability,
    consensus = summary(object$consensus_fit),
    homogeneous = summary(object$homogeneous_fit)
  )
  class(out) <- "summary.angular_two_step"
  out
}

#' @rdname angular_two_step-methods
#' @export
print.summary.angular_two_step <- function(x, ...) {
  cat("Two-step angular workflow summary\n")
  cat("Reference angle:", x$reference, "\n\n")
  if (!is.null(x$reference_stability)) {
    cat(
      "Reference coefficient:",
      format(x$reference_stability$coefficient, digits = 6),
      "| Relative threshold:",
      format(x$reference_stability$threshold, digits = 3),
      "\n\n"
    )
  }
  cat("Consensus step:\n")
  print(x$consensus)
  cat("\nHomogeneous step:\n")
  print(x$homogeneous)
  invisible(x)
}

#' @rdname angular_two_step-methods
#' @export
coef.angular_two_step <- function(
  object,
  step = c("homogeneous", "consensus"),
  ...
) {
  step <- match.arg(step)
  if (step == "consensus") {
    return(coef(object$consensus_fit, ...))
  }
  coef(object$homogeneous_fit, ...)
}

#' @rdname angular_two_step-methods
#' @export
vcov.angular_two_step <- function(object, robust = FALSE,
                                  step = c("homogeneous", "consensus"), ...) {
  step <- match.arg(step)
  fit <- if (step == "consensus") object$consensus_fit else object$homogeneous_fit
  stats::vcov(fit, robust = robust, ...)
}

#' @rdname angular_two_step-methods
#' @export
fitted.angular_two_step <- function(object, ...) {
  fitted(object$homogeneous_fit, ...)
}

#' @rdname angular_two_step-methods
#' @export
residuals.angular_two_step <- function(object, ...) {
  residuals(object$homogeneous_fit, ...)
}

#' @rdname angular_two_step-methods
#' @export
predict.angular_two_step <- function(
  object,
  newdata = NULL,
  type = c("response", "components"),
  ...
) {
  type <- match.arg(type)
  predict(object$homogeneous_fit, newdata = newdata, type = type, ...)
}

#' @rdname angular_two_step-methods
#' @export
plot.angular_two_step <- function(x, ...) {
  plot(x$homogeneous_fit, ...)
  invisible(x)
}

#' @rdname angular_two_step-methods
#' @export
logLik.angular_two_step <- function(object, ...) {
  stats::logLik(object$homogeneous_fit, ...)
}

#' @rdname angular_two_step-methods
#' @export
AIC.angular_two_step <- function(object, ..., k = 2) {
  NextMethod("AIC")
}

#' @rdname angular_two_step-methods
#' @export
BIC.angular_two_step <- function(object, ...) {
  NextMethod("BIC")
}
