#' Two-step angular regression workflow (consensus then homogeneous)
#'
#' Implements the strategy recommended by Rivest et al. (2016): fit the
#' consensus model first, choose the reference angle, then fit the homogeneous
#' model with the reference fixed and initialized from consensus ratios.
#'
#' @param formula Model formula with terms of the form \code{x} or \code{x:z}.
#' @param data Data frame containing model variables.
#' @param weights Optional non-negative weights passed to \code{consensus()}.
#' @param reference Reference-angle rule for the homogeneous step:
#'   \code{"auto"}, \code{"first"}, or \code{c("name", "x_ref")}.
#' @param control_consensus Control list for \code{consensus()}.
#' @param control_angular Control list for \code{angular()}.
#' @return An object of class \code{"angular_two_step"}.
#' @export
angular_two_step <- function(
  formula,
  data,
  weights = NULL,
  reference = c("auto", "first", "name"),
  control_consensus = list(),
  control_angular = list()
) {
  call <- match.call()

  fit_cons <- consensus(
    formula = formula,
    data = data,
    weights = weights,
    control = control_consensus
  )

  ref <- .resolve_reference(
    reference = reference,
    y = fit_cons$y,
    angle_data = fit_cons$angle_data,
    candidates = fit_cons$plain_angles,
    tie_method = "first"
  )
  ref_name <- ref$ref_name
  ref_idx <- .reference_term_index(fit_cons$term_info, ref_name)

  kappa_ref <- fit_cons$kappa[ref_idx]
  if (!is.finite(kappa_ref) || abs(kappa_ref) < .Machine$double.eps) {
    stop(
      "Cannot initialize homogeneous step from consensus: selected reference kappa is zero.",
      call. = FALSE
    )
  }

  initbeta <- as.numeric(fit_cons$kappa[-ref_idx] / kappa_ref)
  names(initbeta) <- fit_cons$paramname[-ref_idx]

  fit_hom <- angular(
    formula = formula,
    data = data,
    reference = c("name", ref_name),
    initbeta = initbeta,
    control = control_angular
  )

  out <- list(
    call = call,
    formula = formula,
    reference = ref_name,
    reference_scores = ref$scores,
    consensus_fit = fit_cons,
    homogeneous_fit = fit_hom,
    initbeta_from_consensus = initbeta
  )
  class(out) <- "angular_two_step"
  out
}

#' @export
print.angular_two_step <- function(x, ...) {
  cat("Two-step angular workflow (Rivest et al., 2016)\n")
  cat("Reference angle:", x$reference, "\n")
  cat(
    "Consensus logLik:",
    format(x$consensus_fit$logLik, digits = 6),
    "| Homogeneous max cosine:",
    format(x$homogeneous_fit$MaxCosine, digits = 6),
    "\n"
  )
  invisible(x)
}

#' @export
summary.angular_two_step <- function(object, ...) {
  out <- list(
    call = object$call,
    reference = object$reference,
    consensus = summary(object$consensus_fit),
    homogeneous = summary(object$homogeneous_fit)
  )
  class(out) <- "summary.angular_two_step"
  out
}

#' @export
print.summary.angular_two_step <- function(x, ...) {
  cat("Two-step angular workflow summary\n")
  cat("Reference angle:", x$reference, "\n\n")
  cat("Consensus step:\n")
  print(x$consensus)
  cat("\nHomogeneous step:\n")
  print(x$homogeneous)
  invisible(x)
}
