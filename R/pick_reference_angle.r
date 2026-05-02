#' Backward-compatible wrapper for reference-angle selection
#'
#' @param formula Model formula with terms \code{x} or \code{x:z}.
#' @param data Data frame containing model variables.
#' @param tie_method Tie-break rule for equal scores.
#' @param ... Unused.
#' @return A \code{select_reference_angle} object, with legacy class
#'   \code{pick_reference_angle} prepended.
#' @export
pick_reference_angle <- function(
  formula,
  data,
  tie_method = c("first", "random"),
  ...
) {
  tie_method <- match.arg(tie_method)
  .Deprecated("select_reference_angle")
  out <- select_reference_angle(formula = formula, data = data, tie_method = tie_method)
  class(out) <- c("pick_reference_angle", class(out))
  out
}

#' Print method for \code{pick_reference_angle}
#'
#' @param x Object returned by \code{pick_reference_angle()}.
#' @param ... Unused.
#' @export
print.pick_reference_angle <- function(x, ...) {
  NextMethod("print", x)
}
