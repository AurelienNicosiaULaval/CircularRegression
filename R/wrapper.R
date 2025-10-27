# Internal helper --------------------------------------------------------------
.safe_name <- function(x) make.names(x, unique = TRUE)

#' Mean Direction Model Wrapper
#'
#' Wraps estimation of the "mean direction model" by automatically creating the
#' required covariates and the appropriate formula.
#'
#' Two columns are added:
#' \itemize{
#'   \item \code{mu0}: 0 if \code{mu0} is \code{NULL} (default), otherwise the provided value;
#'   \item \code{tan_mu}: \eqn{\pi/2} if \code{mu0} is \code{NULL}, otherwise \code{mu0 + pi/2}.
#' }
#' The coefficient of \code{mu0} is fixed at 1 by \code{angular()}, and the
#' estimated coefficient of \code{tan_mu} corresponds to \eqn{\tan(\mu)}.
#'
#' @param data A data frame containing at least the response variable.
#' @param response A string: name of the angular response variable.
#' @param mu0 Optional numeric (scalar or vector). If \code{NULL}, the baseline is
#'   0 and \code{tan_mu = pi/2}; otherwise \code{tan_mu = mu0 + pi/2}.
#' @param ... Additional arguments passed to \code{angular()}.
#'
#' @return A \strong{list} with:
#' \itemize{
#'   \item \code{fit}: the object returned by \code{angular()};
#'   \item \code{data_aug}: the augmented data;
#'   \item \code{formula}: the formula used in the call to \code{angular()}.
#' }
#' @export
meanDirectionModel <- function(data, response, mu0 = NULL, ...) {
  data_aug <- data
  if (is.null(mu0)) {
    data_aug$mu0 <- 0
    data_aug$tan_mu <- pi / 2
  } else {
    data_aug$mu0 <- mu0
    data_aug$tan_mu <- mu0 + pi / 2
  }
  modelFormula <- stats::as.formula(paste(
    .safe_name(response),
    "~ mu0 + tan_mu"
  ))
  out <- list(
    fit = angular(formula = modelFormula, data = data_aug, ...),
    data_aug = data_aug,
    formula = modelFormula
  )
  return(out)
}

#' Decentered Predictor Model Wrapper
#'
#' Wraps the "decentered predictor" model (Rivest, 1997) for a single angular
#' explanatory variable \code{w}. The augmented data include:
#' \itemize{
#'   \item \code{<w>_plus_pi2}: \eqn{w + \pi/2};
#'   \item \code{const0}: 0;
#'   \item \code{constPi2}: \eqn{\pi/2}.
#' }
#' The fitted formula is:
#' \code{response ~ <w> + <w>_plus_pi2 + const0 + constPi2}.
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response A string: name of the angular response variable.
#' @param w A string: name of the angular explanatory variable in \code{data}.
#' @param ... Additional arguments passed to \code{angular()}.
#'
#' @return A \strong{list}: \code{fit}, \code{data_aug}, \code{formula}.
#' @export
decentredPredictorModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The 'w' variable provided (", w, ") is not present in 'data'.")
  }
  w_safe <- .safe_name(w)
  data_aug <- data
  data_aug[[paste0(w_safe, "_plus_pi2")]] <- data_aug[[w]] + pi / 2
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  rhs <- paste(
    w_safe,
    paste0(w_safe, "_plus_pi2"),
    "const0",
    "constPi2",
    sep = " + "
  )
  modelFormula <- stats::as.formula(paste(.safe_name(response), "~", rhs))

  out <- list(
    fit = angular(formula = modelFormula, data = data_aug, ...),
    data_aug = data_aug,
    formula = modelFormula
  )
  return(out)
}

#' Presnell Model Wrapper
#'
#' Wraps the Presnell et al. (1998) model with one non-circular continuous covariate \code{w}.
#' The augmented data include:
#' \itemize{
#'   \item \code{const0}: 0;
#'   \item \code{constPi2}: \eqn{\pi/2}.
#' }
#' The fitted formula is:
#' \code{response ~ const0 + constPi2 + const0:<w> + constPi2:<w>}.
#'
#' @param data A data frame containing \code{response} and \code{w} (continuous).
#' @param response A string: name of the angular response variable.
#' @param w A string: name of the continuous covariate in \code{data}.
#' @param ... Additional arguments passed to \code{angular()}.
#'
#' @return A \strong{list}: \code{fit}, \code{data_aug}, \code{formula}.
#' @export
presnellModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The 'w' variable provided (", w, ") is not present in 'data'.")
  }
  w_safe <- .safe_name(w)
  data_aug <- data
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  # response ~ const0 + constPi2 + const0:w + constPi2:w
  rhs <- paste(
    "const0 + constPi2",
    paste0("const0:", w_safe),
    paste0("constPi2:", w_safe),
    sep = " + "
  )
  modelFormula <- stats::as.formula(paste(.safe_name(response), "~", rhs))

  out <- list(
    fit = angular(formula = modelFormula, data = data_aug, ...),
    data_aug = data_aug,
    formula = modelFormula
  )
  return(out)
}

#' Jammalamadaka–Sengupta Model Wrapper
#'
#' Wraps the Jammalamadaka–Sengupta circular regression with one angular predictor \code{w}.
#' The augmented data include:
#' \itemize{
#'   \item \code{const0}: 0; \code{constPi2}: \eqn{\pi/2}
#'   \item \code{<w>_orig}: \eqn{w}
#'   \item \code{neg_<w>}: \eqn{-w}
#'   \item \code{<w>_pi2}: \eqn{w + \pi/2}
#'   \item \code{neg_<w>_pi2}: \eqn{-w + \pi/2}
#' }
#' The fitted formula is:
#' \code{response ~ const0 + constPi2 + <w>_orig + neg_<w> + <w>_pi2 + neg_<w>_pi2}.
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response A string: name of the angular response variable.
#' @param w A string: name of the angular explanatory variable in \code{data}.
#' @param ... Additional arguments passed to \code{angular()}.
#'
#' @return A \strong{list}: \code{fit}, \code{data_aug}, \code{formula}.
#' @export
jammalamadakaModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The 'w' variable provided (", w, ") is not present in 'data'.")
  }
  w_safe <- .safe_name(w)
  data_aug <- data
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2
  data_aug[[paste0(w_safe, "_orig")]] <- data_aug[[w]]
  data_aug[[paste0("neg_", w_safe)]] <- -data_aug[[w]]
  data_aug[[paste0(w_safe, "_pi2")]] <- data_aug[[w]] + pi / 2
  data_aug[[paste0("neg_", w_safe, "_pi2")]] <- -data_aug[[w]] + pi / 2

  rhs <- paste(
    "const0",
    "constPi2",
    paste0(w_safe, "_orig"),
    paste0("neg_", w_safe),
    paste0(w_safe, "_pi2"),
    paste0("neg_", w_safe, "_pi2"),
    sep = " + "
  )
  modelFormula <- stats::as.formula(paste(.safe_name(response), "~", rhs))

  out <- list(
    fit = angular(formula = modelFormula, data = data_aug, ...),
    data_aug = data_aug,
    formula = modelFormula
  )
  return(out)
}
