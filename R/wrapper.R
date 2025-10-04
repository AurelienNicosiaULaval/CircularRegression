#' Mean Direction Model Wrapper
#'
#' This function wraps the estimation of the "mean direction model" by automatically
#' creating the required covariates and constructing the appropriate formula.
#' It generates two covariates:
#' \itemize{
#'   \item{\code{mu0}: If \code{mu0} is \code{NULL} (default), this is set to 0. Otherwise, it is set to the provided value (or vector).}
#'   \item{\code{tan_mu}: If \code{mu0} is \code{NULL}, this is set to \eqn{\pi/2}; otherwise, it is set to \code{mu0 + pi/2}.}
#' }
#'
#' The model is then specified as:
#' \deqn{y = 1 \times mu0 + tan\_mu,}
#' where the coefficient for \code{mu0} is fixed at 1 (by design in \code{angular}),
#' and the estimated coefficient for \code{tan_mu} corresponds to \eqn{\tan(\mu)}.
#'
#' @param data A data frame containing at least the response variable.
#' @param response A character string representing the name of the response variable (an angular variable).
#' @param mu0 Optional. A numeric vector or value to be used as the baseline covariate.
#'        If \code{NULL} (default), the baseline is set to 0 and \code{tan_mu} is set to pi/2.
#'        Otherwise, \code{tan_mu} is set to \code{mu0 + pi/2}.
#' @param ... Additional arguments to be passed to the \code{angular} function.
#'
#' @return An object representing the fitted model as returned by the \code{angular} function.
#'
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

  modelFormula <- stats::as.formula(paste(response, "~ mu0 + tan_mu"))

  fit <- angular(formula = modelFormula, data = data_aug, ...)

  return(fit)
}

#' Decentered Predictor Model Wrapper
#'
#' This function wraps the estimation of the "decentered predictor" model as described in Rivest (1997),
#' which is based on a single angular explanatory variable. The synthetic model is built using four covariates:
#' \itemize{
#'   \item{\code{w}: the original explanatory angular variable,}
#'   \item{\code{w_plus_pi2}: computed as \eqn{w + \pi/2},}
#'   \item{\code{zero}: a constant covariate with value 0,}
#'   \item{\code{pi2}: a constant covariate with value \eqn{\pi/2}.}
#' }
#'
#' The model is expressed as:
#' \deqn{y = w + \left(w + \frac{\pi}{2}\right) + 0 + \left(\frac{\pi}{2}\right),}
#' which corresponds to the decentered predictor formulation.
#'
#' @param data A data frame containing at least the response variable and the explanatory variable.
#' @param response A character string representing the name of the response variable (an angular variable).
#' @param w A character string representing the name of the explanatory angular variable in \code{data}.
#' @param ... Additional arguments to be passed to the underlying \code{angular} estimation function.
#'
#' @return An object representing the fitted model as returned by the \code{angular} function.
#'
#' @export
decentredPredictorModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The explanatory variable specified in 'w' is not present in the data.")
  }

  data_aug <- data
  data_aug$w_plus_pi2 <- data_aug[[w]] + pi / 2
  data_aug$zero <- 0
  data_aug$pi2 <- pi / 2

  modelFormula <- stats::as.formula(paste(response, "~", w, "+ w_plus_pi2 + zero + pi2"))

  fit <- angular(formula = modelFormula, data = data_aug, ...)

  return(fit)
}


#' Presnell Model Wrapper
#'
#' This function wraps the estimation of the "Presnell model" as described in Presnell et al. (1998),
#' which incorporates one non-circular continuous covariate (w). The synthetic model is specified using
#' the following form:
#' \deqn{y = 0 + \frac{\pi}{2} + (0 : w) + \left(\frac{\pi}{2} : w\right),}
#' where:
#' \itemize{
#'   \item The term 0 is represented by a constant covariate,
#'   \item The term \eqn{\pi/2} is represented by a constant covariate,
#'   \item 0 : w denotes the interaction between the constant 0 and the continuous covariate w,
#'   \item \eqn{(\pi/2) : w} denotes the interaction between the constant \eqn{\pi/2} and w.
#' }
#'
#' In order for the formula to be interpreted correctly by the underlying \code{angular} function,
#' this wrapper first creates two variables in \code{data}: \code{zero} (set to 0) and \code{pi2} (set to pi/2).
#' The model is then specified as:
#' \deqn{y ~ zero + pi2 + zero:w + pi2:w.}
#'
#' @param data A data frame containing at least the response variable and the continuous covariate w.
#' @param response A character string representing the name of the response variable (an angular variable).
#' @param w A character string representing the name of the continuous covariate in \code{data}.
#' @param ... Additional arguments to be passed to the underlying \code{angular} estimation function.
#'
#' @return An object representing the fitted model as returned by the \code{angular} function.
#'
#' @export
presnellModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The variable specified in 'w' is not present in the data.")
  }

  data_aug <- data
  data_aug$zero <- 0
  data_aug$pi2 <- pi / 2

  modelFormula <- stats::as.formula(paste(response, "~ zero + pi2 + zero:", w, " + pi2:", w, sep = ""))

  fit <- angular(formula = modelFormula, data = data_aug, ...)

  return(fit)
}

#' Jammalamadaka-Sengupta Model Wrapper
#'
#' This function wraps the estimation of the circular regression model following the
#' Jammalamadaka-Sengupta formulation (see end of page 447). In this formulation the model is given by:
#' \deqn{y = 0 + \frac{\pi}{2} + w + (-w) + \left(w+\frac{\pi}{2}\right) + \left(-w+\frac{\pi}{2}\right),}
#' where w is an angular explanatory variable. To implement this model, six covariates are created:
#' \itemize{
#'   \item{\code{const0}: a constant set to 0,}
#'   \item{\code{constPi2}: a constant set to \eqn{\pi/2},}
#'   \item{\code{w_orig}: the original value of w,}
#'   \item{\code{neg_w}: the negative of w,}
#'   \item{\code{w_pi2}: w shifted by \eqn{\pi/2} (i.e., w + pi/2),}
#'   \item{\code{neg_w_pi2}: the negative of w shifted by \eqn{\pi/2} (i.e., -w + pi/2).}
#' }
#'
#' The resulting model is estimated using the \code{angular} function with the formula:
#' \deqn{y \sim const0 + constPi2 + w\_orig + neg\_w + w\_pi2 + neg\_w\_pi2.}
#'
#' @param data A data frame containing at least the response variable and the angular explanatory variable.
#' @param response A character string representing the name of the response variable (an angular variable).
#' @param w A character string representing the name of the angular explanatory variable in \code{data}.
#' @param ... Additional arguments to be passed to the underlying \code{angular} estimation function.
#'
#' @return An object representing the fitted model as returned by the \code{angular} function.
#'
#' @export
jammalamadakaModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The explanatory variable specified in 'w' is not present in the data.")
  }

  data_aug <- data
  data_aug$const0    <- 0
  data_aug$constPi2  <- pi / 2
  data_aug$w_orig    <- data_aug[[w]]
  data_aug$neg_w     <- -data_aug[[w]]
  data_aug$w_pi2     <- data_aug[[w]] + pi / 2
  data_aug$neg_w_pi2 <- -data_aug[[w]] + pi / 2

  modelFormula <- stats::as.formula(
    paste(response, "~ const0 + constPi2 + w_orig + neg_w + w_pi2 + neg_w_pi2")
  )

  fit <- angular(formula = modelFormula, data = data_aug, ...)

  return(fit)
}
