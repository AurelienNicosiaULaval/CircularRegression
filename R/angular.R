###############################################################################
# Angular Regression Model
###############################################################################

#' Angular Regression Model
#'
#' Function that fits a regression model for angular variables along with its
#' associated S3 methods.
#'
#' The control argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{pginit}}{The approximate number of points on the grid of possible initial beta values tried when
#'                         \code{initbeta} is not given. The default is 1000, which runs quickly. A large value of
#'                         \code{pginit} makes the function slower.}
#'   \item{\code{maxiter}}{The maximum number of iterations. The default is 1000.}
#'   \item{\code{mindiff}}{The minimum difference between two max cosine values to be reached. It defines the convergence criteria:
#'                          if the difference between the max cosine for the updated parameter values and the previous
#'                          one is below \code{mindiff}, convergence is reached. The default is 0.000001.}
#' }
#'
#' @param formula A formula with the dependent angle on the left of the ~ operator and terms specifying
#'                the explanatory variables on the right. These terms must be written \code{x:z}, where
#'                \code{x} is an explanatory angle whose relative importance might depend on the
#'                positive variable \code{z}. For \code{model = "simplified"}, the first explanatory angle is
#'                the reference direction (any provided \code{z} for this angle is ignored).
#' @param data An optional data frame, list or environment containing the variables in the model formula.
#' @param model A character string, either \code{"complete"} for the complete model with an intercept (the default)
#'              or \code{"simplified"} for the simplified model without an intercept.
#' @param initbeta A numerical vector of initial beta parameter values. The default is to use the best initial
#'                 values found among a grid of possible values.
#' @param control A list of control parameters. See Details.
#'
#' @return An object of class "angular" containing:
#' \describe{
#'   \item{MaxCosine}{the maximum cosine value.}
#'   \item{parameters}{the parameter estimates and their standard errors, z-values and associated p-values.}
#'   \item{varcov0}{the estimated variance-covariance matrix (first definition).}
#'   \item{varcov1}{the estimated variance-covariance matrix (second definition).}
#'   \item{autocorr}{the autocorrelation of the residuals \eqn{\sin(y_i - \mu_i)}.}
#'   \item{long}{the vector of predicted concentrations.}
#'   \item{mui}{the vector of predicted mean angles.}
#'   \item{y}{the response variable.}
#'   \item{iter.detail}{the iteration details.}
#'   \item{call}{the function call.}
#' }
#'
#' @author Sophie Baillargeon, Louis-Paul Rivest, and Aurélien Nicosia
#' @references L.-P. Rivest, T. Duchesne, A. Nicosia & D. Fortin. A general angular regression model for the analysis of data on animal movement in ecology. Journal of the Royal Statistical Society, series C, to appear.
#' @examples
#' \dontrun{
#'   # Example with fictitious data
#'   data(wind)
#'   n <- length(wind)
#'   dat.hom <- data.frame(wind.t = wind[-c(1,2)],
#'                         wind.t_1 = wind[-c(1, n)],
#'                         wind.t_2 = wind[-c(n-1, n)])
#'   an <- angular(wind.t ~ wind.t_1 + wind.t_2, data = dat.hom)
#'   print(an)
#'   summary(an)
#'   plot(an)
#' }
#' @export
angular <- function(
  formula,
  data,
  model = "simplified",
  initbeta = NULL,
  control = list()
) {
  call <- match.call()
  model <- model[1]

  ### Extraction of the model.frame
  mfargs <- match(c("formula", "data"), names(call), 0L)
  call_subset <- call[c(1L, mfargs)]
  call_subset[[1L]] <- as.name("model.frame")
  mf <- eval(call_subset, parent.frame())
  nobs <- nrow(mf)
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  paramname <- nomterms
  nterms <- length(nomterms)
  p <- if (model == "simplified") nterms - 1 else nterms
  nparam <- if (model == "simplified") p else p + 1
  betaname <- if (length(paramname) > 1) paramname[-1] else character(0)

  # Response variable
  y <- mf[, 1]
  # Explanatory variables splitting on ":"
  noms <- strsplit(paramname, split = ":")
  noms <- do.call(rbind, noms)
  if (model == "simplified") {
    x0 <- mf[, noms[1, 1]]
    noms <- noms[-1, , drop = FALSE]
  }
  matx <- as.matrix(mf[, noms[, 1], drop = FALSE])
  if (ncol(noms) == 1) {
    matz <- matrix(1, ncol = ncol(matx), nrow = nrow(matx))
  } else {
    matz <- as.matrix(mf[, noms[, 2], drop = FALSE])
    matz[, noms[, 2] == noms[, 1]] <- 1
  }

  ### Model fitting: define internal function maxcos
  maxcos <- function(beta) {
    sinmu <- sin(if (model == "simplified") x0 else beta[p + 1]) +
      as.vector((matz * sin(matx)) %*% beta[1:p])
    cosmu <- cos(if (model == "simplified") x0 else beta[p + 1]) +
      as.vector((matz * cos(matx)) %*% beta[1:p])
    long <- sqrt(sinmu^2 + cosmu^2)
    mui <- atan2(sinmu, cosmu)
    maxcos_val <- mean(cos(y - mui))
    list(maxcos = maxcos_val, long = long, mui = mui)
  }

  betaUpdate <- function(betak, long, mui) {
    matd <- cbind(
      matz * sin(matx - mui),
      if (model == "simplified") NULL else cos(betak[p + 1] - mui)
    ) /
      long
    res <- sin(y - mui)
    dbeta <- as.vector(solve(t(matd) %*% matd, t(matd) %*% res))
    betak1 <- betak + dbeta
    list(betak1 = betak1, dbeta = dbeta, matd = matd, res = res)
  }

  if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 1000 else control$pginit
    pg <- round(pginit^(1 / nparam))
    possbeta <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]), p)
    if (model == "complete")
      possbeta[[p + 1]] <- seq(0, 2 * pi, length.out = pg + 2)[-c(1, pg + 2)]
    possVal <- cbind(expand.grid(possbeta), NA)
    colnames(possVal) <- c(betaname, "maxcosine")
    maxcos1 <- function(beta) maxcos(beta)$maxcos
    possVal[, nparam + 1] <- apply(
      possVal[, 1:nparam, drop = FALSE],
      1,
      maxcos1
    )
    betak <- unlist(possVal[which.max(possVal[, nparam + 1]), 1:nparam])
  } else {
    if (length(initbeta) != nparam)
      stop("For the requested model, 'initbeta' must be of length ", nparam)
    betak <- initbeta
  }

  calcul <- maxcos(beta = betak)
  maxcosk <- calcul$maxcos
  long <- calcul$long
  mui <- calcul$mui
  iter <- iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 0.000001 else control$mindiff
  conv <- FALSE
  iter.detail <- matrix(NA, nrow = maxiter + 1, ncol = nparam + 3)
  colnames(iter.detail) <- c(betaname, "maxcosine", "iter", "nitersh")
  iter.detail[1, ] <- c(betak, maxcosk, iter, iter.sh)

  while (!conv && iter <= maxiter) {
    maj <- betaUpdate(betak = betak, long = long, mui = mui)
    betak1 <- maj$betak1
    dbeta <- maj$dbeta
    calcul <- maxcos(beta = betak1)
    maxcosk1 <- calcul$maxcos
    long <- calcul$long
    mui <- calcul$mui
    iter.sh <- 0
    while (maxcosk1 < maxcosk) {
      iter.sh <- iter.sh + 1
      betak1 <- betak + dbeta / (2^iter.sh)
      calcul <- maxcos(beta = betak1)
      maxcosk1 <- calcul$maxcos
      long <- calcul$long
      mui <- calcul$mui
      if (iter.sh >= maxiter) break
    }
    if (maxcosk1 < maxcosk) {
      conv <- FALSE
      warning(
        "The algorithm did not converge, it failed to maximize the max-cosine"
      )
      break
    } else {
      conv <- if (maxcosk1 - maxcosk > mindiff) FALSE else TRUE
      betak <- betak1
      maxcosk <- maxcosk1
      iter <- iter + 1
      iter.detail[iter + 1, ] <- c(betak, maxcosk, iter, iter.sh)
    }
  }
  if (iter > maxiter + 1) {
    warning(
      "The algorithm did not converge, the maximum number of iterations was reached"
    )
  } else {
    iter.detail <- iter.detail[1:(iter + 1), , drop = FALSE]
  }

  if (maxcosk == maxcosk1) {
    maj <- betaUpdate(betak = betak, long = long, mui = mui)
  }
  matd <- maj$matd
  res <- maj$res
  Akappa <- mean(maxcosk)
  kappahat <- circular::A1inv(Akappa)

  # Compute variance-covariance matrix v0
  v0 <- solve(t(matd) %*% matd) / (Akappa * kappahat)

  zvalue <- abs(betak) / sqrt(diag(v0))
  pvals <- round(
    2 * stats::pnorm(abs(betak) / sqrt(diag(v0)), lower.tail = FALSE),
    5
  )
  parameters <- cbind(betak, sqrt(diag(v0)), zvalue, pvals)
  colnames(parameters) <- c("estimate", "stderr", "z value", "P(|z|>.)")
  rownames(parameters) <- betaname

  autocorr <- stats::acf(res, plot = FALSE)

  out <- list(
    MaxCosine = maxcosk,
    parameters = parameters,
    kappahat = kappahat,
    varcov0 = v0,
    varcov1 = NULL,
    autocorr = autocorr,
    long = long,
    mui = mui,
    y = y,
    iter.detail = iter.detail,
    call = call
  )
  class(out) <- "angular"
  out
}

###############################################################################
### S3 Methods for Angular Regression Objects
###############################################################################

#' Methods for Angular Regression Objects
#'
#' These functions provide methods for printing, extracting coefficients,
#' summarizing, and retrieving residuals from an angular regression model.
#'
#' @param object An object of class "angular".
#' @param x An object of class "angular".
#' @param ... Additional arguments passed to the corresponding functions.
#'
#' @rdname angular-methods
#' @export
print.angular <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMaximum cosine:", format(x$MaxCosine, digits = 4), "")
  cat("\nConcentration parameter:", format(x$kappahat, digits = 4), "\n\n")
  cat("Parameters:\n")
  coefmat <- x$parameters
  stars <- ifelse(
    coefmat[, "P(|z|>.)"] < 0.001,
    "***",
    ifelse(
      coefmat[, "P(|z|>.)"] < 0.01,
      "**",
      ifelse(
        coefmat[, "P(|z|>.)"] < 0.05,
        "*",
        ifelse(coefmat[, "P(|z|>.)"] < 0.1, ".", " ")
      )
    )
  )

  # Transformation en matrice numérique pour printCoefmat
  mat_numeric <- cbind(
    Estimate = as.numeric(coefmat[, 1]),
    `Robust std` = as.numeric(coefmat[, 2]),
    `z value` = as.numeric(coefmat[, 3]),
    `P(|z|>.)` = as.numeric(coefmat[, 4])
  )

  rownames(mat_numeric) <- rownames(coefmat)

  # Impression avec printCoefmat
  stats::printCoefmat(mat_numeric, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

#' @rdname angular-methods
#' @export
coef.angular <- function(object, ...) {
  coef_vals <- object$parameters[, "estimate"]
  return(coef_vals)
}

#' @rdname angular-methods
#' @export
residuals.angular <- function(object, ...) {
  angle_diff_signed <- function(a, b) {
    diff <- (a - b + pi) %% (2 * pi) - pi
    diff
  }

  angle1 <- circular::circular(object$y)
  angle2 <- circular::circular(object$mui)
  res <- angle_diff_signed(angle1, angle2)
  res
}

#' Summary Method for Angular Regression Objects
#'
#' This function summarizes an object of class "angular" by providing key
#' diagnostic statistics including the model call, maximum cosine, a summary of residuals,
#' the parameter estimates, and the number of observations.
#'
#' @param object An object of class "angular".
#' @param ... Further arguments (currently ignored).
#'
#' @return An object of class "summary.angular", a list containing:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{MaxCosine}{The maximum cosine value.}
#'   \item{parameters}{A matrix with the estimates, standard errors, z-values and associated p-values.}
#'   \item{residuals}{A summary of the residuals (min, 1st Qu., median, 3rd Qu., max).}
#'   \item{nobs}{The number of observations.}
#' }
#'
#' @rdname angular-methods
#' @export
summary.angular <- function(object, ...) {
  res <- residuals.angular(object, ...)
  resid_summary <- summary(res)
  s <- list(
    call = object$call,
    MaxCosine = object$MaxCosine,
    parameters = object$parameters,
    kappahat = object$kappahat,
    residuals = resid_summary,
    nobs = length(object$y)
  )
  class(s) <- "summary.angular"
  s
}

#' Print Method for Summary of Angular Regression Objects
#'
#' This function prints a summary of an angular regression model.
#'
#' @param x An object of class "summary.angular".
#' @param ... Further arguments passed to printing functions.
#'
#' @rdname angular-methods
#' @export
print.summary.angular <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMaximum cosine:", format(x$MaxCosine, digits = 4), "")
  cat("\nConcentration parameter:", format(x$kappahat, digits = 4), "\n\n")
  cat("Residuals:\n")
  print(x$residuals)
  cat("\nParameters:\n")
  stats::printCoefmat(
    x$parameters,
    signif.stars = TRUE,
    signif.legend = TRUE,
    ...
  )
  cat("\nNumber of observations:", x$nobs, "\n")
  invisible(x)
}

#' Diagnostic Plots for Angular Regression Model
#'
#' This function produces a set of diagnostic plots for an object of class "angular"
#' in a layout similar to that of \code{plot.lm}. The diagnostics include:
#' \enumerate{
#'   \item Residuals vs Fitted values.
#'   \item Histogram of residuals with a normal density overlay.
#'   \item Normal Q-Q plot.
#'   \item A text display of a goodness-of-fit measure.
#' }
#'
#' @param x An object of class "angular".
#' @param ... Further arguments passed to plotting functions.
#'
#' @rdname angular-methods
#' @export
plot.angular <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("Package 'gridExtra' is required for arranging plots.")

  # Use ggplot2:: functions directly instead of library()
  p1 <- ggplot2::ggplot(
    data = data.frame(Fitted = x$mui, Residual = residuals.angular(x)),
    ggplot2::aes(x = Fitted, y = Residual)
  ) +
    ggplot2::geom_point(color = "blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed") +
    ggplot2::labs(
      title = "Residuals vs Fitted",
      x = "Fitted Values",
      y = "Residuals"
    )

  res <- as.numeric(residuals.angular(x))
  res_mean <- mean(res)
  res_sd <- stats::sd(res)
  p2 <- ggplot2::ggplot(
    data = data.frame(Residual = res),
    ggplot2::aes(x = Residual)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = after_stat(density)),
      bins = 12,
      color = "black",
      fill = "gray"
    ) +
    ggplot2::stat_function(
      fun = stats::dnorm,
      args = list(mean = res_mean, sd = res_sd),
      color = "red",
      size = 1
    ) +
    ggplot2::labs(
      title = "Histogram of Residuals",
      x = "Residuals",
      y = "Density"
    )

  p3 <- ggplot2::ggplot(
    data = data.frame(Residual = as.numeric(res)),
    ggplot2::aes(sample = Residual)
  ) +
    ggplot2::stat_qq(color = "blue") +
    ggplot2::stat_qq_line(color = "red") +
    ggplot2::labs(
      title = "Normal Q-Q",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    )

  gridExtra::grid.arrange(
    p1,
    p2,
    p3,
    ncol = 2,
    nrow = 2,
    top = "Diagnostic Plots for Angular Regression Model"
  )

  invisible(x)
}
