###############################################################################
# Consensus Model
###############################################################################

#' Consensus Model
#'
#' Function that fits a consensus model for angular variables along with its print,
#' summary, and diagnostic plot methods.
#'
#' @param formula A formula with the dependent angle on the left of the \code{~} operator and terms specifying
#'                the explanatory variables on the right. These terms must be written as \code{x:z}, where
#'                \code{x} is an explanatory angle whose relative importance may depend on the positive variable
#'                \code{z}. For \code{model = "simplified"}, the first explanatory angle listed is the reference
#'                direction (any provided \code{z} for this angle is ignored).
#' @param data An optional data frame, list or environment containing the variables in the model formula.
#'             If not found in \code{data}, the variables are taken from \code{environment(formula)}.
#' @param model A character string, either \code{"complete"} for the complete model with an intercept (default)
#'              or \code{"simplified"} for the simplified model without an intercept.
#' @param initbeta A numerical vector of initial values for the parameters. The default is to use the best
#'                 initial values found among a grid of possible values.
#' @param weights Optional vector of weights to be used in the estimation.
#' @param control A list of control parameters. The following components can be supplied:
#' \describe{
#'   \item{\code{pginit}}{The approximate number of points on the grid of possible initial beta values (default 1000).}
#'   \item{\code{maxiter}}{The maximum number of iterations (default 1000).}
#'   \item{\code{mindiff}}{The minimum difference between two max cosine values to be reached for convergence
#'                          (default 1e-06).}
#' }
#'
#' @return An object of class "consensus" containing:
#' \describe{
#'   \item{MaxLL}{the maximum value of the log likelihood.}
#'   \item{AIC}{the Akaike Information Criterion.}
#'   \item{BIC}{the Bayesian Information Criterion.}
#'   \item{parameters}{the parameter estimates and their standard errors, z-values and associated p-values.}
#'   \item{varcov1}{the estimated variance covariance matrix.}
#'   \item{parambeta}{the beta parameter estimates and their standard errors (obtained by linearization).}
#'   \item{varcovbeta1}{the estimated variance covariance matrix for the beta estimates (by linearization).}
#'   \item{autocorr}{the autocorrelation of the residuals \eqn{\sin(y_i - \mu_i)}.}
#'   \item{iter.detail}{the iteration details.}
#'   \item{call}{the function call.}
#' }
#'
#' @author Sophie Baillargeon, Louis-Paul Rivest, and Aurélien Nicosia
#' @examples
#' \dontrun{
#'   data(bison)
#'   fit <- consensus(y.dir ~ y.prec + y.prec2, data = bison)
#'   print(fit)
#'   summary(fit)
#' }
#' @export
consensus <- function(
  formula,
  data,
  model = "simplified",
  weights = NULL,
  initbeta = NULL,
  control = list()
) {
  call <- mfcall <- match.call()
  model <- model[1]

  ### Extract model.frame using stats::as.formula
  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- as.name("model.frame")
  mf <- eval(mfcall, parent.frame())
  nobs <- nrow(mf)
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  nterms <- length(nomterms)
  p <- if (model == "simplified") nterms - 1 else nterms
  nparam <- if (model == "simplified") p + 1 else p + 2
  paramname <- nomterms
  if (model == "complete") paramname <- c(paramname, paste0("beta", p + 1))

  # Response variable
  y <- as.vector(mf[, 1])
  # Split explanatory variables (assumes format x:z)
  noms <- strsplit(nomterms, split = ":")
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
  weight <- rep(1, nobs) * (is.null(weights)) + (!is.null(weights)) * weights

  ### Define log-likelihood function
  LL <- function(param) {
    angleref <- if (model == "simplified") x0 else rep(param[p + 2], nobs)
    sinmu <- param[1] * sin(angleref) + (matz * sin(matx)) %*% param[2:(p + 1)]
    cosmu <- param[1] * cos(angleref) + (matz * cos(matx)) %*% param[2:(p + 1)]
    long <- as.vector(sqrt(sinmu^2 + cosmu^2))
    mui <- as.vector(atan2(sinmu, cosmu))
    term1 <- param[1] *
      cos(y - angleref) +
      (matz * cos(y - matx)) %*% param[2:(p + 1)]
    LL_val <- sum(term1) - sum(log(besselI(long, 0, expon.scaled = FALSE)))
    list(LL = LL_val, long = long, mui = mui)
  }

  # Parameter update function
  betaUpdate <- function(paramk, long, mui) {
    angleref <- if (model == "simplified") x0 else rep(paramk[p + 2], nobs)
    matx0 <- cbind(angleref, matx)
    matz0 <- cbind(rep(1, nobs), matz)
    Along <- as.vector(
      besselI(long, 1, expon.scaled = FALSE) /
        besselI(long, 0, expon.scaled = FALSE)
    )
    matu <- matz0 * (cos(y - matx0) - cos(matx0 - mui) * Along)
    if (model == "complete")
      matu <- cbind(
        matu,
        paramk[1] * sin(y - angleref) - sin(mui - angleref) * Along
      )
    vecs <- colSums(matu)
    names(vecs) <- paramname
    Xc <- matz0 * cos(matx0 - mui)
    Xs <- matz0 * sin(matx0 - mui)
    if (model == "complete") {
      Xc <- cbind(Xc, paramk[1] * sin(mui - paramk[p + 2]))
      Xs <- cbind(Xs, paramk[1] * cos(mui - paramk[p + 2]))
    }
    Dc <- diag(1 - Along / long - Along^2, nrow = nobs, ncol = nobs)
    Ds <- diag(Along / long, nrow = nobs, ncol = nobs)
    matI <- t(Xc) %*% Dc %*% Xc + t(Xs) %*% Ds %*% Xs
    colnames(matI) <- rownames(matI) <- paramname
    dparam <- as.vector(solve(matI, vecs))
    paramk1 <- paramk + dparam
    list(paramk1 = paramk1, dparam = dparam, matu = matu, matI = matI)
  }

  if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 1000 else control$pginit
    pg <- round(pginit^(1 / nparam))
    possparam <- rep(
      list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]),
      p + 1
    )
    if (model == "complete")
      possparam[[nparam]] <- seq(0, 2 * pi, length.out = pg + 2)[-c(1, pg + 2)]
    possVal <- cbind(expand.grid(possparam), NA)
    colnames(possVal) <- c(paramname, "LL")
    maxLLfun <- function(param) {
      LL(param = param)$LL
    }
    possVal[, nparam + 1] <- apply(
      possVal[, 1:nparam, drop = FALSE],
      1,
      maxLLfun
    )
    paramk <- unlist(possVal[which.max(possVal[, nparam + 1]), 1:nparam])
  } else {
    if (length(initbeta) != nparam)
      stop("For the requested model, 'initbeta' must be of length ", nparam)
    paramk <- initbeta
  }

  calcul <- LL(param = paramk)
  maxLLk <- calcul$LL
  long <- calcul$long
  mui <- calcul$mui
  iter <- iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 1e-06 else control$mindiff
  conv <- FALSE
  iter.detail <- matrix(NA, nrow = maxiter + 1, ncol = nparam + 3)
  colnames(iter.detail) <- c(paramname, "LL", "iter", "nitersh")
  iter.detail[1, ] <- c(paramk, maxLLk, iter, iter.sh)

  while (!conv && iter <= maxiter) {
    maj <- betaUpdate(paramk = paramk, long = long, mui = mui)
    paramk1 <- maj$paramk1
    dparam <- maj$dparam
    calcul <- LL(param = paramk1)
    maxLLk1 <- calcul$LL
    long <- calcul$long
    mui <- calcul$mui
    iter.sh <- 0
    while (maxLLk1 < maxLLk) {
      iter.sh <- iter.sh + 1
      paramk1 <- paramk + dparam / (2^iter.sh)
      calcul <- LL(param = paramk1)
      maxLLk1 <- calcul$LL
      long <- calcul$long
      mui <- calcul$mui
      if (iter.sh >= maxiter) break
    }
    if (maxLLk1 < maxLLk) {
      conv <- FALSE
      warning(
        "The algorithm did not converge, it failed to maximize the log likelihood"
      )
      break
    } else {
      conv <- if (maxLLk1 - maxLLk > mindiff) FALSE else TRUE
      paramk <- paramk1
      maxLLk <- maxLLk1
      iter <- iter + 1
      iter.detail[iter + 1, ] <- c(paramk, maxLLk, iter, iter.sh)
    }
  }
  if (iter > maxiter + 1) {
    warning(
      "The algorithm did not converge, the maximum number of iterations was reached"
    )
  } else {
    iter.detail <- iter.detail[1:(iter + 1), , drop = FALSE]
  }

  if (maxLLk == maxLLk1) {
    maj <- betaUpdate(paramk = paramk, long = long, mui = mui)
  }
  matd <- maj$matu

  v1 <- solve(crossprod(matd))

  paramb <- paramk[2:(p + 1)] / paramk[1]
  matDeriv <- rbind(
    -paramk[2:(p + 1)] / paramk[1]^2,
    diag(1 / paramk[1], nrow = p, ncol = p)
  )
  vb <- t(matDeriv) %*% v1[1:(p + 1), 1:(p + 1)] %*% matDeriv
  names(paramb) <- colnames(vb) <- rownames(vb) <- paramname[-1]

  zvalue <- abs(paramk) / sqrt(diag(v1))
  pvals <- round(
    2 * stats::pnorm(abs(paramk) / sqrt(diag(v1)), lower.tail = FALSE),
    5
  )
  parameters <- cbind(paramk, sqrt(diag(v1)), zvalue, pvals)
  colnames(parameters) <- c("estimate", "Std. Error", "z value", "P(|z|>.)")
  rownames(parameters) <- paramname

  zvaluebeta <- abs(paramb) / sqrt(diag(vb))
  pbeta <- round(
    2 * stats::pnorm(abs(paramb) / sqrt(diag(vb)), lower.tail = FALSE),
    5
  )
  parambeta <- cbind(paramb, sqrt(diag(vb)), zvaluebeta, pbeta)
  colnames(parambeta) <- c("estimate", "Std. Error", "z value", "P(|z|>.)")
  rownames(parambeta) <- names(paramb)

  # Use residuals.consensus() defined in the S3 methods below
  autocorr <- stats::acf(
    residuals.consensus(object = list(y = y, mui = mui, long = long)),
    plot = FALSE
  )

  # Calculate AIC and BIC as in lm/glm
  k <- nparam
  logLik <- maxLLk
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + log(nobs) * k

  out <- list(
    MaxLL = maxLLk,
    AIC = AIC,
    BIC = BIC,
    parameters = parameters,
    varcov1 = v1,
    parambeta = parambeta,
    varcovbeta1 = vb,
    #   autocorr = autocorr,
    matx = matx,
    matz = matz,
    y = y,
    long = long,
    mui = mui,
    iter.detail = iter.detail,
    call = call,
    nobs = nobs,
    k = k,
    logLik = logLik
  )
  class(out) <- "consensus"
  out
}

###############################################################################
### S3 Methods for Consensus Objects
###############################################################################

#' Methods for Consensus Objects
#'
#' These functions provide methods for printing, extracting coefficients,
#' summarizing, and retrieving residuals from a consensus model.
#'
#' @param object An object of class "consensus".
#' @param x An object of class "consensus".
#' @param ... Additional arguments passed to the corresponding functions.
#'
#' @rdname consensus-methods
#' @export
print.consensus <- function(x, ...) {
  cat(
    "Call:
"
  )
  print(x$call)
  cat(
    "
Maximum log-likelihood:",
    format(x$MaxLL, digits = 4),
    "
"
  )
  cat(
    "AIC:",
    format(x$AIC, digits = 4),
    "
"
  )
  cat(
    "BIC:",
    format(x$BIC, digits = 4),
    "

"
  )
  cat(
    "Coefficients (Kappa estimates):
"
  )

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
    `Std. Error` = as.numeric(coefmat[, 2]),
    `z value` = as.numeric(coefmat[, 3]),
    `P(|z|>.)` = as.numeric(coefmat[, 4])
  )

  rownames(mat_numeric) <- rownames(coefmat)

  # Impression avec printCoefmat
  stats::printCoefmat(mat_numeric, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

#' @rdname consensus-methods
#' @export
coef.consensus <- function(object, ...) {
  object$parameters[, "estimate"]
}

#' @rdname consensus-methods
#' @export
residuals.consensus <- function(object, ...) {
  # tout en radians pour éviter les ambiguïtés
  y <- circular::conversion.circular(
    circular::as.circular(object$y),
    units = "radians"
  )
  mu <- circular::conversion.circular(
    circular::as.circular(object$mui),
    units = "radians"
  )

  # différence angulaire signée minimale dans (-pi, pi]
  d <- atan2(sin(y - mu), cos(y - mu))

  # renvoyer un 'circular' en radians (convertis ensuite si tu veux)
  circular::as.circular(d, units = "radians", modulo = "asis")
}

#' Summary Method for Consensus Objects
#'
#' This function summarizes an object of class "consensus" by providing key
#' diagnostic statistics including the model call, maximum log-likelihood, AIC, BIC, a summary of residuals,
#' the coefficient estimates, and the number of observations.
#'
#' @param object An object of class "consensus".
#' @param ... Further arguments (currently ignored).
#'
#' @return An object of class "summary.consensus", a list containing:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{MaxLL}{The maximum log-likelihood value.}
#'   \item{AIC}{The Akaike information criterion.}
#'   \item{BIC}{The Bayesian information criterion.}
#'   \item{coefficients}{A matrix with the estimates, standard errors, z-values and associated p-values.}
#'   \item{residuals}{A summary of the residuals (min, 1st Qu., median, 3rd Qu., max).}
#'   \item{nobs}{The number of observations.}
#' }
#'
#' @rdname consensus-methods
#' @export
summary.consensus <- function(object, ...) {
  res <- residuals.consensus(object, ...)
  resid_summary <- summary(res)
  s <- list(
    call = object$call,
    MaxLL = object$MaxLL,
    AIC = object$AIC,
    BIC = object$BIC,
    coefficients = object$parameters,
    residuals = resid_summary,
    nobs = length(object$y)
  )
  class(s) <- "summary.consensus"
  s
}

#' Print Method for Summary of Consensus Objects
#'
#' This function prints a summary of a consensus model.
#'
#' @param x An object of class "summary.consensus".
#' @param ... Further arguments passed to printing functions.
#'
#' @rdname consensus-methods
#' @export
print.summary.consensus <- function(x, ...) {
  cat(
    "Call:
"
  )
  print(x$call)
  cat(
    "
Maximum log-likelihood:",
    format(x$MaxLL, digits = 4),
    "
"
  )
  cat(
    "AIC:",
    format(x$AIC, digits = 4),
    "
"
  )
  cat(
    "BIC:",
    format(x$BIC, digits = 4),
    "

"
  )
  cat(
    "Residuals:
"
  )
  print(x$residuals)
  cat(
    "
Coefficients:
"
  )
  stats::printCoefmat(
    x$coefficients,
    signif.stars = TRUE,
    signif.legend = TRUE,
    ...
  )
  cat(
    "
Number of observations:",
    x$nobs,
    "
"
  )
  invisible(x)
}

#' Diagnostic Plots for Consensus Model
#'
#' This function produces a set of diagnostic plots for an object of class "consensus"
#' in a layout similar to that of \code{plot.lm}. The diagnostics include:
#' \enumerate{
#'   \item Residuals vs Fitted values.
#'   \item Histogram of scaled residuals with a normal density overlay.
#'   \item Normal Q-Q plot.
#'   \item A text display of a goodness-of-fit measure: the Spearman correlation between
#'         the absolute raw residuals and the inverse of the length parameter.
#' }
#'
#' @param x An object of class "consensus".
#' @param ... Further arguments passed to plotting functions.
#'
#' @rdname consensus-methods
#' @export
plot.consensus <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("Package 'gridExtra' is required for arranging plots.")

  p1 <- ggplot2::ggplot(
    data = data.frame(Fitted = x$mui, Residual = residuals.consensus(x, ...)),
    ggplot2::aes(x = Fitted, y = Residual)
  ) +
    ggplot2::geom_point(color = "blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed") +
    ggplot2::labs(
      title = "Residuals vs Fitted",
      x = "Fitted Values",
      y = "Residuals"
    )

  res <- as.numeric(residuals.consensus(x, ...))
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
      x = "Scaled Residuals",
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

  p4 <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste(
        "Spearman Correlation:\n|y - mui| vs 1/long =",
        round(stats::cor(abs(res), 1 / x$long, method = "spearman"), 4)
      ),
      size = 5,
      hjust = 0.5
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Goodness-of-Fit")

  gridExtra::grid.arrange(
    p1,
    p2,
    p3,
    p4,
    ncol = 2,
    nrow = 2,
    top = "Diagnostic Plots for Consensus Model"
  )

  invisible(x)
}

# Méthodes S3 AIC, BIC, logLik
#' @rdname consensus-methods
#' @param object An object of class "consensus".
#' @param k Numeric. Penalty term used in the definition of AIC. Defaults to 2.
#' @param ... Additional arguments passed to or from other methods.
#' @export
AIC.consensus <- function(object, ..., k = 2) {
  ll <- object$logLik
  nparam <- object$k
  out <- -2 * ll + k * nparam
  attr(out, "df") <- nparam
  attr(out, "nobs") <- object$nobs
  out
}

#' @rdname consensus-methods
#' @param object An object of class "consensus".
#' @param ... Additional arguments passed to or from other methods.
#' @export
BIC.consensus <- function(object, ...) {
  n <- object$nobs
  nparam <- object$k
  ll <- object$logLik
  out <- -2 * ll + log(n) * nparam
  attr(out, "df") <- nparam
  attr(out, "nobs") <- n
  out
}

#' @rdname consensus-methods
#' @param object An object of class "consensus".
#' @param ... Additional arguments passed to or from other methods.
#' @export
logLik.consensus <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$k
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}
