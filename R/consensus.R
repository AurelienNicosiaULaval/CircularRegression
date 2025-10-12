###############################################################################
# Consensus Model
###############################################################################

#' Consensus Model
#'
#' Fit the simplified consensus regression model for angular data as described
#' by Rivest et al. (2016). The first term in the formula supplies the reference
#' direction with coefficient fixed to one. Subsequent terms must be written
#' \code{x:z} where \code{x} is an angular covariate and \code{z} is an optional
#' non-negative scaling variable (a value of one is assumed when \code{z} is
#' omitted or identical to \code{x}).
#'
#' @param formula A formula with the dependent angle on the left of the
#'   \code{~} operator and terms specifying the explanatory variables on the
#'   right. The first term provides the reference direction.
#' @param data An optional data frame, list or environment containing the
#'   variables in the model formula.
#' @param weights Optional non-negative observation weights. If supplied, they
#'   must have the same length as the response.
#' @param initbeta A numeric vector of initial parameter values. When
#'   \code{NULL} (default) a quasi-uniform grid search is used to obtain starting
#'   values. The required length equals the number of \code{x:z} terms plus one.
#' @param control A list of control parameters. The following components can be
#'   supplied:
#'   \describe{
#'     \item{\code{pginit}}{Approximate number of grid points used to obtain
#'       starting values when \code{initbeta} is not provided (default 1000).}
#'     \item{\code{maxiter}}{Maximum number of Gauss--Newton iterations
#'       (default 1000).}
#'     \item{\code{mindiff}}{Convergence tolerance on the increase of the log
#'       likelihood (default 1e-6).}
#'   }
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
  weights = NULL,
  initbeta = NULL,
  control = list()
) {
  call <- mfcall <- match.call()

  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- quote(model.frame)
  mf <- eval(mfcall, parent.frame())
  nobs <- nrow(mf)
  if (nobs == 0) stop("No observations available.")

  y <- as.numeric(mf[[1]])
  term_labels <- attr(attr(mf, "terms"), "term.labels")
  if (length(term_labels) == 0) {
    stop(
      "The model must include at least one term specifying the reference direction."
    )
  }

  split_terms <- strsplit(term_labels, ":", fixed = TRUE)
  term_matrix <- t(vapply(
    split_terms,
    function(x) {
      if (length(x) == 1) {
        c(x, x)
      } else if (length(x) == 2) {
        x
      } else {
        stop("Each term must be of the form 'x' or 'x:z'.", call. = FALSE)
      }
    },
    character(2)
  ))

  ref_name <- term_matrix[1, 1]
  if (!ref_name %in% names(mf)) {
    stop(
      sprintf(
        "Reference direction '%s' not found in the supplied data.",
        ref_name
      ),
      call. = FALSE
    )
  }
  x0 <- mf[[ref_name]]

  paramname <- term_labels
  betaname <- if (length(term_labels) > 1) term_labels[-1] else character(0)
  p <- length(betaname)
  if (p > 0) {
    x_names <- term_matrix[-1, 1]
    z_names <- term_matrix[-1, 2]
    matx <- as.matrix(mf[, x_names, drop = FALSE])
    matz <- matrix(1, nrow = nobs, ncol = p)
    for (j in seq_len(p)) {
      if (!identical(z_names[j], x_names[j])) {
        if (!z_names[j] %in% names(mf)) {
          stop(
            sprintf(
              "Modifier '%s' not found in the supplied data.",
              z_names[j]
            ),
            call. = FALSE
          )
        }
        matz[, j] <- mf[[z_names[j]]]
      }
    }
  } else {
    matx <- matrix(0, nrow = nobs, ncol = 0)
    matz <- matrix(0, nrow = nobs, ncol = 0)
  }

  weight <- if (is.null(weights)) rep(1, nobs) else as.numeric(weights)
  if (length(weight) != nobs) {
    stop("'weights' must have the same length as the response.", call. = FALSE)
  }
  if (any(!is.finite(weight)) || any(weight < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (all(weight == 0)) {
    stop("'weights' cannot be all zero.", call. = FALSE)
  }

  pginit <- if (is.null(control$pginit)) 1000 else control$pginit
  maxiter <- if (is.null(control$maxiter)) 1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 1e-06 else control$mindiff

  compute_components <- function(param) {
    kappa0 <- param[1]
    sinmu <- kappa0 * sin(x0)
    cosmu <- kappa0 * cos(x0)
    if (p > 0) {
      beta <- param[2:(p + 1)]
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta)
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta)
    }
    long <- sqrt(sinmu^2 + cosmu^2)
    mui <- atan2(sinmu, cosmu)
    list(long = long, mui = mui)
  }

  loglik_components <- function(param) {
    comp <- compute_components(param)
    long <- comp$long
    mui <- comp$mui
    term1 <- param[1] * cos(y - x0)
    if (p > 0) {
      term1 <- term1 + as.vector((matz * cos(y - matx)) %*% param[2:(p + 1)])
    }
    LL <- sum(weight * term1) -
      sum(weight * log(besselI(long, 0, expon.scaled = FALSE)))
    list(LL = LL, long = long, mui = mui)
  }

  betaUpdate <- function(paramk, long, mui) {
    matx0 <- if (p > 0) cbind(x0, matx) else matrix(x0, ncol = 1)
    matz0 <- if (p > 0) cbind(rep(1, nobs), matz) else
      matrix(1, ncol = 1, nrow = nobs)
    ratio_num <- besselI(long, 1, expon.scaled = FALSE)
    ratio_den <- besselI(long, 0, expon.scaled = FALSE)
    Along <- ratio_num / ratio_den
    Along[!is.finite(Along)] <- 0
    Along_over_long <- ifelse(long > 0, Along / long, 0)
    matu <- matz0 * (cos(y - matx0) - cos(matx0 - mui) * Along)
    vecs <- colSums(weight * matu)
    Xc <- matz0 * cos(matx0 - mui)
    Xs <- matz0 * sin(matx0 - mui)
    Dc_vec <- 1 - Along_over_long - Along^2
    Ds_vec <- Along_over_long
    sqrt_Dc <- sqrt(pmax(weight * Dc_vec, 0))
    sqrt_Ds <- sqrt(pmax(weight * Ds_vec, 0))
    matI <- crossprod(Xc * sqrt_Dc, Xc) + crossprod(Xs * sqrt_Ds, Xs)
    matI <- (matI + t(matI)) / 2
    qrI <- qr(matI)
    if (qrI$rank < ncol(matI)) {
      stop(
        "Information matrix is singular; cannot update parameters.",
        call. = FALSE
      )
    }
    dparam <- as.vector(qr.coef(qrI, vecs))
    list(paramk1 = paramk + dparam, dparam = dparam, matu = matu, matI = matI)
  }

  nparam <- p + 1
  if (is.null(initbeta)) {
    pg <- max(1L, round(pginit^(1 / nparam)))
    grid_vals <- rep(
      list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]),
      nparam
    )
    possVal <- cbind(expand.grid(grid_vals), NA_real_)
    colnames(possVal) <- c(paramname, "LL")
    possVal[, nparam + 1] <- apply(
      possVal[, seq_len(nparam), drop = FALSE],
      1,
      function(par) {
        loglik_components(as.numeric(par))$LL
      }
    )
    paramk <- as.numeric(possVal[
      which.max(possVal[, nparam + 1]),
      seq_len(nparam),
      drop = TRUE
    ])
  } else {
    if (length(initbeta) != nparam) {
      stop(
        sprintf("'initbeta' must have length %d for this model.", nparam),
        call. = FALSE
      )
    }
    paramk <- initbeta
  }

  calcul <- loglik_components(paramk)
  maxLLk <- calcul$LL
  long <- calcul$long
  mui <- calcul$mui
  iter <- iter.sh <- 0
  conv <- FALSE
  iter.detail <- matrix(NA_real_, nrow = maxiter + 1, ncol = nparam + 3)
  colnames(iter.detail) <- c(paramname, "LL", "iter", "nitersh")
  iter.detail[1, ] <- c(paramk, maxLLk, iter, iter.sh)
  maxLLk1 <- maxLLk

  while (!conv && iter <= maxiter) {
    maj <- betaUpdate(paramk = paramk, long = long, mui = mui)
    paramk1 <- maj$paramk1
    dparam <- maj$dparam
    calcul <- loglik_components(paramk1)
    maxLLk1 <- calcul$LL
    long <- calcul$long
    mui <- calcul$mui
    iter.sh <- 0
    while (maxLLk1 < maxLLk) {
      iter.sh <- iter.sh + 1
      paramk1 <- paramk + dparam / (2^iter.sh)
      calcul <- loglik_components(paramk1)
      maxLLk1 <- calcul$LL
      long <- calcul$long
      mui <- calcul$mui
      if (iter.sh >= maxiter) break
    }
    if (maxLLk1 < maxLLk) {
      warning(
        "The algorithm did not converge; it failed to maximise the log-likelihood."
      )
      conv <- FALSE
      break
    } else {
      conv <- (maxLLk1 - maxLLk) <= mindiff
      paramk <- paramk1
      maxLLk <- maxLLk1
      iter <- iter + 1
      iter.detail[iter + 1, ] <- c(paramk, maxLLk, iter, iter.sh)
    }
  }
  if (iter > maxiter + 1) {
    warning(
      "The algorithm did not converge: maximum number of iterations reached."
    )
  } else {
    iter.detail <- iter.detail[seq_len(iter + 1), , drop = FALSE]
  }

  maj <- betaUpdate(paramk = paramk, long = long, mui = mui)
  matI <- maj$matI

  invert_information <- function(M) {
    M_sym <- (M + t(M)) / 2
    chol_res <- tryCatch(chol(M_sym), error = function(e) NULL)
    if (is.null(chol_res)) {
      warning(
        "Information matrix is not positive definite; returning NA variances."
      )
      matrix(NA_real_, nrow = nrow(M_sym), ncol = ncol(M_sym))
    } else {
      chol2inv(chol_res)
    }
  }

  v1 <- invert_information(matI)

  if (p > 0) {
    paramb <- paramk[2:(p + 1)] / paramk[1]
    matDeriv <- rbind(
      -paramk[2:(p + 1)] / paramk[1]^2,
      diag(1 / paramk[1], nrow = p, ncol = p)
    )
    vb <- t(matDeriv) %*% v1[1:(p + 1), 1:(p + 1), drop = FALSE] %*% matDeriv
    names(paramb) <- colnames(vb) <- rownames(vb) <- betaname
  } else {
    paramb <- numeric(0)
    vb <- matrix(numeric(0), nrow = 0, ncol = 0)
  }

  se_param <- sqrt(diag(v1))
  zvalue <- paramk / se_param
  pvals <- 2 * stats::pnorm(abs(zvalue), lower.tail = FALSE)
  parameters <- cbind(
    estimate = paramk,
    `Std. Error` = se_param,
    `z value` = zvalue,
    `P(|z|>.)` = pvals
  )
  rownames(parameters) <- paramname

  if (p > 0) {
    se_beta <- sqrt(diag(vb))
    zbeta <- paramb / se_beta
    pbeta <- 2 * stats::pnorm(abs(zbeta), lower.tail = FALSE)
    parambeta <- cbind(
      estimate = paramb,
      `Std. Error` = se_beta,
      `z value` = zbeta,
      `P(|z|>.)` = pbeta
    )
  } else {
    parambeta <- matrix(
      numeric(0),
      nrow = 0,
      ncol = 4,
      dimnames = list(NULL, c("estimate", "Std. Error", "z value", "P(|z|>.)"))
    )
  }

  residual_obj <- list(y = y, mui = mui, long = long)
  autocorr <- stats::acf(
    residuals.consensus(object = residual_obj),
    plot = FALSE
  )

  k <- nparam
  logLik <- maxLLk
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + log(sum(weight)) * k

  out <- list(
    MaxLL = maxLLk,
    AIC = AIC,
    BIC = BIC,
    parameters = parameters,
    varcov1 = v1,
    parambeta = parambeta,
    varcovbeta1 = vb,
    autocorr = autocorr,
    matx = matx,
    matz = matz,
    y = y,
    long = long,
    mui = mui,
    weights = weight,
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
    ) +
    ggplot2::theme_minimal()

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
    ) +
    ggplot2::theme_minimal()

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
    ) +
    ggplot2::theme_minimal()

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
