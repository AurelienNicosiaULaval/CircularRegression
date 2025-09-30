###############################################################################
# Angular Regression Model (minimal changes + random intercept κ_a)
###############################################################################

#' Angular Regression Model
#'
#' Function that fits a regression model for angular variables along with its
#' associated S3 methods. Optionally, a random intercept a_i ~ VM(0, kappa_a)
#' may be estimated for clustered data (Rivest et al., 2019).
#'
#' The control argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{pginit}}{Grid size for initial beta search when \code{initbeta} is missing. Default 1000.}
#'   \item{\code{maxiter}}{Maximum number of iterations. Default 1000.}
#'   \item{\code{mindiff}}{Convergence tolerance on max-cosine improvement. Default 1e-6.}
#' }
#'
#' @param formula A formula with the dependent angle on the left of the ~ operator and terms specifying
#'                the explanatory variables on the right. These terms must be written \code{x:z}.
#' @param data An optional data frame, list or environment containing the variables in the model formula.
#' @param model A character string, either \code{"complete"} (with intercept) or \code{"simplified"} (no intercept).
#' @param initbeta Optional numeric vector of initial beta values.
#' @param control A list of control parameters. See Details.
#' @param random Optional one-sided formula of the form \code{~ 1 | id} to enable a VM random intercept per id.
#' @param cluster Optional factor giving the cluster for each row (ignored if \code{random} is supplied).
#' @param alternate Logical; if TRUE (default), alternates once between updating \eqn{\beta} on \eqn{y-\hat a}
#'        and re-estimating \eqn{(\kappa_e,\kappa_a)}.
#'
#' @return An object of class "angular" containing:
#' \describe{
#'   \item{MaxCosine}{the maximum cosine value.}
#'   \item{parameters}{matrix of estimates, std. errors, z-values and p-values for beta.}
#'   \item{varcov0}{model-based variance-covariance for beta.}
#'   \item{varcov1}{robust (sandwich-like) variance-covariance for beta.}
#'   \item{long}{vector of predicted concentrations.}
#'   \item{mui}{vector of predicted mean angles.}
#'   \item{y}{the response variable.}
#'   \item{iter.detail}{iteration details.}
#'   \item{call}{the function call.}
#'   \item{kappa_e}{(if random intercept) estimated residual concentration.}
#'   \item{kappa_a}{(if random intercept) estimated random-intercept concentration.}
#'   \item{a_i}{(if random intercept) empirical Bayes predictors per observation (length n).}
#'   \item{cluster}{(if random intercept) factor of cluster membership.}
#' }
#'
#' @export
angular <- function(
  formula,
  data,
  model = "simplified",
  initbeta = NULL,
  control = list(),
  random = NULL,
  cluster = NULL,
  alternate = TRUE
) {
  call <- match.call()
  model <- model[1]

  # helpers
  A1inv <- circular::A1inv
  A1 <- function(k) ifelse(k == 0, 0, besselI(k, 1) / besselI(k, 0))

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

  # clusters (minimal API change)
  if (!is.null(random)) {
    # expect ~ 1 | id
    rhs <- random[[3L]]
    if (length(rhs) == 3L && as.character(rhs[[1L]]) == "|") {
      idname <- as.character(rhs[[3L]])
      cluster <- as.factor(data[[idname]])
    } else stop("'random' must be of the form ~ 1 | id")
  } else if (!is.null(cluster)) {
    cluster <- as.factor(cluster)
    if (length(cluster) != nobs) stop("'cluster' must have same length as data")
  }
  has_re <- !is.null(cluster)

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
  maxcos <- function(
    beta,
    yi = y,
    x0i = if (model == "simplified") x0 else NULL
  ) {
    sinmu <- sin(if (model == "simplified") x0i else beta[p + 1]) +
      as.vector((matz * sin(matx)) %*% beta[1:p])
    cosmu <- cos(if (model == "simplified") x0i else beta[p + 1]) +
      as.vector((matz * cos(matx)) %*% beta[1:p])
    long <- sqrt(sinmu^2 + cosmu^2)
    mui <- atan2(sinmu, cosmu)
    maxcos_val <- mean(cos(yi - mui))
    list(maxcos = maxcos_val, long = long, mui = mui)
  }

  betaUpdate <- function(betak, long, mui, yi = y) {
    matd <- cbind(
      matz * sin(matx - mui),
      if (model == "simplified") NULL else cos(betak[p + 1] - mui)
    ) /
      long
    res <- sin(yi - mui)
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

  if (exists("maxcosk1") && maxcosk == maxcosk1) {
    maj <- betaUpdate(betak = betak, long = long, mui = mui)
  }
  matd <- maj$matd
  res <- maj$res
  Akappa <- mean(maxcosk)
  kappahat <- A1inv(Akappa)

  # Base SEs as before (for no-RE case and as fallback)
  XtX <- crossprod(matd)
  XtX_inv <- solve(XtX)
  v0 <- XtX_inv / (Akappa * kappahat)

  res_sq <- res^2
  matd_res <- matd * matrix(res_sq, nrow = nrow(matd), ncol = ncol(matd))
  v2 <- XtX_inv %*% crossprod(matd, matd_res) %*% XtX_inv

  # -------- Random intercept branch (minimal, optional) --------
  kappa_e <- kappa_a <- NA_real_
  a_i <- rep(0, nobs)
  if (has_re) {
    eps <- y - mui
    eps <- atan2(sin(eps), cos(eps))
    lev <- levels(cluster)
    n_i <- as.integer(table(cluster)[lev])

    sumC <- tapply(eps, cluster, function(e) sum(cos(e)))
    sumS <- tapply(eps, cluster, function(e) sum(sin(e)))

    loglik_vm_marg <- function(par) {
      ke <- pmax(par[1], 0)
      ka <- pmax(par[2], 0)
      sig <- sqrt((ka + ke * sumC)^2 + (ke * sumS)^2)
      sum(
        log(besselI(sig, 0)) -
          log(besselI(ka, 0)) -
          n_i * log(besselI(ke, 0)) -
          n_i * log(2 * pi)
      )
    }

    # crude inits via moments
    Rbar <- mean(cos(eps))
    ke0 <- max(kappahat, 1e-1)
    ka0 <- max(A1inv(pmin(pmax(Rbar / pmax(A1(ke0), 1e-8), 1e-6), 1 - 1e-6)), 0)

    fitK <- optim(
      c(ke0, ka0),
      fn = function(p) -loglik_vm_marg(p),
      method = "L-BFGS-B",
      lower = c(0, 0)
    )
    kappa_e <- max(fitK$par[1], 0)
    kappa_a <- max(fitK$par[2], 0)

    # EB predictors per cluster, then per obs
    a_hat_lev <- atan2(kappa_e * sumS, kappa_a + kappa_e * sumC)
    a_i <- as.numeric(a_hat_lev[cluster])

    if (isTRUE(alternate)) {
      # one alternation: refit beta on tilted response and re-opt κ
      y_t <- atan2(sin(y - a_i), cos(y - a_i))
      # re-run short maxcos updates
      calcul2 <- maxcos(beta = betak, yi = y_t)
      maxcosk2 <- calcul2$maxcos
      long2 <- calcul2$long
      mui2 <- calcul2$mui
      iter2 <- 0
      repeat {
        maj2 <- betaUpdate(betak = betak, long = long2, mui = mui2, yi = y_t)
        betak_n <- maj2$betak1
        calc_n <- maxcos(beta = betak_n, yi = y_t)
        if ((calc_n$maxcos - maxcosk2) <= mindiff || iter2 >= 50) {
          betak <- betak_n
          mui <- calc_n$mui
          long <- calc_n$long
          break
        } else {
          betak <- betak_n
          mui <- calc_n$mui
          long <- calc_n$long
          maxcosk2 <- calc_n$maxcos
          iter2 <- iter2 + 1
        }
      }
      # update eps and sums, then κ
      eps <- atan2(sin(y - mui), cos(y - mui))
      sumC <- tapply(eps, cluster, function(e) sum(cos(e)))
      sumS <- tapply(eps, cluster, function(e) sum(sin(e)))
      loglik2 <- function(par) {
        ke <- pmax(par[1], 0)
        ka <- pmax(par[2], 0)
        sig <- sqrt((ka + ke * sumC)^2 + (ke * sumS)^2)
        sum(
          log(besselI(sig, 0)) -
            log(besselI(ka, 0)) -
            n_i * log(besselI(ke, 0)) -
            n_i * log(2 * pi)
        )
      }
      fitK2 <- optim(
        c(kappa_e, kappa_a),
        fn = function(p) -loglik2(p),
        method = "L-BFGS-B",
        lower = c(0, 0)
      )
      kappa_e <- max(fitK2$par[1], 0)
      kappa_a <- max(fitK2$par[2], 0)
      a_hat_lev <- atan2(kappa_e * sumS, kappa_a + kappa_e * sumC)
      a_i <- as.numeric(a_hat_lev[cluster])
    }

    # recompute SEs for beta conditionally on κ_e using conditional residuals
    res_cond <- sin(y - mui - a_i)
    XtX <- crossprod(matd)
    XtX_inv <- solve(XtX)
    v0 <- XtX_inv / (pmax(A1(kappa_e), 1e-8) * pmax(kappa_e, 1e-8))
    matd_res <- matd * matrix(res_cond^2, nrow = nrow(matd), ncol = ncol(matd))
    v2 <- XtX_inv %*% crossprod(matd, matd_res) %*% XtX_inv
  }

  zvalue <- abs(betak) / sqrt(diag(v0))
  pvals <- round(
    2 * stats::pnorm(abs(betak) / sqrt(diag(v0)), lower.tail = FALSE),
    5
  )
  parameters <- cbind(betak, sqrt(diag(v0)), zvalue, pvals)
  colnames(parameters) <- c("estimate", "stderr", "z value", "P(|z|>.)")
  rownames(parameters) <- betaname

  out <- list(
    MaxCosine = maxcosk,
    parameters = parameters,
    kappahat = kappahat,
    varcov0 = v0,
    varcov1 = v2,
    long = long,
    mui = mui,
    y = y,
    iter.detail = iter.detail,
    call = call,
    kappa_e = if (has_re) kappa_e else NULL,
    kappa_a = if (has_re) kappa_a else NULL,
    a_i = if (has_re) a_i else NULL,
    cluster = if (has_re) cluster else NULL
  )
  class(out) <- "angular"
  out
}

###############################################################################
### S3 Methods for Angular Regression Objects (minimal edits)
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
  cat("\nConcentration parameter:", format(x$kappahat, digits = 4), "\n")
  if (!is.null(x$kappa_e) && !is.null(x$kappa_a)) {
    cat(
      "Random-intercept concentrations: kappa_e =",
      format(x$kappa_e, digits = 4),
      ", kappa_a =",
      format(x$kappa_a, digits = 4),
      "\n\n"
    )
  } else cat("\n")
  cat("Parameters:\n")
  coefmat <- x$parameters
  mat_numeric <- cbind(
    Estimate = as.numeric(coefmat[, 1]),
    `Robust std` = as.numeric(coefmat[, 2]),
    `z value` = as.numeric(coefmat[, 3]),
    `P(|z|>.)` = as.numeric(coefmat[, 4])
  )
  rownames(mat_numeric) <- rownames(coefmat)
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
  # work entirely in radians and return a circular object (radians)
  as_circ_rad <- function(x)
    circular::conversion.circular(
      circular::as.circular(x),
      units = "radians"
    )

  y <- as_circ_rad(object$y)
  mu <- as_circ_rad(object$mui)

  if (!is.null(object$a_i)) {
    ai <- as_circ_rad(object$a_i)
    delta <- y - (mu + ai)
  } else {
    delta <- y - mu
  }

  d <- atan2(sin(delta), cos(delta))
  circular::as.circular(d, units = "radians", modulo = "asis")
}

#' Summary Method for Angular Regression Objects
#'
#' Summarizes an object of class "angular" by providing key diagnostic statistics.
#'
#' @param object An object of class "angular".
#' @param ... Further arguments (currently ignored).
#'
#' @return An object of class "summary.angular".
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
    kappa_e = object$kappa_e,
    kappa_a = object$kappa_a,
    residuals = resid_summary,
    nobs = length(object$y)
  )
  class(s) <- "summary.angular"
  s
}

#' Print Method for Summary of Angular Regression Objects
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
  cat("\nConcentration parameter:", format(x$kappahat, digits = 4), "\n")
  if (!is.null(x$kappa_e) && !is.null(x$kappa_a)) {
    cat(
      "Random-intercept: kappa_e =",
      format(x$kappa_e, digits = 4),
      ", kappa_a =",
      format(x$kappa_a, digits = 4),
      "\n"
    )
  }
  cat("\nResiduals:\n")
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
#' This function produces diagnostic plots for an object of class "angular".
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

  res <- as.numeric(residuals.angular(x))

  p1 <- ggplot2::ggplot(
    data = data.frame(Fitted = x$mui, Residual = res),
    ggplot2::aes(x = Fitted, y = Residual)
  ) +
    ggplot2::geom_point(color = "blue") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed") +
    ggplot2::labs(
      title = "Residuals vs Fitted",
      x = "Fitted Values",
      y = "Residuals"
    )

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
