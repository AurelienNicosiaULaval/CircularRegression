###############################################################################
# Angular Regression Model
###############################################################################

#' Angular Regression Model
#'
#' Fit the simplified angular regression model of Rivest et al. (2016). The
#' first term in the formula provides the reference direction whose coefficient
#' is fixed to one. Subsequent terms must be written \code{x:z}, where
#' \code{x} is an angular covariate and \code{z} is a non-negative scaling
#' variable. When \code{z} is omitted (or equal to \code{x}), a value of one is
#' assumed. A reference-only model is obtained by including only the reference
#' term on the right-hand side.
#'
#' The control argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{pginit}}{Approximate number of points on the grid of possible
#'     initial beta values when \code{initbeta} is not supplied (default 1000).}
#'   \item{\code{maxiter}}{Maximum number of iterations for the Gauss--Newton
#'     updates (default 1000).}
#'   \item{\code{mindiff}}{Convergence tolerance on the increase of the mean
#'     cosine (default 1e-6).}
#' }
#'
#' @param formula A formula with the dependent angle on the left of the
#'   \code{~} operator and terms specifying the explanatory variables on the
#'   right. The first term provides the reference direction.
#' @param data An optional data frame, list or environment containing the
#'   variables in the model formula.
#' @param initbeta A numeric vector of initial beta parameter values. When
#'   \code{NULL} (default) a quasi-uniform grid search is used to select
#'   starting values. The required length equals the number of \code{x:z} terms
#'   after the reference direction.
#' @param control A list of control parameters. See Details.
#' @return An object of class "angular" containing:
#' \describe{
#'   \item{MaxCosine}{the maximum cosine value.}
#'   \item{parameters}{the parameter estimates and their standard errors, z-values and associated p-values.}
#'   \item{varcov0}{the asymptotic estimated variance-covariance matrix.}
#'   \item{varcov1}{the non-parametric estimated variance-covariance matrix (sandwich).}
#'   \item{long}{the vector of predicted concentrations.}
#'   \item{mui}{the vector of predicted mean angles.}
#'   \item{y}{the response variable.}
#'   \item{iter.detail}{the iteration details.}
#'   \item{call}{the function call.}
#' }
#'
#' @author Sophie Baillargeon, Louis-Paul Rivest, and Aurelien Nicosia
#' @references L.-P. Rivest, T. Duchesne, A. Nicosia & D. Fortin. A general angular regression model for the analysis of data on animal movement in ecology. Journal of the Royal Statistical Society, series C, to appear.
#' @examples
#' \dontrun{
#'   # Example with the bison data set included in the package
#'   data(bison)
#'   an <- angular(y.dir ~ y.prec + y.prec2, data = bison)
#'   print(an)
#'   summary(an)
#'   plot(an)
#' }
#' @export
angular <- function(
  formula,
  data,
  initbeta = NULL,
  control = list()
) {
  call <- match.call()

  mfargs <- match(c("formula", "data"), names(call), 0L)
  call_subset <- call[c(1L, mfargs)]
  call_subset[[1L]] <- quote(model.frame)
  mf <- eval(call_subset, parent.frame())
  nobs <- nrow(mf)
  if (nobs == 0) stop("No observations available.")

  y <- mf[[1]]
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

  compute_components <- function(beta) {
    sinmu <- sin(x0)
    cosmu <- cos(x0)
    if (p > 0) {
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta)
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta)
    }
    long <- sqrt(sinmu^2 + cosmu^2)
    mui <- atan2(sinmu, cosmu)
    list(long = long, mui = mui)
  }

  maxcos <- function(beta) {
    cm <- compute_components(beta)
    maxcos_val <- mean(cos(y - cm$mui))
    list(maxcos = maxcos_val, long = cm$long, mui = cm$mui)
  }

  betaUpdate <- function(betak, long, mui) {
    matd <- (matz * sin(matx - mui)) / as.vector(long)
    res <- sin(y - mui)
    qr_obj <- qr(matd)
    if (qr_obj$rank < ncol(matd)) {
      stop(
        "Design matrix is rank deficient; cannot update coefficients.",
        call. = FALSE
      )
    }
    dbeta <- as.vector(qr.coef(qr_obj, res))
    betak1 <- betak + dbeta
    list(betak1 = betak1, dbeta = dbeta, matd = matd, res = res)
  }

  if (p == 0) {
    betak <- numeric(0)
  } else if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 1000 else control$pginit
    pg <- max(1L, round(pginit^(1 / p)))
    grid_vals <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]), p)
    possVal <- cbind(expand.grid(grid_vals), NA_real_)
    colnames(possVal) <- c(betaname, "maxcosine")
    maxcos1 <- function(beta) maxcos(beta)$maxcos
    possVal[, p + 1] <- apply(possVal[, seq_len(p), drop = FALSE], 1, maxcos1)
    betak <- as.numeric(possVal[
      which.max(possVal[, p + 1]),
      seq_len(p),
      drop = TRUE
    ])
  } else {
    if (length(initbeta) != p) {
      stop(
        sprintf("'initbeta' must have length %d for this model.", p),
        call. = FALSE
      )
    }
    betak <- initbeta
  }

  calcul <- maxcos(beta = betak)
  maxcosk <- calcul$maxcos
  long <- calcul$long
  mui <- calcul$mui
  iter <- iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 1e-06 else control$mindiff
  iter.detail <- matrix(NA_real_, nrow = maxiter + 1, ncol = p + 3)
  colnames(iter.detail) <- c(betaname, "maxcosine", "iter", "nitersh")
  iter.detail[1, ] <- c(betak, maxcosk, iter, iter.sh)
  maxcosk1 <- maxcosk

  if (p > 0) {
    conv <- FALSE
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
        warning(
          "The algorithm did not converge; it failed to maximise the mean cosine."
        )
        conv <- FALSE
        break
      } else {
        conv <- (maxcosk1 - maxcosk) <= mindiff
        betak <- betak1
        maxcosk <- maxcosk1
        iter <- iter + 1
        iter.detail[iter + 1, ] <- c(betak, maxcosk, iter, iter.sh)
      }
    }
    if (iter > maxiter + 1) {
      warning(
        "The algorithm did not converge: maximum number of iterations reached."
      )
    } else {
      iter.detail <- iter.detail[seq_len(iter + 1), , drop = FALSE]
    }
    if (maxcosk == maxcosk1) {
      maj <- betaUpdate(betak = betak, long = long, mui = mui)
    }
    matd <- maj$matd
    res <- maj$res
  } else {
    iter.detail <- iter.detail[1, , drop = FALSE]
    matd <- matrix(0, nrow = nobs, ncol = 0)
    res <- sin(y - mui)
  }

  Akappa <- maxcosk
  kappahat <- circular::A1inv(Akappa)

  if (p > 0) {
    XtX <- crossprod(matd)
    qr_X <- qr(XtX)
    if (qr_X$rank < ncol(XtX)) {
      stop(
        "Information matrix is singular; variance estimation failed.",
        call. = FALSE
      )
    }
    XtX_inv <- qr.solve(qr_X, diag(ncol(XtX)))
    v0 <- XtX_inv / (Akappa * kappahat)

    res_sq <- res^2
    matd_res <- matd * res_sq
    v2 <- XtX_inv %*% crossprod(matd, matd_res) %*% XtX_inv

    zvalue <- betak / sqrt(diag(v0))
    pvals <- 2 * stats::pnorm(abs(zvalue), lower.tail = FALSE)
    parameters <- cbind(
      estimate = betak,
      stderr = sqrt(diag(v0)),
      `z value` = zvalue,
      `P(|z|>.)` = pvals
    )
    rownames(parameters) <- betaname
  } else {
    v0 <- v2 <- matrix(0, nrow = 0, ncol = 0)
    parameters <- matrix(
      numeric(0),
      nrow = 0,
      ncol = 4,
      dimnames = list(NULL, c("estimate", "stderr", "z value", "P(|z|>.)"))
    )
  }

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
  if (nrow(coefmat) == 0) {
    cat("  (no regression coefficients; reference-only model)\n")
  } else {
    mat_numeric <- cbind(
      Estimate = as.numeric(coefmat[, 1]),
      `Std. Error` = as.numeric(coefmat[, 2]),
      `z value` = as.numeric(coefmat[, 3]),
      `P(|z|>.)` = as.numeric(coefmat[, 4])
    )
    rownames(mat_numeric) <- rownames(coefmat)
    stats::printCoefmat(mat_numeric, P.values = TRUE, has.Pvalue = TRUE)
  }
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
  # keep everything in radians to avoid ambiguity
  y <- circular::conversion.circular(
    circular::as.circular(object$y),
    units = "radians"
  )
  mu <- circular::conversion.circular(
    circular::as.circular(object$mui),
    units = "radians"
  )

  # signed angular difference in (-pi, pi]
  d <- atan2(sin(y - mu), cos(y - mu))

  # renvoyer un 'circular' en radians (convertis ensuite si tu veux)
  circular::as.circular(d, units = "radians", modulo = "asis")
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
#' This function produces a set of diagnostic plots for an object of class
#' "angular" in a layout similar to that of \code{plot.lm}. Three panels are
#' displayed:
#' \enumerate{
#'   \item Residuals vs fitted values.
#'   \item Histogram of residuals with a normal density overlay.
#'   \item Normal Q-Q plot.
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
  fitted_vals <- as.numeric(x$mui)
  residuals_vec <- residuals.angular(x)
  residuals_num <- as.numeric(residuals_vec)

  p1 <- ggplot2::ggplot(
    data = data.frame(Fitted = fitted_vals, Residual = residuals_num),
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

  res_mean <- mean(residuals_num)
  res_sd <- stats::sd(residuals_num)
  p2 <- ggplot2::ggplot(
    data = data.frame(Residual = residuals_num),
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
    ) +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(
    data = data.frame(Residual = residuals_num),
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

  layout <- rbind(c(1, 2), c(3, 3))
  gridExtra::grid.arrange(
    p1,
    p2,
    p3,
    layout_matrix = layout,
    top = "Diagnostic Plots for Angular Regression Model"
  )

  invisible(x)
}
