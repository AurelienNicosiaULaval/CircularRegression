###############################################################################
# Angular Regression Model (Homogeneous errors)
###############################################################################

#' Angular Regression Model with Homogeneous von Mises Errors
#'
#' Fits the homogeneous angular regression model of Rivest et al. (2016), where
#' the mean direction is the orientation of a weighted resultant vector and one
#' reference angle has coefficient fixed to either 1 or -1.
#'
#' @param formula A formula with terms of the form \code{x} or \code{x:z}.
#' @param data A data frame, list or environment containing model variables.
#' @param reference Reference-angle rule. Use \code{"auto"} (default, article
#'   criterion), \code{"first"}, or \code{c("name", "x_ref")}. A plain
#'   angle term \code{x_ref} must exist in the formula.
#' @param reference_sign Sign of the fixed reference coefficient. The default is
#'   \code{1}; \code{-1} selects the opposite signed coordinate chart without
#'   changing the scale normalization.
#' @param initbeta Optional numeric vector of initial values for non-reference
#'   coefficients.
#' @param control Control list with optional components \code{pginit},
#'   \code{maxiter}, and \code{mindiff}.
#' @param na.action Function used by \code{\link[stats]{model.frame}} to handle
#'   missing values. The default is \code{\link[stats]{na.omit}}.
#'   With \code{na.exclude}, in-sample fitted values, residuals and predictions
#'   restore omitted rows as \code{NA}. Angular \code{circular} objects must
#'   declare radian units; degree or hour units are rejected before coercion.
#' @return An object of class \code{"angular"}.
#' @details Convergence requires a scaled score norm no larger than
#'   \code{control$mindiff}; a small objective increment alone is insufficient.
#'   The object stores \code{convergence}, \code{score_norm}, and the observed
#'   derivative matrix. A positive-definite derivative matrix is required for
#'   sandwich inference. Both covariance choices assume independent observations;
#'   neither is a serial-dependence correction. Zero resultant vectors are
#'   excluded because they have no mean direction.
#'   The returned \code{fit_status} distinguishes convergence from zero or
#'   infinite concentration and singular information. Wald covariances are
#'   unavailable when convergence is not established.
#' @examples
#' data(bison)
#' d <- na.omit(bison[1:100, c("y.dir", "y.prec", "x.meadow", "z.meadow")])
#' fit <- angular(y.dir ~ y.prec + x.meadow:z.meadow, d, reference = "first")
#' coef(fit)
#' fit$convergence
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @export
angular <- function(
  formula,
  data,
  reference = c("auto", "first", "name"),
  initbeta = NULL,
  control = list(),
  na.action = stats::na.omit,
  reference_sign = 1
) {
  call <- match.call()
  des <- .angular_design(
    formula = formula,
    data = data,
    reference = reference,
    na.action = na.action
  )

  y <- des$y
  x0 <- des$x0
  matx <- des$matx
  matz <- des$matz
  betaname <- des$betaname
  p <- length(betaname)
  nobs <- des$nobs

  if (!is.numeric(reference_sign) || length(reference_sign) != 1L ||
      !is.finite(reference_sign) || !(reference_sign %in% c(-1, 1))) {
    stop("'reference_sign' must be either 1 or -1.", call. = FALSE)
  }
  reference_sign <- as.numeric(reference_sign)

  if (!is.null(initbeta) && (!is.numeric(initbeta) || any(!is.finite(initbeta)) || length(initbeta) != p)) {
    stop(sprintf("'initbeta' must have length %d for this model.", p), call. = FALSE)
  }

  ctrl <- .check_fit_control(
    control = control,
    defaults = list(pginit = 1000, maxiter = 1000, mindiff = 1e-06)
  )
  pginit <- ctrl$pginit
  maxiter <- ctrl$maxiter
  mindiff <- ctrl$mindiff

  compute_components <- function(beta) {
    sinmu <- reference_sign * sin(x0)
    cosmu <- reference_sign * cos(x0)
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
    scale <- 1 + if (p > 0L) as.vector(abs(matz) %*% abs(beta)) else 0
    if (any(cm$long <= sqrt(.Machine$double.eps) * scale)) {
      return(list(maxcos = -Inf, long = cm$long, mui = cm$mui))
    }
    val <- mean(cos(y - cm$mui))
    list(maxcos = val, long = cm$long, mui = cm$mui)
  }

  beta_update <- function(betak, long, mui) {
    deriv <- .angular_derivatives(y, mui, long, matx, matz)
    matd <- deriv$D
    res <- deriv$residual

    qr_obj <- qr(matd)
    if (qr_obj$rank < ncol(matd)) {
      stop(
        "Design matrix is rank deficient; homogeneous model is not identifiable for this reference.",
        call. = FALSE
      )
    }

    dbeta <- as.vector(qr.coef(qr_obj, res))
    # Newton acceleration near a local maximum; otherwise retain an ascent
    # Gauss-Newton step and let the line search control the step length.
    ch <- tryCatch(chol(deriv$information), error = function(e) NULL)
    if (!is.null(ch)) dbeta <- as.vector(chol2inv(ch) %*% deriv$score)
    list(betak1 = betak + dbeta, dbeta = dbeta, matd = matd, res = res,
         score = deriv$score, information = deriv$information,
         score_norm = .score_norm(deriv$score, crossprod(matd), nobs))
  }

  if (p == 0) {
    betak <- numeric(0)
  } else if (is.null(initbeta)) {
    pg <- max(1L, round(pginit^(1 / p)))
    grid_vals <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]), p)
    poss <- cbind(expand.grid(grid_vals), NA_real_)
    colnames(poss) <- c(betaname, "maxcos")
    poss[, p + 1] <- apply(
      poss[, seq_len(p), drop = FALSE],
      1,
      function(b) maxcos(as.numeric(b))$maxcos
    )
    betak <- as.numeric(poss[which.max(poss[, p + 1]), seq_len(p), drop = TRUE])
  } else {
    betak <- as.numeric(initbeta)
  }

  calc <- maxcos(betak)
  if (!is.finite(calc$maxcos)) {
    stop("Initial values produce a zero or numerically singular resultant vector.", call. = FALSE)
  }
  maxcosk <- calc$maxcos
  long <- calc$long
  mui <- calc$mui

  iter <- 0L
  conv <- p == 0L
  iter.detail <- matrix(NA_real_, nrow = maxiter + 1L, ncol = p + 3L)
  colnames(iter.detail) <- c(betaname, "maxcosine", "iter", "nitersh")
  iter.detail[1L, ] <- c(betak, maxcosk, iter, 0)

  if (p > 0 && maxiter > 0) {
    while (!conv && iter < maxiter) {
      step <- beta_update(betak = betak, long = long, mui = mui)
      if (step$score_norm <= mindiff) {
        conv <- TRUE
        break
      }
      betak1 <- step$betak1
      dbeta <- step$dbeta

      calc1 <- maxcos(beta = betak1)
      maxcosk1 <- calc1$maxcos
      long1 <- calc1$long
      mui1 <- calc1$mui

      iter.sh <- 0L
      while (maxcosk1 < maxcosk && iter.sh < maxiter) {
        iter.sh <- iter.sh + 1L
        betak1 <- betak + dbeta / (2^iter.sh)
        calc1 <- maxcos(beta = betak1)
        maxcosk1 <- calc1$maxcos
        long1 <- calc1$long
        mui1 <- calc1$mui
      }

      if (maxcosk1 < maxcosk) {
        warning(
          "The algorithm did not converge; step-halving failed to improve the objective.",
          call. = FALSE
        )
        break
      }

      if (identical(betak1, betak) ||
          (maxcosk1 == maxcosk && iter.sh > 20L)) {
        warning("Angular iteration stagnated with a nonzero score; convergence was not established.", call. = FALSE)
        break
      }
      betak <- betak1
      maxcosk <- maxcosk1
      long <- long1
      mui <- mui1

      iter <- iter + 1L
      iter.detail[iter + 1L, ] <- c(betak, maxcosk, iter, iter.sh)
    }

    if (!conv && iter >= maxiter) {
      warning("Maximum number of iterations reached before convergence.", call. = FALSE)
    }
  }

  iter.detail <- iter.detail[seq_len(iter + 1L), , drop = FALSE]

  # explicit identifiability check (Proposition 1 matrix rank p-1)
  nterm <- length(des$term_labels)
  beta_full <- numeric(nterm)
  beta_full[des$ref_idx] <- reference_sign
  if (p > 0) beta_full[-des$ref_idx] <- betak

  id_mat <- .identifiability_matrix(
    beta_full = beta_full,
    matx_full = des$full_matx,
    matz_full = des$full_matz
  )
  id_rank <- qr(id_mat)$rank
  if (id_rank < (nterm - 1L)) {
    stop(
      sprintf(
        "Model is not identifiable for this reference: rank=%d < %d.",
        id_rank,
        nterm - 1L
      ),
      call. = FALSE
    )
  }

  if (p > 0) {
    last_step <- beta_update(betak = betak, long = long, mui = mui)
    matd <- last_step$matd
    res <- last_step$res
    conv <- conv || (maxiter > 0L && last_step$score_norm <= mindiff)
  } else {
    matd <- matrix(0, nrow = nobs, ncol = 0)
    res <- sin(y - mui)
  }

  Akappa <- as.numeric(maxcosk)
  kappahat <- .A1_inverse(min(max(Akappa, 0), 1))
  if (!is.finite(kappahat)) {
    warning("A unit mean cosine has no finite concentration MLE; Wald inference is unavailable.", call. = FALSE)
  }
  if (!is.finite(kappahat) || kappahat == 0) conv <- FALSE

  if (p > 0) {
    XtX <- crossprod(matd)
    XtX_sym <- (XtX + t(XtX)) / 2
    chol_xtx <- tryCatch(chol(XtX_sym), error = function(e) NULL)
    if (is.null(chol_xtx)) {
      stop("Information matrix is singular; variance estimation failed.", call. = FALSE)
    }
    XtX_inv <- chol2inv(chol_xtx)

    if (Akappa > 0 && kappahat > 0 && is.finite(kappahat)) {
      v0 <- XtX_inv / (Akappa * kappahat)
      se0 <- sqrt(diag(v0))
      zvalue <- betak / se0
      pvals <- 2 * stats::pnorm(abs(zvalue), lower.tail = FALSE)
    } else {
      v0 <- matrix(NA_real_, nrow = p, ncol = p)
      se0 <- rep(NA_real_, p)
      zvalue <- rep(NA_real_, p)
      pvals <- rep(NA_real_, p)
    }

    # Rivest et al. (2016), equations 13-14: the observed derivative
    # includes the second derivative of the mean direction, not just D'D.
    bread <- .information_inverse(last_step$information)
    v2 <- bread %*% crossprod(matd * res) %*% bread
    if (!is.finite(kappahat) || kappahat == 0) v2[,] <- NA_real_

    parameters <- cbind(
      estimate = betak,
      stderr = se0,
      `z value` = zvalue,
      `P(|z|>.)` = pvals
    )
    rownames(parameters) <- betaname
  } else {
    v0 <- matrix(0, nrow = 0, ncol = 0)
    v2 <- matrix(0, nrow = 0, ncol = 0)
    parameters <- matrix(
      numeric(0),
      nrow = 0,
      ncol = 4,
      dimnames = list(NULL, c("estimate", "stderr", "z value", "P(|z|>.)"))
    )
  }

  inference_available <- conv && all(is.finite(v0)) && all(is.finite(v2))
  if (!conv) {
    v0[,] <- NA_real_
    v2[,] <- NA_real_
    if (p > 0L) parameters[, -1L] <- NA_real_
  }
  v0 <- .named_covariance(v0, betaname)
  v2 <- .named_covariance(v2, betaname)
  fit_status <- if (!is.finite(kappahat)) "infinite_concentration" else if (kappahat == 0) {
    "zero_concentration"
  } else if (!conv) "not_converged" else if (!inference_available) {
    "singular_information"
  } else "converged"

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
    convergence = conv,
    fit_status = fit_status,
    inference_available = inference_available,
    convergence_tolerance = mindiff,
    score = if (p > 0L) last_step$score else numeric(0),
    score_norm = if (p > 0L) last_step$score_norm else 0,
    information = if (p > 0L) last_step$information else matrix(numeric(0), 0, 0),
    model_data = des$mf,
    call = call,
    formula = formula,
    nobs = nobs,
    k = p + 1L,
    reference = des$reference$ref_name,
    reference_sign = reference_sign,
    reference_scores = des$reference$scores,
    term_info = des$term_info,
    term_labels = des$term_labels,
    ref_idx = des$ref_idx,
    full_coefficients = beta_full,
    response_name = des$response_name,
    na.action = attr(des$mf, "na.action")
  )
  class(out) <- "angular"
  out
}

###############################################################################
### S3 Methods for Angular Regression Objects
###############################################################################

#' Methods for Angular Regression Objects
#'
#' @details The \code{plot()} method uses wrapped residuals
#' \eqn{e_i = \mathrm{atan2}(\sin(y_i-\hat\mu_i),
#' \cos(y_i-\hat\mu_i))}. It provides a circular rose and density display,
#' a von Mises probability-integral-transform (PIT) Q-Q plot, and residual
#' plots against fitted direction and observation order. For an \code{angular}
#' fit, the rose overlays the fitted \eqn{\mathrm{VM}(0,\hat\kappa)} density
#' and the PIT uses the same common concentration. The Q-Q plot is descriptive:
#' its diagonal is not a formal goodness-of-fit test because the parameters were
#' estimated from the displayed observations.
#'
#' The diagnostic plots require the suggested packages \pkg{ggplot2} and
#' \pkg{ggcircular}; \code{type = "all"} additionally requires
#' \pkg{gridExtra}.
#'
#' @param object An object of class \code{"angular"}.
#' @param x An object of class \code{"angular"}.
#' @param newdata Optional data frame for prediction.
#' @param type For \code{predict()}, use \code{"response"} for fitted angles in
#'   radians or \code{"components"} for cosine and sine components. For
#'   \code{plot()}, select \code{"all"}, \code{"rose"}, \code{"pit"},
#'   \code{"fitted"}, or \code{"index"}.
#' @param bins Number of bins in a circular residual rose plot.
#' @param robust Logical; for \code{vcov()}, return the sandwich covariance
#'   matrix when \code{TRUE}.
#' @param ... Additional arguments.
#' @rdname angular-methods
#' @examples
#' data(bison)
#' fit <- angular(y.dir ~ y.prec, bison[1:100, ])
#' summary(fit)
#' head(predict(fit))
#' head(residuals(fit))
#' vcov(fit, robust = TRUE)
#' @export
print.angular <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  .print_fit_status(x)
  cat("\nReference angle:", x$reference, "\n")
  cat("Reference coefficient:", x$reference_sign %||% 1, "\n")
  cat("Maximum cosine:", format(x$MaxCosine, digits = 4), "\n")
  cat("Concentration parameter:", format(x$kappahat, digits = 4), "\n\n")
  cat("Parameters (non-reference terms):\n")

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
  object$parameters[, "estimate"]
}

#' @rdname angular-methods
#' @export
vcov.angular <- function(object, robust = FALSE, ...) {
  .check_logical_scalar(robust, "robust")
  if (isTRUE(robust)) object$varcov1 else object$varcov0
}

#' @rdname angular-methods
#' @export
fitted.angular <- function(object, ...) {
  .restore_excluded(object, object$mui)
}

#' @rdname angular-methods
#' @export
residuals.angular <- function(object, ...) {
  .restore_excluded(object, .wrap_angle(object$y - object$mui), residual = TRUE)
}

#' @rdname angular-methods
#' @export
predict.angular <- function(
  object,
  newdata = NULL,
  type = c("response", "components"),
  ...
) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    mu <- .restore_excluded(object, object$mui)
  } else {
    design <- .new_term_design(
      term_info = object$term_info,
      formula = object$formula,
      newdata = newdata
    )
    mu <- .mean_direction_from_terms(
      coef = object$full_coefficients,
      matx = design$matx,
      matz = design$matz
    )
  }

  if (type == "components") {
    return(cbind(cos = cos(mu), sin = sin(mu)))
  }
  mu
}

#' @rdname angular-methods
#' @export
summary.angular <- function(object, ...) {
  res <- residuals.angular(object)
  s <- list(
    call = object$call,
    MaxCosine = object$MaxCosine,
    parameters = object$parameters,
    kappahat = object$kappahat,
    residuals = summary(res),
    nobs = object$nobs,
    reference = object$reference,
    reference_sign = object$reference_sign %||% 1
  )
  s$convergence <- object$convergence
  s$fit_status <- object$fit_status
  s$inference_available <- object$inference_available
  class(s) <- "summary.angular"
  s
}

#' @rdname angular-methods
#' @export
print.summary.angular <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  .print_fit_status(x)
  cat("\nReference angle:", x$reference, "\n")
  cat("Reference coefficient:", x$reference_sign %||% 1, "\n")
  cat("Maximum cosine:", format(x$MaxCosine, digits = 4), "\n")
  cat("Concentration parameter:", format(x$kappahat, digits = 4), "\n\n")
  cat("Residuals:\n")
  print(x$residuals)
  cat("\nParameters:\n")
  if (nrow(x$parameters) > 0) {
    stats::printCoefmat(
      x$parameters,
      signif.stars = TRUE,
      signif.legend = TRUE,
      ...
    )
  } else {
    cat("  (no regression coefficients; reference-only model)\n")
  }
  cat("\nNumber of observations:", x$nobs, "\n")
  invisible(x)
}

#' @rdname angular-methods
#' @export
plot.angular <- function(
  x,
  type = c("all", "rose", "pit", "fitted", "index"),
  bins = 24,
  ...
) {
  type <- match.arg(type)
  .plot_fixed_circular_diagnostics(
    x = x,
    residuals = .wrap_angle(x$y - x$mui),
    fitted = x$mui,
    kappa = x$kappahat,
    type = type,
    bins = bins,
    title = sprintf("Homogeneous angular model diagnostics (reference: %s)", x$reference)
  )
}
