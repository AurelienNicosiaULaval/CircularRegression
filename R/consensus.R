###############################################################################
# Consensus Model
###############################################################################

#' Consensus Angular Regression Model
#'
#' Fits the consensus model of Rivest et al. (2016), where each formula term of
#' the form \code{x} or \code{x:z} receives a \eqn{\kappa_j} coefficient in the
#' log-likelihood parameterization.
#'
#' @param formula Model formula containing terms of the form \code{x} or
#'   \code{x:z}.
#' @param data Data frame, list or environment containing model variables.
#' @param weights Optional non-negative numeric case weights. The model
#'   covariance is inverse weighted information; the independent-observation
#'   sandwich uses the squared weights. These are not cluster weights.
#' @param initkappa Optional numeric vector of initial \eqn{\kappa_j} values.
#' @param initbeta Deprecated alias for \code{initkappa}, retained for
#'   backward compatibility.
#' @param control Optional list with \code{pginit}, \code{maxiter},
#'   \code{mindiff}.
#' @param na.action Function used by \code{\link[stats]{model.frame}} to handle
#'   missing values. The default is \code{\link[stats]{na.omit}}.
#'   \code{na.exclude} restores omitted rows in in-sample fitted values,
#'   residuals and predictions. Weights may match the original or retained
#'   rows, including when \code{data} is a list or environment.
#' @return An object of class \code{"consensus"}.
#' @details The object stores \code{convergence}, \code{score_norm}, model
#'   covariance \code{varcov_kappa}, and \code{varcov_robust}. The legacy
#'   \code{varcov1} remains the model covariance. Convergence requires a scaled
#'   score norm at most \code{control$mindiff}. At zero concentration, the
#'   likelihood and its derivatives remain defined but the fitted mean direction
#'   is \code{NA}. The weighted BIC uses the number of retained rows, not the sum
#'   of weights; information criteria under arbitrary case weighting should not
#'   be interpreted as ordinary unweighted likelihood comparisons.
#'
#'   A recession direction aligns every nonzero directional contribution with
#'   its response and has no oppositely aligned contribution. The likelihood
#'   then has no finite maximum, including when some modifiers are zero.
#'   The object records \code{fit_status = "unbounded_likelihood"}, sets
#'   \code{convergence = FALSE}, and returns unavailable Wald covariances.
#'   In this case the coefficients are starting values, not maximum likelihood
#'   estimates. The numerical check ignores zero-weight rows and verifies the
#'   returned direction within floating-point tolerance. It is a boundary
#'   diagnostic, not a global convergence certificate for arbitrary data.
#'   \code{inference_available} indicates whether the fitted covariance is usable.
#' @examples
#' data(bison)
#' fit <- consensus(y.dir ~ y.prec, bison[1:100, ])
#' coef(fit)
#' vcov(fit, robust = TRUE)
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @export
consensus <- function(
  formula,
  data,
  weights = NULL,
  initkappa = NULL,
  initbeta = NULL,
  control = list(),
  na.action = stats::na.omit
) {
  call <- match.call()
  des <- .consensus_design(formula = formula, data = data, na.action = na.action)

  y <- des$y
  matx <- des$matx
  matz <- des$matz
  paramname <- des$term_labels
  term_info <- des$term_info
  nobs <- des$nobs
  nparam <- length(paramname)

  if (nparam == 0) {
    stop("The model must include at least one explanatory term.", call. = FALSE)
  }

  if (is.null(weights)) {
    weight <- rep(1, nobs)
  } else {
    if (!is.numeric(weights)) stop("'weights' must be numeric.", call. = FALSE)
    w_raw <- as.numeric(weights)
    weight <- .align_model_vector(w_raw, des$mf, "weights")
  }

  if (any(!is.finite(weight)) || any(weight < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (all(weight == 0)) {
    stop("'weights' cannot be all zero.", call. = FALSE)
  }
  weight_total <- sum(weight)
  if (!is.finite(weight_total)) stop("The sum of weights must be finite.", call. = FALSE)
  recession_direction <- .consensus_recession(y, matx, matz, weight)
  alignment_boundary <- !is.null(recession_direction)

  if (!is.null(initkappa) && !is.null(initbeta)) {
    stop("Use only one of 'initkappa' or 'initbeta'.", call. = FALSE)
  }
  if (is.null(initkappa) && !is.null(initbeta)) {
    initkappa <- initbeta
  }

  if (!is.null(initkappa) && (!is.numeric(initkappa) || any(!is.finite(initkappa)) || length(initkappa) != nparam)) {
    stop(sprintf("'initkappa' must have length %d.", nparam), call. = FALSE)
  }

  ctrl <- .check_fit_control(
    control = control,
    defaults = list(pginit = 1000, maxiter = 1000, mindiff = 1e-06)
  )
  pginit <- ctrl$pginit
  maxiter <- ctrl$maxiter
  mindiff <- ctrl$mindiff

  cos_y_minus_x <- cos(matrix(y, nrow = nobs, ncol = nparam) - matx)

  compute_components <- function(kappa) {
    sinmu <- as.vector((matz * sin(matx)) %*% kappa)
    cosmu <- as.vector((matz * cos(matx)) %*% kappa)
    long <- sqrt(sinmu^2 + cosmu^2)
    mui <- atan2(sinmu, cosmu)
    list(long = long, mui = mui)
  }

  loglik_components <- function(kappa) {
    comp <- compute_components(kappa)
    term1 <- as.vector((matz * cos_y_minus_x) %*% kappa)
    ll <- sum(weight * term1) -
      sum(weight * .logI0(comp$long)) -
      weight_total * log(2 * pi)
    if (!is.finite(ll)) {
      ll <- -Inf
    }
    list(LL = ll, long = comp$long, mui = comp$mui)
  }

  fisher_step <- function(kappa, long, mui) {
    deriv <- .consensus_derivatives(y, mui, long, matx, matz, weight)
    score <- deriv$score
    I <- deriv$information

    qrI <- qr(I)
    if (qrI$rank < ncol(I)) {
      stop(
        "Consensus information matrix is singular; parameters are not identifiable.",
        call. = FALSE
      )
    }

    dparam <- as.vector(qr.coef(qrI, score))
    list(kappa1 = kappa + dparam, dparam = dparam, I = I,
         scores = deriv$scores, score = score,
         score_norm = .score_norm(score, I, weight_total))
  }

  if (is.null(initkappa)) {
    pg <- max(1L, round(pginit^(1 / nparam)))
    grid_vals <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]), nparam)
    poss <- cbind(expand.grid(grid_vals), NA_real_)
    colnames(poss) <- c(paramname, "LL")
    poss[, nparam + 1] <- apply(
      poss[, seq_len(nparam), drop = FALSE],
      1,
      function(par) loglik_components(as.numeric(par))$LL
    )
    finite_start <- is.finite(poss[, nparam + 1])
    if (!any(finite_start)) {
      stop("Unable to find a finite initial log-likelihood for the consensus model.", call. = FALSE)
    }
    poss[!finite_start, nparam + 1] <- -Inf
    kappak <- as.numeric(poss[which.max(poss[, nparam + 1]), seq_len(nparam), drop = TRUE])
  } else {
    kappak <- as.numeric(initkappa)
  }

  calc <- loglik_components(kappak)
  maxLLk <- calc$LL
  long <- calc$long
  mui <- calc$mui

  if (!is.finite(maxLLk)) {
    stop("Initial values produce a non-finite consensus log-likelihood.", call. = FALSE)
  }

  is_worse_loglik <- function(new, old) {
    !is.finite(new) || new < old
  }

  iter <- 0L
  conv <- FALSE
  iter.detail <- matrix(NA_real_, nrow = maxiter + 1L, ncol = nparam + 3L)
  colnames(iter.detail) <- c(paramname, "LL", "iter", "nitersh")
  iter.detail[1L, ] <- c(kappak, maxLLk, iter, 0)

  if (maxiter > 0 && !alignment_boundary) {
    while (!conv && iter < maxiter) {
      step <- fisher_step(kappak, long = long, mui = mui)
      if (step$score_norm <= mindiff) {
        conv <- TRUE
        break
      }
      kappak1 <- step$kappa1
      dparam <- step$dparam

      calc1 <- loglik_components(kappak1)
      maxLLk1 <- calc1$LL
      long1 <- calc1$long
      mui1 <- calc1$mui

      iter.sh <- 0L
      while (is_worse_loglik(maxLLk1, maxLLk) && iter.sh < maxiter) {
        iter.sh <- iter.sh + 1L
        kappak1 <- kappak + dparam / (2^iter.sh)
        calc1 <- loglik_components(kappak1)
        maxLLk1 <- calc1$LL
        long1 <- calc1$long
        mui1 <- calc1$mui
      }

      if (is_worse_loglik(maxLLk1, maxLLk)) {
        warning(
          "The algorithm did not converge; step-halving failed to improve the log-likelihood.",
          call. = FALSE
        )
        break
      }

      if (identical(kappak1, kappak) ||
          (maxLLk1 == maxLLk && iter.sh > 20L)) {
        warning("Consensus iteration stagnated with a nonzero score; convergence was not established.", call. = FALSE)
        break
      }
      kappak <- kappak1
      maxLLk <- maxLLk1
      long <- long1
      mui <- mui1

      iter <- iter + 1L
      iter.detail[iter + 1L, ] <- c(kappak, maxLLk, iter, iter.sh)
    }

    if (!conv && iter >= maxiter) {
      warning("Maximum number of iterations reached before convergence.", call. = FALSE)
    }
  }

  iter.detail <- iter.detail[seq_len(iter + 1L), , drop = FALSE]

  final_step <- fisher_step(kappak, long = long, mui = mui)
  conv <- !alignment_boundary &&
    (conv || (maxiter > 0L && final_step$score_norm <= mindiff))
  matI <- final_step$I
  if (alignment_boundary) {
    warning("A direction of perfect alignment makes the consensus likelihood unbounded; no finite maximum exists and Wald inference is unavailable.", call. = FALSE)
  }

  invert_information <- function(M) {
    M_sym <- (M + t(M)) / 2
    chol_res <- tryCatch(chol(M_sym), error = function(e) NULL)
    if (is.null(chol_res)) {
      warning(
        "Information matrix is not positive definite; returning NA variances.",
        call. = FALSE
      )
      matrix(NA_real_, nrow = nrow(M_sym), ncol = ncol(M_sym))
    } else {
      chol2inv(chol_res)
    }
  }

  v_kappa <- invert_information(matI)
  inference_available <- conv && all(is.finite(v_kappa))
  if (!inference_available) v_kappa[,] <- NA_real_
  v_robust <- v_kappa %*% crossprod(final_step$scores) %*% v_kappa
  v_kappa <- .named_covariance(v_kappa, paramname)
  v_robust <- .named_covariance(v_robust, paramname)
  fit_status <- if (alignment_boundary) "unbounded_likelihood" else if (!conv) {
    "not_converged"
  } else if (!inference_available) "singular_information" else "converged"
  se_kappa <- sqrt(diag(v_kappa))
  zvalue <- kappak / se_kappa
  pvals <- 2 * stats::pnorm(abs(zvalue), lower.tail = FALSE)

  parameters <- cbind(
    estimate = kappak,
    `Std. Error` = se_kappa,
    `z value` = zvalue,
    `P(|z|>.)` = pvals
  )
  rownames(parameters) <- paramname

  k <- nparam
  logLik <- as.numeric(maxLLk)
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + log(nobs) * k
  if (!conv) AIC <- BIC <- NA_real_

  plain_angles <- des$plain_angles
  reference_scores <- if (length(plain_angles) > 0) {
    .reference_scores(y = y, angle_data = des$angle_data, candidates = plain_angles)
  } else {
    numeric(0)
  }

  out <- list(
    MaxLL = maxLLk,
    AIC = AIC,
    BIC = BIC,
    parameters = parameters,
    kappa = stats::setNames(kappak, paramname),
    varcov_kappa = v_kappa,
    varcov1 = v_kappa,
    varcov_robust = v_robust,
    information = matI,
    score = final_step$score,
    score_norm = final_step$score_norm,
    convergence = conv,
    fit_status = fit_status,
    inference_available = inference_available,
    alignment_boundary = alignment_boundary,
    recession_direction = if (alignment_boundary) stats::setNames(recession_direction, paramname) else NULL,
    convergence_tolerance = mindiff,
    parambeta = NULL,
    varcovbeta1 = NULL,
    matx = matx,
    matz = matz,
    y = y,
    long = long,
    mui = ifelse(long == 0, NA_real_, mui),
    mean_defined = long > 0,
    weights = weight,
    iter.detail = iter.detail,
    call = call,
    formula = formula,
    nobs = nobs,
    k = k,
    logLik = logLik,
    term_info = term_info,
    paramname = paramname,
    angle_data = des$angle_data,
    plain_angles = plain_angles,
    reference_scores = reference_scores,
    autocorr = if (all(long > 0)) stats::acf(residuals.consensus(list(y = y, mui = mui)), plot = FALSE) else NULL,
    na.action = attr(des$mf, "na.action")
  )
  class(out) <- "consensus"
  out
}

###############################################################################
### S3 Methods for Consensus Objects
###############################################################################

#' Methods for Consensus Objects
#'
#' @details The \code{plot()} method uses wrapped residuals
#' \eqn{e_i = \mathrm{atan2}(\sin(y_i-\hat\mu_i),
#' \cos(y_i-\hat\mu_i))}. It provides a circular rose and density display,
#' a von Mises probability-integral-transform (PIT) Q-Q plot, and residual
#' plots against fitted direction and observation order. The consensus model has
#' observation-specific fitted concentrations \eqn{L_i}; the PIT is therefore
#' computed with \eqn{\mathrm{VM}(0,L_i)} for each observation. The Q-Q plot is
#' descriptive rather than a formal goodness-of-fit test because the parameters
#' were estimated from the displayed observations.
#' Observations with an undefined fitted direction have no residual. They are
#' omitted from the rose, fitted-direction and PIT panels, with the matching
#' concentrations omitted from the PIT calculation. The index panel preserves
#' their positions as gaps. A note reports their number. Only requested panels
#' are constructed; no diagnostic is available when every residual is undefined.
#'
#' The diagnostic plots require the suggested packages \pkg{ggplot2} and
#' \pkg{ggcircular}; \code{type = "all"} additionally requires
#' \pkg{gridExtra}.
#'
#' @param object An object of class \code{"consensus"}.
#' @param x An object of class \code{"consensus"}.
#' @param type For \code{coef.consensus()}, either \code{"kappa"} or
#'   \code{"beta"}. The latter returns orientation-preserving ratios normalized
#'   by the absolute selected reference coefficient, with the signed reference
#'   chart and stability diagnostics recorded as attributes.
#'   For \code{plot()}, select \code{"all"}, \code{"rose"}, \code{"pit"},
#'   \code{"fitted"}, or \code{"index"}.
#' @param bins Number of bins in a circular residual rose plot.
#' @param reference Reference selection rule for \code{type = "beta"}.
#' @param newdata Optional data frame for prediction.
#' @param robust Logical; return the independent-observation sandwich based on
#'   weighted scores when \code{TRUE}. Weights are case weights, not cluster
#'   identifiers; the sandwich uses squared weights in its middle matrix.
#' @param ... Additional arguments.
#' @rdname consensus-methods
#' @examples
#' data(bison)
#' fit <- consensus(y.dir ~ y.prec, bison[1:100, ])
#' summary(fit)
#' head(predict(fit))
#' vcov(fit)
#' vcov(fit, robust = TRUE)
#' @export
print.consensus <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  .print_fit_status(x)
  label <- if (isTRUE(x$convergence)) "Maximum log-likelihood:" else "Log-likelihood at returned coefficients:"
  cat("\n", label, format(x$MaxLL, digits = 5), "\n", sep = " ")
  cat("AIC:", format(x$AIC, digits = 5), "\n")
  cat("BIC:", format(x$BIC, digits = 5), "\n\n")
  cat("Kappa coefficients:\n")
  stats::printCoefmat(
    x$parameters,
    signif.stars = TRUE,
    signif.legend = TRUE
  )
  invisible(x)
}

#' @rdname consensus-methods
#' @export
coef.consensus <- function(
  object,
  type = c("kappa", "beta"),
  reference = c("auto", "first", "name"),
  ...
) {
  type <- match.arg(type)

  if (type == "kappa") {
    return(object$kappa)
  }

  ref <- .resolve_reference(
    reference = reference,
    y = object$y,
    angle_data = object$angle_data,
    candidates = object$plain_angles,
    tie_method = "first"
  )
  ref_name <- ref$ref_name
  ref_idx <- .reference_term_index(object$term_info, ref_name)

  reference_stability <- .reference_coefficient_stability(
    kappa = object$kappa,
    ref_idx = ref_idx
  )
  kappa_ref <- reference_stability$coefficient
  if (!isTRUE(reference_stability$stable)) {
    stop(
      "Cannot derive beta coefficients: the selected reference coefficient ",
      "does not exceed its scale-relative stability threshold.",
      call. = FALSE
    )
  }

  beta <- object$kappa[-ref_idx] / abs(kappa_ref)
  names(beta) <- object$paramname[-ref_idx]
  attr(beta, "reference") <- ref_name
  attr(beta, "reference_sign") <- reference_stability$reference_sign
  attr(beta, "reference_stability") <- reference_stability
  attr(beta, "reference_scores") <- ref$scores
  beta
}

#' @rdname consensus-methods
#' @export
vcov.consensus <- function(object, robust = FALSE, ...) {
  .check_logical_scalar(robust, "robust")
  if (isTRUE(robust)) {
    if (is.null(object$varcov_robust)) {
      stop("Refit this object with CircularRegression >= 0.7.1 to obtain a robust covariance.", call. = FALSE)
    }
    return(object$varcov_robust)
  }
  object$varcov_kappa
}

#' @rdname consensus-methods
#' @export
fitted.consensus <- function(object, ...) {
  .restore_excluded(object, object$mui)
}

#' @rdname consensus-methods
#' @export
residuals.consensus <- function(object, ...) {
  .restore_excluded(object, .wrap_angle(object$y - object$mui), residual = TRUE)
}

#' @rdname consensus-methods
#' @export
predict.consensus <- function(
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
      coef = as.numeric(object$kappa),
      matx = design$matx,
      matz = design$matz
    )
  }

  if (type == "components") {
    return(cbind(cos = cos(mu), sin = sin(mu)))
  }
  mu
}

#' @rdname consensus-methods
#' @export
summary.consensus <- function(object, ...) {
  s <- list(
    call = object$call,
    MaxLL = object$MaxLL,
    AIC = object$AIC,
    BIC = object$BIC,
    coefficients = object$parameters,
    convergence = object$convergence,
    fit_status = object$fit_status,
    inference_available = object$inference_available,
    residuals = summary(residuals.consensus(object)),
    nobs = object$nobs
  )
  class(s) <- "summary.consensus"
  s
}

#' @rdname consensus-methods
#' @export
print.summary.consensus <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  .print_fit_status(x)
  label <- if (isTRUE(x$convergence)) "Maximum log-likelihood:" else "Log-likelihood at returned coefficients:"
  cat("\n", label, format(x$MaxLL, digits = 5), "\n", sep = " ")
  cat("AIC:", format(x$AIC, digits = 5), "\n")
  cat("BIC:", format(x$BIC, digits = 5), "\n\n")
  cat("Residuals:\n")
  print(x$residuals)
  cat("\nKappa coefficients:\n")
  stats::printCoefmat(
    x$coefficients,
    signif.stars = TRUE,
    signif.legend = TRUE,
    ...
  )
  cat("\nNumber of observations:", x$nobs, "\n")
  invisible(x)
}

#' @rdname consensus-methods
#' @export
plot.consensus <- function(
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
    kappa = x$long,
    type = type,
    bins = bins,
    title = "Consensus angular model diagnostics"
  )
}

#' @rdname consensus-methods
#' @param k Numeric penalty for AIC.
#' @export
AIC.consensus <- function(object, ..., k = 2) {
  NextMethod("AIC")
}

#' @rdname consensus-methods
#' @export
BIC.consensus <- function(object, ...) {
  NextMethod("BIC")
}

#' @rdname consensus-methods
#' @export
logLik.consensus <- function(object, ...) {
  if (isTRUE(object$alignment_boundary)) {
    stop("The consensus likelihood has no finite maximum.", call. = FALSE)
  }
  if (identical(object$convergence, FALSE)) {
    stop("Log-likelihood comparison requires a converged fit.", call. = FALSE)
  }
  val <- object$logLik
  attr(val, "df") <- object$k
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}
