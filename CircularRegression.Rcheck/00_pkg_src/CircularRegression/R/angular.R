###############################################################################
# Angular Regression Model (Homogeneous errors)
###############################################################################

#' Angular Regression Model with Homogeneous von Mises Errors
#'
#' Fits the homogeneous angular regression model of Rivest et al. (2016), where
#' the mean direction is the orientation of a weighted resultant vector and one
#' reference angle has coefficient fixed to 1.
#'
#' @param formula A formula with terms of the form \code{x} or \code{x:z}.
#' @param data A data frame, list or environment containing model variables.
#' @param reference Reference-angle rule. Use \code{"auto"} (default, article
#'   criterion), \code{"first"}, or \code{c("name", "x_ref")}. A plain
#'   angle term \code{x_ref} must exist in the formula.
#' @param initbeta Optional numeric vector of initial values for non-reference
#'   coefficients.
#' @param control Control list with optional components \code{pginit},
#'   \code{maxiter}, and \code{mindiff}.
#' @return An object of class \code{"angular"}.
#' @export
angular <- function(
  formula,
  data,
  reference = c("auto", "first", "name"),
  initbeta = NULL,
  control = list()
) {
  call <- match.call()
  des <- .angular_design(formula = formula, data = data, reference = reference)

  y <- des$y
  x0 <- des$x0
  matx <- des$matx
  matz <- des$matz
  betaname <- des$betaname
  p <- length(betaname)
  nobs <- des$nobs

  if (!is.null(initbeta) && length(initbeta) != p) {
    stop(sprintf("'initbeta' must have length %d for this model.", p), call. = FALSE)
  }

  pginit <- control$pginit %||% 1000
  maxiter <- control$maxiter %||% 1000
  mindiff <- control$mindiff %||% 1e-06

  if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter < 0) {
    stop("'control$maxiter' must be a non-negative scalar.", call. = FALSE)
  }
  maxiter <- as.integer(maxiter)

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
    val <- mean(cos(y - cm$mui))
    list(maxcos = val, long = cm$long, mui = cm$mui)
  }

  beta_update <- function(betak, long, mui) {
    denom <- pmax(as.vector(long), .Machine$double.eps)
    matd <- (matz * sin(matx - mui)) / denom
    res <- sin(y - mui)

    qr_obj <- qr(matd)
    if (qr_obj$rank < ncol(matd)) {
      stop(
        "Design matrix is rank deficient; homogeneous model is not identifiable for this reference.",
        call. = FALSE
      )
    }

    dbeta <- as.vector(qr.coef(qr_obj, res))
    list(betak1 = betak + dbeta, dbeta = dbeta, matd = matd, res = res)
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
  maxcosk <- calc$maxcos
  long <- calc$long
  mui <- calc$mui

  iter <- 0L
  conv <- FALSE
  iter.detail <- matrix(NA_real_, nrow = maxiter + 1L, ncol = p + 3L)
  colnames(iter.detail) <- c(betaname, "maxcosine", "iter", "nitersh")
  iter.detail[1L, ] <- c(betak, maxcosk, iter, 0)

  if (p > 0 && maxiter > 0) {
    while (!conv && iter < maxiter) {
      step <- beta_update(betak = betak, long = long, mui = mui)
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

      conv <- (maxcosk1 - maxcosk) <= mindiff
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
  beta_full[des$ref_idx] <- 1
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
  } else {
    matd <- matrix(0, nrow = nobs, ncol = 0)
    res <- sin(y - mui)
  }

  Akappa <- as.numeric(maxcosk)
  kappahat <- if (Akappa <= 0) 0 else circular::A1inv(min(Akappa, 0.999999))

  if (p > 0) {
    XtX <- crossprod(matd)
    XtX_sym <- (XtX + t(XtX)) / 2
    chol_xtx <- tryCatch(chol(XtX_sym), error = function(e) NULL)
    if (is.null(chol_xtx)) {
      stop("Information matrix is singular; variance estimation failed.", call. = FALSE)
    }
    XtX_inv <- chol2inv(chol_xtx)

    if (Akappa > 0 && kappahat > 0) {
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

    res_sq <- res^2
    matd_res <- matd * res_sq
    v2 <- XtX_inv %*% crossprod(matd, matd_res) %*% XtX_inv

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
    nobs = nobs,
    k = p + 1L,
    reference = des$reference$ref_name,
    reference_scores = des$reference$scores,
    term_info = des$term_info,
    term_labels = des$term_labels,
    ref_idx = des$ref_idx,
    full_coefficients = beta_full,
    response_name = des$response_name
  )
  class(out) <- "angular"
  out
}

###############################################################################
### S3 Methods for Angular Regression Objects
###############################################################################

#' Methods for Angular Regression Objects
#'
#' @param object An object of class \code{"angular"}.
#' @param x An object of class \code{"angular"}.
#' @param ... Additional arguments.
#' @rdname angular-methods
#' @export
print.angular <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nReference angle:", x$reference, "\n")
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
residuals.angular <- function(object, ...) {
  atan2(sin(object$y - object$mui), cos(object$y - object$mui))
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
    reference = object$reference
  )
  class(s) <- "summary.angular"
  s
}

#' @rdname angular-methods
#' @export
print.summary.angular <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nReference angle:", x$reference, "\n")
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
plot.angular <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for arranging plots.", call. = FALSE)
  }

  fitted_vals <- as.numeric(x$mui)
  residuals_num <- as.numeric(residuals.angular(x))

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
      linewidth = 1
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
    top = sprintf("Diagnostic Plots for Angular Model (ref: %s)", x$reference)
  )

  invisible(x)
}
