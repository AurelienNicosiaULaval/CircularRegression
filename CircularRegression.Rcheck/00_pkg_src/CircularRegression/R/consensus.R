###############################################################################
# Consensus Model
###############################################################################

.A1_ratio <- function(k) {
  num <- besselI(k, 1, expon.scaled = TRUE)
  den <- besselI(k, 0, expon.scaled = TRUE)
  out <- num / den
  out[!is.finite(out)] <- 1
  out[k == 0] <- 0
  out
}

.logI0 <- function(k) {
  log(besselI(k, 0, expon.scaled = TRUE)) + k
}

#' Consensus Angular Regression Model
#'
#' Fits the consensus model of Rivest et al. (2016), where each formula term of
#' the form \code{x} or \code{x:z} receives a \eqn{\kappa_j} coefficient in the
#' log-likelihood parameterization.
#'
#' @param formula Model formula containing terms of the form \code{x} or
#'   \code{x:z}.
#' @param data Data frame, list or environment containing model variables.
#' @param weights Optional non-negative weights.
#' @param initkappa Optional numeric vector of initial \eqn{\kappa_j} values.
#' @param control Optional list with \code{pginit}, \code{maxiter},
#'   \code{mindiff}.
#' @return An object of class \code{"consensus"}.
#' @export
consensus <- function(
  formula,
  data,
  weights = NULL,
  initkappa = NULL,
  control = list()
) {
  call <- match.call()
  des <- .consensus_design(formula = formula, data = data)

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
    w_raw <- as.numeric(weights)
    if (length(w_raw) == nobs) {
      weight <- w_raw
    } else if (is.data.frame(data) && length(w_raw) == nrow(data)) {
      keep_idx <- as.integer(rownames(des$mf))
      weight <- w_raw[keep_idx]
    } else {
      stop("'weights' must have length equal to nrow(data) or to the model-frame size.", call. = FALSE)
    }
  }

  if (any(!is.finite(weight)) || any(weight < 0)) {
    stop("'weights' must be finite and non-negative.", call. = FALSE)
  }
  if (all(weight == 0)) {
    stop("'weights' cannot be all zero.", call. = FALSE)
  }

  if (!is.null(initkappa) && length(initkappa) != nparam) {
    stop(sprintf("'initkappa' must have length %d.", nparam), call. = FALSE)
  }

  pginit <- control$pginit %||% 1000
  maxiter <- control$maxiter %||% 1000
  mindiff <- control$mindiff %||% 1e-06

  if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter < 0) {
    stop("'control$maxiter' must be a non-negative scalar.", call. = FALSE)
  }
  maxiter <- as.integer(maxiter)

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
    ll <- sum(weight * term1) - sum(weight * .logI0(comp$long))
    list(LL = ll, long = comp$long, mui = comp$mui)
  }

  fisher_step <- function(kappa, long, mui) {
    A <- .A1_ratio(long)
    A_over_long <- ifelse(long > .Machine$double.eps, A / long, 0)

    S <- matz * sin(matx - mui)
    C <- matz * cos(matx - mui)

    matu <- matz * (cos_y_minus_x - C * A)
    score <- colSums(weight * matu)

    w1 <- pmax(weight * A_over_long, 0)
    w2 <- pmax(weight * (1 - A_over_long - A^2), 0)

    I1 <- crossprod(S * sqrt(w1), S * sqrt(w1))
    I2 <- crossprod(C * sqrt(w2), C * sqrt(w2))
    I <- (I1 + I2)
    I <- (I + t(I)) / 2

    qrI <- qr(I)
    if (qrI$rank < ncol(I)) {
      stop(
        "Consensus information matrix is singular; parameters are not identifiable.",
        call. = FALSE
      )
    }

    dparam <- as.vector(qr.coef(qrI, score))
    list(kappa1 = kappa + dparam, dparam = dparam, I = I)
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
    kappak <- as.numeric(poss[which.max(poss[, nparam + 1]), seq_len(nparam), drop = TRUE])
  } else {
    kappak <- as.numeric(initkappa)
  }

  calc <- loglik_components(kappak)
  maxLLk <- calc$LL
  long <- calc$long
  mui <- calc$mui

  iter <- 0L
  conv <- FALSE
  iter.detail <- matrix(NA_real_, nrow = maxiter + 1L, ncol = nparam + 3L)
  colnames(iter.detail) <- c(paramname, "LL", "iter", "nitersh")
  iter.detail[1L, ] <- c(kappak, maxLLk, iter, 0)

  if (maxiter > 0) {
    while (!conv && iter < maxiter) {
      step <- fisher_step(kappak, long = long, mui = mui)
      kappak1 <- step$kappa1
      dparam <- step$dparam

      calc1 <- loglik_components(kappak1)
      maxLLk1 <- calc1$LL
      long1 <- calc1$long
      mui1 <- calc1$mui

      iter.sh <- 0L
      while (maxLLk1 < maxLLk && iter.sh < maxiter) {
        iter.sh <- iter.sh + 1L
        kappak1 <- kappak + dparam / (2^iter.sh)
        calc1 <- loglik_components(kappak1)
        maxLLk1 <- calc1$LL
        long1 <- calc1$long
        mui1 <- calc1$mui
      }

      if (maxLLk1 < maxLLk) {
        warning(
          "The algorithm did not converge; step-halving failed to improve the log-likelihood.",
          call. = FALSE
        )
        break
      }

      conv <- (maxLLk1 - maxLLk) <= mindiff
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
  matI <- final_step$I

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
    parambeta = NULL,
    varcovbeta1 = NULL,
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
    logLik = logLik,
    term_info = term_info,
    paramname = paramname,
    angle_data = des$angle_data,
    plain_angles = plain_angles,
    reference_scores = reference_scores,
    autocorr = stats::acf(residuals.consensus(list(y = y, mui = mui)), plot = FALSE)
  )
  class(out) <- "consensus"
  out
}

###############################################################################
### S3 Methods for Consensus Objects
###############################################################################

#' Methods for Consensus Objects
#'
#' @param object An object of class \code{"consensus"}.
#' @param x An object of class \code{"consensus"}.
#' @param type For \code{coef.consensus()}, either \code{"kappa"} or
#'   \code{"beta"}. The latter returns ratios relative to a selected reference.
#' @param reference Reference selection rule for \code{type = "beta"}.
#' @param ... Additional arguments.
#' @rdname consensus-methods
#' @export
print.consensus <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMaximum log-likelihood:", format(x$MaxLL, digits = 5), "\n")
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

  kappa_ref <- object$kappa[ref_idx]
  if (!is.finite(kappa_ref) || abs(kappa_ref) < .Machine$double.eps) {
    stop("Cannot derive beta coefficients: selected reference kappa is zero.", call. = FALSE)
  }

  beta <- object$kappa[-ref_idx] / kappa_ref
  names(beta) <- object$paramname[-ref_idx]
  attr(beta, "reference") <- ref_name
  attr(beta, "reference_scores") <- ref$scores
  beta
}

#' @rdname consensus-methods
#' @export
residuals.consensus <- function(object, ...) {
  atan2(sin(object$y - object$mui), cos(object$y - object$mui))
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
  cat("\nMaximum log-likelihood:", format(x$MaxLL, digits = 5), "\n")
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
plot.consensus <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for arranging plots.", call. = FALSE)
  }

  res <- as.numeric(residuals.consensus(x, ...))

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
    ) +
    ggplot2::theme_minimal()

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
      linewidth = 1
    ) +
    ggplot2::labs(
      title = "Histogram of Residuals",
      x = "Residuals",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(
    data = data.frame(Residual = res),
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
        "Spearman Correlation:\n|residual| vs 1/long =",
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

#' @rdname consensus-methods
#' @param k Numeric penalty for AIC.
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
#' @export
logLik.consensus <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- object$k
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}
