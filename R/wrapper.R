# Internal helper --------------------------------------------------------------
.build_formula <- function(response, rhs_terms) {
  stats::as.formula(
    paste(.btick(response), "~", paste(rhs_terms, collapse = " + ")),
    env = parent.frame()
  )
}

.wrap_to_pi <- function(x) atan2(sin(x), cos(x))

.numeric_jacobian <- function(g, theta, eps = 1e-06) {
  theta <- as.numeric(theta)
  q <- length(as.numeric(g(theta)))
  p <- length(theta)
  J <- matrix(NA_real_, nrow = q, ncol = p)

  if (p == 0) {
    return(J)
  }

  for (j in seq_len(p)) {
    up <- theta
    dn <- theta
    up[j] <- up[j] + eps
    dn[j] <- dn[j] - eps
    J[, j] <- (as.numeric(g(up)) - as.numeric(g(dn))) / (2 * eps)
  }
  J
}

.delta_method_transform <- function(theta_hat, vcov_theta, g, jacobian = NULL, eps = 1e-06) {
  theta_hat <- as.numeric(theta_hat)
  est <- as.numeric(g(theta_hat))
  est_names <- names(g(theta_hat))
  if (is.null(est_names)) {
    est_names <- paste0("param", seq_along(est))
  }
  q <- length(est)
  p <- length(theta_hat)

  if (is.null(jacobian)) {
    J <- .numeric_jacobian(g = g, theta = theta_hat, eps = eps)
  } else {
    J <- as.matrix(jacobian)
  }

  if (!identical(dim(J), c(q, p))) {
    J <- matrix(NA_real_, nrow = q, ncol = p)
  }

  vc <- matrix(NA_real_, nrow = q, ncol = q)
  if (!is.null(vcov_theta)) {
    V <- as.matrix(vcov_theta)
    if (identical(dim(V), c(p, p)) && all(is.finite(V))) {
      vc <- J %*% V %*% t(J)
      vc <- (vc + t(vc)) / 2
    }
  }

  se <- sqrt(pmax(diag(vc), 0))

  names(est) <- est_names
  names(se) <- est_names
  rownames(vc) <- est_names
  colnames(vc) <- est_names
  rownames(J) <- est_names
  colnames(J) <- names(theta_hat)

  list(
    estimate = est,
    vcov = vc,
    se = se,
    jacobian = J
  )
}

.full_coefficients_from_nonref <- function(
  theta,
  ref_idx,
  n_terms,
  reference_sign = 1
) {
  out <- numeric(n_terms)
  out[ref_idx] <- reference_sign
  out[-ref_idx] <- as.numeric(theta)
  out
}

.full_coefficient_jacobian <- function(ref_idx, n_terms) {
  out <- matrix(0, nrow = n_terms, ncol = n_terms - 1L)
  out[-ref_idx, ] <- diag(n_terms - 1L)
  out
}

.natural_parameter_stability <- function(
  theta_hat,
  model_type,
  ref_idx,
  n_terms,
  reference_sign = 1,
  tolerance = sqrt(.Machine$double.eps)
) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("'tolerance' must be a single finite positive number.", call. = FALSE)
  }

  coefficients <- .full_coefficients_from_nonref(
    theta = theta_hat,
    ref_idx = ref_idx,
    n_terms = n_terms,
    reference_sign = reference_sign
  )
  coefficient_scale <- sqrt(sum(coefficients^2))
  threshold <- tolerance * coefficient_scale
  metrics <- switch(
    model_type,
    mean_direction = c(resultant_length = sqrt(sum(coefficients[1:2]^2))),
    decentred = c(
      baseline_resultant_length = sqrt(sum(coefficients[1:2]^2)),
      decentring_resultant_length = sqrt(sum(coefficients[3:4]^2))
    ),
    presnell = c(ratio_denominator = abs(coefficients[[1L]])),
    jammalamadaka_sarma = c(ratio_denominator = abs(coefficients[[1L]])),
    stop("Unknown model type for stability assessment.", call. = FALSE)
  )

  list(
    stable = is.finite(coefficient_scale) && coefficient_scale > 0 &&
      all(is.finite(metrics)) && all(metrics > threshold),
    metrics = metrics,
    coefficient_scale = coefficient_scale,
    threshold = threshold,
    tolerance = tolerance,
    orientation_sign = if (model_type %in% c("presnell", "jammalamadaka_sarma"))
      sign(coefficients[[1L]]) else NA_real_
  )
}

.warn_natural_parameter_singularity <- function(stability, model_type) {
  if (isTRUE(stability$stable)) return(invisible(stability))

  model_label <- switch(
    model_type,
    mean_direction = "Mean-direction",
    decentred = "Decentred",
    presnell = "Presnell",
    jammalamadaka_sarma = "Jammalamadaka-Sarma",
    "Wrapper"
  )
  metric_text <- paste(
    sprintf("%s = %.3g", names(stability$metrics), stability$metrics),
    collapse = ", "
  )
  warning(
    model_label,
    " natural-parameter transformation is numerically singular: ",
    metric_text,
    " does not exceed the relative threshold ",
    sprintf("%.3g", stability$threshold),
    ". Natural parameters and delta-method standard errors may be undefined ",
    "or unstable. Inspect the fitted resultant coefficients and use a ",
    "nonsingular parameterization.",
    call. = FALSE
  )
  invisible(stability)
}

.natural_parameter_mapping <- function(
  model_type,
  ref_idx,
  n_terms,
  reference_sign = 1
) {
  build_full <- function(theta) {
    .full_coefficients_from_nonref(
      theta = theta,
      ref_idx = ref_idx,
      n_terms = n_terms,
      reference_sign = reference_sign
    )
  }
  coefficient_jacobian <- .full_coefficient_jacobian(
    ref_idx = ref_idx,
    n_terms = n_terms
  )

  ratio_map <- function(b, numerator) {
    if (abs(b[1]) <= .Machine$double.eps) {
      return(rep(NA_real_, length(numerator)))
    }
    b[numerator] / b[1]
  }

  ratio_jacobian <- function(b, numerator) {
    J <- matrix(NA_real_, nrow = length(numerator), ncol = length(b))
    if (abs(b[1]) <= .Machine$double.eps) {
      return(J)
    }
    J[] <- 0
    J[, 1] <- -b[numerator] / b[1]^2
    J[cbind(seq_along(numerator), numerator)] <- 1 / b[1]
    J
  }

  mapping <- switch(
    model_type,
    mean_direction = {
      if (n_terms != 2L) {
        stop("Mean direction wrapper expects exactly 2 synthetic terms.", call. = FALSE)
      }
      list(
        g = function(theta) {
          b <- build_full(theta)
          c(mu = .wrap_to_pi(atan2(b[2], b[1])))
        },
        jacobian = function(theta) {
          b <- build_full(theta)
          den <- b[1]^2 + b[2]^2
          J_full <- matrix(NA_real_, nrow = 1L, ncol = n_terms)
          if (den > .Machine$double.eps) {
            J_full[] <- 0
            J_full[1, 1] <- -b[2] / den
            J_full[1, 2] <- b[1] / den
          }
          J_full %*% coefficient_jacobian
        }
      )
    },
    decentred = {
      if (n_terms != 4L) {
        stop("Decentred wrapper expects exactly 4 synthetic terms.", call. = FALSE)
      }
      list(
        g = function(theta) {
          b <- build_full(theta)
          d <- sqrt(b[1]^2 + b[2]^2)
          e <- sqrt(b[3]^2 + b[4]^2)
          if (d <= .Machine$double.eps || e <= .Machine$double.eps) {
            return(c(r_hat = NA_real_, alpha_hat = NA_real_, beta_hat = NA_real_))
          }
          theta_1 <- atan2(b[2], b[1])
          theta_2 <- atan2(b[4], b[3])
          c(
            r_hat = e / d,
            alpha_hat = .wrap_to_pi(theta_2 - theta_1),
            beta_hat = .wrap_to_pi(theta_2)
          )
        },
        jacobian = function(theta) {
          b <- build_full(theta)
          d_sq <- b[1]^2 + b[2]^2
          e_sq <- b[3]^2 + b[4]^2
          d <- sqrt(d_sq)
          e <- sqrt(e_sq)
          J_full <- matrix(NA_real_, nrow = 3L, ncol = n_terms)
          if (d > .Machine$double.eps && e > .Machine$double.eps) {
            J_full[] <- 0
            J_full[1, 1:2] <- -e * b[1:2] / d^3
            J_full[1, 3:4] <- b[3:4] / (e * d)
            J_full[2, 1:2] <- c(b[2], -b[1]) / d_sq
            J_full[2, 3:4] <- c(-b[4], b[3]) / e_sq
            J_full[3, 3:4] <- c(-b[4], b[3]) / e_sq
          }
          J_full %*% coefficient_jacobian
        }
      )
    },
    presnell = {
      if (n_terms != 4L) {
        stop("Presnell wrapper expects exactly 4 synthetic terms.", call. = FALSE)
      }
      list(
        g = function(theta) {
          b <- build_full(theta)
          out <- ratio_map(b, numerator = 2:4)
          stats::setNames(out, c("beta2", "beta3", "beta4"))
        },
        jacobian = function(theta) {
          b <- build_full(theta)
          ratio_jacobian(b, numerator = 2:4) %*% coefficient_jacobian
        }
      )
    },
    jammalamadaka_sarma = {
      if (n_terms != 6L) {
        stop("Jammalamadaka-Sarma wrapper expects exactly 6 synthetic terms.", call. = FALSE)
      }
      list(
        g = function(theta) {
          b <- build_full(theta)
          q <- ratio_map(b, numerator = 2:6)
          if (anyNA(q)) {
            return(c(
              sin_intercept = NA_real_,
              cos_cos = NA_real_,
              cos_sin = NA_real_,
              sin_cos = NA_real_,
              sin_sin = NA_real_
            ))
          }
          c(
            sin_intercept = q[1],
            cos_cos = q[2] + q[3],
            cos_sin = -q[4] + q[5],
            sin_cos = q[4] + q[5],
            sin_sin = q[2] - q[3]
          )
        },
        jacobian = function(theta) {
          b <- build_full(theta)
          L <- rbind(
            c(1, 0, 0, 0, 0),
            c(0, 1, 1, 0, 0),
            c(0, 0, 0, -1, 1),
            c(0, 0, 0, 1, 1),
            c(0, 1, -1, 0, 0)
          )
          L %*% ratio_jacobian(b, numerator = 2:6) %*% coefficient_jacobian
        }
      )
    },
    stop("Unknown model type for natural-parameter computation.", call. = FALSE)
  )

  mapping
}

.run_two_step_with_default_reference <- function(model_formula, data_aug, ...) {
  extra_args <- list(...)
  arg_names <- names(extra_args)
  if (is.null(arg_names) || !("reference" %in% arg_names)) {
    extra_args$reference <- "first"
  }

  do.call(
    angular_two_step,
    c(
      list(formula = model_formula, data = data_aug),
      extra_args
    )
  )
}

.compute_natural_parameters <- function(homogeneous_fit, model_type, eps = 1e-06) {
  if (!inherits(homogeneous_fit, "angular")) {
    stop("'homogeneous_fit' must be an object of class 'angular'.", call. = FALSE)
  }

  theta_hat <- as.numeric(homogeneous_fit$parameters[, "estimate"])
  names(theta_hat) <- rownames(homogeneous_fit$parameters)
  n_terms <- length(homogeneous_fit$term_labels)
  ref_idx <- homogeneous_fit$ref_idx
  reference_sign <- homogeneous_fit$reference_sign %||% 1

  if (length(theta_hat) != (n_terms - 1L)) {
    stop(
      "Inconsistent dimensions in homogeneous fit: expected n_terms - 1 coefficients.",
      call. = FALSE
    )
  }

  mapping <- .natural_parameter_mapping(
    model_type = model_type,
    ref_idx = ref_idx,
    n_terms = n_terms,
    reference_sign = reference_sign
  )
  stability <- .natural_parameter_stability(
    theta_hat = theta_hat,
    model_type = model_type,
    ref_idx = ref_idx,
    n_terms = n_terms,
    reference_sign = reference_sign
  )
  .warn_natural_parameter_singularity(stability, model_type = model_type)

  out_model <- .delta_method_transform(
    theta_hat = theta_hat,
    vcov_theta = homogeneous_fit$varcov0,
    g = mapping$g,
    jacobian = mapping$jacobian(theta_hat),
    eps = eps
  )
  out_robust <- .delta_method_transform(
    theta_hat = theta_hat,
    vcov_theta = homogeneous_fit$varcov1,
    g = mapping$g,
    jacobian = out_model$jacobian,
    eps = eps
  )

  table <- cbind(
    estimate = out_model$estimate,
    se_model = out_model$se,
    se_robust = out_robust$se
  )
  rownames(table) <- names(out_model$estimate)

  list(
    parameters = table,
    vcov_model = out_model$vcov,
    vcov_robust = out_robust$vcov,
    jacobian = out_model$jacobian,
    stability = stability
  )
}

.wrap_two_step_output <- function(pipeline, data_aug, model_formula, model_type) {
  natural <- .compute_natural_parameters(
    homogeneous_fit = pipeline$homogeneous_fit,
    model_type = model_type
  )

  list(
    fit = pipeline$homogeneous_fit,
    consensus_fit = pipeline$consensus_fit,
    homogeneous_fit = pipeline$homogeneous_fit,
    reference = pipeline$reference,
    reference_stability = pipeline$reference_stability,
    data_aug = data_aug,
    formula = model_formula,
    natural_parameters = natural$parameters,
    natural_vcov_model = natural$vcov_model,
    natural_vcov_robust = natural$vcov_robust,
    natural_jacobian = natural$jacobian,
    natural_stability = natural$stability,
    natural_orientation = natural$stability$orientation_sign
  )
}

.validate_wrapper_data <- function(data, response, w = NULL) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!(response %in% names(data))) {
    stop("The response variable provided in 'response' is not present in 'data'.", call. = FALSE)
  }
  .validate_numeric_vector(data[[response]][!is.na(data[[response]])], response)
  if (!is.null(w)) {
    if (!(w %in% names(data))) {
      stop("The 'w' variable provided (", w, ") is not present in 'data'.", call. = FALSE)
    }
    .validate_numeric_vector(data[[w]][!is.na(data[[w]])], w)
  }
  invisible(TRUE)
}

.validate_mu0 <- function(mu0, n) {
  if (is.null(mu0)) {
    return(NULL)
  }
  if (!is.numeric(mu0) || any(!is.finite(mu0)) || !(length(mu0) %in% c(1L, n))) {
    stop("'mu0' must be a finite numeric scalar or have length nrow(data).", call. = FALSE)
  }
  mu0
}

#' Mean Direction Model Wrapper (two-step workflow)
#'
#' Fits the mean-direction special case through the same two-step workflow used
#' by \code{\link{angular_two_step}()}. The wrapper creates synthetic reference
#' directions internally and returns natural-parameter summaries obtained by a
#' delta-method transformation of the homogeneous fit. With a varying
#' \code{mu0}, the returned \code{mu} parameter is a rotation relative to that
#' baseline rather than an unconditional constant mean direction.
#'
#' @param data A data frame containing at least the response variable.
#' @param response Name of the angular response variable.
#' @param mu0 Optional finite numeric scalar or vector of length
#'   \code{nrow(data)} used as a baseline direction.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{reference_stability}, \code{data_aug},
#'   \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}, and
#'   \code{natural_stability}.
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @examples
#' data(bison)
#' d <- na.omit(bison[seq_len(80), "y.dir", drop = FALSE])
#' fit <- meanDirectionModel(d, response = "y.dir")
#' fit$natural_parameters
#' @export
meanDirectionModel <- function(data, response, mu0 = NULL, ...) {
  .validate_wrapper_data(data = data, response = response)
  if (response %in% c("mu0", "tan_mu")) stop("The response name conflicts with a synthetic direction.", call. = FALSE)
  mu0 <- .validate_mu0(mu0, nrow(data))
  data_aug <- data
  if (is.null(mu0)) {
    data_aug$mu0 <- 0
    data_aug$tan_mu <- pi / 2
  } else {
    data_aug$mu0 <- mu0
    data_aug$tan_mu <- mu0 + pi / 2
  }

  model_formula <- .build_formula(response = response, rhs_terms = c("mu0", "tan_mu"))

  pipeline <- .run_two_step_with_default_reference(
    model_formula = model_formula,
    data_aug = data_aug,
    ...
  )
  .wrap_two_step_output(
    pipeline,
    data_aug = data_aug,
    model_formula = model_formula,
    model_type = "mean_direction"
  )
}

#' Decentered Predictor Model Wrapper (two-step workflow)
#'
#' Fits the decentered-predictor special case by augmenting the data with
#' synthetic orthogonal directions and applying \code{\link{angular_two_step}()}.
#' Natural parameters are returned with model-based and robust delta-method
#' standard errors. The two angular natural parameters are undefined when the
#' fitted decentring component has zero length.
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response Name of the angular response variable.
#' @param w Name of the angular explanatory variable.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{reference_stability}, \code{data_aug},
#'   \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}, and
#'   \code{natural_stability}.
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @examples
#' data(bison)
#' d <- na.omit(bison[seq_len(80), c("y.dir", "y.prec")])
#' fit <- decentredPredictorModel(d, response = "y.dir", w = "y.prec")
#' fit$natural_parameters
#' @export
decentredPredictorModel <- function(data, response, w, ...) {
  .validate_wrapper_data(data = data, response = response, w = w)
  if (response %in% c("const0", "constPi2", paste0(w, "_plus_pi2")) ||
      w %in% c("const0", "constPi2")) stop("A column name conflicts with a synthetic direction.", call. = FALSE)

  data_aug <- data
  w_pi2 <- paste0(w, "_plus_pi2")
  data_aug[[w_pi2]] <- data_aug[[w]] + pi / 2
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  rhs <- c(.btick(w), .btick(w_pi2), "const0", "constPi2")
  model_formula <- .build_formula(response = response, rhs_terms = rhs)

  pipeline <- .run_two_step_with_default_reference(
    model_formula = model_formula,
    data_aug = data_aug,
    ...
  )
  .wrap_two_step_output(
    pipeline,
    data_aug = data_aug,
    model_formula = model_formula,
    model_type = "decentred"
  )
}

#' Affine Vector Mean Curve Wrapper (two-step workflow)
#'
#' Fits an affine vector mean curve by building the linear and orthogonal
#' synthetic terms required by the homogeneous angular formulation. The returned
#' natural parameters are ratios computed from the homogeneous coefficients.
#' The curve is related to the mean construction in the projected multivariate
#' linear model of Presnell et al. (1998), but this function fits the homogeneous
#' von Mises angular model rather than a projected-normal likelihood.
#'
#' The natural ratios divide by the first synthetic coefficient. The package
#' also returns \code{natural_orientation}, its sign. Multiply both reconstructed
#' vector components by this sign before applying \code{atan2()}; otherwise a
#' negative denominator reverses the fitted direction by pi radians.
#' The package
#' compares its absolute value with
#' \eqn{\sqrt{\epsilon_{\mathrm{mach}}}\lVert b\rVert_2} and emits a warning
#' when this ratio chart is numerically singular. In that case, the underlying
#' fitted resultant may remain usable, but the returned ratios and their
#' delta-method standard errors are unstable or undefined.
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response Name of the angular response variable.
#' @param w Name of the continuous covariate.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{reference_stability}, \code{data_aug},
#'   \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}, and
#'   \code{natural_stability}.
#' @references Presnell, B., Morrison, S. P., and Littell, R. C. (1998).
#'   Projected multivariate linear models for directional data. \emph{Journal
#'   of the American Statistical Association}, 93(443), 1068-1077.
#'
#'   Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @examples
#' data(bison)
#' d <- na.omit(bison[500:579, c("y.dir", "z.gap")])
#' fit <- presnellModel(d, response = "y.dir", w = "z.gap")
#' fit$natural_parameters
#' @export
presnellModel <- function(data, response, w, ...) {
  .validate_wrapper_data(data = data, response = response, w = w)
  reserved <- c("const0", "constPi2", "affine_direction", "affine_orthogonal", "affine_magnitude")
  if (any(c(response, w) %in% reserved)) stop("A column name conflicts with a synthetic direction.", call. = FALSE)

  data_aug <- data
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  # Encode a signed real covariate as a non-negative length and a direction.
  # This preserves the general resultant-vector parser's modifier contract.
  data_aug$affine_direction <- ifelse(data_aug[[w]] < 0, pi, 0)
  data_aug$affine_orthogonal <- data_aug$affine_direction + pi / 2
  data_aug$affine_magnitude <- abs(data_aug[[w]])
  rhs <- c("const0", "constPi2", "affine_direction:affine_magnitude",
           "affine_orthogonal:affine_magnitude")
  model_formula <- .build_formula(response = response, rhs_terms = rhs)

  pipeline <- .run_two_step_with_default_reference(
    model_formula = model_formula,
    data_aug = data_aug,
    ...
  )
  .wrap_two_step_output(
    pipeline,
    data_aug = data_aug,
    model_formula = model_formula,
    model_type = "presnell"
  )
}

#' Jammalamadaka-Sarma First-Harmonic Wrapper (two-step workflow)
#'
#' Fits the first-harmonic circular-circular mean-direction curve associated
#' with Sarma and Jammalamadaka (1993) by creating the required sine/cosine
#' angular terms. The returned natural parameters are the coefficients of the
#' cosine and sine components of that curve, with delta-method uncertainty.
#' Multiply both components by \code{natural_orientation} before applying
#' \code{atan2()} to preserve the orientation when the ratio denominator is negative.
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response Name of the angular response variable.
#' @param w Name of the angular explanatory variable.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{reference_stability}, \code{data_aug},
#'   \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}, and
#'   \code{natural_stability}.
#' @references Sarma, Y. R., and Jammalamadaka, S. R. (1993). Circular
#'   regression. In \emph{Statistical Sciences and Data Analysis}, 109-128.
#'
#'   Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#' @examples
#' data(bison)
#' d <- na.omit(bison[seq_len(80), c("y.dir", "y.prec")])
#' fit <- jammalamadakaModel(d, response = "y.dir", w = "y.prec")
#' fit$natural_parameters
#' @export
jammalamadakaModel <- function(data, response, w, ...) {
  .validate_wrapper_data(data = data, response = response, w = w)
  reserved <- c("const0", "constPi2", paste0(w, "_orig"), paste0("neg_", w),
                paste0(w, "_pi2"), paste0("neg_", w, "_pi2"))
  if (any(c(response, w) %in% reserved)) stop("A column name conflicts with a synthetic direction.", call. = FALSE)

  data_aug <- data
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  w_orig <- paste0(w, "_orig")
  w_neg <- paste0("neg_", w)
  w_pi2 <- paste0(w, "_pi2")
  w_neg_pi2 <- paste0("neg_", w, "_pi2")

  data_aug[[w_orig]] <- data_aug[[w]]
  data_aug[[w_neg]] <- -data_aug[[w]]
  data_aug[[w_pi2]] <- data_aug[[w]] + pi / 2
  data_aug[[w_neg_pi2]] <- -data_aug[[w]] + pi / 2

  rhs <- c("const0", "constPi2", .btick(c(w_orig, w_neg, w_pi2, w_neg_pi2)))
  model_formula <- .build_formula(response = response, rhs_terms = rhs)

  pipeline <- .run_two_step_with_default_reference(
    model_formula = model_formula,
    data_aug = data_aug,
    ...
  )
  .wrap_two_step_output(
    pipeline,
    data_aug = data_aug,
    model_formula = model_formula,
    model_type = "jammalamadaka_sarma"
  )
}
