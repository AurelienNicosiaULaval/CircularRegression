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

.full_coefficients_from_nonref <- function(theta, ref_idx, n_terms) {
  out <- numeric(n_terms)
  out[ref_idx] <- 1
  out[-ref_idx] <- as.numeric(theta)
  out
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

  if (length(theta_hat) != (n_terms - 1L)) {
    stop(
      "Inconsistent dimensions in homogeneous fit: expected n_terms - 1 coefficients.",
      call. = FALSE
    )
  }

  build_full <- function(theta) .full_coefficients_from_nonref(
    theta = theta,
    ref_idx = ref_idx,
    n_terms = n_terms
  )

  g <- switch(
    model_type,
    mean_direction = {
      if (n_terms != 2L) {
        stop("Mean direction wrapper expects exactly 2 synthetic terms.", call. = FALSE)
      }
      function(theta) {
        b <- build_full(theta)
        c(mu = .wrap_to_pi(atan2(b[2], b[1])))
      }
    },
    decentred = {
      if (n_terms != 4L) {
        stop("Decentred wrapper expects exactly 4 synthetic terms.", call. = FALSE)
      }
      function(theta) {
        b <- build_full(theta)
        den <- sqrt(b[1]^2 + b[2]^2)
        r_hat <- if (den <= .Machine$double.eps) NA_real_ else sqrt(b[3]^2 + b[4]^2) / den
        theta1 <- atan2(b[2], b[1])
        theta2 <- atan2(b[4], b[3])
        c(
          r_hat = r_hat,
          alpha_hat = .wrap_to_pi(theta2 - theta1),
          beta_hat = .wrap_to_pi(theta2)
        )
      }
    },
    presnell = {
      if (n_terms != 4L) {
        stop("Presnell wrapper expects exactly 4 synthetic terms.", call. = FALSE)
      }
      function(theta) {
        b <- build_full(theta)
        if (abs(b[1]) <= .Machine$double.eps) {
          return(c(beta2 = NA_real_, beta3 = NA_real_, beta4 = NA_real_))
        }
        c(
          beta2 = b[2] / b[1],
          beta3 = b[3] / b[1],
          beta4 = b[4] / b[1]
        )
      }
    },
    jammalamadaka = {
      if (n_terms != 6L) {
        stop("Jammalamadaka wrapper expects exactly 6 synthetic terms.", call. = FALSE)
      }
      function(theta) {
        b <- build_full(theta)
        if (abs(b[1]) <= .Machine$double.eps) {
          return(
            c(
              gamma2 = NA_real_,
              gamma3 = NA_real_,
              gamma4 = NA_real_,
              gamma5 = NA_real_,
              gamma6 = NA_real_
            )
          )
        }
        b2 <- b[2] / b[1]
        b3 <- b[3] / b[1]
        b4 <- b[4] / b[1]
        b5 <- b[5] / b[1]
        b6 <- b[6] / b[1]
        c(
          gamma2 = b2,
          gamma3 = b3 - b5,
          gamma4 = b4 + b6,
          gamma5 = b3 + b5,
          gamma6 = -b4 + b6
        )
      }
    },
    stop("Unknown model type for natural-parameter computation.", call. = FALSE)
  )

  out_model <- .delta_method_transform(
    theta_hat = theta_hat,
    vcov_theta = homogeneous_fit$varcov0,
    g = g,
    eps = eps
  )
  out_robust <- .delta_method_transform(
    theta_hat = theta_hat,
    vcov_theta = homogeneous_fit$varcov1,
    g = g,
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
    jacobian = out_model$jacobian
  )
}

.wrap_two_step_output <- function(pipeline, data_aug, model_formula, model_type) {
  natural <- .compute_natural_parameters(
    homogeneous_fit = pipeline$homogeneous_fit,
    model_type = model_type
  )

  list(
    consensus_fit = pipeline$consensus_fit,
    homogeneous_fit = pipeline$homogeneous_fit,
    reference = pipeline$reference,
    data_aug = data_aug,
    formula = model_formula,
    natural_parameters = natural$parameters,
    natural_vcov_model = natural$vcov_model,
    natural_vcov_robust = natural$vcov_robust,
    natural_jacobian = natural$jacobian
  )
}

#' Mean Direction Model Wrapper (two-step workflow)
#'
#' @param data A data frame containing at least the response variable.
#' @param response Name of the angular response variable.
#' @param mu0 Optional numeric scalar/vector.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{data_aug}, \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}.
#' @export
meanDirectionModel <- function(data, response, mu0 = NULL, ...) {
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
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response Name of the angular response variable.
#' @param w Name of the angular explanatory variable.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{data_aug}, \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}.
#' @export
decentredPredictorModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The 'w' variable provided (", w, ") is not present in 'data'.", call. = FALSE)
  }

  data_aug <- data
  w_pi2 <- paste0(w, "_plus_pi2")
  data_aug[[w_pi2]] <- data_aug[[w]] + pi / 2
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  rhs <- c(w, w_pi2, "const0", "constPi2")
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

#' Presnell Model Wrapper (two-step workflow)
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response Name of the angular response variable.
#' @param w Name of the continuous covariate.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{data_aug}, \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}.
#' @export
presnellModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The 'w' variable provided (", w, ") is not present in 'data'.", call. = FALSE)
  }

  data_aug <- data
  data_aug$const0 <- 0
  data_aug$constPi2 <- pi / 2

  rhs <- c("const0", "constPi2", paste0("const0:", w), paste0("constPi2:", w))
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

#' Jammalamadaka-Sengupta Model Wrapper (two-step workflow)
#'
#' @param data A data frame containing \code{response} and \code{w}.
#' @param response Name of the angular response variable.
#' @param w Name of the angular explanatory variable.
#' @param ... Additional arguments passed to \code{angular_two_step()}.
#' @return A list with \code{consensus_fit}, \code{homogeneous_fit},
#'   \code{reference}, \code{data_aug}, \code{formula}, and natural-parameter
#'   outputs \code{natural_parameters}, \code{natural_vcov_model},
#'   \code{natural_vcov_robust}, \code{natural_jacobian}.
#' @export
jammalamadakaModel <- function(data, response, w, ...) {
  if (!(w %in% names(data))) {
    stop("The 'w' variable provided (", w, ") is not present in 'data'.", call. = FALSE)
  }

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

  rhs <- c("const0", "constPi2", w_orig, w_neg, w_pi2, w_neg_pi2)
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
    model_type = "jammalamadaka"
  )
}
