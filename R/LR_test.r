# ------------------------------------------------------------------------------
# Likelihood-ratio-style test for homogeneous angular models
# ------------------------------------------------------------------------------

#' Likelihood-ratio test (T2) between two homogeneous angular models
#'
#' Computes the statistic from Rivest et al. (2016), Section 3.6:
#' \deqn{T_2 = 2 n \hat\kappa_1 (\bar C_1 - \bar C_0),}
#' where model 1 is the full model and model 0 the nested reduced model.
#'
#' @param full Object of class \code{"angular"} (full model).
#' @param reduced Object of class \code{"angular"} (nested reduced model).
#' @param check Logical; check that samples, shared variables, term nesting,
#'   reference directions and signs agree.
#'   Setting this to \code{FALSE} bypasses these structural checks only and
#'   places responsibility for the test assumptions on the caller.
#'   Both fits must have converged with finite positive concentrations.
#' @return An object of class \code{"htest"}.
#' @details The chi-squared calibration assumes regular, identifiable nested
#'   homogeneous models with independent errors and finite positive concentration.
#'   This is not a serially robust test or a test at singular independence submodels.
#' @examples
#' data(bison)
#' d <- bison[1:150, ]
#' full <- angular(y.dir ~ y.prec + x.meadow:z.meadow, d, reference = "first")
#' reduced <- angular(y.dir ~ y.prec, d, reference = "first")
#' angular_lrtest(full, reduced)
#' @export
angular_lrtest <- function(full, reduced, check = TRUE) {
  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("'check' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!inherits(full, "angular") || !inherits(reduced, "angular")) {
    stop("Both 'full' and 'reduced' must be objects of class 'angular'.", call. = FALSE)
  }

  nm_full <- deparse(substitute(full))
  nm_reduced <- deparse(substitute(reduced))

  n_f <- length(full$y)
  n_r <- length(reduced$y)

  if (check) {
    if (n_f != n_r) {
      stop("Models have different sample sizes.", call. = FALSE)
    }
    if (!identical(full$y, reduced$y) || is.null(full$model_data) ||
        is.null(reduced$model_data) ||
        !identical(rownames(full$model_data), rownames(reduced$model_data))) {
      stop("Models must use identical observations in the same order; refit older objects.", call. = FALSE)
    }
    if (!all(reduced$term_labels %in% full$term_labels)) {
      stop("The reduced model must be nested in the full model (a subset of terms).", call. = FALSE)
    }
    if (is.null(full$term_info) || is.null(reduced$term_info)) {
      stop("Refit older objects before checking directional term nesting.", call. = FALSE)
    }
    term_index <- match(reduced$term_labels, full$term_labels)
    roles <- c("angle", "modifier", "has_modifier")
    if (!isTRUE(all.equal(
      unname(as.list(full$term_info[term_index, roles])),
      unname(as.list(reduced$term_info[, roles])),
      check.attributes = FALSE
    ))) {
      stop("Models are not nested: shared terms have different direction/modifier roles.", call. = FALSE)
    }
    shared <- names(reduced$model_data)
    if (!all(shared %in% names(full$model_data)) ||
        !isTRUE(all.equal(as.list(full$model_data[shared]), as.list(reduced$model_data),
                          check.attributes = FALSE))) {
      stop("Models must use identical values for all shared variables.", call. = FALSE)
    }
    if (!identical(full$reference, reduced$reference) ||
        !identical(full$reference_sign, reduced$reference_sign)) {
      stop("Models must use compatible reference directions and signs.", call. = FALSE)
    }
  }
  if (!isTRUE(full$convergence) || !isTRUE(reduced$convergence)) {
    stop("Both models must have converged before a T2 comparison.", call. = FALSE)
  }
  if (identical(full$fit_status, "singular_information") ||
      identical(reduced$fit_status, "singular_information")) {
    stop("T2 requires regular fitted models with nonsingular information.", call. = FALSE)
  }

  n <- n_f
  Cbar_full <- as.numeric(full$MaxCosine)
  Cbar_red <- as.numeric(reduced$MaxCosine)
  kappa_full <- as.numeric(full$kappahat)
  kappa_red <- as.numeric(reduced$kappahat)
  if (!is.finite(kappa_full) || kappa_full <= 0 || !is.finite(Cbar_full) ||
      !is.finite(kappa_red) || kappa_red <= 0 || !is.finite(Cbar_red)) {
    stop("T2 requires a finite positive concentration and regular fitted models.", call. = FALSE)
  }
  if (Cbar_full < Cbar_red - 1e-10) {
    stop("The full model fits worse than the reduced model; check optimization.", call. = FALSE)
  }

  T2 <- 2 * n * kappa_full * max(Cbar_full - Cbar_red, 0)

  df_full <- nrow(full$parameters)
  df_red <- nrow(reduced$parameters)
  df <- df_full - df_red
  if (df <= 0) {
    stop("'full' must have strictly more beta parameters than 'reduced'.", call. = FALSE)
  }

  pval <- stats::pchisq(T2, df = df, lower.tail = FALSE)

  LRT_ll <- tryCatch(
    {
      2 * (as.numeric(stats::logLik(full)) - as.numeric(stats::logLik(reduced)))
    },
    error = function(e) NA_real_
  )

  out <- list(
    statistic = c(T2 = unname(T2)),
    parameter = c(df = unname(df)),
    p.value = unname(pval),
    method = "T2 test (homogeneous angular model; Rivest et al. 2016)",
    alternative = "full model improves fit over reduced",
    data.name = sprintf("%s vs %s", nm_full, nm_reduced),
    details = list(
      n = n,
      kappa_full = kappa_full,
      Cbar_full = Cbar_full,
      Cbar_reduced = Cbar_red,
      T2 = T2,
      LRT_logLik = LRT_ll
    )
  )
  class(out) <- "htest"
  out
}

#' @rdname angular-methods
#' @param test For \code{anova()}, either \code{"none"} for a model table or
#'   \code{"T2"} for the nested homogeneous comparison.
#' @export
anova.angular <- function(object, ..., test = c("none", "T2")) {
  test <- match.arg(test)
  others <- list(...)

  if (length(others) == 0) {
    tab <- data.frame(
      Df = nrow(object$parameters),
      MaxCosine = object$MaxCosine,
      row.names = deparse(substitute(object))
    )
    class(tab) <- c("anova", "data.frame")
    return(tab)
  }

  if (length(others) > 1) {
    stop("Provide exactly two models: anova(full, reduced, test = 'T2').", call. = FALSE)
  }

  reduced <- others[[1L]]
  if (test == "none") {
    tab <- data.frame(
      Df = c(nrow(object$parameters), nrow(reduced$parameters)),
      MaxCosine = c(object$MaxCosine, reduced$MaxCosine),
      row.names = c(deparse(substitute(object)), deparse(substitute(reduced)))
    )
    class(tab) <- c("anova", "data.frame")
    return(tab)
  }

  angular_lrtest(full = object, reduced = reduced)
}

#' @rdname angular-methods
#' @export
logLik.angular <- function(object, ...) {
  n <- length(object$y)
  kappa <- as.numeric(object$kappahat)
  Cbar <- as.numeric(object$MaxCosine)
  if (!is.finite(kappa)) {
    stop("The homogeneous likelihood has no finite concentration maximum.", call. = FALSE)
  }
  if (identical(object$convergence, FALSE)) {
    stop("Log-likelihood comparison requires a converged regular fit.", call. = FALSE)
  }

  log_i0 <- .logI0(kappa)
  ll <- n * (kappa * Cbar - log(2 * pi) - log_i0)

  attr(ll, "df") <- nrow(object$parameters) + 1L
  attr(ll, "nobs") <- n
  class(ll) <- "logLik"
  ll
}
