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
#' @param check Logical; perform consistency checks.
#' @return An object of class \code{"htest"}.
#' @export
angular_lrtest <- function(full, reduced, check = TRUE) {
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
    if (!isTRUE(all.equal(full$y, reduced$y))) {
      warning(
        "Response vectors differ; assuming models were fitted on comparable data.",
        call. = FALSE
      )
    }
  }

  n <- n_f
  Cbar_full <- as.numeric(full$MaxCosine)
  Cbar_red <- as.numeric(reduced$MaxCosine)
  kappa_full <- as.numeric(full$kappahat)

  T2 <- 2 * n * kappa_full * (Cbar_full - Cbar_red)

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

#' @export
logLik.angular <- function(object, ...) {
  n <- length(object$y)
  kappa <- as.numeric(object$kappahat)
  Cbar <- as.numeric(object$MaxCosine)

  log_i0 <- log(besselI(kappa, nu = 0, expon.scaled = TRUE)) + kappa
  ll <- n * (kappa * Cbar - log(2 * pi) - log_i0)

  attr(ll, "df") <- nrow(object$parameters) + 1L
  attr(ll, "nobs") <- n
  class(ll) <- "logLik"
  ll
}
