# ------------------------------------------------------------------------------
# Likelihood-ratio-style test for homogeneous angular models
# T2 = 2 * n * kappa_full * (Cbar_full - Cbar_reduced),
# where Cbar = mean cos(y - mu_hat). In your fit this is 'MaxCosine'.
# Returns an 'htest' object like base R tests (t.test, etc.).
# ------------------------------------------------------------------------------

#' Likelihood-ratio test between two homogeneous angular models
#'
#' @param full    An object of class "angular" (the larger model).
#' @param reduced An object of class "angular" (the nested model).
#' @param check   Logical; run basic sanity checks (default TRUE).
#' @return An object of class \code{"htest"} with fields:
#' \itemize{
#'   \item \code{statistic} The test statistic \(T_2\).
#'   \item \code{parameter} The degrees of freedom (difference in number of \eqn{\beta}).
#'   \item \code{p.value}   The \eqn{\chi^2} p-value.
#'   \item \code{method}    A short description.
#'   \item \code{data.name} Model names compared.
#'   \item \code{details}   List with \code{n}, \code{kappa_full}, \code{Cbar_full}, \code{Cbar_reduced},
#'                          and \code{LRT_logLik} (= 2*Delta_logLik, if \code{logLik.angular} is available).
#' }
#' @export
angular_lrtest <- function(full, reduced, check = TRUE) {
  if (!inherits(full, "angular") || !inherits(reduced, "angular"))
    stop("Both 'full' and 'reduced' must be objects of class 'angular'.")

  nm_full <- deparse(substitute(full))
  nm_reduced <- deparse(substitute(reduced))

  # Pull required components (all exist in your fit object)
  y_f <- full$y
  y_r <- reduced$y
  n_f <- length(y_f)
  n_r <- length(y_r)
  if (check) {
    if (n_f != n_r) stop("Models have different sample sizes.")
    # Same data in same order is strongly recommended
    if (!isTRUE(all.equal(y_f, y_r))) {
      warning(
        "Response vectors differ; assuming models use the same data/order."
      )
    }
  }
  n <- n_f

  Cbar_full <- as.numeric(full$MaxCosine)
  Cbar_red <- as.numeric(reduced$MaxCosine)
  kappa_full <- as.numeric(full$kappahat)

  # T2 per homogeneous-error theory
  T2 <- 2 * n * kappa_full * (Cbar_full - Cbar_red)

  # df = difference in number of beta parameters (rows of 'parameters')
  df_full <- nrow(full$parameters)
  df_red <- nrow(reduced$parameters)
  df <- df_full - df_red
  if (df <= 0)
    stop("'full' must have strictly more beta parameters than 'reduced'.")

  pval <- stats::pchisq(T2, df = df, lower.tail = FALSE)

  # If logLik.angular is available, also report the plain 2*Delta_logLik
  LRT_ll <- tryCatch(
    {
      2 * (as.numeric(stats::logLik(full)) - as.numeric(stats::logLik(reduced)))
    },
    error = function(e) NA_real_
  )

  out <- list(
    statistic = unname(T2),
    parameter = unname(df),
    p.value = unname(pval),
    method = "LRT (homogeneous von Mises errors) via mean residual cosines",
    alternative = "full model improves fit over reduced",
    data.name = sprintf("%s vs %s", nm_full, nm_reduced),
    details = list(
      n = n,
      kappa_full = kappa_full,
      Cbar_full = Cbar_full,
      Cbar_reduced = Cbar_red,
      LRT_logLik = LRT_ll
    )
  )
  class(out) <- "htest"
  out
}

# ------------------------------------------------------------------------------
# anova() convenience, mirroring lm/glm usage: anova(full, reduced, test="LRT")
# ------------------------------------------------------------------------------

#' @export
anova.angular <- function(object, ..., test = c("none", "LRT")) {
  test <- match.arg(test)
  others <- list(...)
  if (length(others) == 0 || test != "LRT") {
    # no full anova table implemented; fall back to default behavior
    return(NextMethod())
  }
  if (length(others) > 1)
    stop("Provide exactly two models: anova(full, reduced, test = 'LRT').")
  reduced <- others[[1L]]
  angular_lrtest(object, reduced)
}

# ------------------------------------------------------------------------------
# Optional (recommended): logLik method for 'angular'
# Lets you check 2*Delta_logLik and use AIC/BIC if desired.
# Log-likelihood at kappa_hat: sum[kappa_hat * cos(res_i) - log(2*pi*I0(kappa_hat))]
# where res_i = y_i - mu_i. Your MaxCosine = mean cos(res_i) is used for speed.
# ------------------------------------------------------------------------------

#' @export
logLik.angular <- function(object, ...) {
  n <- length(object$y)
  kappa <- as.numeric(object$kappahat)
  # mean residual cosine (already stored)
  Cbar <- as.numeric(object$MaxCosine)
  # von Mises normalizing const: 2*pi*I0(kappa)
  # Use Bessel I0 from base R via besselI(x, 0, expon.scaled = FALSE)
  log_c <- log(2 * pi) + log(besselI(kappa, nu = 0, expon.scaled = FALSE))
  ll <- n * (kappa * Cbar - log_c)
  # number of *beta* parameters (kappa is treated as nuisance here; for AIC/BIC you can decide whether to add +1)
  attr(ll, "df") <- nrow(object$parameters) + 1L # count beta + kappa for AIC/BIC
  class(ll) <- "logLik"
  ll
}
