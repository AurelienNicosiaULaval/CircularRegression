# ------- Helpers -------

.extract_angles_from_formula <- function(formula) {
  # Extract angle names (left part of "x:z" terms) in order
  tl <- attr(terms(formula), "term.labels")
  if (length(tl) == 0) return(character(0))
  vapply(
    strsplit(tl, ":", fixed = TRUE),
    function(p) p[1],
    FUN.VALUE = character(1)
  )
}

.make_angular_formula_with_ref <- function(formula, ref_name) {
  # Replace "x:z" with "ref(x)" for the reference angle
  y <- deparse(formula[[2]])
  tl <- attr(terms(formula), "term.labels")
  new_terms <- vapply(
    strsplit(tl, ":", fixed = TRUE),
    function(p) {
      x <- p[1]
      z <- if (length(p) >= 2) p[2] else NULL
      if (identical(x, ref_name)) {
        sprintf("ref(%s)", x) # ignore z for reference angle
      } else if (!is.null(z)) {
        sprintf("%s:%s", x, z)
      } else {
        x
      }
    },
    FUN.VALUE = character(1)
  )
  as.formula(paste(y, "~", paste(new_terms, collapse = " + ")))
}

.extract_beta_consensus <- function(fit) {
  # Try coef() first
  b <- try(coef(fit), silent = TRUE)
  if (!inherits(b, "try-error")) {
    if (is.list(b) && !is.null(b$beta) && is.numeric(b$beta)) return(b$beta)
    if (is.numeric(b)) return(b)
  }
  # Fallbacks on internal slots
  for (cand in c("beta", "par", "coefficients", "pars")) {
    if (!is.null(fit[[cand]])) {
      obj <- fit[[cand]]
      if (is.numeric(obj)) return(obj)
      if (is.list(obj) && !is.null(obj$beta) && is.numeric(obj$beta))
        return(obj$beta)
    }
  }
  stop("Could not extract beta coefficients from consensus fit object.")
}

# ------- Main function -------

#' Pick reference angle from a consensus model (CircularRegression)
#'
#' @param formula A formula like y ~ x1:z1 + x2:z2 + ...
#' @param data A data.frame with response y, angular predictors x_j and positive variables z_j
#' @param tie_method How to break ties if multiple |beta| are equal: "first" or "random"
#' @param build_angular_formula If TRUE, also returns a formula ready for angular() with ref(angle)
#' @param ... Additional arguments passed to consensus()
#'
#' @return A list with:
#'   - ref_index: index of the chosen reference angle
#'   - ref_name: name of the chosen reference angle
#'   - beta_hat: raw consensus beta estimates (scale-free)
#'   - beta_ref1: beta estimates rescaled so that beta[ref] = 1
#'   - angular_formula (optional): formula with ref(angle_ref) ready for angular()
#'
#' @export
pick_reference_angle <- function(
  formula,
  data,
  tie_method = c("first", "random"),
  build_angular_formula = TRUE,
  ...
) {
  tie_method <- match.arg(tie_method)

  # 1) Fit consensus model
  fit <- consensus(formula = formula, data = data, ...)

  # 2) Extract beta estimates
  beta_hat <- .extract_beta_consensus(fit)
  if (!is.numeric(beta_hat)) stop("Beta coefficients are not numeric.")

  # 3) Map betas to angle names (from formula or names(beta_hat))
  angle_names <- .extract_angles_from_formula(formula)
  if (length(angle_names) && length(angle_names) != length(beta_hat)) {
    warning(
      "Number of angles from formula does not match length of beta_hat. ",
      "Using names(beta_hat) if available."
    )
    angle_names <- if (!is.null(names(beta_hat))) names(beta_hat) else
      rep(NA_character_, length(beta_hat))
  } else if (!length(angle_names) && !is.null(names(beta_hat))) {
    angle_names <- names(beta_hat)
  }

  # 4) Pick reference index (max |beta|)
  absb <- abs(beta_hat)
  if (tie_method == "first") {
    ref_index <- which.max(absb)
  } else {
    ref_index <- sample(which(absb == max(absb)), 1L)
  }
  ref_name <- if (length(angle_names)) angle_names[ref_index] else NA_character_

  # 5) Rescale so that reference beta = 1
  beta_ref1 <- beta_hat / beta_hat[ref_index]

  out <- list(
    ref_index = ref_index,
    ref_name = ref_name,
    beta_hat = beta_hat,
    beta_ref1 = beta_ref1
  )

  # 6) Build angular formula if requested
  if (isTRUE(build_angular_formula) && !is.na(ref_name)) {
    out$angular_formula <- .make_angular_formula_with_ref(formula, ref_name)
  }

  out
}
