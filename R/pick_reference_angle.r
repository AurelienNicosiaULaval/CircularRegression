# ======================================================================
# Utilities
# ======================================================================

#' Safe infix "a %||% b": returns a if not NULL, otherwise b
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ======================================================================
# Helpers: formula parsing & ref() builder
# ======================================================================

#' Extract unique angle names (left of "x:z"), preserving first appearance order
#' @noRd
#'
#' Example:
#'   y ~ x.meadow + x.meadow:z.meadow + x.gap + x.gap:z.gap
#' -> c("x.meadow","x.gap")
.extract_angles_from_formula <- function(formula) {
  tl <- attr(terms(formula), "term.labels")
  if (length(tl) == 0) return(character(0))
  ang <- vapply(
    strsplit(tl, ":", fixed = TRUE),
    function(p) p[1],
    FUN.VALUE = character(1)
  )
  ang[!duplicated(ang)]
}

#' Build an angular formula with ref(angle_ref) for the chosen reference angle
#' @noRd
#'
#' Only the left symbol of any "x:z" pair is considered an angle. The reference
#' term gets rewritten as ref(x).
.make_angular_formula_with_ref <- function(formula, ref_name) {
  y <- deparse(formula[[2]])
  tl <- attr(terms(formula), "term.labels")
  new_terms <- vapply(
    strsplit(tl, ":", fixed = TRUE),
    function(p) {
      x <- p[1]
      z <- if (length(p) >= 2) p[2] else NULL
      if (identical(x, ref_name)) {
        sprintf("ref(%s)", x) # ignore z for ref
      } else if (!is.null(z)) {
        sprintf("%s:%s", x, z)
      } else x
    },
    FUN.VALUE = character(1)
  )
  as.formula(paste(y, "~", paste(new_terms, collapse = " + ")))
}


# ======================================================================
# Helpers: extraction from consensus() fit
# ======================================================================

#' Extract beta (scale-free) from a consensus() fit object, robustly
#' @noRd
.extract_beta_consensus <- function(fit) {
  # try coef() first
  b <- try(coef(fit), silent = TRUE)
  if (!inherits(b, "try-error")) {
    if (is.list(b) && !is.null(b$beta) && is.numeric(b$beta)) return(b$beta)
    if (is.numeric(b)) return(b)
  }
  # common internal slots fallbacks
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

#' Extract kappa table (estimate, sd) from a consensus() fit object
#'
#' Returns a data.frame with columns "estimate" and "sd". Row names are made
#' unique to avoid duplicate 'row.names' errors.
#' @noRd
.extract_kappa_table <- function(fit, angle_names = NULL) {
  # Preferred path: fit$parameters already shaped like a table
  if (!is.null(fit$parameters)) {
    tab <- as.data.frame(fit$parameters, stringsAsFactors = FALSE)

    # normalize column names
    cn <- tolower(names(tab))
    if (!("estimate" %in% cn)) names(tab)[1] <- "estimate"
    if (!("sd" %in% cn)) {
      if (ncol(tab) >= 2) names(tab)[2] <- "sd" else tab$sd <- NA_real_
    }

    # choose rownames: prefer existing rownames from consensus(); else angle_names
    rn <- rownames(tab)
    if (is.null(rn)) rn <- angle_names
    # only apply angle_names if it matches the number of rows
    if (!is.null(angle_names) && length(angle_names) == nrow(tab))
      rn <- angle_names
    # fallback generic names if still missing or length mismatch
    if (is.null(rn) || length(rn) != nrow(tab))
      rn <- paste0("x", seq_len(nrow(tab)))

    # ensure uniqueness to avoid "duplicate 'row.names' are not allowed"
    rownames(tab) <- make.unique(rn)

    return(tab[, c("estimate", "sd"), drop = FALSE])
  }

  # Fallback path: attempt to find a numeric vector of kappa
  kap <- NULL
  for (cand in c("kappa", "kappas", "par_kappa", "pars_kappa", "par", "pars")) {
    if (!is.null(fit[[cand]])) {
      obj <- fit[[cand]]
      if (is.numeric(obj)) {
        kap <- obj
        break
      }
      if (is.list(obj) && !is.null(obj$kappa) && is.numeric(obj$kappa)) {
        kap <- obj$kappa
        break
      }
    }
  }
  if (is.null(kap)) stop("Could not extract the kappa values from the consensus() fit.")

  rn <- angle_names %||% names(kap) %||% paste0("x", seq_along(kap))
  data.frame(
    estimate = as.numeric(kap),
    sd = NA_real_,
    row.names = make.unique(rn),
    stringsAsFactors = FALSE
  )
}

#' Add the beta_ref = 1 row (sd = 0) on top of a 2-column beta table
#' @noRd
.add_beta1_row <- function(beta_tab, ref_name, put_first = TRUE) {
  beta_tab <- beta_tab[, c("estimate", "sd"), drop = FALSE]
  b1 <- data.frame(estimate = 1, sd = 0, row.names = ref_name %||% "ref")
  out <- rbind(b1, beta_tab)
  # guard against duplicated names
  rownames(out) <- make.unique(rownames(out))
  if (put_first) return(out)
  out
}


# ======================================================================
# Main: pick_reference_angle()
# ======================================================================

#' Pick reference angle from a consensus model and assemble kappa/beta tables
#'
#' This function:
#'   (1) fits `consensus(formula, data, ...)`,
#'   (2) picks the reference angle as argmax |beta|,
#'   (3) rescales beta so that beta_ref = 1,
#'   (4) returns kappa (estimate, sd) as `parameters` and beta (estimate, sd) as `parambeta`,
#'       with the first beta_ref row added (estimate=1, sd=0),
#'   (5) optionally returns an `angular_formula` where the reference angle is written `ref(angle)`.
#'
#' @param formula A formula like y ~ x1 + x1:z1 + x2 + x2:z2 + ...
#'                The left symbol of any "x:z" term is treated as an angular predictor.
#' @param data A data.frame with response y, angular predictors x_j and positive variables z_j
#' @param tie_method How to break ties if multiple |beta| are equal: "first" or "random"
#' @param build_angular_formula If TRUE, also returns a formula ready for angular() with ref(angle)
#' @param ... Additional arguments passed to consensus()
#'
#' @return A list with:
#'   - ref_index: index of the chosen reference angle in beta
#'   - ref_name:  name of the chosen reference angle
#'   - beta_hat:  raw consensus beta estimates (scale-free)
#'   - beta_ref1: beta rescaled so that beta[ref] = 1
#'   - parameters: data.frame (kappa estimate, sd)
#'   - parambeta:  data.frame (beta estimate, sd) with beta_ref added on the first row (sd = 0)
#'   - angular_formula (optional): formula with ref(angle_ref) ready for angular()
#' @export
pick_reference_angle <- function(
  formula,
  data,
  tie_method = c("first", "random"),
  build_angular_formula = TRUE,
  ...
) {
  tie_method <- match.arg(tie_method)

  # 1) Fit consensus model (kappa-parameterization internally)
  fit <- consensus(formula = formula, data = data, ...)

  # 2) Extract beta (scale-free) and basic angle names
  beta_hat <- .extract_beta_consensus(fit)
  if (!is.numeric(beta_hat)) stop("Beta coefficients are not numeric.")

  angle_names <- .extract_angles_from_formula(formula)

  # If counts mismatch, prefer names(beta_hat) or stay unnamed; warn once.
  if (length(angle_names) && length(angle_names) != length(beta_hat)) {
    warning(
      "Nombre d'angles (formula) != longueur(beta_hat). On tente names(beta_hat)."
    )
    angle_names <- if (!is.null(names(beta_hat))) names(beta_hat) else
      rep(NA_character_, length(beta_hat))
  } else if (!length(angle_names) && !is.null(names(beta_hat))) {
    angle_names <- names(beta_hat)
  }

  # 3) Choose reference index (max |beta|)
  absb <- abs(beta_hat)
  ref_index <- if (tie_method == "first") which.max(absb) else
    sample(which(absb == max(absb)), 1L)
  ref_name <- if (length(angle_names)) angle_names[ref_index] else NA_character_

  # 4) Rescale beta so that beta_ref = 1
  beta_ref1 <- beta_hat / beta_hat[ref_index]
  names(beta_hat) <- angle_names %||% names(beta_hat)
  names(beta_ref1) <- names(beta_hat)

  # 5) Kappa table (parameters)
  kappa_tab <- .extract_kappa_table(fit, angle_names = angle_names)

  # 6) Beta table (parambeta): prefer fit$parambeta if available, then add beta_ref row
  if (!is.null(fit$parambeta)) {
    pb <- as.data.frame(fit$parambeta, stringsAsFactors = FALSE)

    # normalize to two columns: estimate, sd
    cn <- tolower(names(pb))
    if (!("estimate" %in% cn)) names(pb)[1] <- "estimate"
    if (!("sd" %in% cn)) {
      if (ncol(pb) >= 2) names(pb)[2] <- "sd" else pb$sd <- NA_real_
    }
    pb <- pb[, c("estimate", "sd"), drop = FALSE]

    # keep consensus-provided rownames as-is if present; otherwise, try to infer
    if (is.null(rownames(pb)) || anyNA(rownames(pb))) {
      # remove the (unknown) reference from angle_names
      others <- setdiff(angle_names, ref_name)
      if (length(others) == nrow(pb)) rownames(pb) <- others else
        rownames(pb) <- paste0("beta", seq_len(nrow(pb)))
    }

    parambeta <- .add_beta1_row(pb, ref_name = ref_name, put_first = TRUE)
  } else {
    # fallback: construct from beta_ref1 (sd unknown -> NA), then set sd_ref = 0
    pb <- data.frame(
      estimate = as.numeric(beta_ref1),
      sd = NA_real_,
      row.names = names(beta_ref1),
      stringsAsFactors = FALSE
    )
    pb[ref_name, "sd"] <- 0
    parambeta <- pb
  }

  # 7) Assemble output
  out <- list(
    ref_index = ref_index,
    ref_name = ref_name,
    beta_hat = beta_hat,
    beta_ref1 = beta_ref1,
    parameters = kappa_tab, # kappa (estimate, sd)
    parambeta = parambeta # beta (estimate, sd) with beta_ref row added
  )

  # 8) Optional angular formula with ref()
  if (isTRUE(build_angular_formula) && !is.na(ref_name)) {
    out$angular_formula <- .make_angular_formula_with_ref(formula, ref_name)
  }

  class(out) <- "pick_reference_angle" # <--- ajoute cette ligne ici
  out
}


# ======================================================================
#' Print method for pick_reference_angle objects
#'
#' @param x A list returned by pick_reference_angle()
#' @param digits Number of digits to print (default 4)
#' @param ... Unused
#'
#' @export
print.pick_reference_angle <- function(x, digits = 4, ...) {
  cat("\n--- Picked reference angle ---\n")
  cat("Reference index:", x$ref_index, "\n")
  cat("Reference name : ", x$ref_name, "\n\n")

  # Beta summary
  cat("Rescaled beta coefficients (beta_ref = 1):\n")
  print(round(x$beta_ref1, digits))
  cat("\n")

  # Kappa table
  if (!is.null(x$parameters)) {
    cat("Kappa parameters (kappa):\n")
    df <- x$parameters
    df[] <- lapply(
      df,
      function(col) if (is.numeric(col)) round(col, digits) else col
    )
    print(df)
    cat("\n")
  }

  # Beta table
  if (!is.null(x$parambeta)) {
    cat("Beta parameters (beta):\n")
    dfb <- x$parambeta
    dfb[] <- lapply(
      dfb,
      function(col) if (is.numeric(col)) round(col, digits) else col
    )
    print(dfb)
    cat("\n")
  }

  invisible(x)
}
