# Internal helpers for circular model parsing and reference handling -----------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

.btick <- function(x) {
  x_chr <- as.character(x)
  paste0("`", gsub("`", "\\\\`", x_chr, fixed = TRUE), "`")
}

.make_model_frame <- function(formula, data) {
  mf <- stats::model.frame(formula = formula, data = data)
  if (nrow(mf) == 0) {
    stop("No observations available.", call. = FALSE)
  }
  mf
}

.parse_term_labels <- function(term_labels) {
  if (length(term_labels) == 0) {
    stop(
      "The model must include at least one term specifying explanatory direction(s).",
      call. = FALSE
    )
  }

  parts <- strsplit(term_labels, ":", fixed = TRUE)
  term_info <- lapply(seq_along(parts), function(i) {
    pi <- trimws(parts[[i]])
    if (length(pi) == 1) {
      list(
        label = term_labels[i],
        angle = pi[1],
        modifier = NA_character_,
        has_modifier = FALSE
      )
    } else if (length(pi) == 2) {
      list(
        label = term_labels[i],
        angle = pi[1],
        modifier = pi[2],
        has_modifier = TRUE
      )
    } else {
      stop("Each term must be of the form 'x' or 'x:z'.", call. = FALSE)
    }
  })

  do.call(
    rbind,
    lapply(term_info, function(x) {
      data.frame(
        label = x$label,
        angle = x$angle,
        modifier = x$modifier,
        has_modifier = x$has_modifier,
        stringsAsFactors = FALSE
      )
    })
  )
}

.circular_design <- function(formula, data) {
  mf <- .make_model_frame(formula = formula, data = data)
  y <- as.numeric(mf[[1]])
  term_labels <- attr(attr(mf, "terms"), "term.labels")
  term_info <- .parse_term_labels(term_labels)

  if (any(!term_info$angle %in% names(mf))) {
    miss <- unique(term_info$angle[!term_info$angle %in% names(mf)])
    stop(
      sprintf("Angle variable(s) not found in data: %s", paste(miss, collapse = ", ")),
      call. = FALSE
    )
  }

  nobs <- nrow(mf)
  p <- nrow(term_info)
  matx <- matrix(0, nrow = nobs, ncol = p)
  matz <- matrix(1, nrow = nobs, ncol = p)
  colnames(matx) <- term_info$label
  colnames(matz) <- term_info$label

  for (j in seq_len(p)) {
    ang_name <- term_info$angle[j]
    matx[, j] <- as.numeric(mf[[ang_name]])

    if (isTRUE(term_info$has_modifier[j])) {
      mod_name <- term_info$modifier[j]
      if (!mod_name %in% names(mf)) {
        stop(
          sprintf("Modifier '%s' not found in data.", mod_name),
          call. = FALSE
        )
      }
      z <- as.numeric(mf[[mod_name]])
      if (any(!is.finite(z)) || any(z < 0)) {
        stop(
          sprintf("Modifier '%s' must be finite and non-negative.", mod_name),
          call. = FALSE
        )
      }
      matz[, j] <- z
    }
  }

  angle_names <- unique(term_info$angle)
  angle_data <- stats::setNames(
    lapply(angle_names, function(nm) as.numeric(mf[[nm]])),
    angle_names
  )

  plain_angles <- unique(term_info$angle[!term_info$has_modifier])

  list(
    mf = mf,
    y = y,
    nobs = nobs,
    term_labels = term_labels,
    term_info = term_info,
    matx = matx,
    matz = matz,
    angle_names = angle_names,
    angle_data = angle_data,
    plain_angles = plain_angles,
    response_name = names(mf)[1]
  )
}

.reference_scores <- function(y, angle_data, candidates) {
  stats::setNames(
    vapply(candidates, function(nm) mean(cos(y - angle_data[[nm]])), numeric(1)),
    candidates
  )
}

.resolve_reference <- function(reference, y, angle_data, candidates, tie_method = "first") {
  tie_method <- match.arg(tie_method, c("first", "random"))

  if (length(candidates) == 0) {
    stop(
      "No candidate reference angle found. Add at least one plain angle term 'x'.",
      call. = FALSE
    )
  }

  ref_input <- reference %||% "auto"

  if (length(ref_input) == 1 && !(ref_input %in% c("auto", "first", "name"))) {
    mode <- "name"
    explicit <- as.character(ref_input)
  } else {
    mode <- match.arg(ref_input[1], c("auto", "first", "name"))
    explicit <- if (length(ref_input) >= 2) as.character(ref_input[2]) else NULL
  }

  scores <- .reference_scores(y = y, angle_data = angle_data, candidates = candidates)

  if (mode == "auto") {
    best <- which(scores == max(scores))
    pick <- if (tie_method == "random") sample(best, 1L) else best[1]
    ref_name <- names(scores)[pick]
  } else if (mode == "first") {
    ref_name <- candidates[1]
  } else {
    if (is.null(explicit) || !nzchar(explicit)) {
      stop(
        "When reference = 'name', provide the angle name as second element, e.g. c('name', 'x1').",
        call. = FALSE
      )
    }
    if (!explicit %in% candidates) {
      stop(
        sprintf(
          "Reference '%s' is invalid. Allowed reference angles: %s",
          explicit,
          paste(candidates, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    ref_name <- explicit
  }

  list(
    mode = mode,
    tie_method = tie_method,
    candidates = candidates,
    scores = scores,
    ref_name = ref_name
  )
}

.reference_term_index <- function(term_info, ref_name) {
  idx <- which(term_info$angle == ref_name & !term_info$has_modifier)
  if (length(idx) == 0) {
    stop(
      sprintf(
        "Reference angle '%s' must appear as a plain term 'x' (without ':z') in the formula.",
        ref_name
      ),
      call. = FALSE
    )
  }
  idx[1]
}

.angular_design <- function(formula, data, reference = "auto", tie_method = "first") {
  des <- .circular_design(formula = formula, data = data)

  ref <- .resolve_reference(
    reference = reference,
    y = des$y,
    angle_data = des$angle_data,
    candidates = des$plain_angles,
    tie_method = tie_method
  )

  ref_idx <- .reference_term_index(des$term_info, ref$ref_name)
  keep <- setdiff(seq_len(ncol(des$matx)), ref_idx)

  list(
    y = des$y,
    nobs = des$nobs,
    x0 = as.numeric(des$matx[, ref_idx]),
    matx = des$matx[, keep, drop = FALSE],
    matz = des$matz[, keep, drop = FALSE],
    betaname = des$term_labels[keep],
    term_info = des$term_info,
    term_labels = des$term_labels,
    full_matx = des$matx,
    full_matz = des$matz,
    ref_idx = ref_idx,
    reference = ref,
    response_name = des$response_name,
    angle_data = des$angle_data,
    plain_angles = des$plain_angles
  )
}

.consensus_design <- function(formula, data) {
  .circular_design(formula = formula, data = data)
}

.identifiability_matrix <- function(beta_full, matx_full, matz_full) {
  n <- nrow(matx_full)
  p <- ncol(matx_full)
  out <- matrix(0, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      out[i, j] <- matz_full[i, j] *
        sum(beta_full * matz_full[i, ] * sin(matx_full[i, j] - matx_full[i, ]))
    }
  }
  out
}

#' Select Reference Angle for Homogeneous Angular Model
#'
#' Selects the reference direction using the criterion of Rivest et al. (2016):
#' the angle variable maximizing the empirical mean cosine
#' \eqn{\overline{\cos(y-x_j)}} among plain angle terms \code{x} present in the
#' formula.
#'
#' @param formula Model formula containing terms of the form \code{x} or
#'   \code{x:z}.
#' @param data Data frame used to evaluate the formula.
#' @param tie_method Tie-break rule when multiple candidates share the same
#'   score: \code{"first"} (default) or \code{"random"}.
#' @return An object of class \code{"select_reference_angle"}.
#' @export
select_reference_angle <- function(formula, data, tie_method = c("first", "random")) {
  tie_method <- match.arg(tie_method)
  call <- match.call()
  des <- .circular_design(formula = formula, data = data)
  ref <- .resolve_reference(
    reference = "auto",
    y = des$y,
    angle_data = des$angle_data,
    candidates = des$plain_angles,
    tie_method = tie_method
  )

  out <- list(
    call = call,
    formula = formula,
    ref_name = ref$ref_name,
    candidates = ref$candidates,
    scores = ref$scores,
    tie_method = tie_method
  )
  class(out) <- "select_reference_angle"
  out
}

#' @export
print.select_reference_angle <- function(x, digits = 4, ...) {
  cat("Reference angle selection (Rivest et al., 2016 criterion)\n")
  cat("Selected reference:", x$ref_name, "\n\n")
  tab <- data.frame(
    angle = names(x$scores),
    mean_cos = as.numeric(x$scores),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  tab$mean_cos <- round(tab$mean_cos, digits)
  print(tab, row.names = FALSE)
  invisible(x)
}
