# A consensus recession direction d satisfies B d = 0 and A d >= 0,
# where B_ij = z_ij sin(x_ij - y_i) and A_ij = z_ij cos(x_ij - y_i).
# If at least one A d is positive, increasing its magnitude increases the
# likelihood without a finite maximum. Rows of zero case weight impose no
# constraint. A returned direction is checked against both conditions.
.consensus_recession <- function(y, matx, matz, weights) {
  active <- weights > 0
  x <- matx[active, , drop = FALSE]
  z <- matz[active, , drop = FALSE]
  y <- y[active]
  scale <- apply(abs(z), 2L, max)
  scale[scale == 0] <- 1
  z <- sweep(z, 2L, scale, "/")
  A <- z * cos(x - y)
  B <- z * sin(x - y)
  p <- ncol(B)
  decomposition <- svd(B, nu = 0, nv = p)
  rank_tolerance <- 10 * max(dim(B)) * .Machine$double.eps *
    max(1, decomposition$d)
  rank <- sum(decomposition$d > rank_tolerance)
  if (rank == p) return(NULL)

  null_basis <- decomposition$v[, seq.int(rank + 1L, p), drop = FALSE]
  design <- A %*% null_basis
  average <- colMeans(design)
  if (sum(average^2) <= .Machine$double.eps^2) return(NULL)

  # The squared-hinge criterion is convex. Its zero set intersects mean(A d)=1
  # exactly when the null space contains a nontrivial recession direction.
  objective <- function(v) {
    projection <- as.vector(design %*% v)
    mean(pmin(projection, 0)^2) + (mean(projection) - 1)^2
  }
  gradient <- function(v) {
    projection <- as.vector(design %*% v)
    as.vector(2 * crossprod(design, pmin(projection, 0)) / nrow(design) +
      2 * (mean(projection) - 1) * average)
  }
  fit <- stats::optim(
    average / sum(average^2), objective, gradient,
    method = "BFGS", control = list(reltol = 1e-13, maxit = 1000)
  )
  direction <- as.vector(null_basis %*% fit$par)
  direction <- direction / sqrt(sum(direction^2))
  aligned <- as.vector(A %*% direction)
  perpendicular <- as.vector(B %*% direction)
  tolerance <- 100 * rank_tolerance
  if (max(abs(perpendicular)) > tolerance ||
      min(aligned) < -tolerance || max(aligned) <= tolerance) {
    return(NULL)
  }
  direction / scale
}

.named_covariance <- function(value, parameter_names) {
  dimnames(value) <- list(parameter_names, parameter_names)
  value
}

.print_fit_status <- function(x) {
  cat("Converged:", isTRUE(x$convergence), "\n")
  if (!is.null(x$fit_status)) cat("Fit status:", x$fit_status, "\n")
  if (!is.null(x$inference_available)) {
    cat("Wald inference:", if (isTRUE(x$inference_available)) "available" else "unavailable", "\n")
  }
  invisible(x)
}
