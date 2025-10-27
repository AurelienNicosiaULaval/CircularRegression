#' Angular mixed-effects regression with von Mises random intercepts
#'
#' Fits the homogeneous angular regression of Rivest et al. (2016) augmented
#' with a cluster-specific random intercept \eqn{a_i \sim \mathrm{VM}(0, \kappa_a)}
#' and unit-errors \eqn{e_{ij} \sim \mathrm{VM}(0, \kappa_e)}, following Rivest
#' & Kato (2019). The fixed-effects component shares the simplified syntax of
#' \code{angular()}, where the first term provides the reference direction and
#' subsequent terms are written \code{x:z}.
#'
#' @param formula Model formula of the form \code{y ~ x1:z1 + x2:z2 + ...}. The
#'   first angular covariate supplies the reference direction whose coefficient
#'   is fixed to one.
#' @param data A data frame containing the variables used in \code{formula}.
#' @param cluster Factor or vector of cluster identifiers (must have the same
#'   length as the response).
#' @param init Optional list with components \code{beta}, \code{kappa_e}, and
#'   \code{kappa_a} providing starting values. Missing entries are filled with
#'   moment-based estimates.
#' @param control Optional list controlling the optimisation with components
#'   \code{maxit} (default 2000), \code{reltol} (1e-8), \code{trace} (0/1),
#'   \code{hessian} (\code{TRUE}), and \code{start_from_angular} (\code{TRUE})
#'   indicating whether \code{angular()} should be used to initialise the fixed
#'   effects.
#' @return object of class \code{"angular_re"}
#' @importFrom stats model.frame optim pnorm
#' @examples
#' # See README example below with Sandhopper data
#' @export
angular_re <- function(
  formula,
  data,
  cluster,
  init = list(),
  control = list()
) {
  stopifnot(!missing(cluster))
  call <- match.call()

  ## ---- Parse model frame exactly like angular() ----
  mfargs <- match(c("formula", "data"), names(call), 0L)
  call_subset <- call[c(1L, mfargs)]
  call_subset[[1L]] <- as.name("model.frame")
  mf <- eval(call_subset, parent.frame())
  y <- as.numeric(mf[[1]])
  n <- length(y)

  if (length(cluster) != n)
    stop("cluster must have length equal to the response length.")
  id <- as.factor(cluster)
  lev <- levels(id)
  m <- length(lev)
  cluster_index <- split(seq_len(n), id)

  term_labels <- attr(attr(mf, "terms"), "term.labels")
  if (length(term_labels) == 0) {
    stop("The model must include at least one term specifying the reference direction.")
  }
  split_terms <- strsplit(term_labels, ":", fixed = TRUE)
  term_matrix <- t(vapply(split_terms, function(x) {
    if (length(x) == 1) {
      c(x, x)
    } else if (length(x) == 2) {
      x
    } else {
      stop("Each term must be of the form 'x' or 'x:z'.", call. = FALSE)
    }
  }, character(2)))

  ref_name <- term_matrix[1, 1]
  if (!ref_name %in% names(mf)) {
    stop(
      sprintf("Reference direction '%s' not found in the supplied data.", ref_name),
      call. = FALSE
    )
  }
  x0 <- mf[[ref_name]]

  betaname <- if (length(term_labels) > 1) term_labels[-1] else character(0)
  p <- length(betaname)
  if (p > 0) {
    x_names <- term_matrix[-1, 1]
    z_names <- term_matrix[-1, 2]
    matx <- as.matrix(mf[, x_names, drop = FALSE])
    matz <- matrix(1, nrow = n, ncol = p)
    for (j in seq_len(p)) {
      if (!identical(z_names[j], x_names[j])) {
        if (!z_names[j] %in% names(mf)) {
          stop(
            sprintf("Modifier '%s' not found in the supplied data.", z_names[j]),
            call. = FALSE
          )
        }
        matz[, j] <- mf[[z_names[j]]]
      }
    }
  } else {
    matx <- matrix(0, nrow = n, ncol = 0)
    matz <- matrix(0, nrow = n, ncol = 0)
  }
  nparam <- p

  ## ---- helpers ----
  wrap <- function(a) atan2(sin(a), cos(a))
  A1 <- function(k) {
    # safe ratio I1/I0 using scaled Bessel to avoid overflow
    if (any(k < 0)) stop("kappa must be >= 0")
    num <- besselI(k, 1, expon.scaled = TRUE)
    den <- besselI(k, 0, expon.scaled = TRUE)
    r <- num / den
    r[k == 0] <- 0
    r
  }
  A1inv <- function(r) {
    # use circular::A1inv if available, else Newton
    if (requireNamespace("circular", quietly = TRUE)) {
      return(circular::A1inv(pmin(pmax(r, 0), 0.999999)))
    } else {
      # crude Newton for robustness
      r <- pmin(pmax(r, 0), 0.999999)
      k <- ifelse(
        r < 0.53,
        2 * r + r^3 + (5 * r^5) / 6,
        ifelse(
          r < 0.85,
          -0.4 + 1.39 * r + 0.43 / (1 - r),
          1 / (r^3 - 4 * r^2 + 3 * r)
        )
      )
      for (iter in 1:10) {
        f <- A1(k) - r
        df <- 1 - A1(k) / k - A1(k)^2
        k <- pmax(0, k - f / df)
      }
      return(k)
    }
  }
  logI0 <- function(k) log(besselI(k, 0, expon.scaled = TRUE)) + k

  compute_mu <- function(beta) {
    sinmu <- sin(x0)
    cosmu <- cos(x0)
    if (p > 0) {
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta)
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta)
    }
    long <- sqrt(sinmu^2 + cosmu^2)
    mu <- atan2(sinmu, cosmu)
    list(mu = mu, long = long)
  }

  dmu_dbeta <- function(beta, mu, long) {
    if (p == 0) {
      return(matrix(0, nrow = length(mu), ncol = 0))
    }
    denom <- ifelse(long > 0, long, 1)
    (matz * sin(matx - mu)) / as.vector(denom)
  }

  ## ---- initial values ----
  ctrl <- modifyList(
    list(
      maxit = 2000,
      reltol = 1e-8,
      trace = 0,
      hessian = TRUE,
      start_from_angular = TRUE
    ),
    control
  )

  coalesce_null <- function(x, default) if (is.null(x)) default else x

  if (p == 0) {
    beta0 <- numeric(0)
  } else if (isTRUE(ctrl$start_from_angular) && is.null(init$beta)) {
    if (!exists("angular", mode = "function")) {
      warning("angular() not found; starting beta at zeros.")
      beta0 <- rep(0, p)
    } else {
      an0 <- angular(formula, data = data)
      beta0 <- as.numeric(an0$parameters[betaname, "estimate"])
    }
  } else {
    beta0 <- coalesce_null(init$beta, rep(0, p))
    if (length(beta0) != p) {
      stop(sprintf("'init$beta' must have length %d", p), call. = FALSE)
    }
  }

  mu0 <- compute_mu(beta0)$mu
  r0 <- wrap(y - mu0)

  # moment estimators for kappas (Sec 3.1 RK2019)
  pair_sum <- 0
  pair_count <- 0
  for (idx in cluster_index) {
    ni <- length(idx)
    if (ni >= 2) {
      rr <- r0[idx]
      cos_sum <- sum(cos(rr))
      sin_sum <- sum(sin(rr))
      pair_sum <- pair_sum + (cos_sum^2 + sin_sum^2 - ni)
      pair_count <- pair_count + ni * (ni - 1)
    }
  }
  mean_cosdiff <- if (pair_count > 0) pair_sum / pair_count else 0
  A1_ke <- sqrt(pmax(mean_cosdiff, 0))
  kappa_e0 <- coalesce_null(init$kappa_e, A1inv(A1_ke))
  A1_ka <- if (A1_ke < 1e-8) mean(cos(r0), na.rm = TRUE) else mean(cos(r0), na.rm = TRUE) / A1_ke
  A1_ka <- pmin(pmax(A1_ka, 0), 0.999999)
  kappa_a0 <- coalesce_null(init$kappa_a, if (A1_ka >= 0.999999) 50 else A1inv(A1_ka))

  par0 <- c(beta0, kappa_e0, kappa_a0)
  names(par0) <- c(betaname, "kappa_e", "kappa_a")

  ## ---- loglik + score (by cluster) ----
  pack_par <- function(par) {
    list(
      beta = par[seq_len(nparam)],
      kappa_e = par[nparam + 1],
      kappa_a = par[nparam + 2]
    )
  }

  invert_hessian <- function(M) {
    M_sym <- (M + t(M)) / 2
    chol_res <- tryCatch(chol(M_sym), error = function(e) NULL)
    if (is.null(chol_res)) {
      warning("Observed information matrix is not positive definite; returning NA variances.")
      matrix(NA_real_, nrow = nrow(M_sym), ncol = ncol(M_sym))
    } else {
      chol2inv(chol_res)
    }
  }

  objective <- function(par) {
    pa <- pack_par(par)
    cm <- compute_mu(pa$beta)
    mu <- cm$mu
    long <- cm$long
    r <- wrap(y - mu)

    # cluster loop
    val <- 0
    for (idx in cluster_index) {
      ni <- length(idx)
      if (ni == 0) next
      cg <- sum(cos(r[idx]))
      sg <- sum(sin(r[idx]))
      sigma <- sqrt((pa$kappa_a + pa$kappa_e * cg)^2 + (pa$kappa_e * sg)^2)
      val <- val + (logI0(sigma) - ni * logI0(pa$kappa_e) - logI0(pa$kappa_a))
    }
    as.numeric(val)
  }

  gradient <- function(par) {
    pa <- pack_par(par)
    cm <- compute_mu(pa$beta)
    mu <- cm$mu
    long <- cm$long
    r <- wrap(y - mu)
    D <- dmu_dbeta(pa$beta, mu, long) # n x nparam

    s <- rep(0, nparam + 2)
    # store cluster scores for sandwich
    S_by_cluster <- matrix(0, nrow = m, ncol = nparam + 2)

    for (g in seq_along(cluster_index)) {
      idx <- cluster_index[[g]]
      ni <- length(idx)
      rr <- r[idx]
      Dg <- D[idx, , drop = FALSE]
      cg <- sum(cos(rr))
      sg <- sum(sin(rr))
      sigma <- sqrt((pa$kappa_a + pa$kappa_e * cg)^2 + (pa$kappa_e * sg)^2)
      A1sig <- A1(sigma)
      a_tilde <- atan2(pa$kappa_e * sg, pa$kappa_a + pa$kappa_e * cg)

      # components
      s_beta <- if (p > 0) {
        pa$kappa_e * A1sig * colSums(Dg * (sin(rr - a_tilde)))
      } else {
        numeric(0)
      }
      s_ke <- A1sig * sum(cos(rr - a_tilde)) - ni * A1(pa$kappa_e)
      s_ka <- A1sig * cos(a_tilde) - A1(pa$kappa_a)

      sg_vec <- c(s_beta, s_ke, s_ka)
      s <- s + sg_vec
      S_by_cluster[g, ] <- sg_vec
    }

    attr(s, "cluster_scores") <- S_by_cluster
    s
  }

  ## ---- maximize ----
  opt <- optim(
    par0,
    fn = function(p) -objective(p),
    gr = function(p) -gradient(p),
    method = "BFGS",
    control = list(
      maxit = ctrl$maxit,
      reltol = ctrl$reltol,
      trace = ctrl$trace
    ),
    hessian = isTRUE(ctrl$hessian)
  )

  par_hat <- opt$par
  val <- -opt$value
  conv <- (opt$convergence == 0)

  # observed Fisher (negative Hessian)
  if (isTRUE(ctrl$hessian)) {
    H <- opt$hessian
    V_model <- invert_hessian(H)
  } else {
    H <- matrix(NA_real_, length(par_hat), length(par_hat))
    V_model <- H
  }

  # robust sandwich (cluster scores)
  sc <- attr(gradient(par_hat), "cluster_scores")
  if (isTRUE(ctrl$hessian)) {
    bread <- invert_hessian(H)
    if (all(is.finite(bread))) {
      meat <- crossprod(sc)
      V_robust <- bread %*% meat %*% bread
    } else {
      V_robust <- matrix(NA_real_, nrow(H), ncol(H))
    }
  } else {
    V_robust <- matrix(NA_real_, nrow(H), ncol(H))
  }

  # extract components and diagnostics
  pa_hat <- pack_par(par_hat)
  names(pa_hat$beta) <- betaname
  cmh <- compute_mu(pa_hat$beta)
  mu_hat <- cmh$mu
  r_fix <- wrap(y - mu_hat)

  # random effects predictors per cluster
  a_hat <- numeric(m)
  for (g in seq_along(cluster_index)) {
    idx <- cluster_index[[g]]
    rr <- r_fix[idx]
    cg <- sum(cos(rr))
    sg <- sum(sin(rr))
    a_hat[g] <- atan2(pa_hat$kappa_e * sg, pa_hat$kappa_a + pa_hat$kappa_e * cg)
  }
  names(a_hat) <- lev

  # conditional residuals
  r_cond <- r_fix
  for (g in seq_along(cluster_index)) {
    idx <- cluster_index[[g]]
    r_cond[idx] <- wrap(r_fix[idx] - a_hat[g])
  }

  # table of parameters with both SEs
  se_model <- sqrt(diag(V_model))
  se_robust <- sqrt(diag(V_robust))
  z_model <- par_hat / se_model
  p_model <- 2 * pnorm(abs(z_model), lower.tail = FALSE)

  tab <- cbind(
    estimate = par_hat,
    se_model = se_model,
    se_robust = se_robust,
    z_model = z_model,
    `Pr(>|z|)` = p_model
  )
  rownames(tab) <- c(betaname, "kappa_e", "kappa_a")

  # intra-cluster sine-sine correlation (rho_SS)
  A2 <- function(k) {
    # I2/I0, safer near 0
    besselI(k, 2, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
  }
  rho_SS <- (A1(pa_hat$kappa_e)^2 * (1 - A2(pa_hat$kappa_a))) /
    (1 - A2(pa_hat$kappa_a) * A2(pa_hat$kappa_e))

  out <- list(
    call = call,
    logLik = val,
    convergence = conv,
    coefficients = tab,
    vcov_model = V_model,
    vcov_robust = V_robust,
    beta = pa_hat$beta,
    betaname = betaname,
    kappa_e = pa_hat$kappa_e,
    kappa_a = pa_hat$kappa_a,
    rho_SS = as.numeric(rho_SS),
    mu = mu_hat,
    resid_fixed = r_fix,
    resid_conditional = r_cond,
    cluster = id,
    ranef = a_hat,
    control = ctrl
  )
  class(out) <- "angular_re"
  out
}

#' @export
print.angular_re <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = 5),
    "\nConverged:",
    x$convergence,
    "\n"
  )
  cat("\nEstimates (model SE / robust SE):\n")
  co <- x$coefficients
  show <- cbind(
    co[, "estimate"],
    co[, "se_model"],
    co[, "se_robust"],
    co[, "z_model"],
    co[, "Pr(>|z|)"]
  )
  colnames(show) <- c("Estimate", "SE(model)", "SE(robust)", "z", "Pr(>|z|)")
  print(show, digits = 4)
  cat(
    "\nIntra-cluster sine-sine correlation (rho_SS):",
    round(x$rho_SS, 3),
    "\n"
  )
  invisible(x)
}

#' @export
coef.angular_re <- function(object, ...) {
  object$coefficients[, "estimate"]
}

#' @export
vcov.angular_re <- function(object, robust = FALSE, ...) {
  if (robust) object$vcov_robust else object$vcov_model
}

#' @export
residuals.angular_re <- function(
  object,
  type = c("fixed", "conditional"),
  ...
) {
  type <- match.arg(type)
  if (type == "fixed") object$resid_fixed else object$resid_conditional
}

#' Extract estimated random intercepts from an angular mixed model
#'
#' Returns the cluster-specific random intercepts \eqn{\hat a_i} estimated by
#' \code{angular_re()}.
#'
#' @param object A fitted \code{angular_re} model.
#' @param ... Unused; included for compatibility with the generic.
#' @return A named numeric vector of random intercept estimates, one per cluster.
#' @export
ranef.angular_re <- function(object, ...) object$ranef

#' @export
fitted.angular_re <- function(object, ...) object$mu


#' Plot method for angular_re objects
#'
#' Draws rose/circular plots of residuals: fixed (y - mu) and conditional
#' (y - mu - a_hat[cluster]) as in Rivest & Kato (2019), Fig. 1-2.
#' If the 'circular' package is available, uses it; otherwise falls back
#' to a simple polar scatter with density ticks.
#'
#' @param x an object of class \code{angular_re}
#' @param which one of "both", "fixed", "conditional"
#' @param breaks number of bins for rose histogram when using \pkg{circular}
#' @param points logical; add individual points on the circle
#' @param main base title (panel-specific suffixes are added)
#' @param ... unused
#' @export
plot.angular_re <- function(
  x,
  which = c("both", "fixed", "conditional"),
  breaks = 24,
  points = TRUE,
  main = "Residuals",
  ...
) {
  which <- match.arg(which)
  r_fix <- x$resid_fixed
  r_con <- x$resid_conditional

  has_circ <- requireNamespace("circular", quietly = TRUE)

  draw_panel_circular <- function(angles, ttl) {
    th <- circular::circular(
      angles,
      type = "angles",
      units = "radians",
      modulo = "2pi",
      template = "none"
    )
    op <- par(mar = c(1.8, 1.8, 2.2, 1.8))
    on.exit(par(op))
    circular::rose.diag(
      th,
      bins = breaks,
      main = ttl,
      col = "grey85",
      border = NA
    )
    if (points) circular::points.circular(th, pch = 16, cex = 0.5)
    # resultant vector
    R <- abs(mean(exp(1i * angles)))
    phi <- Arg(mean(exp(1i * angles)))
    circular::arrows.circular(
      circular::circular(phi),
      y = 0,
      lwd = 2,
      length = 0.1
    )
    mtext(sprintf("R = %.3f", R), line = -1.2, cex = 0.9)
  }

  draw_panel_base <- function(angles, ttl) {
    th <- angles
    op <- par(mar = c(1.8, 1.8, 2.2, 1.8))
    on.exit(par(op))
    plot.new()
    plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)
    # circle
    t <- seq(0, 2 * pi, length.out = 360)
    lines(cos(t), sin(t), col = "grey60")
    # axes
    segments(-1, 0, 1, 0, col = "grey80")
    segments(0, -1, 0, 1, col = "grey80")
    title(ttl)
    if (points) points(cos(th), sin(th), pch = 16, cex = 0.5)
    # resultant vector
    v <- mean(exp(1i * th))
    arrows(0, 0, Re(v), Im(v), lwd = 2, length = 0.08)
    mtext(sprintf("R = %.3f", Mod(v)), line = -1.2, cex = 0.9)
  }

  if (which == "both") {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op), add = TRUE)
    if (has_circ) {
      draw_panel_circular(r_fix, sprintf("%s - fixed", main))
      draw_panel_circular(r_con, sprintf("%s - conditional", main))
    } else {
      draw_panel_base(r_fix, sprintf("%s - fixed", main))
      draw_panel_base(r_con, sprintf("%s - conditional", main))
    }
  } else {
    ttl <- sprintf("%s - %s", main, if (which == "fixed") "fixed" else "conditional")
    if (has_circ)
      draw_panel_circular(if (which == "fixed") r_fix else r_con, ttl) else
      draw_panel_base(if (which == "fixed") r_fix else r_con, ttl)
  }

  invisible(x)
}


#' Predict method for angular_re objects
#'
#' @param object an \code{angular_re} fit.
#' @param newdata data.frame; if NULL, uses model frame from the original call.
#' @param cluster optional vector/factor of cluster IDs for newdata.
#' @param type One of \code{"auto"} (default), \code{"marginal"} (mean
#'        direction) or \code{"conditional"} (mean direction plus the predicted
#'        random intercept). With \code{"auto"} the conditional mean is returned
#'        whenever cluster information is provided; otherwise the marginal mean
#'        is used.
#' @param a_hat optional numeric vector of random-intercept predictions to use
#'        for \code{cluster} (same length as newdata). Ignored unless \code{type="conditional"}.
#' @param se.fit not implemented (reserved).
#' @param ... unused
#' @return numeric vector (angles in radians).
#' @export
predict.angular_re <- function(
  object,
  newdata = NULL,
  cluster = NULL,
  type = c("auto", "marginal", "conditional"),
  a_hat = NULL,
  se.fit = FALSE,
  ...
) {
  type <- match.arg(type)
  if (isTRUE(se.fit)) stop("se.fit not implemented for angular_re.")

  form <- eval(object$call$formula)

  parse_design <- function(formula, data) {
    mf <- model.frame(formula, data = data)
    term_labels <- attr(attr(mf, "terms"), "term.labels")
    if (length(term_labels) == 0) {
      stop("The formula must include at least one term for the reference direction.")
    }
    split_terms <- strsplit(term_labels, ":", fixed = TRUE)
    term_matrix <- t(vapply(split_terms, function(x) {
      if (length(x) == 1) {
        c(x, x)
      } else if (length(x) == 2) {
        x
      } else {
        stop("Each term must be of the form 'x' or 'x:z'.", call. = FALSE)
      }
    }, character(2)))
    ref_name <- term_matrix[1, 1]
    x0 <- mf[[ref_name]]
    betaname <- if (length(term_labels) > 1) term_labels[-1] else character(0)
    p <- length(betaname)
    if (p > 0) {
      x_names <- term_matrix[-1, 1]
      z_names <- term_matrix[-1, 2]
      matx <- as.matrix(mf[, x_names, drop = FALSE])
      matz <- matrix(1, nrow = nrow(mf), ncol = p)
      for (j in seq_len(p)) {
        if (!identical(z_names[j], x_names[j])) {
          if (!z_names[j] %in% names(mf)) {
            stop(
              sprintf("Modifier '%s' not found in the supplied data.", z_names[j]),
              call. = FALSE
            )
          }
          matz[, j] <- mf[[z_names[j]]]
        }
      }
    } else {
      matx <- matrix(0, nrow = nrow(mf), ncol = 0)
      matz <- matrix(0, nrow = nrow(mf), ncol = 0)
    }
    list(x0 = x0, matx = matx, matz = matz, betaname = betaname, nobs = nrow(mf))
  }

  data_source <- if (is.null(newdata)) {
    if ("data" %in% names(object$call)) {
      eval(object$call$data, parent.frame())
    } else {
      NULL
    }
  } else {
    newdata
  }

  design <- parse_design(form, data_source)
  beta <- object$beta
  if (length(beta) != ncol(design$matz)) {
    stop("The supplied data are incompatible with the fitted fixed effects.")
  }

  compute_mu_new <- function(beta, x0, matx, matz) {
    sinmu <- sin(x0)
    cosmu <- cos(x0)
    if (length(beta) > 0) {
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta)
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta)
    }
    atan2(sinmu, cosmu)
  }

  mu <- compute_mu_new(beta, design$x0, design$matx, design$matz)
  wrap <- function(a) atan2(sin(a), cos(a))

  if (type == "auto") {
    type <- if (!is.null(cluster) || !is.null(a_hat)) "conditional" else "marginal"
  }

  if (type == "marginal") {
    return(mu)
  }

  if (is.null(cluster) && is.null(a_hat)) {
    if (is.null(newdata)) {
      cluster <- object$cluster
    } else {
      stop("For conditional predictions, provide 'cluster' or 'a_hat'.")
    }
  }

  n_pred <- length(mu)
  a_used <- rep(0, n_pred)
  if (!is.null(a_hat)) {
    if (length(a_hat) == 1) {
      a_used[] <- a_hat
    } else if (length(a_hat) == n_pred) {
      a_used <- a_hat
    } else {
      stop("Length of 'a_hat' must be 1 or equal to the number of predictions.")
    }
  } else {
    cl <- cluster
    if (length(cl) != n_pred) {
      stop("Length of 'cluster' does not match the number of rows in the data.")
    }
    cl_fac <- as.factor(cl)
    known <- intersect(levels(cl_fac), names(object$ranef))
    if (length(known) == 0) {
      warning("No matching clusters found; using random intercept equal to zero.")
    } else {
      a_map <- object$ranef[known]
      a_used <- a_map[as.character(cl_fac)]
      a_used[is.na(a_used)] <- 0
    }
  }

  wrap(mu + a_used)
}
#' Plot random-intercept directions for an angular_re fit
#'
#' Draws one arrow per cluster pointing in the direction of the estimated
#' random intercept \eqn{\hat a_i}. Arrow length is either A1(kappa_a) or 1.
#' If the 'circular' package is available, uses it for nicer axes; otherwise
#' falls back to base polar plotting.
#'
#' @param x an object of class \code{angular_re}
#' @param scale character, one of \code{"A1_kappa_a"} (default) or \code{"unit"}.
#'        Controls the common arrow length.
#' @param labels logical; draw cluster labels at arrow tips.
#' @param cex.labels numeric expansion for labels.
#' @param main plot title.
#' @param ... unused
#' @export
plot_ranef.angular_re <- function(
  x,
  scale = c("A1_kappa_a", "unit"),
  labels = FALSE,
  cex.labels = 0.7,
  main = "Random intercepts (a_hat[i])",
  ...
) {
  scale <- match.arg(scale)
  a_hat <- as.numeric(x$ranef)
  names(a_hat) <- names(x$ranef)

  A1 <- function(k) {
    num <- besselI(k, 1, expon.scaled = TRUE)
    den <- besselI(k, 0, expon.scaled = TRUE)
    r <- num / den
    r[k == 0] <- 0
    r
  }
  L <- if (scale == "unit") 1 else A1(x$kappa_a)
  L <- as.numeric(L)

  has_circ <- requireNamespace("circular", quietly = TRUE)

  if (has_circ) {
    th <- circular::circular(
      a_hat,
      type = "angles",
      units = "radians",
      modulo = "2pi",
      template = "none"
    )
    op <- par(mar = c(1.8, 1.8, 2.2, 1.8))
    on.exit(par(op))
    circular::plot.circular(
      th,
      stack = FALSE,
      shrink = 1,
      zero = 0,
      rotation = "counter",
      cex = 0,
      axes = TRUE,
      main = main
    )
    # draw unit circle
    t <- seq(0, 2 * pi, length.out = 360)
    lines(cos(t), sin(t), col = "grey75")
    # draw arrows
    xend <- L * cos(a_hat)
    yend <- L * sin(a_hat)
    arrows(0, 0, xend, yend, length = 0.08, lwd = 1.8)
    # annotate common arrow length
    mtext(
      sprintf(
        "Arrow length: %s = %.3f",
        if (scale == "unit") "1" else "A1(kappa_a)",
        L
      ),
      line = -1.2,
      cex = 0.9
    )
    if (labels) {
      text(xend * 1.04, yend * 1.04, labels = names(a_hat), cex = cex.labels)
    }
  } else {
    # base plotting fallback
    op <- par(mar = c(1.8, 1.8, 2.2, 1.8))
    on.exit(par(op))
    plot.new()
    plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)
    title(main)
    t <- seq(0, 2 * pi, length.out = 360)
    lines(cos(t), sin(t), col = "grey75")
    segments(-1, 0, 1, 0, col = "grey85")
    segments(0, -1, 0, 1, col = "grey85")
    xend <- L * cos(a_hat)
    yend <- L * sin(a_hat)
    arrows(0, 0, xend, yend, length = 0.08, lwd = 1.8)
    mtext(
      sprintf(
        "Arrow length: %s = %.3f",
        if (scale == "unit") "1" else "A1(kappa_a)",
        L
      ),
      line = -1.2,
      cex = 0.9
    )
    if (labels) {
      text(xend * 1.04, yend * 1.04, labels = names(a_hat), cex = cex.labels)
    }
  }

  invisible(x)
}
