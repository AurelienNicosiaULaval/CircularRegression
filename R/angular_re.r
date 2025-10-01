#' Angular mixed-effects regression with von Mises random intercepts
#'
#' Fits the homogeneous angular regression of Rivest et al. (2016)
#' augmented with a cluster-specific random intercept a_i ~ VM(0, kappa_a)
#' and unit-errors e_ij ~ VM(0, kappa_e), following Rivest & Kato (2019).
#'
#' @param formula Same syntax as \code{angular()}: y ~ x1:z1 + x2:z2 + ...
#'        The first angular covariate is the reference if \code{model="simplified"}.
#' @param data data.frame
#' @param cluster factor or vector of cluster IDs (same length as response)
#' @param model "simplified" (default) or "complete" (as in \code{angular()})
#' @param init list with optional entries: beta (numeric), kappa_e, kappa_a
#' @param control list: maxit (default 2000), reltol (1e-8), trace (0/1), hessian (TRUE),
#'        start_from_angular (TRUE) uses \code{angular()} for beta init.
#' @return object of class \code{"angular_re"}
#' @importFrom stats model.frame optim pnorm
#' @examples
#' # See README example below with Sandhopper data
#' @export
angular_re <- function(
  formula,
  data,
  cluster,
  model = "simplified",
  init = list(),
  control = list()
) {
  stopifnot(!missing(cluster))
  call <- match.call()
  model <- model[1]

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

  nomterms <- attr(attr(mf, "terms"), "term.labels")
  nterms <- length(nomterms)
  p <- if (model == "simplified") nterms - 1 else nterms
  nparam <- if (model == "simplified") p else p + 1
  betaname <- if (nterms > 1) nomterms[-1] else character(0)

  # Split "x:z" pairs
  sp <- strsplit(nomterms, split = ":")
  sp <- do.call(rbind, sp)
  if (model == "simplified") {
    x0 <- mf[[sp[1, 1]]]
    sp <- sp[-1, , drop = FALSE]
  }
  matx <- as.matrix(mf[, sp[, 1], drop = FALSE])
  if (ncol(sp) == 1) {
    matz <- matrix(1, nrow = nrow(matx), ncol = ncol(matx))
  } else {
    matz <- as.matrix(mf[, sp[, 2], drop = FALSE])
    # if someone wrote x:x, its z is forced to 1 as in angular()
    matz[, sp[, 2] == sp[, 1]] <- 1
  }

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
    # as in angular(): builds mean direction and its length
    sinmu <- if (model == "simplified") sin(x0) else sin(beta[p + 1])
    cosmu <- if (model == "simplified") cos(x0) else cos(beta[p + 1])
    if (p > 0) {
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta[1:p])
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta[1:p])
    }
    long <- sqrt(sinmu^2 + cosmu^2)
    mu <- atan2(sinmu, cosmu)
    list(mu = mu, long = long)
  }

  dmu_dbeta <- function(beta, mu, long) {
    # derivative matrix per obs (n x nparam) as in angular()
    D <- cbind(
      matz * sin(matx - mu),
      if (model == "simplified") NULL else cos(beta[p + 1] - mu)
    )
    D / as.vector(long)
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

  if (isTRUE(ctrl$start_from_angular) && is.null(init$beta)) {
    if (!exists("angular", mode = "function")) {
      warning("angular() not found; starting beta at zeros.")
      beta0 <- rep(0, nparam)
      if (model == "complete") beta0[p + 1] <- 0
    } else {
      # reuse angular() to get preliminary beta (iid von Mises)
      an0 <- angular(formula, data = data, model = model)
      beta0 <- as.numeric(an0$parameters[, "estimate"])
    }
  } else {
    beta0 <- init$beta %||% rep(0, nparam)
  }

  mu0 <- compute_mu(beta0)$mu
  r0 <- wrap(y - mu0)

  # moment estimators for kappas (Sec 3.1 RK2019)
  # pairwise cos over all clusters, variable ni allowed
  sum_cosdiff <- 0
  n_pairs <- 0
  for (g in seq_len(m)) {
    idx <- which(id == lev[g])
    if (length(idx) >= 2) {
      rr <- r0[idx]
      # all unordered pairs
      k <- combn(rr, 2, FUN = function(v) cos(v[1] - v[2]))
      sum_cosdiff <- sum_cosdiff + sum(k)
      n_pairs <- n_pairs + length(k)
    }
  }
  mean_cosdiff <- if (n_pairs > 0) sum_cosdiff / n_pairs else 0
  A1_ke <- sqrt(pmax(mean_cosdiff, 0))
  kappa_e0 <- init$kappa_e %||% A1inv(A1_ke)
  A1_ka <- mean(cos(r0), na.rm = TRUE) / pmax(A1_ke, 1e-8)
  kappa_a0 <- init$kappa_a %||% if (A1_ka >= 0.999999) 50 else A1inv(A1_ka)

  par0 <- c(beta0, kappa_e0, kappa_a0)
  names(par0) <- c(
    betaname,
    if (model == "complete") "intercept",
    "kappa_e",
    "kappa_a"
  )

  ## ---- loglik + score (by cluster) ----
  pack_par <- function(par) {
    list(
      beta = par[seq_len(nparam)],
      kappa_e = par[nparam + 1],
      kappa_a = par[nparam + 2]
    )
  }

  objective <- function(par) {
    pa <- pack_par(par)
    cm <- compute_mu(pa$beta)
    mu <- cm$mu
    long <- cm$long
    r <- wrap(y - mu)

    # cluster loop
    val <- 0
    for (g in seq_len(m)) {
      idx <- which(id == lev[g])
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

    for (g in seq_len(m)) {
      idx <- which(id == lev[g])
      ni <- length(idx)
      rr <- r[idx]
      Dg <- D[idx, , drop = FALSE]
      cg <- sum(cos(rr))
      sg <- sum(sin(rr))
      sigma <- sqrt((pa$kappa_a + pa$kappa_e * cg)^2 + (pa$kappa_e * sg)^2)
      A1sig <- A1(sigma)
      a_tilde <- atan2(pa$kappa_e * sg, pa$kappa_a + pa$kappa_e * cg)

      # components
      s_beta <- pa$kappa_e * A1sig * colSums(Dg * (sin(rr - a_tilde)))
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
    V_model <- tryCatch(solve(H), error = function(e) {
      matrix(NA, nrow(H), ncol(H))
    })
  } else {
    H <- matrix(NA, length(par_hat), length(par_hat))
    V_model <- H
  }

  # robust sandwich (cluster scores)
  # recompute to get scores at optimum
  sc <- attr(gradient(par_hat), "cluster_scores")
  bread <- tryCatch(solve(H), error = function(e) {
    matrix(NA, nrow(H), ncol(H))
  })
  meat <- crossprod(sc) # sum_i s_i s_i^T
  V_robust <- tryCatch(bread %*% meat %*% bread, error = function(e) {
    matrix(NA, nrow(H), ncol(H))
  })

  # extract components and diagnostics
  pa_hat <- pack_par(par_hat)
  cmh <- compute_mu(pa_hat$beta)
  mu_hat <- cmh$mu
  r_fix <- wrap(y - mu_hat)

  # random effects predictors per cluster
  a_hat <- numeric(m)
  for (g in seq_len(m)) {
    idx <- which(id == lev[g])
    rr <- r_fix[idx]
    cg <- sum(cos(rr))
    sg <- sum(sin(rr))
    a_hat[g] <- atan2(pa_hat$kappa_e * sg, pa_hat$kappa_a + pa_hat$kappa_e * cg)
  }
  names(a_hat) <- lev

  # conditional residuals
  r_cond <- r_fix
  for (g in seq_len(m))
    r_cond[id == lev[g]] <- wrap(r_fix[id == lev[g]] - a_hat[lev[g]])

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
  rownames(tab) <- c(
    betaname,
    if (model == "complete") "intercept",
    "kappa_e",
    "kappa_a"
  )

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

#' @export
ranef.angular_re <- function(object, ...) object$ranef

#' @export
fitted.angular_re <- function(object, ...) object$mu


#' Plot method for angular_re objects
#'
#' Draws rose/circular plots of residuals: fixed (y - mu) and conditional
#' (y - mu - a_hat[cluster]) as in Rivest & Kato (2019), Fig. 1–2.
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
      draw_panel_circular(r_fix, paste0(main, " — fixed"))
      draw_panel_circular(r_con, paste0(main, " — conditional"))
    } else {
      draw_panel_base(r_fix, paste0(main, " — fixed"))
      draw_panel_base(r_con, paste0(main, " — conditional"))
    }
  } else {
    ttl <- paste0(main, if (which == "fixed") " — fixed" else " — conditional")
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
#' @param type "marginal" (mu) or "conditional" (mu + a_hat). "auto" chooses
#'        "conditional" if \code{cluster} is supplied and matches known clusters, else "marginal".
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
  type = c("marginal", "conditional", "auto"),
  a_hat = NULL,
  se.fit = FALSE,
  ...
) {
  type <- match.arg(type)
  if (isTRUE(se.fit)) stop("se.fit not implemented for angular_re.")

  # retrieve model flavor used at fit time
  mdl <- tryCatch(
    as.character(object$call$model),
    error = function(e) "simplified"
  )
  if (length(mdl) == 0 || is.na(mdl)) mdl <- "simplified"
  mdl <- tolower(mdl)

  # rebuild model.frame like in angular_re()
  form <- eval(object$call$formula)
  if (is.null(newdata)) {
    # reconstruct from the original call environment
    env <- parent.frame()
    mf <- model.frame(form, data = eval(object$call$data, env))
  } else {
    mf <- model.frame(form, data = newdata)
  }

  y_dummy <- as.numeric(mf[[1]]) # not used but keeps indexing
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  sp <- do.call(rbind, strsplit(nomterms, split = ":"))

  if (mdl == "simplified") {
    x0 <- mf[[sp[1, 1]]]
    sp <- sp[-1, , drop = FALSE]
  }

  matx <- as.matrix(mf[, sp[, 1], drop = FALSE])
  if (ncol(sp) == 1) {
    matz <- matrix(1, nrow = nrow(matx), ncol = ncol(matx))
  } else {
    matz <- as.matrix(mf[, sp[, 2], drop = FALSE])
    matz[, sp[, 2] == sp[, 1]] <- 1
  }

  beta <- object$beta
  wrap <- function(a) atan2(sin(a), cos(a))

  compute_mu <- function(beta) {
    sinmu <- if (mdl == "simplified") sin(x0) else sin(beta[length(beta)])
    cosmu <- if (mdl == "simplified") cos(x0) else cos(beta[length(beta)])
    p <- if (mdl == "simplified") length(beta) else length(beta) - 1
    if (p > 0) {
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta[1:p])
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta[1:p])
    }
    atan2(sinmu, cosmu)
  }

  mu <- compute_mu(beta)

  # Decide conditional vs marginal
  if (type == "auto") {
    type <- if (!is.null(cluster)) "conditional" else "marginal"
  }

  if (type == "marginal") {
    return(mu)
  }

  # conditional
  if (is.null(cluster) && is.null(a_hat)) {
    stop("For conditional predictions, provide 'cluster' or 'a_hat'.")
  }

  # build a_used for each row
  a_used <- rep(0, length(mu))
  if (!is.null(a_hat)) {
    if (length(a_hat) == 1) a_used[] <- a_hat else if (
      length(a_hat) == length(mu)
    )
      a_used <- a_hat else
      stop("Length of 'a_hat' must be 1 or equal to nrow(newdata).")
  } else {
    # use object's ranef if cluster matches
    cl <- as.factor(cluster)
    known <- intersect(levels(cl), names(object$ranef))
    if (length(known) == 0) {
      warning("No matching clusters found in the fit; using a = 0 for all.")
    } else {
      a_map <- object$ranef[known]
      a_used <- a_map[as.character(cl)]
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
  main = "Random intercepts (â_i)",
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
    # cercle unité
    t <- seq(0, 2 * pi, length.out = 360)
    lines(cos(t), sin(t), col = "grey75")
    # flèches
    xend <- L * cos(a_hat)
    yend <- L * sin(a_hat)
    arrows(0, 0, xend, yend, length = 0.08, lwd = 1.8)
    # étiquette de la longueur commune
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
