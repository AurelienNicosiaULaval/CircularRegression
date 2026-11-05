.angular_re_formula <- function(formula) {
  group <- NULL
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    return(list(fixed = formula, group = group))
  }

  fixed_terms <- function(term) {
    if (!is.call(term)) return(term)
    operator <- term[[1L]]
    if (identical(operator, as.name("("))) return(fixed_terms(term[[2L]]))
    if (identical(operator, as.name("+")) && length(term) == 3L) {
      lhs <- fixed_terms(term[[2L]])
      rhs <- fixed_terms(term[[3L]])
      if (is.null(lhs)) return(rhs)
      if (is.null(rhs)) return(lhs)
      return(call("+", lhs, rhs))
    }
    if (identical(operator, as.name("|")) && length(term) == 3L) {
      if (!is.null(group)) {
        stop("Supply only one grouping term in 'formula'.", call. = FALSE)
      }
      intercept_only <- identical(term[[2L]], 1) || identical(term[[2L]], 1L)
      if (!intercept_only || !is.symbol(term[[3L]])) {
        stop("The grouping term must have the form '(1 | group)'.", call. = FALSE)
      }
      group <<- as.character(term[[3L]])
      return(NULL)
    }
    if (contains_grouping_call(term)) {
      stop("Add the grouping term as '+ (1 | group)'.", call. = FALSE)
    }
    term
  }

  contains_grouping_call <- function(term) {
    is.call(term) && (
      identical(term[[1L]], as.name("|")) || identical(term[[1L]], as.name("||")) ||
        any(vapply(as.list(term)[-1L], contains_grouping_call, logical(1)))
    )
  }

  rhs <- fixed_terms(formula[[3L]])
  if (is.null(group)) return(list(fixed = formula, group = group))
  if (is.null(rhs) || "." %in% all.vars(rhs)) {
    stop("Write the fixed direction terms explicitly with '(1 | group)'.", call. = FALSE)
  }
  fixed <- stats::as.formula(call("~", formula[[2L]], rhs), env = environment(formula))
  list(fixed = fixed, group = group)
}

#' Angular mixed-effects regression with von Mises random intercepts
#'
#' Fits the homogeneous angular regression of Rivest et al. (2016) augmented
#' with a cluster-specific random intercept \eqn{a_i \sim \mathrm{VM}(0, \kappa_a)}
#' and unit-errors \eqn{e_{ij} \sim \mathrm{VM}(0, \kappa_e)}, following Rivest
#' & Kato (2019). The fixed-effects component shares the simplified syntax of
#' \code{angular()}, where the first term provides the reference direction and
#' subsequent terms are written \code{x:z}.
#'
#' @param formula Model formula with a circular response and terms of the form
#'   \code{x} or \code{x:z}. The first angular covariate supplies the reference
#'   direction whose coefficient is fixed to one. Add a single \code{(1 | group)}
#'   term for one random intercept per level of the grouping variable.
#' @param data A data frame containing the variables used in \code{formula}.
#' @param cluster Optional factor or vector of cluster identifiers, supplied
#'   in place of a grouping term in \code{formula}. It must have length equal
#'   to \code{nrow(data)} or to the model-frame size after missing-value handling.
#' @param init Optional list with components \code{beta}, \code{kappa_e}, and
#'   \code{kappa_a} providing starting values. Missing entries are filled with
#'   moment-based estimates.
#' @param control Optional list controlling the optimisation with components
#'   \code{maxit} (default 2000), \code{reltol} (1e-8), \code{trace} (0/1),
#'   \code{hessian} (\code{TRUE}), and \code{start_from_angular} (\code{TRUE})
#'   indicating whether \code{angular()} should be used to initialise the fixed
#'   effects.
#' @param na.action Function used by \code{\link[stats]{model.frame}} to handle
#'   missing values, including missing group identifiers in the formula.
#'   The default is \code{\link[stats]{na.omit}}.
#' @return An object of class \code{"angular_re"} containing fixed-effect
#'   estimates, model-based and cluster-robust covariance matrices, concentration
#'   estimates, fitted marginal mean directions, fixed and conditional
#'   residuals, and estimated random intercepts.
#' @details The reported log-likelihood includes all normalizing constants.
#'   The object records the projected \code{score_norm}, \code{convergence},
#'   \code{boundary}, observed \code{information}, and \code{cluster_scores}.
#'   Regular Wald covariances are unavailable at a zero-concentration boundary.
#'   Cluster-robust inference assumes independent clusters and sufficiently many
#'   clusters; it has no small-sample correction. Singleton-only groups cannot
#'   distinguish the two concentrations and are rejected. The independent
#'   limit as the random-intercept concentration tends to infinity is
#'   evaluated at the fitted fixed effects and error concentration.
#'   If it improves the likelihood, \code{fit_status = "independence_limit"}
#'   and the finite fit is retained for diagnosis with Wald inference unavailable.
#'   Perfect within-group alignment gives an unbounded likelihood as the error
#'   concentration increases. Such fits have
#'   \code{fit_status = "unbounded_likelihood"}, \code{convergence = FALSE}, and
#'   unavailable covariance matrices. If alignment is detected at the initial
#'   coefficients, optimization is not run and those coefficients are returned
#'   for diagnosis only. \code{inference_available} records availability of
#'   regular Wald inference. After an objective-based optimizer stop, the fit
#'   is refined when necessary to meet its reported score tolerance; this is
#'   recorded in \code{optimization_polished}.
#' @references Rivest, L.-P., Duchesne, T., Nicosia, A., and Fortin, D. (2016).
#'   A general angular regression model for the analysis of data on animal
#'   movement in ecology. \emph{Journal of the Royal Statistical Society:
#'   Series C (Applied Statistics)}, 65(3), 445-463.
#'
#'   Rivest, L.-P., and Kato, S. (2019). A random-effects model for clustered
#'   circular data. \emph{Canadian Journal of Statistics}, 47(4), 712-728.
#' @importFrom stats model.frame optim pnorm
#' @examples
#' data(sandhopper_trials)
#'
#' fit <- angular_re(
#'   escape ~ sun + wind:wind_speed + (1 | Anim),
#'   data = sandhopper_trials
#' )
#' coef(fit)
#' predict(fit, newdata = head(sandhopper_trials), type = "conditional")
#' @export
angular_re <- function(
  formula,
  data,
  cluster = NULL,
  init = list(),
  control = list(),
  na.action = stats::na.omit
) {
  call <- match.call()
  model_formula <- .angular_re_formula(formula)
  fixed_formula <- model_formula$fixed
  group <- model_formula$group
  if (is.null(group) && is.null(cluster)) {
    stop("Add '(1 | group)' to 'formula' or supply 'cluster'.", call. = FALSE)
  }
  if (!is.null(group)) {
    if (!is.null(cluster)) {
      stop("Supply a grouping term in 'formula' or 'cluster', but not both.", call. = FALSE)
    }
    if (is.null(data[[group]])) {
      stop(sprintf("Grouping variable '%s' was not found in 'data'.", group), call. = FALSE)
    }
  }

  ## ---- Shared direction/modifier parsing and complete-case alignment ----
  des <- .angular_design(fixed_formula, data, reference = "first",
                         na.action = na.action, group = group)
  mf <- des$mf
  y <- des$y
  n <- des$nobs
  cluster_mf <- if (is.null(group)) {
    .align_model_vector(as.vector(cluster), mf, "cluster")
  } else {
    mf[[ncol(mf)]]
  }
  if (!is.null(dim(cluster_mf))) {
    stop("Grouping variable must be a vector or factor.", call. = FALSE)
  }
  if (anyNA(cluster_mf)) {
    stop("cluster contains missing values after applying model-frame NA handling.", call. = FALSE)
  }

  id <- as.factor(as.vector(cluster_mf))
  lev <- levels(id)
  m <- length(lev)
  cluster_index <- split(seq_len(n), id)
  if (m < 2L || all(lengths(cluster_index) == 1L)) {
    stop("Random-intercept estimation requires at least two clusters and repeated observations within clusters.", call. = FALSE)
  }

  if (des$ref_idx != 1L) {
    stop("The first term must be a plain reference direction without a modifier.", call. = FALSE)
  }
  x0 <- des$x0
  matx <- des$matx
  matz <- des$matz
  betaname <- des$betaname
  ref_name <- des$reference$ref_name
  term_labels <- des$term_labels
  p <- length(betaname)
  nparam <- p

  ## ---- helpers ----
  wrap <- function(a) atan2(sin(a), cos(a))
  # Smooth signed extension is needed for numerical Hessians near the lower
  # bound. The statistical fit itself restricts both concentrations to >= 0.
  A1 <- function(k) sign(k) * .A1_ratio(abs(k))
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
        df <- 1 - A1(k) / pmax(k, .Machine$double.eps) - A1(k)^2
        k <- pmax(0, k - f / df)
      }
      return(k)
    }
  }
  logI0 <- function(k) .logI0(abs(k))

  compute_mu <- function(beta) {
    sinmu <- sin(x0)
    cosmu <- cos(x0)
    if (p > 0) {
      sinmu <- sinmu + as.vector((matz * sin(matx)) %*% beta)
      cosmu <- cosmu + as.vector((matz * cos(matx)) %*% beta)
    }
    long <- sqrt(sinmu^2 + cosmu^2)
    scale <- 1 + if (p > 0L) as.vector(abs(matz) %*% abs(beta)) else 0
    if (any(long <= sqrt(.Machine$double.eps) * scale)) {
      stop("Zero or numerically singular resultant in the grouped mean model.", call. = FALSE)
    }
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
  if (!is.list(control) || !is.list(init)) stop("'control' and 'init' must be lists.", call. = FALSE)
  valid_control <- c("maxit", "reltol", "trace", "hessian", "start_from_angular")
  valid_init <- c("beta", "kappa_e", "kappa_a")
  for (input in c("control", "init")) {
    value <- if (input == "control") control else init
    allowed <- if (input == "control") valid_control else valid_init
    if (length(value) && (is.null(names(value)) || anyDuplicated(names(value)) ||
        any(!names(value) %in% allowed) || any(vapply(value, is.null, logical(1))))) {
      stop(sprintf("Unknown, duplicated, unnamed or NULL component in '%s'.", input), call. = FALSE)
    }
  }
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
  if (!is.numeric(ctrl$maxit) || length(ctrl$maxit) != 1L || !is.finite(ctrl$maxit) || ctrl$maxit < 0) {
    stop("'control$maxit' must be a finite non-negative scalar.", call. = FALSE)
  }
  if (ctrl$maxit != floor(ctrl$maxit) || ctrl$maxit >= .Machine$integer.max) {
    stop("'control$maxit' must be an integer below .Machine$integer.max.", call. = FALSE)
  }
  if (!is.numeric(ctrl$reltol) || length(ctrl$reltol) != 1L || !is.finite(ctrl$reltol) || ctrl$reltol <= 0) {
    stop("'control$reltol' must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.numeric(ctrl$trace) || length(ctrl$trace) != 1L || !is.finite(ctrl$trace) || ctrl$trace < 0) {
    stop("'control$trace' must be a finite non-negative scalar.", call. = FALSE)
  }
  ctrl$trace <- as.integer(ctrl$trace)
  if (!is.logical(ctrl$hessian) || length(ctrl$hessian) != 1L || is.na(ctrl$hessian)) {
    stop("'control$hessian' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(ctrl$start_from_angular) || length(ctrl$start_from_angular) != 1L || is.na(ctrl$start_from_angular)) {
    stop("'control$start_from_angular' must be TRUE or FALSE.", call. = FALSE)
  }

  coalesce_null <- function(x, default) if (is.null(x)) default else x
  if (!is.null(init$beta) && (!is.numeric(init$beta) || length(init$beta) != p ||
      any(!is.finite(init$beta)))) {
    stop(sprintf("'init$beta' must have length %d and finite numeric values.", p), call. = FALSE)
  }

  if (p == 0) {
    beta0 <- numeric(0)
  } else if (isTRUE(ctrl$start_from_angular) && is.null(init$beta)) {
    if (!exists("angular", mode = "function")) {
      warning("angular() not found; starting beta at zeros.")
      beta0 <- rep(0, p)
    } else {
      initial_na_action <- na.action
      if (!is.null(group)) {
        # Use the same rows, including omissions due to missing groups.
        initial_na_action <- function(frame) {
          frame[match(rownames(mf), rownames(frame)), , drop = FALSE]
        }
      }
      an0 <- angular(fixed_formula, data = data, reference = "first",
                     na.action = initial_na_action)
      beta0 <- as.numeric(an0$parameters[betaname, "estimate"])
      beta0[is.na(beta0)] <- 0
    }
  } else {
    beta0 <- coalesce_null(init$beta, rep(0, p))
    if (!is.numeric(beta0) || any(!is.finite(beta0)) || length(beta0) != p) {
      stop(sprintf("'init$beta' must have length %d", p), call. = FALSE)
    }
  }

  mu0 <- compute_mu(beta0)$mu
  r0 <- wrap(y - mu0)
  perfect_within <- function(residual) {
    all(vapply(cluster_index, function(idx) {
      max(abs(wrap(residual[idx] - residual[idx[1L]]))) <=
        100 * .Machine$double.eps
    }, logical(1)))
  }
  initial_alignment_boundary <- perfect_within(r0)

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
  if (!is.null(init$kappa_e) && (!.is_scalar_number(kappa_e0) || kappa_e0 <= 0)) {
    stop("'init$kappa_e' must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.null(init$kappa_a) && (!.is_scalar_number(kappa_a0) || kappa_a0 <= 0)) {
    stop("'init$kappa_a' must be a finite positive scalar.", call. = FALSE)
  }

  kappa_floor <- sqrt(.Machine$double.eps)
  par0 <- c(beta0, max(as.numeric(kappa_e0), kappa_floor), max(as.numeric(kappa_a0), kappa_floor))
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
    val <- -n * log(2 * pi)
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
  optimize_group <- function(start, factr, pgtol = 0) optim(
    start,
    fn = function(p) -objective(p),
    gr = function(p) -gradient(p),
    method = "L-BFGS-B",
    lower = c(rep(-Inf, nparam), kappa_floor, kappa_floor),
    control = list(
      maxit = ctrl$maxit,
      factr = factr,
      pgtol = pgtol,
      trace = ctrl$trace
    ),
    hessian = isTRUE(ctrl$hessian)
  )

  if (initial_alignment_boundary) {
    opt <- list(
      par = par0, value = -objective(par0), convergence = 1L,
      message = "Optimization not run: within-group alignment gives an unbounded likelihood.",
      hessian = matrix(NA_real_, length(par0), length(par0))
    )
  } else {
    opt <- optimize_group(par0, ctrl$reltol / .Machine$double.eps)
  }

  projected_gradient <- function(par) {
    score <- as.numeric(gradient(par))
    active <- seq_along(par) > nparam & par <= 10 * kappa_floor
    score[active] <- pmax(score[active], 0)
    score
  }
  polished <- FALSE
  if (!initial_alignment_boundary && ctrl$maxit > 0L &&
      opt$convergence == 0L &&
      max(abs(projected_gradient(opt$par))) / n > sqrt(ctrl$reltol)) {
    # An objective-based stop can precede the advertised score tolerance.
    refined <- optimize_group(opt$par, factr = 10, pgtol = n * sqrt(ctrl$reltol))
    if (is.finite(refined$value) && refined$value <= opt$value) {
      opt <- refined
      polished <- TRUE
    }
  }

  par_hat <- opt$par
  val <- -opt$value
  score_hat <- gradient(par_hat)
  projected_score <- projected_gradient(par_hat)
  lower_active <- seq_along(par_hat) > nparam & par_hat <= 10 * kappa_floor
  projected_score[lower_active] <- pmax(projected_score[lower_active], 0)
  score_norm <- max(abs(projected_score)) / n
  conv <- (opt$convergence == 0 && score_norm <= sqrt(ctrl$reltol))
  pa_hat <- pack_par(par_hat)
  mu_hat <- compute_mu(pa_hat$beta)$mu
  r_fix <- wrap(y - mu_hat)
  alignment_boundary <- initial_alignment_boundary || perfect_within(r_fix)
  independence_logLik <- pa_hat$kappa_e * sum(cos(r_fix)) -
    n * (log(2 * pi) + logI0(pa_hat$kappa_e))
  independence_boundary <- conv && !alignment_boundary && !any(lower_active) &&
    independence_logLik > val + ctrl$reltol * (1 + abs(val))
  boundary <- any(lower_active) || alignment_boundary || independence_boundary
  if (alignment_boundary) {
    conv <- FALSE
    warning("Perfect within-group alignment gives an unbounded likelihood; no finite maximum exists and Wald inference is unavailable.", call. = FALSE)
  } else if (independence_boundary) {
    conv <- FALSE
    warning("The independent limit has higher likelihood than the finite grouped fit. Wald inference is unavailable.", call. = FALSE)
  } else if (!conv && ctrl$maxit > 0L) {
    warning("Grouped model convergence was not established; inspect score_norm and optim_message.", call. = FALSE)
  }
  if (any(lower_active)) warning("A concentration is on the zero boundary; regular Wald inference is unavailable.", call. = FALSE)

  # observed Fisher (negative Hessian)
  if (isTRUE(ctrl$hessian) && !alignment_boundary) {
    H <- opt$hessian
    V_model <- invert_hessian(H)
  } else {
    H <- matrix(NA_real_, length(par_hat), length(par_hat))
    V_model <- H
  }

  inference_available <- conv && !boundary && all(is.finite(V_model))
  if (!inference_available) V_model[,] <- NA_real_

  # robust sandwich (cluster scores)
  sc <- attr(gradient(par_hat), "cluster_scores")
  if (isTRUE(ctrl$hessian)) {
    bread <- V_model
    if (all(is.finite(bread))) {
      meat <- crossprod(sc)
      V_robust <- bread %*% meat %*% bread
    } else {
      V_robust <- matrix(NA_real_, nrow(H), ncol(H))
    }
  } else {
    V_robust <- matrix(NA_real_, nrow(H), ncol(H))
  }
  parameter_names <- c(betaname, "kappa_e", "kappa_a")
  V_model <- .named_covariance(V_model, parameter_names)
  V_robust <- .named_covariance(V_robust, parameter_names)
  fit_status <- if (alignment_boundary) "unbounded_likelihood" else if (independence_boundary) {
    "independence_limit"
  } else if (any(lower_active)) {
    "zero_concentration"
  } else if (!conv) "not_converged" else if (isTRUE(ctrl$hessian) && !inference_available) {
    "singular_information"
  } else "converged"

  # extract components and diagnostics
  names(pa_hat$beta) <- betaname

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
  if (alignment_boundary) rho_SS <- NA_real_

  out <- list(
    call = call,
    formula = formula,
    fixed_formula = fixed_formula,
    group_name = group,
    model_data = mf,
    logLik = val,
    independence_logLik = independence_logLik,
    convergence = conv,
    fit_status = fit_status,
    inference_available = inference_available,
    alignment_boundary = alignment_boundary,
    optimization_polished = polished,
    score = as.numeric(score_hat),
    score_norm = score_norm,
    cluster_scores = sc,
    information = H,
    boundary = boundary,
    optim_message = opt$message,
    na.action = attr(mf, "na.action"),
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
    control = ctrl,
    ref_name = ref_name,
    term_labels = term_labels
  )
  class(out) <- "angular_re"
  out
}

#' Methods for Angular Random-Effects Fits
#'
#' @param x,object An object of class \code{"angular_re"}.
#' @param robust Logical. For \code{vcov()}, return the cluster-robust
#'   covariance matrix when \code{TRUE}; otherwise return the model-based
#'   covariance matrix.
#' @param type Residual type. \code{"fixed"} returns residuals from the marginal
#'   fixed-effect mean direction. \code{"conditional"} subtracts the estimated
#'   cluster random intercepts.
#' @param ... Additional arguments, currently unused.
#' @return \code{coef()} returns the parameter estimates. \code{vcov()} returns
#'   the requested covariance matrix. \code{fitted()} returns marginal fitted
#'   mean directions. \code{residuals()} returns wrapped circular residuals.
#'   \code{summary()} returns a list with convergence information, likelihood,
#'   sample sizes, concentration estimates, residual summaries and coefficient
#'   table. \code{print()} invisibly returns the input object.
#' @rdname angular_re-methods
#' @examples
#' data(Sandhopper)
#' sh <- Sandhopper[1:24, ]
#' d <- sh[rep(1:24, each = 5), ]
#' d$y <- as.vector(t(as.matrix(sh[paste0("LN", 1:5)]))) * pi / 180
#' d$ref <- d$Azimuth * pi / 180
#' fit <- angular_re(y ~ ref + (1 | Anim), data = d)
#' summary(fit)
#' vcov(fit, robust = TRUE)
#' @export
print.angular_re <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  .print_fit_status(x)
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

#' @rdname angular_re-methods
#' @export
summary.angular_re <- function(object, ...) {
  coef_tab <- object$coefficients
  model_se <- coef_tab[, "se_model"]
  robust_se <- coef_tab[, "se_robust"]
  model_ok <- all(is.finite(model_se))
  robust_ok <- all(is.finite(robust_se))
  se_type <- if (model_ok && robust_ok) {
    "model-based and cluster-robust"
  } else if (model_ok) {
    "model-based only; cluster-robust SE contain NA"
  } else if (robust_ok) {
    "cluster-robust only; model-based SE contain NA"
  } else {
    "not available; covariance matrices contain NA"
  }

  out <- list(
    call = object$call,
    logLik = object$logLik,
    convergence = object$convergence,
    fit_status = object$fit_status,
    inference_available = object$inference_available,
    nobs = length(object$mu),
    n_clusters = length(object$ranef),
    coefficients = coef_tab,
    kappa_e = object$kappa_e,
    kappa_a = object$kappa_a,
    rho_SS = object$rho_SS,
    se_type = se_type,
    se_model_available = model_ok,
    se_robust_available = robust_ok,
    residuals_fixed = summary(object$resid_fixed),
    residuals_conditional = summary(object$resid_conditional)
  )
  class(out) <- "summary.angular_re"
  out
}

#' @rdname angular_re-methods
#' @export
print.summary.angular_re <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  .print_fit_status(x)
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = 5),
    "\nConverged:",
    x$convergence,
    "\nNumber of observations:",
    x$nobs,
    "\nNumber of clusters:",
    x$n_clusters,
    "\n"
  )
  cat(
    "\nConcentration estimates:\n",
    "kappa_e:",
    format(x$kappa_e, digits = 5),
    "\n kappa_a:",
    format(x$kappa_a, digits = 5),
    "\n rho_SS:",
    format(x$rho_SS, digits = 5),
    "\n",
    sep = ""
  )
  cat("\nStandard errors:", x$se_type, "\n")
  if (!isTRUE(x$se_robust_available)) {
    cat("Cluster-robust SE contain NA; inspect vcov(object, robust = TRUE).\n")
  }
  if (!isTRUE(x$se_model_available)) {
    cat("Model-based SE contain NA; inspect vcov(object).\n")
  }
  cat("\nFixed-effect and concentration parameters:\n")
  stats::printCoefmat(
    x$coefficients,
    P.values = TRUE,
    has.Pvalue = TRUE,
    ...
  )
  cat("\nFixed residuals:\n")
  print(x$residuals_fixed)
  cat("\nConditional residuals:\n")
  print(x$residuals_conditional)
  invisible(x)
}

#' @rdname angular_re-methods
#' @export
coef.angular_re <- function(object, ...) {
  object$coefficients[, "estimate"]
}

#' @rdname angular_re-methods
#' @export
vcov.angular_re <- function(object, robust = FALSE, ...) {
  .check_logical_scalar(robust, "robust")
  if (robust) object$vcov_robust else object$vcov_model
}

#' @rdname angular_re-methods
#' @export
residuals.angular_re <- function(
  object,
  type = c("fixed", "conditional"),
  ...
) {
  type <- match.arg(type)
  value <- if (type == "fixed") object$resid_fixed else object$resid_conditional
  .restore_excluded(object, value, residual = TRUE)
}

#' Extract estimated random intercepts from an angular mixed model
#'
#' Returns the cluster-specific random intercepts \eqn{\hat a_i} estimated by
#' \code{angular_re()}.
#'
#' @param object A fitted \code{angular_re} model.
#' @param ... Unused; included for compatibility with the generic.
#' @return A named numeric vector of random intercept estimates, one per cluster.
#' @examples
#' data(Sandhopper)
#' sh <- Sandhopper[1:24, ]
#' d <- sh[rep(1:24, each = 5), ]
#' d$y <- as.vector(t(as.matrix(sh[paste0("LN", 1:5)]))) * pi / 180
#' d$ref <- d$Azimuth * pi / 180
#' fit <- angular_re(y ~ ref + (1 | Anim), data = d)
#' head(ranef.angular_re(fit))
#' @export
ranef.angular_re <- function(object, ...) object$ranef

#' @rdname angular_re-methods
#' @export
fitted.angular_re <- function(object, ...) .restore_excluded(object, object$mu)

.angular_re_residual_plots <- function(
  x,
  breaks,
  bandwidth,
  points,
  main
) {
  if (!is.numeric(breaks) || length(breaks) != 1L || !is.finite(breaks) ||
      breaks < 2 || breaks != as.integer(breaks)) {
    stop("'breaks' must be a single integer greater than one.", call. = FALSE)
  }
  breaks <- as.integer(breaks)
  bandwidth <- .shared_circular_bandwidth(
    c(x$resid_fixed, x$resid_conditional),
    bandwidth = bandwidth
  )

  fixed_plot <- .residual_rose_plot(
    residuals = x$resid_fixed,
    bins = breaks,
    bandwidth = bandwidth,
    rug = points,
    title = sprintf("%s: fixed residuals", main),
    subtitle = "Includes between-cluster variation"
  )
  conditional_plot <- .residual_rose_plot(
    residuals = x$resid_conditional,
    kappa = x$kappa_e,
    bins = breaks,
    bandwidth = bandwidth,
    rug = points,
    title = sprintf("%s: conditional residuals", main),
    subtitle = "Fitted von Mises bin expectation"
  )
  radial_limit <- .shared_rose_radial_limit(list(fixed_plot, conditional_plot))

  list(
    fixed = .apply_rose_radial_limit(fixed_plot, radial_limit),
    conditional = .apply_rose_radial_limit(conditional_plot, radial_limit),
    settings = list(
      bins = breaks,
      bandwidth = bandwidth,
      radial_limit = radial_limit
    )
  )
}


#' Plot method for angular_re objects
#'
#' Draws circular residual diagnostics for a random-intercept angular model.
#'
#' The fixed residuals, \eqn{y - \hat\mu}, retain between-cluster variation.
#' The conditional residuals, \eqn{y - \hat\mu - \hat a_{g}}, are compared
#' with the fitted von Mises bin expectation. The paired rose plots follow the
#' diagnostic distinction illustrated by Rivest and Kato (2019, Figures 1--2).
#' When \code{bandwidth = NULL}, one kernel-density bandwidth is selected from
#' the pooled fixed and conditional residuals. Both panels use that bandwidth,
#' identical angular bins, and identical radial limits, including when they are
#' requested in separate calls.
#'
#' @param x an object of class \code{angular_re}
#' @param which one of "both", "fixed", "conditional"
#' @param breaks number of bins for the rose histogram.
#' @param bandwidth optional finite positive kernel-density bandwidth. If
#'   \code{NULL}, a common bandwidth is selected from the pooled residuals.
#' @param points logical; add individual observations as circular rug marks.
#' @param main base title (panel-specific suffixes are added)
#' @param ... unused
#' @examples
#' data(Sandhopper)
#' sh <- Sandhopper[1:24, ]
#' d <- sh[rep(1:24, each = 5), ]
#' d$y <- as.vector(t(as.matrix(sh[paste0("LN", 1:5)]))) * pi / 180
#' d$ref <- d$Azimuth * pi / 180
#' fit <- angular_re(y ~ ref + (1 | Anim), data = d)
#' if (all(vapply(c("ggplot2", "ggcircular", "gridExtra"),
#'                requireNamespace, logical(1), quietly = TRUE))) plot(fit)
#' @export
plot.angular_re <- function(
  x,
  which = c("both", "fixed", "conditional"),
  breaks = 24,
  bandwidth = NULL,
  points = TRUE,
  main = "Residuals",
  ...
) {
  which <- match.arg(which)
  .require_circular_diagnostic_packages(arrange = identical(which, "both"))
  panels <- .angular_re_residual_plots(
    x = x,
    breaks = breaks,
    bandwidth = bandwidth,
    points = points,
    main = main
  )

  if (which == "both") {
    comparison_note <- sprintf(
      "Common comparison settings: %d bins; bandwidth %.3f; radial limit %.3f",
      panels$settings$bins,
      panels$settings$bandwidth,
      panels$settings$radial_limit
    )
    gridExtra::grid.arrange(
      panels$fixed,
      panels$conditional,
      ncol = 2,
      top = grid::textGrob(comparison_note, gp = grid::gpar(fontsize = 9))
    )
  } else {
    print(panels[[which]])
  }

  invisible(x)
}


#' Predict method for angular_re objects
#'
#' @param object an \code{angular_re} fit.
#' @param newdata Data frame, list or environment containing all fixed-effect
#'   predictors. If \code{NULL}, use the original model frame.
#' @param cluster Optional vector or factor of cluster identifiers. For
#'   conditional prediction, the grouping variable in \code{newdata} is used
#'   by default if it was specified in the fitted formula.
#' @param type One of \code{"auto"} (default), \code{"marginal"} (mean
#'        direction) or \code{"conditional"} (mean direction plus the predicted
#'        random intercept). With \code{"auto"} the conditional mean is returned
#'        when \code{cluster} or \code{a_hat} is supplied. Otherwise the
#'        marginal mean is used.
#' @param a_hat Optional numeric vector of random intercepts for conditional
#'        prediction, of length one or one value per prediction.
#' @param se.fit not implemented (reserved).
#' @param ... unused
#' @return numeric vector (angles in radians).
#' @details New-data rows with missing predictors are rejected rather than
#'   silently dropped. New cluster identifiers use a zero random intercept with
#'   a warning; this is a population-mean fallback, not an estimated individual
#'   effect. At zero marginal concentration the fixed component is an orientation
#'   parameter only, not an identified marginal mean direction.
#' @examples
#' data(Sandhopper)
#' sh <- Sandhopper[1:24, ]
#' d <- sh[rep(1:24, each = 5), ]
#' d$y <- as.vector(t(as.matrix(sh[paste0("LN", 1:5)]))) * pi / 180
#' d$ref <- d$Azimuth * pi / 180
#' fit <- angular_re(y ~ ref + (1 | Anim), data = d)
#' predict(fit, newdata = head(d), type = "conditional")
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

  form <- object$fixed_formula %||% object$formula

  parse_design <- function(formula, data) {
    info <- .parse_term_labels(object$term_labels, formula)
    design <- .new_term_design(info, data)
    list(x0 = design$matx[, 1L],
         matx = design$matx[, -1L, drop = FALSE],
         matz = design$matz[, -1L, drop = FALSE],
         betaname = object$betaname, nobs = nrow(design$matx))
  }

  data_source <- if (is.null(newdata)) {
    object$model_data
  } else {
    newdata
  }
  if (is.null(data_source)) {
    stop("No data available for prediction.", call. = FALSE)
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
    .mean_direction_from_terms(c(1, beta), cbind(x0, matx), cbind(1, matz))
  }

  mu <- compute_mu_new(beta, design$x0, design$matx, design$matz)
  wrap <- function(a) atan2(sin(a), cos(a))

  if (type == "auto") {
    type <- if (!is.null(cluster) || !is.null(a_hat)) "conditional" else "marginal"
  }

  if (type == "marginal") {
    return(if (is.null(newdata)) .restore_excluded(object, mu) else mu)
  }

  if (is.null(cluster) && is.null(a_hat)) {
    if (is.null(newdata)) {
      cluster <- object$cluster
    } else {
      group <- object$group_name
      if (!is.null(group)) cluster <- newdata[[group]]
      if (is.null(cluster)) {
        message <- if (is.null(group)) {
          "For conditional predictions, provide 'cluster' or 'a_hat'."
        } else {
          sprintf("Include '%s' in 'newdata' or supply 'cluster' or 'a_hat'.", group)
        }
        stop(message, call. = FALSE)
      }
    }
  }

  n_pred <- length(mu)
  a_used <- rep(0, n_pred)
  if (!is.null(a_hat)) {
    if (!is.numeric(a_hat) || any(!is.finite(a_hat))) stop("'a_hat' must be numeric and finite.", call. = FALSE)
    if (length(a_hat) == 1) {
      a_used[] <- a_hat
    } else if (length(a_hat) == n_pred) {
      a_used <- a_hat
    } else {
      stop("Length of 'a_hat' must be 1 or equal to the number of predictions.")
    }
  } else {
    cl <- cluster
    if (anyNA(cl)) stop("Prediction cluster identifiers must not be missing.", call. = FALSE)
    if (length(cl) != n_pred) {
      stop("Length of 'cluster' does not match the number of rows in the data.")
    }
    cl_fac <- as.factor(cl)
    known <- intersect(levels(cl_fac), names(object$ranef))
    if (length(known) == 0) {
      warning("No matching clusters found; using random intercept equal to zero.")
    } else {
      if (any(!as.character(cl_fac) %in% known)) {
        warning("Some clusters are new; using random intercept equal to zero for them.", call. = FALSE)
      }
      a_map <- object$ranef[known]
      a_used <- a_map[as.character(cl_fac)]
      a_used[is.na(a_used)] <- 0
    }
  }

  result <- as.numeric(wrap(mu + a_used))
  if (is.null(newdata)) .restore_excluded(object, result) else result
}
#' Plot random-intercept directions for an angular_re fit
#'
#' Draws one arrow per cluster pointing in the direction of the estimated
#' random intercept \eqn{\hat a_i}. Arrow length is either A1(kappa_a) or 1.
#' Uses \pkg{ggcircular} when available, with a \pkg{circular} fallback.
#'
#' @param x an object of class \code{angular_re}
#' @param scale character, one of \code{"A1_kappa_a"} (default) or \code{"unit"}.
#'        Controls the common arrow length.
#' @param labels logical; draw cluster labels at arrow tips.
#' @param cex.labels numeric expansion for labels.
#' @param main plot title.
#' @param ... unused
#' @examples
#' data(Sandhopper)
#' sh <- Sandhopper[1:24, ]
#' d <- sh[rep(1:24, each = 5), ]
#' d$y <- as.vector(t(as.matrix(sh[paste0("LN", 1:5)]))) * pi / 180
#' d$ref <- d$Azimuth * pi / 180
#' fit <- angular_re(y ~ ref + (1 | Anim), data = d)
#' plot_ranef.angular_re(fit)
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

  L <- if (scale == "unit") 1 else .A1_ratio(x$kappa_a)

  if (requireNamespace("ggcircular", quietly = TRUE)) {
    arrows_data <- data.frame(angle = a_hat, group = names(a_hat),
                              tip_x = 1.04 * L * cos(a_hat), tip_y = 1.04 * L * sin(a_hat))
    angle <- seq(0, 2 * pi, length.out = 361)
    circle <- data.frame(x = cos(angle), y = sin(angle))
    ticks <- data.frame(x = c(1.12, 0, -1.12, 0), y = c(0, 1.12, 0, -1.12),
                        label = c("0", "pi/2", "pi", "-pi/2"))
    plot <- ggplot2::ggplot() +
      ggplot2::geom_path(data = circle, ggplot2::aes(x = .data$x, y = .data$y),
                         colour = "grey65") +
      ggcircular::geom_direction_arrow(
        data = arrows_data, ggplot2::aes(x = 0, y = 0, angle = .data$angle),
        length = L, alpha = 0.55, colour = "#1F5A7A"
      ) +
      ggplot2::geom_text(data = ticks, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                         size = 3) +
      ggplot2::coord_equal(xlim = c(-1.22, 1.22), ylim = c(-1.22, 1.22)) +
      ggplot2::labs(title = main, subtitle = sprintf("Common arrow length = %.3f", L)) +
      ggplot2::theme_void()
    if (labels) {
      plot <- plot + ggplot2::geom_text(
        data = arrows_data, ggplot2::aes(x = .data$tip_x, y = .data$tip_y, label = .data$group),
        size = 3.5 * cex.labels
      )
    }
    print(plot)
    return(invisible(x))
  }

  has_circ <- requireNamespace("circular", quietly = TRUE)

  if (has_circ) {
    th <- circular::circular(
      a_hat,
      type = "angles",
      units = "radians",
      modulo = "2pi",
      template = "none"
    )
    op <- par(mar = c(1.8, 1.8, 4.2, 1.8))
    on.exit(par(op))
    circular::plot.circular(
      th,
      stack = FALSE,
      shrink = 1,
      zero = 0,
      rotation = "counter",
      cex = 0,
      axes = TRUE,
      main = ""
    )
    mtext(main, line = 3, font = 2, cex = 1.1)
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
      line = 1,
      cex = 0.9
    )
    if (labels) {
      text(xend * 1.04, yend * 1.04, labels = names(a_hat), cex = cex.labels)
    }
  } else {
    # base plotting fallback
    op <- par(mar = c(1.8, 1.8, 4.2, 1.8))
    on.exit(par(op))
    plot.new()
    plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)
    title(main, line = 3)
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
      line = 1,
      cex = 0.9
    )
    if (labels) {
      text(xend * 1.04, yend * 1.04, labels = names(a_hat), cex = cex.labels)
    }
  }

  invisible(x)
}
