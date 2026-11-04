# Independent deterministic checks: no Monte Carlo experiment is used here.
library(CircularRegression)

fd_gradient <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    delta <- rep(0, length(x)); delta[j] <- h
    (f(x + delta) - f(x - delta)) / (2 * h)
  }, numeric(1))
}

fd_hessian <- function(f, x, h = 2e-4) {
  p <- length(x)
  H <- matrix(0, p, p)
  for (j in seq_len(p)) for (k in seq_len(p)) {
    a <- b <- rep(0, p); a[j] <- h; b[k] <- h
    H[j, k] <- (f(x + a + b) - f(x + a - b) -
                 f(x - a + b) + f(x - a - b)) / (4 * h^2)
  }
  H
}

test_that("modified consensus terms reach the independently computed optimum", {
  d <- data.frame(y = rep(c(-pi / 3, pi / 3), 20), x = 0, z = 2)
  fit <- consensus(y ~ x:z, d, initkappa = 0.1, control = list(mindiff = 1e-10))
  root <- uniroot(function(k) besselI(k, 1) / besselI(k, 0) - 0.5,
                  c(0, 10), tol = 1e-12)$root / 2
  expect_equal(unname(coef(fit)), root, tolerance = 1e-9)
  expect_true(fit$convergence)
  expect_lte(fit$score_norm, 1e-10)
  ll <- function(k) sum(2 * k * cos(d$y) - log(2 * pi * besselI(2 * k, 0)))
  expect_lt(abs(fd_gradient(ll, coef(fit))), 1e-7)
  expect_equal(as.numeric(logLik(fit)), ll(root), tolerance = 1e-10)
})

test_that("a shared modifier is not confused with an angle by terms reordering", {
  d <- data.frame(y = c(-0.2, 0.3, 0.5), x1 = c(0, 1, 2),
                  x2 = c(-0.5, -1, -2), z = c(0.1, 0.2, 0.3))
  des <- CircularRegression:::.circular_design(y ~ x1:z + x2:z, d)
  expect_equal(unname(des$matx), unname(as.matrix(d[c("x1", "x2")])))
  expect_equal(unname(des$matz), cbind(d$z, d$z))
})

test_that("consensus gradients, information and weighted sandwich are correct", {
  t <- seq(-2, 2, length.out = 31)
  X <- cbind(t, cos(t)); Z <- cbind(1, 1.3 + sin(t)); y <- sin(2 * t)
  w <- seq(0, 2, length.out = length(t)); k <- c(0.4, -0.3)
  ll <- function(b) {
    u <- rowSums(sweep(Z * cos(X), 2, b, "*"))
    v <- rowSums(sweep(Z * sin(X), 2, b, "*"))
    sum(w * (u * cos(y) + v * sin(y) - log(2 * pi * besselI(sqrt(u^2 + v^2), 0))))
  }
  u <- (Z * cos(X)) %*% k; v <- (Z * sin(X)) %*% k
  deriv <- CircularRegression:::.consensus_derivatives(y, as.vector(atan2(v, u)),
    as.vector(sqrt(u^2 + v^2)), X, Z, w)
  expect_equal(unname(deriv$score), fd_gradient(ll, k), tolerance = 1e-7)
  expect_equal(unname(deriv$information), -fd_hessian(ll, k), tolerance = 2e-6)
  zero <- CircularRegression:::.consensus_derivatives(y, rep(0, length(t)), rep(0, length(t)), X, Z, w)
  exact <- (crossprod(Z * cos(X) * sqrt(w)) + crossprod(Z * sin(X) * sqrt(w))) / 2
  expect_equal(zero$information, exact, tolerance = 1e-12)
  fit <- consensus(y ~ x + q:z, data.frame(y, x = X[, 1], q = X[, 2], z = Z[, 2]), weights = w)
  deriv <- CircularRegression:::.consensus_derivatives(fit$y, fit$mui, fit$long, fit$matx, fit$matz, w)
  bread <- solve(deriv$information)
  expected <- bread %*% crossprod(deriv$scores) %*% bread
  expect_equal(unname(vcov(fit, robust = TRUE)), unname(expected), tolerance = 1e-10)
  expect_gt(max(abs(vcov(fit) - expected)), 1e-4)
  scaled <- consensus(y ~ x + q:z, data.frame(y, x = X[, 1], q = X[, 2], z = Z[, 2]), weights = 3 * w)
  expect_equal(coef(scaled), coef(fit), tolerance = 1e-6)
  expect_equal(vcov(scaled, robust = TRUE), vcov(fit, robust = TRUE), tolerance = 1e-6)
})

test_that("homogeneous sandwich uses the full observed derivative", {
  d <- data.frame(y = rep(c(-pi / 3, pi / 3), 20), x0 = 0, x1 = pi / 2)
  fit <- angular(y ~ x0 + x1, d, reference = "first", initbeta = 0)
  expect_equal(unname(vcov(fit, robust = TRUE)[1, 1]), 0.075, tolerance = 1e-12)
  expect_equal(besselI(fit$kappahat, 1) / besselI(fit$kappahat, 0), 0.5, tolerance = 1e-10)
  t <- seq(-2, 2, length.out = 41)
  X <- cbind(t, cos(t)); Z <- cbind(1.5 + sin(t), 1)
  y <- sin(t); b <- c(0.2, -0.4)
  components <- function(b) {
    u <- 1 + as.vector((Z * cos(X)) %*% b)
    v <- as.vector((Z * sin(X)) %*% b)
    list(mu = atan2(v, u), L = sqrt(u^2 + v^2))
  }
  objective <- function(b) sum(cos(y - components(b)$mu))
  cm <- components(b)
  deriv <- CircularRegression:::.angular_derivatives(y, cm$mu, cm$L, X, Z)
  expect_equal(unname(deriv$score), fd_gradient(objective, b), tolerance = 1e-7)
  expect_equal(unname(deriv$information), -fd_hessian(objective, b), tolerance = 3e-6)
})

test_that("all signed natural charts have correct Jacobians and orientations", {
  models <- c(mean_direction = 2L, decentred = 4L, presnell = 4L, jammalamadaka_sarma = 6L)
  for (model in names(models)) for (ref in seq_len(models[[model]])) for (sgn in c(-1, 1)) {
    n <- models[[model]]
    b <- seq(0.3, 0.7, length.out = n - 1L)
    mapping <- CircularRegression:::.natural_parameter_mapping(model, ref, n, sgn)
    # Wrapped angular outputs are locally smooth at these strictly nonzero fixtures.
    J <- t(vapply(seq_along(mapping$g(b)), function(k) {
      fd_gradient(function(x) unname(mapping$g(x)[k]), b)
    }, numeric(length(b))))
    expect_equal(unname(mapping$jacobian(b)), unname(J), tolerance = 1e-7)
    if (model %in% c("presnell", "jammalamadaka_sarma")) {
      full <- CircularRegression:::.full_coefficients_from_nonref(b, ref, n, sgn)
      q <- unname(mapping$g(b)); w <- seq(-2, 2, length.out = 19)
      if (model == "presnell") {
        original <- atan2(full[2] + full[4] * w, full[1] + full[3] * w)
        u <- 1 + q[2] * w; v <- q[1] + q[3] * w
      } else {
        u0 <- full[1] + (full[3] + full[4]) * cos(w) + (-full[5] + full[6]) * sin(w)
        v0 <- full[2] + (full[5] + full[6]) * cos(w) + (full[3] - full[4]) * sin(w)
        original <- atan2(v0, u0)
        u <- 1 + q[2] * cos(w) + q[3] * sin(w)
        v <- q[1] + q[4] * cos(w) + q[5] * sin(w)
      }
      orientation <- sign(full[1])
      reconstructed <- atan2(orientation * v, orientation * u)
      expect_lt(max(abs(atan2(sin(original - reconstructed), cos(original - reconstructed)))), 1e-12)
    }
  }
})

test_that("T2 refuses incompatible samples, references and unconverged fits", {
  t <- seq(-2.5, 2.5, length.out = 60)
  d <- data.frame(y = 0.2 * sin(t) + 0.4 * cos(4 * t), ref = 0, x = t, q = cos(t))
  full <- angular(y ~ ref + x, d, reference = "first")
  reduced <- angular(y ~ ref, d, reference = "first")
  expect_s3_class(angular_lrtest(full, reduced), "htest")
  bad <- reduced; bad$y[1] <- bad$y[1] + 0.1
  expect_error(angular_lrtest(full, bad), "identical observations")
  bad <- reduced; bad$model_data$ref[1] <- 0.1
  expect_error(angular_lrtest(full, bad), "identical values")
  bad <- reduced; bad$reference_sign <- -1
  expect_error(angular_lrtest(full, bad), "compatible reference")
  bad <- full; bad$convergence <- FALSE
  expect_error(angular_lrtest(bad, reduced), "converged")
  bad <- reduced; bad$term_labels <- "q"
  expect_error(angular_lrtest(full, bad), "nested")
})

test_that("uniform consensus and input boundary handling are explicit", {
  d <- data.frame(y = c(0, pi, -pi / 2, pi / 2), x = 0)
  fit <- consensus(y ~ x, d, initkappa = 0)
  expect_true(fit$convergence)
  expect_equal(vcov(fit)[1, 1], 0.5, tolerance = 1e-12)
  expect_false(consensus(y ~ x, d, initkappa = 0.3, control = list(maxiter = 0))$convergence)
  d$x <- factor(d$x)
  expect_error(consensus(y ~ x, d), "numeric")
  expect_error(autoregressivedata(data.frame(y = 1:2), "y", 1.5), "integer")
  expect_error(autoregressivedata(data.frame(y = 1:2), "y", Inf), "integer")
  lagged <- autoregressivedata(data.frame(y = 1:3, id = c(1, 1, 2)), "y", 4, "id")
  expect_equal(lagged$y_tm1, c(NA_integer_, 1L, NA_integer_))
  expect_true(all(is.na(lagged$y_tm4)))
  expect_warning(mu <- CircularRegression:::.mean_direction_from_terms(c(1, -1), matrix(0, 1, 2), matrix(1, 1, 2)), "singular")
  expect_true(is.na(mu))
})

test_that("group marginal likelihood agrees with direct integration and derivatives", {
  d <- data.frame(y = c(-0.4, 0.2, 0.7, -0.2, 0.5, 0.1), ref = 0,
                  x = c(0.1, 1.2, 2.3, -0.8, -1.5, 0.7), id = rep(1:3, each = 2))
  loglik <- function(par) {
    mu <- atan2(par[1] * sin(d$x), 1 + par[1] * cos(d$x))
    sum(vapply(split(seq_len(nrow(d)), d$id), function(idx) {
      log(integrate(function(a) vapply(a, function(ai) {
        exp(par[3] * cos(ai) + sum(par[2] * cos(d$y[idx] - mu[idx] - ai))) /
          ((2 * pi)^(length(idx) + 1) * besselI(par[3], 0) * besselI(par[2], 0)^length(idx))
      }, numeric(1)), -pi, pi, rel.tol = 1e-11)$value)
    }, numeric(1)))
  }
  fit <- suppressWarnings(angular_re(y ~ ref + x, d, d$id,
    init = list(beta = 0.2, kappa_e = 2, kappa_a = 1.5),
    control = list(maxit = 0, start_from_angular = FALSE)))
  par <- unname(coef(fit))
  expect_equal(fit$logLik, loglik(par), tolerance = 1e-9)
  expect_equal(fit$score, fd_gradient(loglik, par), tolerance = 1e-7)
  expect_equal(unname(fit$information), -fd_hessian(loglik, par), tolerance = 2e-5)
  expect_equal(colSums(fit$cluster_scores), fit$score, tolerance = 1e-12)
  if (all(is.finite(vcov(fit)))) {
    expected <- vcov(fit) %*% crossprod(fit$cluster_scores) %*% vcov(fit)
    expect_equal(vcov(fit, robust = TRUE), expected, tolerance = 1e-10)
  }
  expect_error(angular_re(y ~ ref, d, seq_len(nrow(d))), "repeated observations")
  expect_warning(pred <- predict(fit, d[1:2, ], cluster = c("1", "new")), "new")
  expect_equal(pred[2], predict(fit, d[2, ], type = "marginal"))
  expect_equal(predict(fit, type = "conditional"),
    atan2(sin(fitted(fit) + fit$ranef[as.character(fit$cluster)]),
          cos(fitted(fit) + fit$ranef[as.character(fit$cluster)])), ignore_attr = TRUE)
})

test_that("PIT endpoints and high-concentration CDF retain von Mises probabilities", {
  skip_if_not_installed("ggcircular")
  cdf <- CircularRegression:::.vonmises_cdf(c(0, pi, 2 * pi), 200)
  expect_equal(cdf, c(0, 0.5, 1), tolerance = 1e-12)
  a <- pi + 0.07
  expected <- integrate(function(x) exp(200 * (cos(x - pi) - 1)) /
      (2 * pi * besselI(200, 0, expon.scaled = TRUE)), 0, a, rel.tol = 1e-10)$value
  expect_equal(CircularRegression:::.vonmises_cdf(a, 200), expected, tolerance = 1e-9)
})

test_that("a negative change of reference preserves the curve and both covariances", {
  t <- seq(-2.9, 2.9, length.out = 100)
  d <- data.frame(y = atan2(-0.3 * sin(t), 1 - 0.3 * cos(t)) + 0.2 * cos(5 * t),
                  x0 = 0, x1 = t)
  a <- angular(y ~ x0 + x1, d, reference = "x0", initbeta = -0.3,
               control = list(mindiff = 1e-9))
  b <- angular(y ~ x0 + x1, d, reference = "x1", reference_sign = -1,
               initbeta = -1 / unname(coef(a)), control = list(mindiff = 1e-9))
  expect_true(a$convergence && b$convergence)
  expect_equal(fitted(a), fitted(b), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8)
  J <- 1 / unname(coef(a))^2
  expect_equal(unname(vcov(b)), J^2 * unname(vcov(a)), tolerance = 1e-6)
  expect_equal(unname(vcov(b, robust = TRUE)), J^2 * unname(vcov(a, robust = TRUE)), tolerance = 1e-6)
})

test_that("affine wrappers accept negative covariates and propagate orientation", {
  w <- seq(-2, 2, length.out = 81)
  d <- data.frame(y = atan2(0.4 - 0.2 * w, -1 + 0.3 * w) + 0.1 * cos(7 * w), w)
  fit <- presnellModel(d, "y", "w")
  q <- fit$natural_parameters[, "estimate"]
  s <- fit$natural_orientation
  mu <- atan2(s * (q["beta2"] + q["beta4"] * w), s * (1 + q["beta3"] * w))
  expect_equal(unname(mu), fitted(fit$homogeneous_fit), tolerance = 1e-10)
  expect_equal(s, -1)
  J <- fit$natural_jacobian
  expect_equal(fit$natural_vcov_robust,
    J %*% vcov(fit$homogeneous_fit, robust = TRUE) %*% t(J), tolerance = 1e-10)
  d$y[3] <- NA
  missing <- presnellModel(d, "y", "w")
  explicit <- presnellModel(d[-3, ], "y", "w")
  expect_equal(coef(missing$fit), coef(explicit$fit), tolerance = 1e-10)
  names(d)[1] <- "const0"
  expect_error(presnellModel(d, "const0", "w"), "conflicts")
})

test_that("missing data, case weights and public controls are checked", {
  d <- data.frame(y = c(0.2, 0.4, NA, -0.2, -0.4), ref = 0)
  a <- consensus(y ~ ref, d, weights = c(1, 2, 100, 3, 4))
  b <- consensus(y ~ ref, d[-3, ], weights = c(1, 2, 3, 4))
  expect_equal(coef(a), coef(b), tolerance = 1e-12)
  expect_equal(vcov(a, robust = TRUE), vcov(b, robust = TRUE), tolerance = 1e-12)
  expect_error(predict(a, data.frame(ref = c(0, NA))), "finite")
  expect_error(consensus(y ~ ref, d, weights = factor(rep(1, 5))), "numeric")
  expect_error(angular(y ~ ref, d, control = list(maxiter = 0.5)), "integer")
  expect_error(angular(y ~ ref, d, control = list(unknown = 1)), "Unknown")
  expect_error(circular_regression(y ~ ref, d, method = "homogeneous", weights = rep(1, 5)), "weights")
  expect_error(circular_regression(y ~ ref, d, control = list(maxiter = 1)), "control_consensus")
})

test_that("group covariance uses independent group scores at an interior fit", {
  data(Sandhopper)
  sh <- Sandhopper[1:24, ]
  d <- sh[rep(1:24, each = 5), ]
  d$y <- as.vector(t(as.matrix(sh[paste0("LN", 1:5)]))) * pi / 180
  d$ref <- d$Azimuth * pi / 180
  fit <- angular_re(y ~ ref, d, d$Anim)
  expect_true(fit$convergence)
  expect_false(fit$boundary)
  expect_true(all(is.finite(vcov(fit))))
  bread <- solve(fit$information)
  expect_equal(unname(vcov(fit, robust = TRUE)),
               unname(bread %*% crossprod(fit$cluster_scores) %*% bread), tolerance = 1e-9)
  expect_error(predict(fit, data.frame(ref = c(0, NA))), "finite")
})

test_that("concentration boundaries are not reported as regular Wald fits", {
  d <- data.frame(y = rep(0, 10), ref = 0)
  expect_warning(hom <- angular(y ~ ref, d), "finite concentration")
  expect_true(is.infinite(hom$kappahat))
  expect_warning(con <- consensus(y ~ ref, d, initkappa = 1), "perfect alignment")
  expect_false(con$convergence)
  expect_true(con$alignment_boundary)
  expect_true(all(is.na(vcov(con))))
  uniform <- data.frame(y = rep(c(0, pi / 2, pi, -pi / 2), 4), ref = 0,
                       id = rep(1:4, each = 4))
  fit <- suppressWarnings(angular_re(y ~ ref, uniform, uniform$id,
                                     init = list(kappa_e = 0.5, kappa_a = 2)))
  expect_true(fit$boundary)
  expect_true(all(is.na(vcov(fit))))
  expect_true(all(is.na(vcov(fit, robust = TRUE))))
})
