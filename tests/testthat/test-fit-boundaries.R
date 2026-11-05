test_that("consensus detects an unbounded likelihood with zero modifiers", {
  d <- data.frame(
    y = c(rep(0, 10), rep(c(-1, 1), 5)), x = 0,
    z = rep(c(1, 0), each = 10)
  )
  expect_warning(fit <- consensus(y ~ x:z, d, initkappa = 1), "no finite maximum")
  expect_false(fit$convergence)
  expect_identical(fit$fit_status, "unbounded_likelihood")
  expect_false(fit$inference_available)
  expect_true(all(is.na(vcov(fit))))
  expect_true(all(is.na(vcov(fit, robust = TRUE))))
  expect_true(is.na(fit$AIC))
  expect_error(logLik(fit), "no finite maximum")
})

test_that("recession detection permits zero and opposite signed components", {
  detect <- CircularRegression:::.consensus_recession
  expect_lt(detect(pi, matrix(0, 1, 1), matrix(1, 1, 1), 1), 0)
  expect_gt(detect(c(0, 1), matrix(0, 2, 1), matrix(1, 2, 1), c(1, 0)), 0)

  # Neither coordinate alone is feasible: d1 - 2 d2 >= 0 and -d1 + 3 d2 >= 0.
  x <- matrix(c(0, pi, pi, 0), 2, 2)
  z <- matrix(c(1, 1, 2, 3), 2, 2)
  direction <- detect(c(0, 0), x, z, c(1, 1))
  expect_false(is.null(direction))
  expect_gt(min(as.vector((z * cos(x)) %*% direction)), 0)
  expect_lt(max(abs(as.vector((z * sin(x)) %*% direction))), 1e-12)

  expect_null(detect(c(-1e-7, 1e-7), matrix(0, 2, 1), matrix(1, 2, 1), c(1, 1)))
  expect_null(detect(c(-pi / 3, pi / 3), matrix(0, 2, 1), matrix(1, 2, 1), c(1, 1)))
})

test_that("perfect within-group alignment is not a regular converged fit", {
  d <- data.frame(y = 0, ref = 0, id = rep(1:8, each = 5))
  expect_warning(fit <- angular_re(y ~ ref, d, d$id), "no finite maximum")
  expect_false(fit$convergence)
  expect_true(fit$boundary)
  expect_true(fit$alignment_boundary)
  expect_identical(fit$fit_status, "unbounded_likelihood")
  expect_false(fit$inference_available)
  expect_true(all(is.na(vcov(fit))))
  expect_true(all(is.na(vcov(fit, robust = TRUE))))
  expect_true(is.na(fit$rho_SS))
})

test_that("unconverged fixed fits do not report Wald standard errors", {
  d <- data.frame(y = c(-.3, .1, .4, -.2), ref = 0, x = c(-1, .2, 1, 2))
  homogeneous <- angular(y ~ ref + x, d, reference = "first", initbeta = 0,
                         control = list(maxiter = 0))
  combined <- consensus(y ~ ref + x, d, initkappa = c(1, .2),
                        control = list(maxiter = 0))
  for (fit in list(homogeneous, combined)) {
    expect_false(fit$convergence)
    expect_false(fit$inference_available)
    expect_true(all(is.na(vcov(fit))))
    expect_true(all(is.na(vcov(fit, robust = TRUE))))
    expect_output(print(fit), "Converged: FALSE")
    expect_error(logLik(fit), "converged")
  }
})

test_that("the two-step workflow stops when initialization is unbounded", {
  d <- data.frame(y = rep(0, 10), ref = 0)
  expect_warning(
    expect_error(angular_two_step(y ~ ref, d), "use angular\\(\\) directly"),
    "no finite maximum"
  )
})
