library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

test_that("angular handles reference-only model", {
  set.seed(123)
  x <- runif(20, -pi, pi)
  y <- wrap_angle(x + rnorm(20, sd = 0.05))
  df <- data.frame(y = y, ref = x)

  fit <- angular(y ~ ref, data = df, reference = "auto")
  expect_s3_class(fit, "angular")
  expect_equal(nrow(fit$parameters), 0)
  expect_equal(fit$reference, "ref")
})

test_that("angular supports maxiter=0 without crashing", {
  set.seed(11)
  n <- 40
  df <- data.frame(
    y = runif(n, -pi, pi),
    x1 = runif(n, -pi, pi),
    x2 = runif(n, -pi, pi),
    z2 = runif(n, 0.1, 2)
  )

  fit <- angular(y ~ x1 + x2:z2, data = df, control = list(maxiter = 0))
  expect_s3_class(fit, "angular")
  expect_true(is.matrix(fit$iter.detail))
})

test_that("automatic reference selection follows mean-cos criterion", {
  set.seed(42)
  n <- 80
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  y <- wrap_angle(x2 + rnorm(n, sd = 0.01))
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- angular(y ~ x1 + x2, data = df, reference = "auto")
  expect_equal(fit$reference, "x2")
})

test_that("angular detects non-identifiable fixed-effects design", {
  set.seed(5)
  n <- 80
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  y <- runif(n, -pi, pi)
  df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x2)

  expect_error(
    angular(y ~ x1 + x2 + x3, data = df, reference = "first"),
    regexp = "rank deficient|identifiable"
  )
})

test_that("logLik/angular model-comparison attributes are present", {
  set.seed(8)
  n <- 100
  df <- data.frame(
    y = runif(n, -pi, pi),
    x1 = runif(n, -pi, pi),
    x2 = runif(n, -pi, pi)
  )

  fit <- angular(y ~ x1 + x2, data = df)
  ll <- logLik(fit)

  expect_true(inherits(ll, "logLik"))
  expect_equal(attr(ll, "nobs"), n)
  expect_true(is.finite(BIC(fit)))
})

test_that("T2 statistic is available through anova and angular_lrtest", {
  set.seed(77)
  n <- 120
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  y <- wrap_angle(x1 + 0.3 * sin(x2) + rnorm(n, sd = 0.2))
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  full <- angular(y ~ x1 + x2, data = df)
  reduced <- angular(y ~ x1, data = df)

  tst <- angular_lrtest(full, reduced)
  expect_true(is.list(tst))
  expect_true(is.finite(tst$statistic))

  tst2 <- anova(full, reduced, test = "T2")
  expect_true(is.list(tst2))
  expect_true(is.finite(tst2$statistic))
})
