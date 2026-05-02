library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

test_that("consensus rejects invalid weights", {
  df <- data.frame(
    y = wrap_angle(runif(20, -pi, pi)),
    x1 = runif(20, -pi, pi),
    x2 = runif(20, -pi, pi),
    z2 = runif(20, 0.5, 1.5)
  )

  expect_error(consensus(y ~ x1 + x2:z2, data = df, weights = c(-1, rep(1, 19))))
  expect_error(consensus(y ~ x1 + x2:z2, data = df, weights = rep(0, 20)))
})

test_that("consensus supports maxiter=0 without crashing", {
  set.seed(21)
  n <- 40
  df <- data.frame(
    y = runif(n, -pi, pi),
    x1 = runif(n, -pi, pi),
    x2 = runif(n, -pi, pi),
    z2 = runif(n, 0.1, 2)
  )

  fit <- consensus(y ~ x1 + x2:z2, data = df, control = list(maxiter = 0))
  expect_s3_class(fit, "consensus")
  expect_true(is.matrix(fit$iter.detail))
})

test_that("consensus BIC in object and method are coherent", {
  set.seed(22)
  n <- 60
  df <- data.frame(
    y = runif(n, -pi, pi),
    x1 = runif(n, -pi, pi),
    x2 = runif(n, -pi, pi),
    z2 = runif(n, 0.1, 2)
  )
  w <- runif(n, 0.2, 2)

  fit <- consensus(y ~ x1 + x2:z2, data = df, weights = w)
  expect_equal(as.numeric(BIC(fit)), as.numeric(fit$BIC), tolerance = 1e-8)
})

test_that("coef.consensus supports both kappa and beta outputs", {
  set.seed(23)
  n <- 80
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  y <- wrap_angle(x1 + 0.2 * sin(x2) + rnorm(n, sd = 0.2))
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- consensus(y ~ x1 + x2, data = df)
  kappa <- coef(fit, type = "kappa")
  beta <- coef(fit, type = "beta", reference = c("name", "x1"))

  expect_length(kappa, 2)
  expect_length(beta, 1)
  expect_equal(attr(beta, "reference"), "x1")
})
