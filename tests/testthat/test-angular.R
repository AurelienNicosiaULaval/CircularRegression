library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

test_that("angular handles reference-only model", {
  set.seed(123)
  x <- runif(20, -pi, pi)
  y <- wrap_angle(x + rnorm(20, sd = 0.05))
  df <- data.frame(y = y, ref = x)
  fit <- angular(y ~ ref, data = df)
  expect_s3_class(fit, "angular")
  expect_equal(nrow(fit$parameters), 0)
  expect_true(is.matrix(fit$varcov0))
})
