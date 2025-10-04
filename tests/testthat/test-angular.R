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

test_that("angular detects singular designs", {
  set.seed(321)
  x <- runif(15, -pi, pi)
  df <- data.frame(y = wrap_angle(x + rnorm(15, sd = 0.1)), ref = x, dup = x)
  expect_error(
    angular(y ~ ref + dup:dup, data = df),
    "rank deficient"
  )
})

test_that("angular predictions are invariant to angle wrap-around", {
  set.seed(456)
  x <- runif(30, -pi, pi)
  z <- runif(30, 0.5, 1.5)
  y <- wrap_angle(x + 0.4 * sin(x) * z + rnorm(30, sd = 0.05))
  df <- data.frame(y = y, ref = x, x1 = x, z1 = z)
  fit1 <- angular(y ~ ref + x1:z1, data = df)
  df2 <- df
  df2$y <- wrap_angle(df$y + 2 * pi)
  fit2 <- angular(y ~ ref + x1:z1, data = df2)
  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-6)
})
