library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

test_that("consensus fits with and without weights", {
  set.seed(99)
  n <- 25
  ref <- runif(n, -pi, pi)
  x1 <- runif(n, -pi, pi)
  z1 <- runif(n, 0.2, 1.5)
  y <- wrap_angle(ref + 0.6 * sin(x1) * z1 + rnorm(n, sd = 0.1))
  w <- runif(n, 0.5, 1.5)
  df <- data.frame(y = y, ref = ref, x1 = x1, z1 = z1)
  fit_w <- consensus(y ~ ref + x1:z1, data = df, weights = w, control = list(maxiter = 50))
  fit_u <- consensus(y ~ ref + x1:z1, data = df, control = list(maxiter = 50))
  expect_s3_class(fit_w, "consensus")
  expect_equal(nrow(fit_w$parameters), 2)
  expect_equal(dim(fit_w$varcov1), c(2, 2))
  expect_false(isTRUE(all.equal(coef(fit_w), coef(fit_u))))
})

test_that("consensus rejects invalid weights", {
  df <- data.frame(
    y = wrap_angle(runif(5, -pi, pi)),
    ref = runif(5, -pi, pi),
    x1 = runif(5, -pi, pi),
    z1 = runif(5, 0.5, 1.5)
  )
  expect_error(consensus(y ~ ref + x1:z1, data = df, weights = c(-1, rep(1, 4))))
  expect_error(consensus(y ~ ref + x1:z1, data = df, weights = rep(0, 5)))
})
