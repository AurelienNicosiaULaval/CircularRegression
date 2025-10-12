library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))


test_that("consensus rejects invalid weights", {
  df <- data.frame(
    y = wrap_angle(runif(5, -pi, pi)),
    ref = runif(5, -pi, pi),
    x1 = runif(5, -pi, pi),
    z1 = runif(5, 0.5, 1.5)
  )
  expect_error(consensus(
    y ~ ref + x1:z1,
    data = df,
    weights = c(-1, rep(1, 4))
  ))
  expect_error(consensus(y ~ ref + x1:z1, data = df, weights = rep(0, 5)))
})
