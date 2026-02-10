library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

test_that("angular_re fits and predicts", {
  set.seed(222)
  n_cluster <- 4
  per_cluster <- 8
  cluster <- rep(letters[1:n_cluster], each = per_cluster)
  n <- length(cluster)

  ref <- runif(n, -pi, pi)
  cov <- runif(n, -pi, pi)
  scale <- runif(n, 0.8, 1.5)
  a_true <- rep(rnorm(n_cluster, sd = 0.4), each = per_cluster)
  mu_true <- wrap_angle(ref + 0.3 * sin(cov) * scale + a_true)
  y <- wrap_angle(mu_true + rnorm(n, sd = 0.15))
  df <- data.frame(y = y, ref = ref, cov = cov, scale = scale)

  fit <- angular_re(
    y ~ ref + cov:scale,
    data = df,
    cluster = cluster,
    control = list(maxit = 200, reltol = 1e-6, trace = 0)
  )

  preds <- predict(fit, type = "auto", cluster = cluster)
  expect_length(preds, n)
  expect_true(all(is.finite(preds)))
})

test_that("predict.angular_re no longer depends on calling environment", {
  set.seed(333)
  n <- 30
  cluster <- rep(letters[1:5], each = 6)
  df <- data.frame(
    y = runif(n, -pi, pi),
    ref = runif(n, -pi, pi),
    cov = runif(n, -pi, pi),
    scale = runif(n, 0.8, 1.5)
  )

  f <- y ~ ref + cov:scale
  fit <- angular_re(f, data = df, cluster = cluster, control = list(maxit = 50, reltol = 1e-6))

  rm(f)
  rm(df)

  preds <- predict(fit)
  expect_length(preds, length(cluster))
})

test_that("angular_re aligns cluster with model-frame NA handling", {
  set.seed(444)
  n <- 40
  cluster <- rep(letters[1:5], each = 8)
  df <- data.frame(
    y = runif(n, -pi, pi),
    ref = runif(n, -pi, pi),
    cov = runif(n, -pi, pi),
    scale = runif(n, 0.8, 1.5)
  )
  df$y[1] <- NA

  fit <- angular_re(
    y ~ ref + cov:scale,
    data = df,
    cluster = cluster,
    control = list(maxit = 50, reltol = 1e-6)
  )

  expect_equal(length(fit$cluster), n - 1)
})
