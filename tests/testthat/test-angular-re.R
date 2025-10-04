library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

test_that("angular_re fits and predicts", {
  set.seed(222)
  n_cluster <- 4
  per_cluster <- 5
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
    control = list(maxit = 100, reltol = 1e-6, trace = 0)
  )
  expect_true(fit$convergence)
  preds <- predict(fit, type = "auto", cluster = cluster)
  expect_length(preds, n)
  expect_true(all(is.finite(preds)))
})

test_that("predict.angular_re handles newdata and a_hat", {
  set.seed(333)
  cluster <- rep(letters[1:3], each = 4)
  n <- length(cluster)
  ref <- runif(n, -pi, pi)
  cov <- runif(n, -pi, pi)
  scale <- runif(n, 0.5, 1.2)
  y <- wrap_angle(ref + 0.2 * sin(cov) * scale + rep(rnorm(3, sd = 0.3), each = 4) + rnorm(n, sd = 0.1))
  df <- data.frame(y = y, ref = ref, cov = cov, scale = scale)
  fit <- angular_re(y ~ ref + cov:scale, data = df, cluster = cluster, control = list(maxit = 80))

  new_df <- df[1:3, ]
  pred_marg <- predict(fit, newdata = new_df, type = "marginal")
  expect_length(pred_marg, 3)
  pred_cond <- predict(fit, newdata = new_df, type = "conditional", a_hat = 0.1)
  expect_length(pred_cond, 3)
})
