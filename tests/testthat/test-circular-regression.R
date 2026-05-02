library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

simulate_angular_data <- function(n = 300, beta = 0.35, sd = 0.06) {
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  z2 <- runif(n, 0.3, 1.8)
  mu <- atan2(
    sin(x1) + beta * z2 * sin(x2),
    cos(x1) + beta * z2 * cos(x2)
  )
  y <- wrap_angle(mu + rnorm(n, sd = sd))
  data.frame(y = y, x1 = x1, x2 = x2, z2 = z2, mu = mu)
}

test_that("angular recovers a simple simulated coefficient", {
  set.seed(1001)
  df <- simulate_angular_data(n = 350, beta = 0.4, sd = 0.04)

  fit <- angular(y ~ x1 + x2:z2, data = df, reference = c("name", "x1"))

  expect_s3_class(fit, "angular")
  expect_equal(unname(coef(fit)), 0.4, tolerance = 0.08)
})

test_that("predict methods match stored fitted values and manual directions", {
  set.seed(1002)
  df <- simulate_angular_data(n = 120, beta = 0.3, sd = 0.05)
  newdata <- df[1:5, c("x1", "x2", "z2")]

  fit_ang <- angular(y ~ x1 + x2:z2, data = df, reference = c("name", "x1"))
  expect_equal(predict(fit_ang), fitted(fit_ang))
  manual_ang <- atan2(
    sin(newdata$x1) + coef(fit_ang)[1] * newdata$z2 * sin(newdata$x2),
    cos(newdata$x1) + coef(fit_ang)[1] * newdata$z2 * cos(newdata$x2)
  )
  expect_equal(predict(fit_ang, newdata = newdata), manual_ang, tolerance = 1e-10)

  fit_cons <- consensus(y ~ x1 + x2:z2, data = df)
  expect_equal(predict(fit_cons), fitted(fit_cons))
  kap <- coef(fit_cons)
  manual_cons <- atan2(
    kap[1] * sin(newdata$x1) + kap[2] * newdata$z2 * sin(newdata$x2),
    kap[1] * cos(newdata$x1) + kap[2] * newdata$z2 * cos(newdata$x2)
  )
  expect_equal(predict(fit_cons, newdata = newdata), unname(manual_cons), tolerance = 1e-10)
})

test_that("circular_regression exposes coherent glm-like S3 methods", {
  set.seed(1003)
  df <- simulate_angular_data(n = 130, beta = 0.25, sd = 0.08)
  newdata <- df[1:4, c("x1", "x2", "z2")]

  fit <- circular_regression(
    y ~ x1 + x2:z2,
    data = df,
    method = "homogeneous",
    reference = c("name", "x1")
  )

  expect_s3_class(fit, "circular_regression")
  expect_s3_class(summary(fit), "summary.circular_regression")
  expect_length(coef(fit), 1)
  expect_equal(fitted(fit), predict(fit))
  expect_length(residuals(fit), fit$nobs)
  expect_length(predict(fit, newdata = newdata), nrow(newdata))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("circular_regression default two-step fit works on bison subset", {
  data(bison)
  d <- bison[seq_len(100), ]

  fit <- circular_regression(y.dir ~ y.prec + x.meadow:z.meadow, data = d)

  expect_s3_class(fit, "circular_regression")
  expect_equal(fit$method, "two_step")
  expect_true(all(is.finite(coef(fit))))
  expect_length(predict(fit), fit$nobs)
})

test_that("NA handling aligns weights with the model frame", {
  set.seed(1004)
  df <- simulate_angular_data(n = 70, beta = 0.2, sd = 0.07)
  rownames(df) <- paste0("row", seq_len(nrow(df)))
  df$y[c(2, 5)] <- NA
  weights <- seq_len(nrow(df))

  fit <- consensus(y ~ x1 + x2:z2, data = df, weights = weights)

  expect_equal(fit$nobs, nrow(df) - 2)
  expect_equal(fit$weights, weights[-c(2, 5)])
})

test_that("input validation catches invalid modifiers, weights, and controls", {
  set.seed(1005)
  df <- simulate_angular_data(n = 40)

  df_bad_z <- df
  df_bad_z$z2[1] <- -0.1
  expect_error(
    angular(y ~ x1 + x2:z2, data = df_bad_z),
    regexp = "non-negative"
  )

  expect_error(
    consensus(y ~ x1 + x2:z2, data = df, weights = rep(1, 5)),
    regexp = "weights"
  )

  expect_error(
    angular(y ~ x1 + x2:z2, data = df, control = list(maxiter = -1)),
    regexp = "maxiter"
  )
})

test_that("angle wrapping is invariant to integer multiples of 2*pi", {
  set.seed(1006)
  df <- simulate_angular_data(n = 110, beta = 0.3, sd = 0.05)
  df_shift <- df
  df_shift$y <- df_shift$y + 4 * pi
  df_shift$x1 <- df_shift$x1 - 2 * pi
  df_shift$x2 <- df_shift$x2 + 6 * pi

  fit <- angular(y ~ x1 + x2:z2, data = df, reference = c("name", "x1"))
  fit_shift <- angular(y ~ x1 + x2:z2, data = df_shift, reference = c("name", "x1"))

  expect_equal(coef(fit), coef(fit_shift), tolerance = 1e-8)
  expect_equal(
    wrap_angle(predict(fit)),
    wrap_angle(predict(fit_shift)),
    tolerance = 1e-8
  )
})

test_that("small reference-only samples fit and predict", {
  set.seed(1007)
  n <- 6
  x <- runif(n, -pi, pi)
  y <- wrap_angle(x + rnorm(n, sd = 0.05))
  df <- data.frame(y = y, x = x)

  fit <- circular_regression(y ~ x, data = df, method = "homogeneous")

  expect_length(coef(fit), 0)
  expect_length(predict(fit), n)
  expect_true(all(is.finite(residuals(fit))))
})
