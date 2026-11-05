test_that("homogeneous inference is withheld when the reference limit improves the fit", {
  d <- data.frame(ref = 0, x = pi / 2, y = rep(pi / 2 + c(-0.3, 0.3), 15))
  expect_warning(fit <- angular(y ~ ref + x, d, reference = "first"),
                 "zero reference coefficient")
  expect_identical(fit$fit_status, "reference_boundary")
  expect_equal(fit$reference_limit_cosine, cos(0.3))
  expect_gt(fit$reference_limit_cosine, fit$MaxCosine)
  expect_false(fit$convergence)
  expect_false(fit$inference_available)
  expect_true(all(is.na(vcov(fit))))
  other <- angular(y ~ x + ref, d, reference = "first", initbeta = 0)
  expect_true(other$convergence)
  expect_equal(other$MaxCosine, cos(0.3))
})

test_that("grouped inference is withheld when the independent limit fits better", {
  d <- data.frame(ref = 0, y = rep(c(-0.3, 0.3), 15), id = rep(1:15, each = 2))
  expect_warning(fit <- angular_re(y ~ ref + (1 | id), d), "independent limit")
  expect_identical(fit$fit_status, "independence_limit")
  expect_gt(fit$independence_logLik, fit$logLik)
  expect_false(fit$convergence)
  expect_false(fit$inference_available)
  expect_true(all(is.na(vcov(fit))))
})

test_that("new data supply fitted predictors including resolved dot terms", {
  data(sandhopper_trials)
  d <- stats::na.omit(sandhopper_trials)
  wind <- 1.2
  fits <- list(
    angular(escape ~ sun + wind:wind_speed, d, reference = "first"),
    consensus(escape ~ sun + wind:wind_speed, d),
    angular_re(escape ~ sun + wind:wind_speed + (1 | Anim), d)
  )
  new <- d[1, c("sun", "wind_speed", "Anim")]
  for (fit in fits) {
    expect_error(predict(fit, new), "not found in 'newdata': wind", fixed = TRUE)
  }

  dot_data <- d[c("escape", "sun", "wind")]
  fit <- angular(escape ~ ., dot_data, reference = "first")
  new <- dot_data[1:5, c("sun", "wind")]
  expect_equal(predict(fit, new), fitted(fit)[1:5])
  expect_equal(predict(fit, list2env(as.list(new))), fitted(fit)[1:5])
})

test_that("T2 always rejects incomplete optimization", {
  data(bison)
  d <- bison[1:120, ]
  reduced <- angular(y.dir ~ y.prec, d, reference = "first")
  unfinished <- suppressWarnings(angular(y.dir ~ y.prec + x.meadow:z.meadow, d,
    reference = "first", initbeta = 0.1, control = list(maxiter = 0)))
  expect_false(unfinished$convergence)
  expect_error(angular_lrtest(unfinished, reduced, check = FALSE), "must have converged")
})

test_that("lag construction preserves pre-existing columns", {
  d <- data.frame(y = c(0, 0.1, 0.2), y_tm2 = 1:3)
  expect_error(autoregressivedata(d, "y", 2), "Lag column(s) already present: y_tm2",
               fixed = TRUE)
  expect_equal(autoregressivedata(d, "y", 1)$y_tm2, 1:3)
})
