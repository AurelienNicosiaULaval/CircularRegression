test_that("AIC and BIC compare every supplied fitted model", {
  data(bison)
  d <- na.omit(bison[1:150, c("y.dir", "y.prec", "x.meadow", "z.meadow")])
  a <- consensus(y.dir ~ y.prec, d)
  b <- consensus(y.dir ~ y.prec + x.meadow:z.meadow, d)
  outer <- circular_regression(y.dir ~ y.prec, d, method = "consensus")
  two <- angular_two_step(y.dir ~ y.prec, d)
  for (criterion in list(AIC, BIC)) {
    compared <- criterion(a, b, outer, two)
    expect_s3_class(compared, "data.frame")
    expect_equal(nrow(compared), 4L)
    expected <- vapply(list(a, b, outer, two), criterion, numeric(1))
    expect_equal(compared[[2L]], unname(expected))
    expect_equal(nrow(criterion(outer, a)), 2L)
    expect_equal(nrow(criterion(two, a)), 2L)
  }
  expect_equal(AIC(a, k = 3), -2 * as.numeric(logLik(a)) + 3 * length(coef(a)))
})

test_that("high-level covariance methods use the selected fitted scale", {
  data(bison)
  d <- bison[1:100, ]
  for (method in c("homogeneous", "consensus", "two_step")) {
    fit <- circular_regression(y.dir ~ y.prec, d, method = method)
    core <- CircularRegression:::.primary_circular_fit(fit)
    expect_equal(vcov(fit), vcov(core))
    expect_equal(vcov(fit, robust = TRUE), vcov(core, robust = TRUE))
  }
  fit <- angular_two_step(y.dir ~ y.prec, d)
  expect_equal(vcov(fit), vcov(fit$homogeneous_fit))
  expect_equal(vcov(fit, step = "consensus"), vcov(fit$consensus_fit))
  expect_error(vcov(fit, robust = NA), "TRUE or FALSE")
})

test_that("the public homogeneous positional arguments retain their order", {
  d <- data.frame(y = c(-.3, .1, .4, -.2), ref = 0, x = c(-1, .2, 1, 2))
  positional <- angular(y ~ ref + x, d, "first", .2, list(mindiff = 1e-8))
  named <- angular(y ~ ref + x, d, reference = "first", initbeta = .2,
                   control = list(mindiff = 1e-8))
  expect_equal(coef(positional), coef(named))
  expect_equal(vcov(positional), vcov(named))
})

test_that("model-frame alignment works for list and environment data", {
  d <- data.frame(y = c(.1, .2, NA, -.1, -.2), ref = 0)
  weights <- c(1, 2, 100, 3, 4)
  expected <- consensus(y ~ ref, d, weights = weights)
  for (data in list(as.list(d), list2env(as.list(d)))) {
    fit <- consensus(y ~ ref, data, weights = weights)
    expect_equal(coef(fit), coef(expected))
    expect_equal(vcov(fit, robust = TRUE), vcov(expected, robust = TRUE))
  }
})

test_that("circular units and unsupported response shapes are checked", {
  y <- circular::circular(c(-10, 10, -15, 15), units = "degrees")
  d <- data.frame(y = y, ref = 0)
  expect_error(angular(y ~ ref, d), "radians")
  d$y <- as.numeric(y) * pi / 180
  d$ref <- circular::circular(rep(0, 4), units = "degrees")
  expect_error(consensus(y ~ ref, d), "radians")
  d$ref <- 0
  fit <- angular(y ~ ref, d)
  expect_error(predict(fit, data.frame(ref = circular::circular(0, units = "degrees"))), "radians")
  expect_error(angular(~ ref, d), "response")
  expect_error(angular(y ~ ref + offset(y), d), "Offset")
  d$x <- I(matrix(1:8, 4, 2))
  expect_error(consensus(y ~ x, d), "numeric and finite")
})

test_that("quoted column names keep colons inside the variable name", {
  d <- data.frame(y = c(-.2, .1, .3, -.1), check.names = FALSE)
  d[["reference:angle"]] <- 0
  fit <- angular(y ~ `reference:angle`, d)
  expect_equal(fitted(fit), rep(0, 4))
  expect_error(consensus(y ~ `reference:angle`:y + y:`reference:angle`, d), "opposite")
})

test_that("na.exclude restores fitted values, residuals and in-sample predictions", {
  d <- data.frame(y = c(.1, .2, NA, -.1, -.2), ref = 0)
  for (fit in list(angular(y ~ ref, d, na.action = na.exclude),
                   consensus(y ~ ref, d, na.action = "na.exclude"),
                   circular_regression(y ~ ref, d, na.action = na.exclude))) {
    expect_length(fitted(fit), nrow(d))
    expect_length(residuals(fit), nrow(d))
    expect_length(predict(fit), nrow(d))
    expect_true(is.na(fitted(fit)[3L]))
    expect_true(is.na(residuals(fit)[3L]))
    expect_equal(nrow(predict(fit, type = "components")), nrow(d))
    expect_length(predict(fit, newdata = d[-3L, ]), 4L)
  }
})

test_that("grouped exclusions retain original group alignment and output rows", {
  data(Sandhopper)
  wide <- Sandhopper[1:24, ]
  d <- wide[rep(seq_len(nrow(wide)), each = 5), ]
  d$y <- as.vector(t(as.matrix(wide[paste0("LN", 1:5)]))) * pi / 180
  d$ref <- d$Azimuth * pi / 180
  d$y[3] <- NA
  fit <- angular_re(y ~ ref, as.list(d), d$Anim, na.action = na.exclude)
  expect_length(fitted(fit), nrow(d))
  expect_length(residuals(fit), nrow(d))
  expect_length(predict(fit, type = "conditional"), nrow(d))
  expect_true(is.na(predict(fit)[3]))
  expect_equal(as.character(fit$cluster), as.character(d$Anim[-3]))
})

test_that("control components are checked before optimization", {
  d <- data.frame(y = c(.1, .2, -.1, -.2), ref = 0, id = rep(1:2, each = 2))
  expect_error(angular(y ~ ref, d, control = list(maxiter = NULL)), "NULL")
  expect_error(consensus(y ~ ref, d, control = list(maxiter = 5, maxiter = 8)), "duplicated")
  expect_error(angular_re(y ~ ref, d, d$id, control = list(5)), "unnamed")
  expect_error(angular_re(y ~ ref, d, d$id, init = list(kappa_unknown = 1)), "Unknown")
})
