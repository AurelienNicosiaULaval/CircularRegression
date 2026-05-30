library(CircularRegression)

test_that("package datasets support fixed-effect workflows", {
  data(bison)
  d <- bison[seq_len(180), c("y.dir", "y.prec", "x.meadow", "z.meadow")]
  d <- na.omit(d)

  sel <- select_reference_angle(y.dir ~ y.prec + x.meadow:z.meadow, data = d)
  expect_s3_class(sel, "select_reference_angle")
  expect_true(sel$ref_name %in% c("y.prec"))

  fit_hom <- angular(
    y.dir ~ y.prec + x.meadow:z.meadow,
    data = d,
    reference = c("name", sel$ref_name)
  )
  fit_cons <- consensus(y.dir ~ y.prec + x.meadow:z.meadow, data = d)
  fit_two <- angular_two_step(y.dir ~ y.prec + x.meadow:z.meadow, data = d)

  expect_s3_class(fit_hom, "angular")
  expect_s3_class(fit_cons, "consensus")
  expect_s3_class(fit_two, "angular_two_step")
  expect_equal(length(predict(fit_hom, newdata = d[1:5, ])), 5)
  expect_equal(length(predict(fit_cons, newdata = d[1:5, ])), 5)
  expect_equal(length(residuals(fit_two)), fit_two$nobs)
})

test_that("noshiro data fit a reference-only homogeneous model", {
  data(noshiro)
  d <- noshiro[seq_len(120), ]

  fit <- angular(DIRMV ~ DIRDSC, data = d, reference = "first")

  expect_s3_class(fit, "angular")
  expect_equal(length(coef(fit)), 0)
  expect_equal(length(predict(fit)), nrow(d))
  expect_true(all(is.finite(residuals(fit))))
})

test_that("special-case wrappers expose natural-parameter outputs on package data", {
  data(bison)
  d <- na.omit(bison[seq_len(80), c("y.dir", "y.prec", "z.gap")])

  mean_fit <- meanDirectionModel(d["y.dir"], response = "y.dir")
  dec_fit <- decentredPredictorModel(d[c("y.dir", "y.prec")], response = "y.dir", w = "y.prec")
  d_pres <- na.omit(bison[500:579, c("y.dir", "z.gap")])
  pres_fit <- presnellModel(d_pres, response = "y.dir", w = "z.gap")
  jam_fit <- jammalamadakaModel(d[c("y.dir", "y.prec")], response = "y.dir", w = "y.prec")

  for (fit in list(mean_fit, dec_fit, pres_fit, jam_fit)) {
    expect_true(all(c(
      "natural_parameters",
      "natural_vcov_model",
      "natural_vcov_robust",
      "natural_jacobian"
    ) %in% names(fit)))
    expect_true(all(c("estimate", "se_model", "se_robust") %in% colnames(fit$natural_parameters)))
  }

  expect_error(meanDirectionModel(d["y.dir"], response = "missing"), "response")
  expect_error(meanDirectionModel(d["y.dir"], response = "y.dir", mu0 = c(1, 2)), "mu0")
})

test_that("angular_re methods are coherent on a small Sandhopper example", {
  data(Sandhopper)
  sh <- Sandhopper[seq_len(24), ]
  sh$y <- sh$LN1 * pi / 180
  sh$ref <- sh$Azimuth * pi / 180
  sh$wind <- sh$DirW * pi / 180
  sh$wind_speed <- sh$SpeedW / max(sh$SpeedW, na.rm = TRUE)

  fit <- angular_re(
    y ~ ref + wind:wind_speed,
    data = sh,
    cluster = sh$Anim,
    control = list(maxit = 50, reltol = 1e-8)
  )

  expect_s3_class(fit, "angular_re")
  expect_true(isTRUE(fit$convergence))
  expect_named(coef(fit), c("wind:wind_speed", "kappa_e", "kappa_a"))
  expect_equal(dim(vcov(fit)), c(3L, 3L))
  expect_equal(dim(vcov(fit, robust = TRUE)), c(3L, 3L))
  expect_equal(length(fitted(fit)), nrow(sh))
  expect_equal(length(residuals(fit, type = "fixed")), nrow(sh))
  expect_equal(length(residuals(fit, type = "conditional")), nrow(sh))
  expect_equal(length(predict(fit)), nrow(sh))
  expect_equal(length(ranef.angular_re(fit)), length(unique(sh$Anim)))

  bad_new <- sh[1:3, c("ref", "wind", "wind_speed")]
  bad_new$wind_speed[1] <- -1
  expect_error(predict(fit, newdata = bad_new), "non-negative")

  expect_error(
    angular_re(y ~ ref + wind:wind_speed, data = sh, cluster = sh$Anim, control = list(hessian = NA)),
    "hessian"
  )
})
