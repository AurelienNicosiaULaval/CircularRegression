library(CircularRegression)

wrap_angle <- function(a) atan2(sin(a), cos(a))

.full_coef <- function(fit_hom) {
  beta <- as.numeric(fit_hom$parameters[, "estimate"])
  out <- numeric(length(fit_hom$term_labels))
  out[fit_hom$ref_idx] <- 1
  out[-fit_hom$ref_idx] <- beta
  out
}

test_that("select_reference_angle applies article criterion", {
  set.seed(900)
  n <- 50
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  y <- wrap_angle(x2 + rnorm(n, sd = 0.03))
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  sel <- select_reference_angle(y ~ x1 + x2, data = df)
  expect_s3_class(sel, "select_reference_angle")
  expect_equal(sel$ref_name, "x2")
})

test_that("angular_two_step returns coherent pipeline object", {
  set.seed(901)
  n <- 90
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  z2 <- runif(n, 0.2, 2)
  y <- wrap_angle(x1 + 0.25 * sin(x2) * z2 + rnorm(n, sd = 0.2))
  df <- data.frame(y = y, x1 = x1, x2 = x2, z2 = z2)

  fit <- angular_two_step(y ~ x1 + x2:z2, data = df)
  expect_s3_class(fit, "angular_two_step")
  expect_s3_class(fit$consensus_fit, "consensus")
  expect_s3_class(fit$homogeneous_fit, "angular")
  expect_true(is.character(fit$reference))
})

test_that("legacy pick_reference_angle remains available as deprecated wrapper", {
  set.seed(902)
  n <- 40
  df <- data.frame(y = runif(n, -pi, pi), x1 = runif(n, -pi, pi), x2 = runif(n, -pi, pi))

  expect_warning(
    out <- pick_reference_angle(y ~ x1 + x2, data = df),
    regexp = "deprecated",
    ignore.case = TRUE
  )
  expect_true(inherits(out, "pick_reference_angle"))
  expect_true(inherits(out, "select_reference_angle"))
})

test_that("special-case wrapper returns two-step components", {
  set.seed(903)
  n <- 80
  w <- runif(n, -pi, pi)
  y <- wrap_angle(w + rnorm(n, sd = 0.25))
  df <- data.frame(y = y, w = w)

  out <- decentredPredictorModel(df, response = "y", w = "w")
  expect_true(all(c(
    "consensus_fit",
    "homogeneous_fit",
    "reference",
    "data_aug",
    "formula",
    "natural_parameters",
    "natural_vcov_model",
    "natural_vcov_robust",
    "natural_jacobian"
  ) %in% names(out)))
  expect_s3_class(out$consensus_fit, "consensus")
  expect_s3_class(out$homogeneous_fit, "angular")
})

test_that("wrapper natural parameters are returned and match homogeneous coefficients", {
  set.seed(904)
  n <- 120

  # Mean direction
  y_mean <- wrap_angle(0.35 + rnorm(n, sd = 0.25))
  fit_mean <- meanDirectionModel(data.frame(y = y_mean), response = "y")
  b_mean <- .full_coef(fit_mean$homogeneous_fit)
  mu_manual <- wrap_angle(atan2(b_mean[2], b_mean[1]))
  expect_equal(rownames(fit_mean$natural_parameters), "mu")
  expect_equal(as.numeric(fit_mean$natural_parameters["mu", "estimate"]), mu_manual, tolerance = 1e-8)
  expect_true(all(c("estimate", "se_model", "se_robust") %in% colnames(fit_mean$natural_parameters)))

  # Decentred
  w_dec <- runif(n, -pi, pi)
  y_dec <- wrap_angle(w_dec + 0.2 * sin(w_dec) + rnorm(n, sd = 0.25))
  fit_dec <- decentredPredictorModel(data.frame(y = y_dec, w = w_dec), response = "y", w = "w")
  b_dec <- .full_coef(fit_dec$homogeneous_fit)
  theta1 <- atan2(b_dec[2], b_dec[1])
  theta2 <- atan2(b_dec[4], b_dec[3])
  dec_manual <- c(
    r_hat = sqrt(b_dec[3]^2 + b_dec[4]^2) / sqrt(b_dec[1]^2 + b_dec[2]^2),
    alpha_hat = wrap_angle(theta2 - theta1),
    beta_hat = wrap_angle(theta2)
  )
  expect_equal(
    as.numeric(fit_dec$natural_parameters[, "estimate"]),
    as.numeric(dec_manual),
    tolerance = 1e-8
  )

  # Presnell
  z <- runif(n, 0, 2)
  mu_pres <- atan2(0.3 + 0.2 * z, 1 + 0.4 * z)
  y_pres <- wrap_angle(mu_pres + rnorm(n, sd = 0.2))
  fit_pres <- presnellModel(data.frame(y = y_pres, z = z), response = "y", w = "z")
  b_pres <- .full_coef(fit_pres$homogeneous_fit)
  pres_manual <- c(
    beta2 = b_pres[2] / b_pres[1],
    beta3 = b_pres[3] / b_pres[1],
    beta4 = b_pres[4] / b_pres[1]
  )
  expect_equal(
    as.numeric(fit_pres$natural_parameters[, "estimate"]),
    as.numeric(pres_manual),
    tolerance = 1e-8
  )

  # JamSen linear map
  w_js <- runif(n, -pi, pi)
  mu_js <- atan2(
    0.2 + 0.45 * sin(w_js) + 0.15 * cos(w_js),
    1 + 0.35 * cos(w_js) - 0.1 * sin(w_js)
  )
  y_js <- wrap_angle(mu_js + rnorm(n, sd = 0.25))
  fit_js <- jammalamadakaModel(data.frame(y = y_js, w = w_js), response = "y", w = "w")
  b_js <- .full_coef(fit_js$homogeneous_fit)
  b2 <- b_js[2] / b_js[1]
  b3 <- b_js[3] / b_js[1]
  b4 <- b_js[4] / b_js[1]
  b5 <- b_js[5] / b_js[1]
  b6 <- b_js[6] / b_js[1]
  js_manual <- c(
    gamma2 = b2,
    gamma3 = b3 - b5,
    gamma4 = b4 + b6,
    gamma5 = b3 + b5,
    gamma6 = -b4 + b6
  )
  expect_equal(
    as.numeric(fit_js$natural_parameters[, "estimate"]),
    as.numeric(js_manual),
    tolerance = 1e-8
  )
})
