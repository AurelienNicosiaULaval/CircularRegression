library(CircularRegression)

wrap_angle <- function(angle) atan2(sin(angle), cos(angle))

test_that("von Mises PIT values are finite and strictly interior", {
  skip_if_not_installed("ggcircular")

  residuals <- c(-pi, -1, 0, 0.5, pi)
  pit <- CircularRegression:::.vonmises_pit(residuals, kappa = 2)

  expect_length(pit, length(residuals))
  expect_true(all(is.finite(pit)))
  expect_true(all(pit > 0 & pit < 1))
})

test_that("angular and consensus diagnostics render through their S3 plot methods", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggcircular")
  skip_if_not_installed("gridExtra")

  set.seed(20260828)
  n <- 60
  ref <- runif(n, -pi, pi)
  data <- data.frame(
    y = wrap_angle(ref + stats::rnorm(n, sd = 0.25)),
    ref = ref
  )
  angular_fit <- angular(y ~ ref, data = data)
  consensus_fit <- consensus(y ~ ref, data = data)

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 9, height = 7)
  on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)

  expect_silent(plot(angular_fit, type = "all", bins = 12))
  expect_silent(plot(consensus_fit, type = "all", bins = 12))
  grDevices::dev.off()
  expect_gt(file.info(output)$size, 0)
})

test_that("angular_re residual roses retain fixed and conditional diagnostics", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggcircular")
  skip_if_not_installed("gridExtra")

  data(Sandhopper, package = "CircularRegression")
  wide <- Sandhopper[stats::complete.cases(Sandhopper[c("LN1", "LN2", "LN3", "LN4", "LN5", "Azimuth", "DirW", "SpeedW")]), ]
  wide <- wide[seq_len(24), ]
  sandhopper <- wide[rep(seq_len(nrow(wide)), each = 5L), ]
  sandhopper$y <- as.vector(t(as.matrix(wide[paste0("LN", 1:5)]))) * pi / 180
  sandhopper$ref <- sandhopper$Azimuth * pi / 180
  sandhopper$wind <- sandhopper$DirW * pi / 180
  sandhopper$wind_speed <- sandhopper$SpeedW / max(sandhopper$SpeedW, na.rm = TRUE)
  fit <- angular_re(
    y ~ ref + wind:wind_speed,
    data = sandhopper,
    cluster = sandhopper$Anim,
    control = list(maxit = 50, reltol = 1e-8)
  )

  panels <- CircularRegression:::.angular_re_residual_plots(
    x = fit,
    breaks = 12,
    bandwidth = NULL,
    points = TRUE,
    main = "Residuals"
  )
  expect_equal(panels$settings$bins, 12L)
  expect_true(is.finite(panels$settings$bandwidth))
  expect_gt(panels$settings$bandwidth, 0)
  expect_equal(
    panels$fixed$layers[[2]]$stat_params$bw,
    panels$settings$bandwidth
  )
  expect_equal(
    panels$conditional$layers[[2]]$stat_params$bw,
    panels$settings$bandwidth
  )
  expect_equal(
    panels$fixed$scales$get_scales("y")$limits,
    panels$conditional$scales$get_scales("y")$limits
  )
  expect_equal(
    panels$fixed$scales$get_scales("y")$limits[[2]],
    panels$settings$radial_limit
  )
  expect_error(
    CircularRegression:::.angular_re_residual_plots(
      x = fit,
      breaks = 12,
      bandwidth = 0,
      points = TRUE,
      main = "Residuals"
    ),
    "finite positive"
  )

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 9, height = 4.5)
  on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)

  expect_silent(plot(fit, which = "both", breaks = 12, bandwidth = 0.4))
  grDevices::dev.off()
  expect_gt(file.info(output)$size, 0)
})
