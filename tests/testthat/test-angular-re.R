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

test_that("summary.angular_re is available and coherent", {
  data(Sandhopper)
  wide <- Sandhopper[stats::complete.cases(Sandhopper[c("LN1", "LN2", "LN3", "LN4", "LN5", "Azimuth", "DirW", "SpeedW")]), ]
  wide <- wide[seq_len(min(24L, nrow(wide))), ]
  sh <- wide[rep(seq_len(nrow(wide)), each = 5L), ]
  sh$y <- as.vector(t(as.matrix(wide[paste0("LN", 1:5)]))) * pi / 180
  sh$ref <- sh$Azimuth * pi / 180
  sh$wind <- sh$DirW * pi / 180
  sh$wind_speed <- sh$SpeedW / max(sh$SpeedW, na.rm = TRUE)

  fit <- angular_re(
    y ~ ref + wind:wind_speed,
    data = sh,
    cluster = sh$Anim,
    control = list(maxit = 50, reltol = 1e-8)
  )

  s <- summary(fit)
  expect_s3_class(s, "summary.angular_re")
  expect_true(is.finite(s$logLik))
  expect_true("coefficients" %in% names(s))
  expect_true("n_clusters" %in% names(s))
  expect_equal(s$nobs, nrow(sh))
  expect_equal(s$n_clusters, length(unique(sh$Anim)))
})

test_that("grouping formulas match vector fits after missing-value handling", {
  data(Sandhopper)
  trials <- paste0("LN", 1:5)
  columns <- c(trials, "Azimuth", "DirW", "SpeedW")
  wide <- head(Sandhopper[stats::complete.cases(Sandhopper[columns]), ], 24)
  d <- wide[rep(seq_len(nrow(wide)), each = length(trials)), ]
  d$y <- as.vector(t(as.matrix(wide[trials]))) * pi / 180
  d$ref <- d$Azimuth * pi / 180
  d$wind <- d$DirW * pi / 180
  d$speed <- d$SpeedW / max(d$SpeedW)
  d[["animal id"]] <- factor(d$Anim, levels = c(rev(unique(d$Anim)), "unused"))
  d[["animal id"]][1] <- NA
  d$y[6] <- NA
  kept <- d[-c(1, 6), ]

  grouped <- angular_re(y ~ (1 | `animal id`) + ref + wind:speed, data = d,
                        na.action = stats::na.exclude)
  vector_fit <- angular_re(y ~ ref + wind:speed, data = kept,
                           cluster = kept[["animal id"]])
  expect_true(grouped$convergence)
  expect_equal(coef(grouped), coef(vector_fit))
  expect_equal(grouped$logLik, vector_fit$logLik)
  expect_identical(grouped$cluster, vector_fit$cluster)

  conditional <- predict(grouped, type = "conditional")
  expect_identical(which(is.na(conditional)), c(1L, 6L))
  expect_equal(
    predict(grouped, newdata = kept, type = "conditional"),
    predict(vector_fit, newdata = kept, cluster = kept[["animal id"]])
  )
  expect_equal(predict(grouped, newdata = kept[c("ref", "wind", "speed")]),
               predict(vector_fit))
  vector_fit$fixed_formula <- NULL
  expect_equal(predict(vector_fit, type = "conditional"), conditional[-c(1, 6)])
  expect_error(angular_re(y ~ ref + (1 | `animal id`), d, na.action = stats::na.fail),
               "missing values")
  expect_error(angular_re(y ~ ref + (wind | Anim), d),
               "must have the form", fixed = TRUE)
  expect_error(angular_re(y ~ ref + (1 | Anim) + (1 | `animal id`), d),
               "only one grouping term", fixed = TRUE)
  expect_error(angular_re(y ~ ref + (1 | Anim), d, cluster = d$Anim),
               "but not both", fixed = TRUE)
})
