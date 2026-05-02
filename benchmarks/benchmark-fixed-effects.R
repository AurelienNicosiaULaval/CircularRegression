# Lightweight timing script for fixed-effect circular regression workflows.
# This script is not run during R CMD check.

library(CircularRegression)

wrap_angle <- function(x) atan2(sin(x), cos(x))

simulate_data <- function(n, beta = 0.35, sd = 0.08) {
  x1 <- runif(n, -pi, pi)
  x2 <- runif(n, -pi, pi)
  z2 <- runif(n, 0.2, 1.8)
  mu <- atan2(
    sin(x1) + beta * z2 * sin(x2),
    cos(x1) + beta * z2 * cos(x2)
  )
  y <- wrap_angle(mu + rnorm(n, sd = sd))
  data.frame(y = y, x1 = x1, x2 = x2, z2 = z2)
}

set.seed(20260502)
sizes <- c(100, 500, 1000)

results <- do.call(rbind, lapply(sizes, function(n) {
  dat <- simulate_data(n)
  formula <- y ~ x1 + x2:z2

  t_angular <- system.time(
    angular(formula, data = dat, reference = c("name", "x1"))
  )[["elapsed"]]

  t_consensus <- system.time(
    consensus(formula, data = dat)
  )[["elapsed"]]

  t_two_step <- system.time(
    circular_regression(formula, data = dat)
  )[["elapsed"]]

  data.frame(
    n = n,
    angular_seconds = t_angular,
    consensus_seconds = t_consensus,
    two_step_seconds = t_two_step
  )
}))

print(results)
