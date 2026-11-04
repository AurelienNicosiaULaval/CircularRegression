test_that("T2 checks ordered direction/modifier roles after R reorders terms", {
  set.seed(731)
  n <- 160
  d <- data.frame(ref = 0, x = runif(n, .1, 2), z = runif(n, .1, 2))
  d$y <- atan2(.4 * sin(d$x) + .7 * d$x * sin(d$z),
               1 + .4 * cos(d$x) + .7 * d$x * cos(d$z)) + .12 * cos(seq_len(n))
  full <- angular(y ~ ref + x + z:x, d, reference = "first")
  wrong <- angular(y ~ ref + x:z, d, reference = "first")
  nested <- angular(y ~ ref + x, d, reference = "first")
  expect_true(all(wrong$term_labels %in% full$term_labels))
  expect_error(angular_lrtest(full, wrong), "direction/modifier roles")
  expect_error(anova(full, wrong, test = "T2"), "direction/modifier roles")
  expect_s3_class(angular_lrtest(full, nested), "htest")
})
