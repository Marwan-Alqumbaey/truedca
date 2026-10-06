test_that("logistic working model keeps unbiasedness and full review is exact", {
  set.seed(11)
  n <- 3000
  risk <- runif(n, 0, 0.5); truth <- rbinom(n, 1, risk); y <- truth * rbinom(n, 1, 0.6)
  full <- twophase_dca(risk, y, truth, rep(1, n), c(0.1, 0.2), B = 20)
  tru <- nb_curve(risk, truth, c(0.1, 0.2))
  expect_equal(full$nb, tru$nb)
  p <- 0.15; tru1 <- nb_curve(risk, truth, p)$delta_all
  prob <- rep(0.2, n)
  est <- replicate(150, {
    seen <- ifelse(rbinom(n, 1, prob) == 1, truth, NA)
    twophase_dca(risk, y, seen, prob, p, B = 2)$delta_all
  })
  expect_lt(abs(mean(est) - tru1), 3 * sd(est) / sqrt(150))
})

test_that("design_gain equals one for constant v * h and exceeds one otherwise", {
  set.seed(12)
  risk <- runif(2000, 0, 0.04); y <- rbinom(2000, 1, 0.02)
  g1 <- design_gain(risk, y, 200, c(0.05, 0.3), "all", v = rep(0.1, 2000))
  expect_equal(unname(g1["gain"]), 1, tolerance = 1e-6)
  risk2 <- runif(2000, 0, 0.5)
  g2 <- design_gain(risk2, y, 200, c(0.05, 0.3), "all", v = runif(2000, 0.01, 0.25))
  expect_gt(g2[["gain"]], 1)
  expect_gt(g2[["bound"]], 1)
})

test_that("bias analysis recovers truth when priors are tight and correct", {
  set.seed(13)
  n <- 2e5
  risk <- plogis(rnorm(n, -2.2, 1.1)); truth <- rbinom(n, 1, risk)
  y <- ifelse(truth == 1, rbinom(n, 1, 0.8), rbinom(n, 1, 0.02))
  b <- nb_bias_analysis(risk, y, c(0.1, 0.2), se0 = c(8000, 2000), sp0 = c(9800, 200), draws = 400)
  tru <- nb_curve(risk, truth, c(0.1, 0.2))
  expect_true(all(b$delta_lo <= tru$delta_all & tru$delta_all <= b$delta_hi))
  expect_true(all(b$kept > 0.9))
})

test_that("confidence limits widen the bounds", {
  set.seed(14)
  risk <- runif(4000, 0, 0.5); y <- rbinom(4000, 1, risk * 0.7)
  b <- nb_bounds(risk, y, c(0.1, 0.2), se0 = c(0.6, 1), sp0 = c(0.98, 1))
  expect_true(all(b$delta_lower_ci <= b$delta_lower & b$delta_upper_ci >= b$delta_upper))
  expect_true(all(b$se_tip_upper >= b$se_tip))
})

test_that("pilot_v is finite and small where records are certain", {
  set.seed(9)
  n <- 5000; risk <- runif(n, 0, 0.5); truth <- rbinom(n, 1, risk)
  y <- truth * rbinom(n, 1, 0.7)                      # every record is true
  pil <- ifelse(seq_len(n) %in% sample.int(n, 150), truth, NA)
  v <- pilot_v(risk, y, pil)
  expect_true(all(is.finite(v)) && all(v >= 0 & v <= 0.25))
  expect_lt(mean(v[y == 1]), mean(v[y == 0]))
  expect_error(pilot_v(risk, y, rep(NA, n)))
})
