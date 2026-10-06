sim <- function(n = 4e5, seed = 2) {
  set.seed(seed)
  risk <- plogis(rnorm(n, -2.2, 1.1))
  truth <- rbinom(n, 1, risk)
  list(risk = risk, truth = truth)
}

test_that("region decomposition is an exact identity", {
  d <- sim(5e4)
  flag <- d$risk >= 0.1
  # differential error: worse recording in unflagged patients
  se_i <- ifelse(flag, 0.95, 0.5)
  y <- ifelse(d$truth == 1, rbinom(length(flag), 1, se_i), rbinom(length(flag), 1, 0.01))
  se0 <- mean(y[!flag & d$truth == 1]); sp0 <- mean(1 - y[!flag & d$truth == 0])
  se1 <- mean(y[flag & d$truth == 1]); sp1 <- mean(1 - y[flag & d$truth == 0])
  adj <- nb_adjust(d$risk, y, 0.1, se0 = se0, sp0 = sp0, se1 = se1, sp1 = sp1)
  tru <- nb_curve(d$risk, d$truth, 0.1)
  expect_equal(adj$nb, tru$nb, tolerance = 1e-10)
  expect_equal(adj$delta_all, tru$delta_all, tolerance = 1e-10)
})

test_that("correction recovers true net benefit under non-differential error", {
  d <- sim()
  y <- ifelse(d$truth == 1, rbinom(length(d$truth), 1, 0.8), rbinom(length(d$truth), 1, 0.03))
  p <- c(0.05, 0.1, 0.2)
  adj <- nb_adjust(d$risk, y, p, se = 0.8, sp = 0.97)
  tru <- nb_curve(d$risk, d$truth, p)
  expect_lt(max(abs(adj$nb - tru$nb)), 0.003)
  expect_lt(max(abs(adj$delta_all - tru$delta_all)), 0.003)
})

test_that("bounds contain truth, collapse when rates are known, and tipping points agree", {
  d <- sim(1e5)
  y <- d$truth * rbinom(length(d$truth), 1, 0.7)
  p <- c(0.05, 0.1, 0.2)
  b <- nb_bounds(d$risk, y, p, se0 = c(0.6, 0.9), sp0 = c(1, 1))
  tru <- nb_curve(d$risk, d$truth, p)
  expect_true(all(b$delta_lower <= tru$delta_all + 0.003 & tru$delta_all <= b$delta_upper + 0.003))
  k <- nb_bounds(d$risk, y, p, se0 = c(0.7, 0.7), sp0 = c(1, 1))
  expect_equal(k$delta_lower, k$delta_upper)
  obs <- nb_curve(d$risk, y, p)
  expect_equal(b$se_tip, obs$y0 / p)
  expect_equal(b$beats_all, 0.6 > b$se_tip)
})

test_that("flagged-group sensitivity tipping point is exact", {
  set.seed(21)
  n <- 50000; risk <- runif(n); truth <- rbinom(n, 1, risk)
  y <- truth * rbinom(n, 1, 0.3)                      # Se 0.3, Sp 1 everywhere
  b <- nb_bounds(risk, y, c(0.2, 0.5), se0 = c(0.2, 1), sp0 = c(0.99, 1))
  tr <- nb_curve(risk, truth, c(0.2, 0.5))
  # true NB > 0 exactly when the actual sensitivity (0.3) is below se1_tip
  expect_equal(tr$nb > 0, 0.3 < b$se1_tip)
  expect_true(all(b$se1_tip_lower <= b$se1_tip & b$se1_tip <= b$se1_tip_upper))
})
