test_that("full review reproduces the true curve", {
  set.seed(3)
  n <- 3000
  risk <- runif(n, 0, 0.5); truth <- rbinom(n, 1, risk); y <- truth * rbinom(n, 1, 0.6)
  fit <- twophase_dca(risk, y, truth, rep(1, n), c(0.1, 0.2), B = 100)
  tru <- nb_curve(risk, truth, c(0.1, 0.2))
  expect_equal(fit$nb, tru$nb)
  expect_equal(fit$delta_all, tru$delta_all)
})

test_that("two-phase estimator is approximately unbiased with valid intervals", {
  set.seed(4)
  n <- 4000; p <- 0.15
  risk <- runif(n, 0, 0.5); truth <- rbinom(n, 1, risk); y <- truth * rbinom(n, 1, 0.6)
  tru <- nb_curve(risk, truth, p)$delta_all
  prob <- ifelse(risk < p, 0.25, 0.05)
  est <- cover <- numeric(300)
  for (i in 1:300) {
    seen <- ifelse(rbinom(n, 1, prob) == 1, truth, NA)
    f <- twophase_dca(risk, y, seen, prob, p, strata = make_strata(risk, y, breaks = p), B = 50)
    est[i] <- f$delta_all
    cover[i] <- f$delta_lo <= tru && tru <= f$delta_hi
  }
  expect_lt(abs(mean(est) - tru), 3 * sd(est) / sqrt(300))
  expect_gt(mean(cover), 0.9)
})

test_that("design probabilities respect the budget and the floor", {
  set.seed(5)
  risk <- runif(1000, 0, 0.6); y <- rbinom(1000, 1, risk)
  v <- rep(0.1, 1000)
  pr <- twophase_design(risk, y, 150, c(0.05, 0.3), "all", v = v, floor = 0.02, defensive = 0)
  expect_equal(sum(pr), 150, tolerance = 1e-6)
  expect_true(all(pr >= 0.02 & pr <= 1))
  expect_true(all(pr[risk >= 0.3] == 0.02))
  expect_true(mean(pr[risk < 0.05]) > mean(pr[risk > 0.2 & risk < 0.3]))
  pd <- twophase_design(risk, y, 150, c(0.05, 0.3), "all", v = v, floor = 0.02, defensive = 0.5)
  expect_equal(sum(pd), 150, tolerance = 1e-6)
  expect_true(all(pd >= 0.5 * 150 / 1000))
})

test_that("split interval keeps coverage when hidden events are rare", {
  set.seed(11)
  n <- 3000; risk <- rbeta(n, 1, 8); truth <- rbinom(n, 1, risk)
  y <- truth * rbinom(n, 1, ifelse(risk > 0.1, 0.9, 0.7))   # missed cases only
  tr <- nb_curve(risk, truth, c(0.05, 0.1))$delta_all
  hit <- replicate(100, {
    prob <- rep(0.03, n); tt <- ifelse(rbinom(n, 1, prob) == 1, truth, NA)
    f <- suppressWarnings(twophase_dca(risk, y, tt, prob, c(0.05, 0.1), B = 2))
    f$delta_lo <= tr & tr <= f$delta_hi
  })
  expect_true(all(rowMeans(hit) >= 0.9))
})

test_that("split interval includes phase-one error (full review)", {
  set.seed(12)
  big <- rbeta(1e6, 2, 8); tr <- nb_curve(big, rbinom(1e6, 1, big), 0.2)
  hit <- replicate(100, {
    n <- 1500; r <- rbeta(n, 2, 8); t <- rbinom(n, 1, r)
    y <- ifelse(t == 1, rbinom(n, 1, 0.9), rbinom(n, 1, 0.01))
    f <- suppressWarnings(twophase_dca(r, y, t, rep(1, n), 0.2, B = 2))
    c(f$nb_lo <= tr$nb & tr$nb <= f$nb_hi, f$delta_lo <= tr$delta_all & tr$delta_all <= f$delta_hi)
  })
  expect_true(all(rowMeans(hit) >= 0.88))
})
