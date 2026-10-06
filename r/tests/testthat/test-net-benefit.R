test_that("nb_curve matches hand calculation", {
  risk <- c(0.1, 0.3, 0.6, 0.8)
  y <- c(0, 1, 0, 1)
  out <- nb_curve(risk, y, 0.5)
  expect_equal(out$nb, (1 * 0.5 + (0 - 0.5)) / 4 / 0.5)
  expect_equal(out$nb_all, (0.5 - 0.5) / 0.5)
  expect_equal(out$nb - out$nb_all, out$delta_all)
})

test_that("window and warp behave as stated", {
  w <- dca_window(0.9, 0.97)
  expect_equal(unname(w["floor"]), 0.03)
  expect_equal(unname(w["ceiling"]), 0.9)
  expect_equal(warp_threshold(w[["fixed_point"]], 0.9, 0.97), w[["fixed_point"]])
  expect_error(dca_window(0.4, 0.5))
})

test_that("no model beats treat all below 1 - sp on recorded labels", {
  set.seed(1)
  n <- 2e5
  truth <- rbinom(n, 1, 0.1)
  y <- ifelse(truth == 1, rbinom(n, 1, 0.9), rbinom(n, 1, 0.05))
  oracle <- nb_curve(truth, y, c(0.01, 0.03, 0.049))
  expect_true(all(oracle$delta_all <= 0))
})

test_that("input checks work", {
  expect_error(nb_curve(c(0.2, 1.2), c(0, 1), 0.1))
  expect_error(nb_curve(c(0.2, 0.3), c(0, 2), 0.1))
  expect_error(nb_curve(c(0.2, 0.3), c(0, 1), 1))
})
