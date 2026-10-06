test_that("Corollary 1 identity holds exactly for the package quantities", {
  set.seed(11)
  n <- 5000
  risk <- runif(n)
  y <- rbinom(n, 1, 0.2 + 0.6 * risk)
  se <- 0.9; sp <- 0.95
  p <- c(0.1, 0.25, 0.5)
  obs <- nb_curve(risk, y, p)
  adj <- suppressWarnings(nb_adjust(risk, y, p, se = se, sp = sp))
  d <- warp_threshold(p, se, sp)
  cp <- (se - p) / (1 - p)
  # NB* = c(p) NB_t(delta), D* = c(p) D_t(delta), with NB_t, D_t the true
  # (corrected) curves of the same classifier scored at threshold delta
  expect_equal(obs$nb, cp * adj$q * (adj$pi1 - d) / (1 - d))
  expect_equal(obs$delta_all, cp * (1 - adj$q) * (d - adj$pi0) / (1 - d))
})

test_that("ties, logical outcomes, unsorted and extreme thresholds", {
  risk <- c(0, 0.2, 0.2, 0.5, 1)
  y <- c(FALSE, TRUE, FALSE, TRUE, TRUE)
  a <- nb_curve(risk, y, c(0.2, 0.9))
  b <- nb_curve(risk, as.numeric(y), c(0.9, 0.2))
  expect_equal(a$q[1], 4 / 5)                      # risk == p counts as flagged
  expect_equal(as.data.frame(a), as.data.frame(b[2:1, ]), ignore_attr = TRUE)
  # q = 1 and q = 0
  all_flag <- nb_curve(rep(0.6, 4), c(0, 1, 1, 0), 0.5)
  expect_equal(all_flag$delta_all, 0)
  expect_true(is.na(all_flag$y0))
  none_flag <- nb_curve(rep(0.1, 4), c(0, 1, 1, 0), 0.5)
  expect_equal(none_flag$nb, 0)
  adj <- nb_adjust(rep(0.1, 4), c(0, 1, 1, 0), 0.5, se = 0.9, sp = 0.9)
  expect_equal(adj$nb, 0)
  expect_equal(adj$nb_all, -adj$delta_all)
  expect_error(nb_curve(risk, factor(y), 0.2))
  expect_error(nb_curve(risk, y, c(0.2, NA)))
  expect_error(nb_curve(risk, y, numeric(0)))
})

test_that("nb_bounds: tipping points match corrections, conflicts give NA", {
  set.seed(12)
  n <- 20000
  risk <- runif(n)
  y <- rbinom(n, 1, 0.6 * risk)
  p <- c(0.1, 0.3)
  b <- nb_bounds(risk, y, p, se0 = c(0.7, 1), sp0 = c(0.97, 1), tip_sp0 = 0.98, tip_se1 = 0.8)
  # at the tipping points the corrected comparison is exactly zero
  for (k in seq_along(p)) {
    a0 <- suppressWarnings(nb_adjust(risk, y, p[k], se0 = b$se_tip[k], sp0 = 0.98, se1 = 0.9, sp1 = 0.9))
    expect_equal(a0$delta_all, 0, tolerance = 1e-10)
    a1 <- suppressWarnings(nb_adjust(risk, y, p[k], se1 = 0.8, sp1 = b$sp_tip[k], se0 = 0.9, sp0 = 0.9))
    expect_equal(a1$nb, 0, tolerance = 1e-10)
  }
  # degenerate intervals reproduce nb_adjust
  e <- nb_bounds(risk, y, p, se0 = c(0.8, 0.8), sp0 = c(0.99, 0.99))
  a <- nb_adjust(risk, y, p, se = 0.8, sp = 0.99)
  expect_equal(e$nb_lower, a$nb)
  expect_equal(e$delta_upper, a$delta_all)
  # recorded rate among flagged exceeds the largest admissible sensitivity
  c1 <- nb_bounds(risk, y, 0.5, se0 = c(0.6, 1), sp0 = c(0.95, 1),
                  se1 = c(0.2, 0.3), sp1 = c(0.95, 1))
  expect_true(c1$conflict)
  expect_true(is.na(c1$nb_lower) && !c1$beats_none)
  expect_false(is.na(c1$delta_lower))
})

test_that("twophase_dca input checks and empty-stratum fallback", {
  set.seed(13)
  n <- 600
  risk <- runif(n, 0, 0.5); truth <- rbinom(n, 1, risk); y <- truth * rbinom(n, 1, 0.7)
  prob <- rep(0.5, n)
  seen <- ifelse(rbinom(n, 1, prob) == 1, truth, NA)
  expect_error(twophase_dca(risk, y, seen, replace(prob, which(is.na(seen))[1], 1), 0.2, B = 20))
  expect_error(twophase_dca(risk, y, replace(seen, which(!is.na(seen))[1], 2), prob, 0.2, B = 20))
  expect_error(twophase_dca(risk, y, seen, prob, 0.2, strata = replace(y, 1, NA), B = 20))
  st <- factor(ifelse(risk > 0.45, "top", "rest"))
  seen2 <- replace(seen, st == "top", NA)
  expect_warning(f <- twophase_dca(risk, y, seen2, prob, c(0.2, 0.3), strata = st, B = 50, working = "strata"),
                 "no reviewed")
  expect_true(all(is.finite(f$nb) & is.finite(f$nb_se)))
  # single threshold; simultaneous band equals pointwise band in distribution
  g <- twophase_dca(risk, y, seen, prob, 0.2, B = 400)
  expect_equal(nrow(g), 1)
  expect_lt(abs(attr(g, "crit")[["nb"]] - qnorm(0.975)), 0.25)
  # full review: influence-function SE is the usual SE of a mean (up to the n_h/(n_h - 1) correction)
  h <- twophase_dca(risk, y, truth, rep(1, n), 0.2, B = 20)
  a <- (risk >= 0.2) * (truth - 0.2) / 0.8
  expect_equal(h$nb_se, sqrt(mean((a - mean(a))^2) / n), tolerance = 0.01)
})

test_that("twophase_design matches h from numerical integration and handles big budgets", {
  pl <- 0.05; pu <- 0.3
  h_num <- function(r, side) {
    lim <- if (side == "all") c(max(r, pl), pu) else c(pl, min(r, pu))
    if (lim[1] >= lim[2]) return(0)
    stats::integrate(function(p) 1 / (1 - p)^2 / (pu - pl), lim[1], lim[2])$value
  }
  risk <- c(0.01, 0.1, 0.2, 0.29, 0.31, 0.6)
  y <- rep(0, 6)
  for (side in c("all", "none")) {
    h <- vapply(risk, h_num, numeric(1), side = side)
    pr <- twophase_design(risk, y, 0.5, c(pl, pu), side, v = rep(1, 6),
                          floor = 1e-6, defensive = 0)
    ok <- h > 0
    expect_equal(pr[ok] / sum(pr[ok]), sqrt(h[ok]) / sum(sqrt(h[ok])), tolerance = 1e-5)
  }
  # budget larger than the number of informative patients: no uniroot error
  big <- twophase_design(risk, y, 5, c(pl, pu), "all", v = rep(1, 6), defensive = 0)
  expect_equal(sum(big), 5)
  expect_true(all(big[risk < pu] == 1))
  expect_error(twophase_design(risk, y, 2, c(pl, pu), v = rep(1, 6), floor = -0.1))
})
