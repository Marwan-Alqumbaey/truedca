#' Strata for chart review designs
#'
#' Crosses risk groups with the recorded outcome.
#'
#' @param risk Predicted risks.
#' @param y Recorded outcome.
#' @param breaks Cut points for risk, or a number of quantile groups.
#' @return A factor.
#' @examples
#' make_strata(c(0.01, 0.2, 0.5, 0.9), c(0, 1, 0, 1), breaks = c(0.1, 0.4))
#' @export
make_strata <- function(risk, y, breaks = 5) {
  if (length(breaks) == 1) {
    cuts <- unique(stats::quantile(risk, probs = seq(0, 1, length.out = breaks + 1)))
    g <- if (length(cuts) < 2) factor(rep("all", length(risk)))
         else cut(risk, cuts, include.lowest = TRUE)
  } else {
    g <- cut(risk, c(-Inf, sort(breaks), Inf), right = FALSE)
  }
  interaction(g, factor(as.numeric(y), levels = c(0, 1)), drop = TRUE)
}

#' Decision curve from a two-phase chart review design
#'
#' Estimates true net benefit when the recorded outcome is available for
#' everyone and the true outcome only for a reviewed subsample drawn with
#' known probabilities. Uses an augmented inverse-probability-weighted
#' estimator with a working model for the true outcome (see `working`),
#' pointwise confidence intervals (see `interval`), and simultaneous bands
#' over the thresholds from a Gaussian multiplier bootstrap. The
#' cross-fitting folds and the bootstrap use the random number generator;
#' call [set.seed()] for reproducible results. Patients with
#' `risk >= threshold` count as flagged.
#'
#' @param risk Predicted risks for all patients.
#' @param y Recorded outcome for all patients.
#' @param truth True outcome (0/1 or logical), `NA` for patients not reviewed.
#' @param prob Probability that each patient was selected for review, in
#'   (0, 1], known by design. Patients with `prob = 1` must be reviewed.
#' @inheritParams nb_curve
#' @param strata Factor used for the working model `E[truth | strata]`;
#'   defaults to [make_strata()]. Use the design strata when available.
#' @param level Confidence level.
#' @param B Number of multiplier bootstrap draws for the simultaneous band.
#' @param working Working model for `P(truth = 1 | risk, y)`: `"logistic"`
#'   (default) fits a weighted logistic regression on logit risk, the recorded
#'   outcome and their interaction, cross-fitted over `folds`; `"strata"` uses
#'   weighted stratum means. Estimates stay unbiased either way; a better
#'   working model gives smaller standard errors. The logistic model falls
#'   back to stratum means when fewer than 40 charts or 10 events or
#'   non-events were reviewed.
#' @param folds Number of cross-fitting folds for the logistic working model.
#' @param interval `"split"` (default) splits flagged and unflagged patients
#'   by their recorded outcome, builds a Wilson score interval for the true
#'   event rate in each part (with an effective sample size capped at the
#'   number of reviewed charts in that part), and combines the two with the
#'   method of variance estimates recovery (MOVER). This keeps coverage when
#'   few or no hidden events (missed or false records) are reviewed.
#'   `"wilson"` applies one Wilson interval to each whole group; `"wald"`
#'   gives symmetric estimate plus-or-minus intervals.
#' @return A data frame with estimates, standard errors, pointwise limits
#'   (`*_lo`, `*_hi`) and simultaneous limits (`*_slo`, `*_shi`) for net
#'   benefit (`nb`) and the gain over treat all (`delta_all`). The
#'   simultaneous critical values are stored in attribute `crit`.
#' @examples
#' set.seed(4)
#' n <- 4000
#' risk <- runif(n, 0, 0.6)
#' truth <- rbinom(n, 1, risk)
#' y <- truth * rbinom(n, 1, 0.7)
#' prob <- ifelse(risk < 0.2, 0.3, 0.1)
#' seen <- ifelse(rbinom(n, 1, prob) == 1, truth, NA)
#' twophase_dca(risk, y, seen, prob, thresholds = c(0.1, 0.2, 0.3), B = 200)
#' @export
twophase_dca <- function(risk, y, truth, prob,
                         thresholds = seq(0.01, 0.5, by = 0.01),
                         strata = make_strata(risk, y), level = 0.95, B = 1000,
                         interval = c("split", "wilson", "wald"), working = c("logistic", "strata"),
                         folds = 5) {
  interval <- match.arg(interval)
  working <- match.arg(working)
  check_inputs(risk, y, thresholds)
  n <- length(risk)
  if (length(truth) != n || length(prob) != n || length(strata) != n)
    stop("`truth`, `prob` and `strata` must have the same length as `risk`.", call. = FALSE)
  if (!is.numeric(prob) || anyNA(prob) || any(prob <= 0 | prob > 1))
    stop("`prob` must lie in (0, 1].", call. = FALSE)
  if (is.factor(truth) || !all(suppressWarnings(as.numeric(truth[!is.na(truth)])) %in% c(0, 1)))
    stop("`truth` must be 0/1 or logical, with `NA` for patients not reviewed.", call. = FALSE)
  if (any(prob == 1 & is.na(truth)))
    stop("Patients selected with `prob` = 1 must have `truth` recorded.", call. = FALSE)
  if (anyNA(strata)) stop("`strata` must not contain missing values.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  r <- !is.na(truth)
  tt <- ifelse(r, as.numeric(truth), 0)
  strata <- droplevels(as.factor(strata))
  wsum <- tapply(r / prob, strata, sum)
  mu_h <- tapply(r * tt / prob, strata, sum) / wsum
  empty <- any(wsum == 0)
  mu_h[wsum == 0] <- sum(r * tt / prob) / sum(r / prob)
  mu <- as.numeric(mu_h[as.character(strata)])
  # small-sample correction n_h / (n_h - 1) for the within-stratum residuals (variance only)
  nh <- tapply(r, strata, sum)
  cf <- as.numeric(ifelse(nh > 1, nh / (nh - 1), 1)[as.character(strata)])
  fit_mu <- if (working == "logistic") logistic_mu(risk, y, r, tt, prob, folds) else NULL
  if (!is.null(fit_mu)) {
    mu <- fit_mu
    cf <- rep(sum(r) / max(sum(r) - 4, 1), n)
  }
  if (empty && is.null(fit_mu))
    warning("Some strata have no reviewed patients; their working mean uses all reviewed patients. ",
            "Estimates stay approximately unbiased but may be less precise.", call. = FALSE)
  tstar <- mu + r * (tt - mu) / prob
  tvar <- mu + r * (tt - mu) / prob * sqrt(cf)
  K <- length(thresholds)
  phi_nb <- phi_d <- matrix(0, n, K)
  est_nb <- est_d <- q <- n1 <- n0 <- numeric(K)
  yy <- as.numeric(y)
  # per threshold, group (flagged = 1, unflagged = 2) and recorded outcome (0, 1):
  # share of the group, estimated true event rate, its standard error, charts reviewed
  sub <- array(0, c(K, 2, 2, 4))
  for (k in seq_len(K)) {
    p <- thresholds[k]
    flag <- risk >= p
    for (g in 1:2) for (j in 0:1) {
      idx <- (if (g == 1) flag else !flag) & yy == j
      ng <- sum(if (g == 1) flag else !flag)
      if (any(idx)) {
        tv <- tvar[idx]
        sub[k, g, j + 1, ] <- c(sum(idx) / ng, mean(tstar[idx]), sqrt(sum((tv - mean(tv))^2)) / sum(idx), sum(r[idx]))
      }
    }
    a <- flag * (tstar - p) / (1 - p)
    b <- (1 - flag) * (p - tstar) / (1 - p)
    est_nb[k] <- mean(a); est_d[k] <- mean(b)
    av <- flag * (tvar - p) / (1 - p)
    bv <- (1 - flag) * (p - tvar) / (1 - p)
    phi_nb[, k] <- av - mean(av); phi_d[, k] <- bv - mean(bv)
    q[k] <- mean(flag); n1[k] <- sum(r & flag); n0[k] <- sum(r & !flag)
  }
  se_nb <- sqrt(colSums(phi_nb^2)) / n
  se_d <- sqrt(colSums(phi_d^2)) / n
  z <- stats::qnorm(1 - (1 - level) / 2)
  crit <- c(nb = sup_crit(phi_nb, se_nb, level, B),
            delta_all = sup_crit(phi_d, se_d, level, B))
  p <- thresholds
  if (interval == "split") {
    # MOVER over three independent sources of error in the event rate of group g:
    # the two parts (recorded outcome 0 or 1) and the phase-one share of each part
    lims <- function(g, zz) {
      cc0 <- sub[, g, 1, 1]; cc1 <- sub[, g, 2, 1]
      m <- cc0 * sub[, g, 1, 2] + cc1 * sub[, g, 2, 2]
      ng <- if (g == 1) q * n else (1 - q) * n
      share <- zz^2 * cc1 * (1 - cc1) * (sub[, g, 2, 2] - sub[, g, 1, 2])^2 / pmax(ng, 1)
      lo <- hi <- share
      for (j in 1:2) {
        cc <- sub[, g, j, 1]; pp <- sub[, g, j, 2]; ss <- sub[, g, j, 3]; nn <- sub[, g, j, 4]
        w <- wilson(pp, ss, nn, zz, cap = TRUE)
        lo <- lo + (cc * (pp - w[, 1]))^2; hi <- hi + (cc * (w[, 2] - pp))^2
      }
      cbind(m, pmax(m - sqrt(lo), 0), pmin(m + sqrt(hi), 1))
    }
    # map to the net benefit scale, adding the phase-one error in the share flagged q
    to_scale <- function(l, zz, flagged) {
      k <- if (flagged) q / (1 - p) else (1 - q) / (1 - p)
      sgn <- if (flagged) 1 else -1
      est <- k * sgn * (l[, 1] - p)
      a <- k * sgn * (l[, 2] - p); b <- k * sgn * (l[, 3] - p)
      cond <- cbind(pmin(a, b), pmax(a, b))
      vq <- zz^2 * q * (1 - q) / n * ((l[, 1] - p) / (1 - p))^2
      cbind(est - sqrt((est - cond[, 1])^2 + vq), est + sqrt((cond[, 2] - est)^2 + vq))
    }
    nb_pw <- to_scale(lims(1, z), z, TRUE); nb_sb <- to_scale(lims(1, crit[["nb"]]), crit[["nb"]], TRUE)
    d_pw <- to_scale(lims(2, z), z, FALSE); d_sb <- to_scale(lims(2, crit[["delta_all"]]), crit[["delta_all"]], FALSE)
    nb_pw[q == 0, ] <- nb_sb[q == 0, ] <- 0
    d_pw[q == 1, ] <- d_sb[q == 1, ] <- 0
  } else if (interval == "wald") {
    lim <- function(est, se, zz) cbind(est - zz * se, est + zz * se)
    nb_pw <- lim(est_nb, se_nb, z); nb_sb <- lim(est_nb, se_nb, crit[["nb"]])
    d_pw <- lim(est_d, se_d, z); d_sb <- lim(est_d, se_d, crit[["delta_all"]])
  } else {
    # event rate among flagged (m1) and unflagged (m0) patients, then map back
    m1 <- p + est_nb * (1 - p) / q; s1 <- se_nb * (1 - p) / q
    m0 <- p - est_d * (1 - p) / (1 - q); s0 <- se_d * (1 - p) / (1 - q)
    to_nb <- function(m) q * (m - p) / (1 - p)
    to_d <- function(m) (1 - q) * (p - m) / (1 - p)
    nb_pw <- to_nb(wilson(m1, s1, n1, z)); nb_sb <- to_nb(wilson(m1, s1, n1, crit[["nb"]]))
    d_pw <- to_d(wilson(m0, s0, n0, z))[, 2:1, drop = FALSE]
    d_sb <- to_d(wilson(m0, s0, n0, crit[["delta_all"]]))[, 2:1, drop = FALSE]
    nb_pw[q == 0, ] <- nb_sb[q == 0, ] <- 0
    d_pw[q == 1, ] <- d_sb[q == 1, ] <- 0
  }
  out <- data.frame(
    threshold = thresholds,
    nb = est_nb, nb_se = se_nb,
    nb_lo = nb_pw[, 1], nb_hi = nb_pw[, 2], nb_slo = nb_sb[, 1], nb_shi = nb_sb[, 2],
    delta_all = est_d, delta_se = se_d,
    delta_lo = d_pw[, 1], delta_hi = d_pw[, 2], delta_slo = d_sb[, 1], delta_shi = d_sb[, 2]
  )
  attr(out, "crit") <- crit
  out
}

# Cross-fitted logistic working model for P(truth = 1 | risk, y), fitted with
# design weights on reviewed patients from the other folds. Returns NULL when
# there are too few reviewed patients or events, so stratum means are used.
logistic_mu <- function(risk, y, r, tt, prob, folds) {
  if (sum(r) < 40 || sum(tt[r]) < 10 || sum(1 - tt[r]) < 10) return(NULL)
  lr <- stats::qlogis(pmin(pmax(risk, 1e-6), 1 - 1e-6))
  yy <- as.numeric(y)
  fold <- sample(rep(seq_len(folds), length.out = length(risk)))
  mu <- numeric(length(risk))
  for (f in seq_len(folds)) {
    tr <- r & fold != f
    d <- data.frame(t = tt[tr], lr = lr[tr], y = yy[tr], w = 1 / prob[tr])
    form <- if (length(unique(d$y)) > 1) t ~ lr * y else t ~ lr
    fit <- tryCatch(suppressWarnings(stats::glm(form, family = stats::quasibinomial(),
                                                data = d, weights = d$w)), error = function(e) NULL)
    if (is.null(fit) || !fit$converged) return(NULL)
    mu[fold == f] <- stats::predict(fit, data.frame(lr = lr[fold == f], y = yy[fold == f]), type = "response")
  }
  mu
}

# Wilson score interval for a weighted proportion m with standard error s,
# using the effective sample size m(1 - m)/s^2 (fallback: n reviewed). With
# cap = TRUE the effective sample size is not allowed to exceed the number of
# reviewed charts, so a part with few or no observed events is not treated as
# known.
wilson <- function(m, s, nrev, z, cap = FALSE) {
  m <- pmin(pmax(m, 0), 1)
  ne <- ifelse(m > 0 & m < 1 & s > 0, m * (1 - m) / s^2, pmax(nrev, 1))
  if (cap) ne <- pmax(pmin(ne, nrev), 1)
  den <- 1 + z^2 / ne
  ctr <- (m + z^2 / (2 * ne)) / den
  half <- z * sqrt(m * (1 - m) / ne + z^2 / (4 * ne^2)) / den
  cbind(pmax(ctr - half, 0), pmin(ctr + half, 1))
}

sup_crit <- function(phi, se, level, B) {
  n <- nrow(phi)
  chunk <- max(1, min(100, floor(1e7 / n)))
  keep <- se > 0
  if (!any(keep)) return(stats::qnorm(1 - (1 - level) / 2))
  phi <- phi[, keep, drop = FALSE]
  sc <- se[keep] * n
  draws <- numeric(0)
  while (length(draws) < B) {
    m <- min(chunk, B - length(draws))
    xi <- matrix(stats::rnorm(m * n), m, n)
    zs <- abs(sweep(xi %*% phi, 2, sc, "/"))
    draws <- c(draws, apply(zs, 1, max))
  }
  as.numeric(stats::quantile(draws, level))
}

#' Optimal chart review probabilities for decision curves
#'
#' Returns review probabilities that minimise the chart review contribution
#' to the variance of the estimated decision curve, averaged uniformly over a
#' range of thresholds, for a fixed expected number of reviews. The optimum
#' is proportional to `sqrt(v * h)`, where `v` is the variance of the true
#' outcome given risk and recorded outcome, and `h` weights each patient by
#' how many thresholds in the range use their outcome. For the comparison
#' with treat all, only unflagged patients enter `h`. Because `v` is only
#' guessed (for example from a small pilot review), a share of the budget is
#' spread evenly (`defensive`).
#'
#' @param risk Predicted risks.
#' @param y Recorded outcome (used only for checks; it enters through `v`).
#' @param n_review Expected number of charts to review.
#' @param range Threshold range c(lower, upper) the decision curve must cover.
#' @param comparator `"all"` (gain over treat all), `"none"` (net benefit
#'   itself) or `"both"`.
#' @param v Guess of `Var(truth | risk, y)` for each patient, for example
#'   `m * (1 - m)` with `m` from a pilot review.
#' @param floor Smallest probability allowed, so that every patient can be
#'   selected.
#' @param defensive Share of the review budget spread evenly over all
#'   patients (0 to 1). Mixing in simple random sampling protects against a
#'   poor guess of `v`: with review fraction `k = n_review / n`, the chart
#'   review part of the variance is at most
#'   `(1 - defensive * k) / (defensive * (1 - k))` times that of simple random
#'   sampling (about `1/defensive` when `k` is small).
#' @return A numeric vector of review probabilities summing to `n_review`.
#'   If the budget covers every patient with `v * h > 0`, those patients get
#'   probability 1 and the remainder is spread evenly over the others.
#' @examples
#' set.seed(5)
#' risk <- runif(2000, 0, 0.5)
#' y <- rbinom(2000, 1, risk * 0.7)
#' m <- ifelse(y == 1, 0.95, risk * 0.3 / (1 - risk * 0.7))
#' p <- twophase_design(risk, y, 200, range = c(0.05, 0.3), v = m * (1 - m))
#' sum(p)
#' @export
twophase_design <- function(risk, y, n_review, range = c(0.05, 0.3),
                            comparator = c("all", "none", "both"), v,
                            floor = 0.01, defensive = 0.3) {
  comparator <- match.arg(comparator)
  n <- length(risk)
  if (!is.numeric(risk) || anyNA(risk) || any(risk < 0 | risk > 1) || length(y) != n)
    stop("`risk` must be in [0, 1] and `y` the same length.", call. = FALSE)
  if (!is.numeric(floor) || floor < 0 || floor >= 1) stop("`floor` must lie in [0, 1).", call. = FALSE)
  if (!is.numeric(defensive) || defensive < 0 || defensive > 1)
    stop("`defensive` must lie in [0, 1].", call. = FALSE)
  if (length(v) != n || any(v < 0, na.rm = TRUE) || anyNA(v))
    stop("`v` must be a non-negative vector, one value per patient.", call. = FALSE)
  if (n_review <= floor * n || n_review > n)
    stop("`n_review` must exceed floor * n and be at most n.", call. = FALSE)
  if (!(0 < range[1] && range[1] < range[2] && range[2] < 1)) stop("`range` must satisfy 0 < lower < upper < 1.", call. = FALSE)
  h <- threshold_weight(risk, range, comparator)
  a <- sqrt(v * h)
  if (all(a == 0)) return(rep(n_review / n, n))
  pos <- a > 0
  if (n_review >= sum(pos) + floor * sum(!pos)) {
    # budget covers every informative patient; spread the rest evenly
    opt <- ifelse(pos, 1, (n_review - sum(pos)) / sum(!pos))
  } else {
    f <- function(lambda) sum(pmin(1, pmax(floor, lambda * a))) - n_review
    lambda <- stats::uniroot(f, c(0, 1 / min(a[pos])), tol = 1e-12)$root
    opt <- pmin(1, pmax(floor, lambda * a))
  }
  (1 - defensive) * opt + defensive * n_review / n
}

#' Plot a decision curve
#'
#' @param x Output of [nb_curve()], [nb_adjust()] or [twophase_dca()].
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' risk <- runif(1000); y <- rbinom(1000, 1, risk)
#' plot_dca(nb_curve(risk, y, seq(0.05, 0.6, 0.05)))
#' @export
plot_dca <- function(x, ...) {
  x0 <- x
  x <- x[order(x$threshold), , drop = FALSE]
  nb_all <- if ("nb_all" %in% names(x)) x$nb_all else x$nb - x$delta_all
  rng <- range(c(x$nb, nb_all, 0, x$nb_slo, x$nb_shi), na.rm = TRUE)
  graphics::plot(x$threshold, x$nb, type = "l", lwd = 2, ylim = rng,
                 xlab = "Threshold probability", ylab = "Net benefit", ...)
  graphics::lines(x$threshold, nb_all, lty = 2)
  graphics::abline(h = 0, lty = 3)
  if (!is.null(x$nb_slo))
    graphics::polygon(c(x$threshold, rev(x$threshold)), c(x$nb_slo, rev(x$nb_shi)),
                      col = grDevices::adjustcolor("grey", 0.4), border = NA)
  invisible(x0)
}
