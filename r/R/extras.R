#' How much an optimal chart review design can gain
#'
#' Compares the chart review part of the variance of the estimated decision
#' curve under simple random sampling with that under the optimal design of
#' [twophase_design()] for the same number of reviews. Without the caps at the
#' floor and at 1, the gain is \eqn{(1 - f) E[a] / \{(E\sqrt{a})^2 - f E[a]\}},
#' with \eqn{a = v h} and \eqn{f} the sampling fraction; for small \eqn{f} it is
#' at least \eqn{1 + CV^2(\sqrt{a})}, one plus the squared coefficient of
#' variation of \eqn{\sqrt{vh}}, with equality as \eqn{f \to 0}. A value close
#' to 1 means simple random sampling is as good. The gain refers to the optimal
#' design without the defensive share; with a guessed `v` it is an estimate.
#'
#' @inheritParams twophase_design
#' @return A named vector: `gain` (exact, with the caps), `bound` (the small
#'   sampling fraction limit \eqn{E[a]/(E\sqrt{a})^2}, a lower bound for the
#'   uncapped gain) and `fraction`.
#'   Both refer to the chart review part of the variance only; the gain in
#'   total variance is smaller.
#' @examples
#' set.seed(6)
#' risk <- runif(5000, 0, 0.5)
#' y <- rbinom(5000, 1, risk * 0.6)
#' m <- ifelse(y == 1, 0.95, risk * 0.4 / (1 - risk * 0.6))
#' design_gain(risk, y, 400, range = c(0.05, 0.3), comparator = "all", v = m * (1 - m))
#' @export
design_gain <- function(risk, y, n_review, range = c(0.05, 0.3),
                        comparator = c("all", "none", "both"), v, floor = 0.001) {
  comparator <- match.arg(comparator)
  n <- length(risk); f <- n_review / n
  h <- threshold_weight(risk, range, comparator)
  a <- v * h
  srs <- mean(a) * (1 / f - 1)
  opt <- twophase_design(risk, y, n_review, range, comparator, v = v, floor = floor, defensive = 0)
  c(gain = srs / mean(a * (1 / opt - 1)), bound = mean(a) / mean(sqrt(a))^2, fraction = f)
}

threshold_weight <- function(risk, range, comparator) {
  pl <- range[1]; pu <- range[2]
  h_all <- ifelse(risk < pu, (1 / (1 - pu) - 1 / (1 - pmax(risk, pl))) / (pu - pl), 0)
  h_none <- ifelse(risk >= pl, (1 / (1 - pmin(risk, pu)) - 1 / (1 - pl)) / (pu - pl), 0)
  switch(comparator, all = h_all, none = h_none, both = h_all + h_none)
}

#' Probabilistic bias analysis for a decision curve
#'
#' When no chart review is possible, uncertainty about the accuracy of the
#' recorded outcome can be expressed as Beta distributions for sensitivity
#' and specificity. Each draw takes error rates from these distributions and
#' recorded event rates from their Jeffreys posteriors, corrects the event
#' rates among flagged and unflagged patients, and recomputes net benefit.
#' Draws that give event rates outside \[0, 1\] are dropped, because those
#' error rates are incompatible with the data.
#'
#' @inheritParams nb_curve
#' @param se0,sp0 Beta shape parameters `c(a, b)` for sensitivity and
#'   specificity among unflagged patients.
#' @param se1,sp1 The same among flagged patients. By default the same draw
#'   is used in both groups (errors unrelated to risk).
#' @param draws Number of Monte Carlo draws.
#' @param level Width of the reported intervals.
#' @return A data frame with the median and interval for net benefit and for
#'   the gain over treat all, the probability that the model beats treat none
#'   (`p_beats_none`) and treat all (`p_beats_all`), and the share of draws
#'   kept.
#' @examples
#' set.seed(7)
#' risk <- runif(3000, 0, 0.5)
#' y <- rbinom(3000, 1, risk * 0.7)
#' nb_bias_analysis(risk, y, c(0.1, 0.2), se0 = c(70, 30), sp0 = c(995, 5), draws = 500)
#' @export
nb_bias_analysis <- function(risk, y, thresholds = seq(0.01, 0.5, by = 0.01),
                             se0 = c(80, 20), sp0 = c(990, 10), se1 = NULL, sp1 = NULL,
                             draws = 2000, level = 0.95) {
  check_inputs(risk, y, thresholds)
  y <- as.numeric(y); n <- length(y)
  shared <- is.null(se1) && is.null(sp1)
  if (is.null(se1)) se1 <- se0
  if (is.null(sp1)) sp1 <- sp0
  s0 <- stats::rbeta(draws, se0[1], se0[2]); c0 <- stats::rbeta(draws, sp0[1], sp0[2])
  if (shared) { s1 <- s0; c1 <- c0 } else {
    s1 <- stats::rbeta(draws, se1[1], se1[2]); c1 <- stats::rbeta(draws, sp1[1], sp1[2])
  }
  al <- (1 - level) / 2
  out <- lapply(thresholds, function(p) {
    flag <- risk >= p; n1 <- sum(flag); n0 <- n - n1; x1 <- sum(y[flag]); x0 <- sum(y[!flag])
    y1 <- stats::rbeta(draws, x1 + 0.5, n1 - x1 + 0.5); y0 <- stats::rbeta(draws, x0 + 0.5, n0 - x0 + 0.5)
    p1 <- (y1 - 1 + c1) / (s1 + c1 - 1); p0 <- (y0 - 1 + c0) / (s0 + c0 - 1)
    ok <- p1 >= 0 & p1 <= 1 & p0 >= 0 & p0 <= 1 & (s1 + c1 > 1) & (s0 + c0 > 1)
    q <- n1 / n
    nb <- (q * (p1 - p) / (1 - p))[ok]; d <- ((1 - q) * (p - p0) / (1 - p))[ok]
    qq <- function(x) if (length(x)) stats::quantile(x, c(0.5, al, 1 - al), names = FALSE) else rep(NA_real_, 3)
    a <- qq(nb); b <- qq(d)
    data.frame(threshold = p, nb = a[1], nb_lo = a[2], nb_hi = a[3],
               delta_all = b[1], delta_lo = b[2], delta_hi = b[3],
               p_beats_none = if (length(nb)) mean(nb > 0) else NA_real_,
               p_beats_all = if (length(d)) mean(d > 0) else NA_real_, kept = mean(ok))
  })
  do.call(rbind, out)
}

#' Guess the outcome uncertainty from a pilot chart review
#'
#' Fits the logistic working model of [twophase_dca()] (logit risk, recorded
#' outcome and their interaction) to the reviewed charts of a pilot sample and
#' returns \eqn{v = \mu(1 - \mu)} for every patient, the input `v` of
#' [twophase_design()] and [design_gain()]. Half an event and half a non-event
#' are added at the median risk of each recorded-outcome group, which keeps the
#' fit finite when, for example, every reviewed record is confirmed. A model
#' with four parameters learns more from a pilot of about 100 charts than
#' separate stratum means do.
#'
#' @inheritParams twophase_dca
#' @param truth True outcome for pilot patients, `NA` for the rest.
#' @param prob Pilot selection probabilities (default: equal).
#' @return Numeric vector of guessed variances, one per patient.
#' @examples
#' set.seed(8)
#' n <- 5000; risk <- runif(n, 0, 0.5); truth <- rbinom(n, 1, risk)
#' y <- truth * rbinom(n, 1, 0.7)
#' pilot <- ifelse(seq_len(n) %in% sample.int(n, 100), truth, NA)
#' summary(pilot_v(risk, y, pilot))
#' @export
pilot_v <- function(risk, y, truth, prob = rep(1, length(risk))) {
  if (length(y) != length(risk) || length(truth) != length(risk) || length(prob) != length(risk))
    stop("`risk`, `y`, `truth` and `prob` must have the same length.", call. = FALSE)
  ok <- !is.na(truth)
  if (sum(ok) < 10) stop("Need at least 10 reviewed pilot charts.", call. = FALSE)
  lr <- stats::qlogis(pmin(pmax(risk, 1e-6), 1 - 1e-6)); yy <- as.numeric(y)
  d <- data.frame(t = as.numeric(truth[ok]), lr = lr[ok], y = yy[ok], w = 1 / prob[ok])
  for (j in unique(yy))
    d <- rbind(d, data.frame(t = c(0, 1), lr = stats::median(lr[yy == j]), y = j, w = 0.5))
  form <- if (length(unique(yy)) > 1) t ~ lr * y else t ~ lr
  fit <- suppressWarnings(stats::glm(form, family = stats::quasibinomial(), data = d, weights = d$w))
  if (anyNA(stats::coef(fit)))   # interaction not estimable (e.g. too few reviewed records): drop it
    fit <- suppressWarnings(stats::glm(t ~ lr + y, family = stats::quasibinomial(), data = d, weights = d$w))
  m <- stats::predict(fit, data.frame(lr = lr, y = yy), type = "response")
  as.numeric(m * (1 - m))
}
