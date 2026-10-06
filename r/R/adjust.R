#' Net benefit corrected for known outcome error
#'
#' Corrects a decision curve when the sensitivity and specificity of the
#' recorded outcome are known. Rates may differ between patients the model
#' flags (`se1`, `sp1`) and those it does not (`se0`, `sp0`). The comparison
#' with treat all depends only on `se0` and `sp0`; the comparison with treat
#' none depends only on `se1` and `sp1`.
#'
#' @inheritParams nb_curve
#' @param se,sp Outcome sensitivity and specificity used for both groups
#'   unless overridden.
#' @param se0,sp0 Rates among unflagged patients (risk below threshold).
#' @param se1,sp1 Rates among flagged patients.
#' @return A `truedca_curve` data frame with corrected `nb`, `nb_all`,
#'   `delta_all`, and the corrected event rates `pi1`, `pi0`. A warning is
#'   given if a corrected rate falls outside \[0, 1\], which means the
#'   assumed error rates do not fit the data.
#' @examples
#' set.seed(2)
#' risk <- runif(5000)
#' truth <- rbinom(5000, 1, risk)
#' y <- ifelse(truth == 1, rbinom(5000, 1, 0.8), rbinom(5000, 1, 0.05))
#' nb_adjust(risk, y, c(0.1, 0.3), se = 0.8, sp = 0.95)
#' @export
nb_adjust <- function(risk, y, thresholds = seq(0.01, 0.5, by = 0.01),
                      se = 1, sp = 1, se0 = se, sp0 = sp, se1 = se, sp1 = sp) {
  check_rates(c(se0, se1), c(sp0, sp1))
  obs <- nb_curve(risk, y, thresholds)
  p <- obs$threshold
  pi1 <- (obs$y1 - 1 + sp1) / (se1 + sp1 - 1)
  pi0 <- (obs$y0 - 1 + sp0) / (se0 + sp0 - 1)
  if (any(c(pi1, pi0) < 0 | c(pi1, pi0) > 1, na.rm = TRUE))
    warning("Corrected event rates fall outside [0, 1]; the assumed error rates conflict with the data.",
            call. = FALSE)
  nb <- ifelse(obs$q > 0, obs$q * (pi1 - p) / (1 - p), 0)
  delta <- ifelse(obs$q < 1, (1 - obs$q) * (p - pi0) / (1 - p), 0)
  out <- data.frame(threshold = p, nb = nb, nb_all = nb - delta,
                    delta_all = delta, q = obs$q, pi1 = pi1, pi0 = pi0)
  structure(out, class = c("truedca_curve", "data.frame"))
}

#' Sharp bounds and tipping points when outcome error is partly known
#'
#' When outcome sensitivity and specificity are only known to lie in
#' intervals, true net benefit is not identified. This function returns the
#' sharp (pointwise) bounds on true net benefit and on its gain over treat
#' all, together with tipping points: the smallest outcome sensitivity among
#' unflagged patients at which the model still beats treat all, and the
#' smallest specificity among flagged patients at which it still beats treat
#' none.
#'
#' @inheritParams nb_curve
#' @param se0,sp0 Length-2 intervals for sensitivity and specificity among
#'   unflagged patients.
#' @param se1,sp1 Intervals among flagged patients (default: same as unflagged).
#' @param level Confidence level for confidence limits on the bounds and
#'   tipping points (`NULL` to skip). The limits plug Wilson limits for the
#'   recorded event rates into the bounds, so that the interval
#'   (`*_lower_ci`, `*_upper_ci`) covers the whole set of values consistent
#'   with the assumptions with probability at least `level`; `se_tip_upper`
#'   and `sp_tip_upper` are the matching upper confidence limits of the
#'   tipping points, and `se1_tip_lower`, `se1_tip_upper` the confidence
#'   limits of `se1_tip`.
#' @param tip_sp0 Specificity among unflagged patients assumed when computing
#'   the sensitivity tipping point `se_tip` (default 1, no false records).
#' @param tip_se1 Sensitivity among flagged patients assumed when computing
#'   the specificity tipping point `sp_tip` (default 1, no missed cases).
#' @param tip_sp1 Specificity among flagged patients assumed when computing
#'   the sensitivity tipping point `se1_tip` (default 1, no false records).
#' @return A data frame with one row per threshold: `nb_lower`, `nb_upper`
#'   (bounds on net benefit), `delta_lower`, `delta_upper` (bounds on the gain
#'   over treat all), `se_tip` (the model beats treat all if and only if
#'   sensitivity among unflagged patients exceeds it; a value above 1 means no
#'   sensitivity is enough), `sp_tip` (the same for specificity among flagged
#'   patients and treat none), `se1_tip` (the model beats treat none if and
#'   only if sensitivity among flagged patients is below it; useful when missed
#'   cases among flagged patients may hide a real benefit), `conflict` (the assumed intervals contradict the
#'   recorded event rates in at least one group; the bounds for that group are
#'   then `NA`), and logical `beats_none` and `beats_all` (true when the model
#'   wins for every admissible error rate). The bounds themselves are plug-in
#'   estimates; use the confidence limits for sampling error.
#' @examples
#' set.seed(3)
#' risk <- runif(5000)
#' y <- rbinom(5000, 1, risk * 0.7)
#' nb_bounds(risk, y, c(0.1, 0.2), se0 = c(0.6, 1), sp0 = c(0.98, 1))
#' @export
nb_bounds <- function(risk, y, thresholds = seq(0.01, 0.5, by = 0.01),
                      se0 = c(0.5, 1), sp0 = c(0.95, 1),
                      se1 = se0, sp1 = sp0, tip_sp0 = 1, tip_se1 = 1, tip_sp1 = 1, level = 0.95) {
  for (b in list(se0, sp0, se1, sp1))
    if (length(b) != 2 || b[1] > b[2] || any(b < 0 | b > 1))
      stop("Each interval must be c(lower, upper) within [0, 1].", call. = FALSE)
  if (!all(c(tip_sp0, tip_se1, tip_sp1) >= 0 & c(tip_sp0, tip_se1, tip_sp1) <= 1))
    stop("`tip_sp0`, `tip_se1` and `tip_sp1` must lie in [0, 1].", call. = FALSE)
  if (se0[1] + sp0[1] <= 1 || se1[1] + sp1[1] <= 1)
    stop("Requires lower(se) + lower(sp) > 1.", call. = FALSE)
  obs <- nb_curve(risk, y, thresholds)
  p <- obs$threshold
  lo <- function(yk, se, sp) pmax(0, (yk - 1 + sp[1]) / (se[2] + sp[1] - 1))
  hi <- function(yk, se, sp) pmin(1, (yk - 1 + sp[2]) / (se[1] + sp[2] - 1))
  feasible <- function(yk, se, sp) yk >= 1 - sp[2] & yk <= se[2]
  L1 <- lo(obs$y1, se1, sp1); U1 <- hi(obs$y1, se1, sp1)
  L0 <- lo(obs$y0, se0, sp0); U0 <- hi(obs$y0, se0, sp0)
  bad1 <- !(feasible(obs$y1, se1, sp1) | is.na(obs$y1))
  bad0 <- !(feasible(obs$y0, se0, sp0) | is.na(obs$y0))
  k1 <- obs$q / (1 - p); k0 <- (1 - obs$q) / (1 - p)
  na_if <- function(x, bad) ifelse(bad, NA_real_, x)
  out <- data.frame(
    threshold = p,
    nb_lower = na_if(ifelse(obs$q > 0, k1 * (L1 - p), 0), bad1),
    nb_upper = na_if(ifelse(obs$q > 0, k1 * (U1 - p), 0), bad1),
    delta_lower = na_if(ifelse(obs$q < 1, k0 * (p - U0), 0), bad0),
    delta_upper = na_if(ifelse(obs$q < 1, k0 * (p - L0), 0), bad0),
    se_tip = 1 - tip_sp0 + (obs$y0 - 1 + tip_sp0) / p,
    sp_tip = (1 - obs$y1 - p * (1 - tip_se1)) / (1 - p),
    se1_tip = 1 - tip_sp1 + (obs$y1 - 1 + tip_sp1) / p,
    conflict = bad1 | bad0
  )
  out$beats_none <- !is.na(out$nb_lower) & out$nb_lower > 0
  out$beats_all <- !is.na(out$delta_lower) & out$delta_lower > 0
  # confidence limits: plug Wilson limits of the recorded event rates into the
  # (monotone) bounds, which covers the whole identified set with prob >= level
  if (!is.null(level)) {
    z <- stats::qnorm(1 - (1 - level) / 2)
    n1 <- vapply(p, function(pp) sum(risk >= pp), numeric(1)); n0 <- length(risk) - n1
    ci <- function(m, nn) wilson(ifelse(is.na(m), 0, m), sqrt(pmax(m * (1 - m), 0) / pmax(nn, 1)), nn, z)
    c1 <- ci(obs$y1, n1); c0 <- ci(obs$y0, n0)
    out$nb_lower_ci <- ifelse(obs$q > 0, k1 * (lo(c1[, 1], se1, sp1) - p), 0)
    out$nb_upper_ci <- ifelse(obs$q > 0, k1 * (hi(c1[, 2], se1, sp1) - p), 0)
    out$delta_lower_ci <- ifelse(obs$q < 1, k0 * (p - hi(c0[, 2], se0, sp0)), 0)
    out$delta_upper_ci <- ifelse(obs$q < 1, k0 * (p - lo(c0[, 1], se0, sp0)), 0)
    out$se_tip_upper <- 1 - tip_sp0 + (c0[, 2] - 1 + tip_sp0) / p
    out$sp_tip_upper <- (1 - c1[, 1] - p * (1 - tip_se1)) / (1 - p)
    out$se1_tip_lower <- 1 - tip_sp1 + (c1[, 1] - 1 + tip_sp1) / p
    out$se1_tip_upper <- 1 - tip_sp1 + (c1[, 2] - 1 + tip_sp1) / p
  }
  out
}
