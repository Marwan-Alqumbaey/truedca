#' Net benefit of a prediction model on recorded outcomes
#'
#' Computes the usual decision curve from the outcomes as recorded. Net benefit
#' is \eqn{NB(p) = E[t_p (Y - p)]/(1 - p)} with \eqn{t_p = 1\{risk \ge p\}}.
#'
#' @param risk Numeric vector of predicted risks in \[0, 1\].
#' @param y Recorded binary outcome (0/1 or logical).
#' @param thresholds Threshold probabilities in (0, 1).
#' @return A data frame of class `truedca_curve` with one row per threshold:
#'   `threshold`, `nb` (model), `nb_all` (treat all), `delta_all`
#'   (model minus treat all), `q` (share flagged), `y1` and `y0` (recorded
#'   event rates among flagged and unflagged patients).
#' @examples
#' set.seed(1)
#' risk <- runif(1000)
#' y <- rbinom(1000, 1, risk)
#' nb_curve(risk, y, c(0.1, 0.2, 0.3))
#' @export
nb_curve <- function(risk, y, thresholds = seq(0.01, 0.5, by = 0.01)) {
  check_inputs(risk, y, thresholds)
  y <- as.numeric(y)
  out <- lapply(thresholds, function(p) {
    t <- risk >= p
    q <- mean(t)
    data.frame(
      threshold = p,
      nb = mean(t * (y - p)) / (1 - p),
      nb_all = (mean(y) - p) / (1 - p),
      delta_all = mean((1 - t) * (p - y)) / (1 - p),
      q = q,
      y1 = if (q > 0) mean(y[t]) else NA_real_,
      y0 = if (q < 1) mean(y[!t]) else NA_real_
    )
  })
  structure(do.call(rbind, out), class = c("truedca_curve", "data.frame"))
}

#' Thresholds at which a decision curve can be informative
#'
#' With outcome sensitivity `se` and specificity `sp` (equal among flagged and
#' unflagged patients), no model can show a gain over treat all for thresholds
#' at or below `1 - sp`, and none can show a gain over treat none at or above
#' `se`. Between these limits, recorded net benefit is a rescaled true net
#' benefit judged at a shifted threshold (see [warp_threshold()]).
#'
#' @param se,sp Outcome sensitivity and specificity, with `se + sp > 1`.
#' @return A named numeric vector: `floor`, `ceiling`, and `fixed_point`, the
#'   threshold at which the shift is zero. Below the fixed point the treat-all
#'   comparison errs towards "no gain" and the treat-none comparison towards
#'   "gain"; above it the directions swap.
#' @examples
#' dca_window(se = 0.9, sp = 0.97)
#' @export
dca_window <- function(se, sp) {
  check_rates(se, sp)
  fp <- if (se + sp < 2) (1 - sp) / ((1 - se) + (1 - sp)) else NA_real_
  c(floor = 1 - sp, ceiling = se, fixed_point = fp)
}

#' Threshold at which recorded outcomes judge a decision
#'
#' Returns \eqn{\delta(p) = (p - (1 - sp))/(se + sp - 1)}. With
#' non-differential outcome error, recorded net benefit at threshold `p`
#' equals \eqn{(se - p)/(1 - p)} times the true net benefit of the same
#' decision rule judged at threshold \eqn{\delta(p)} (Natarajan et al., 2018).
#'
#' @param p Threshold probabilities.
#' @inheritParams dca_window
#' @return Numeric vector of shifted thresholds (outside (0, 1) when `p` is
#'   outside the informative window).
#' @examples
#' warp_threshold(c(0.05, 0.1, 0.2), se = 0.9, sp = 0.97)
#' @export
warp_threshold <- function(p, se, sp) {
  check_rates(se, sp)
  (p - (1 - sp)) / (se + sp - 1)
}

check_inputs <- function(risk, y, thresholds) {
  if (!is.numeric(risk) || length(risk) == 0 || anyNA(risk) || any(risk < 0 | risk > 1))
    stop("`risk` must be numeric in [0, 1] without missing values.", call. = FALSE)
  if (is.factor(y) || length(y) != length(risk) || anyNA(y) ||
      !all(suppressWarnings(as.numeric(y)) %in% c(0, 1)))
    stop("`y` must be 0/1 or logical (not a factor), the same length as `risk`, without missing values.",
         call. = FALSE)
  if (!is.numeric(thresholds) || length(thresholds) == 0 || anyNA(thresholds) ||
      any(thresholds <= 0 | thresholds >= 1))
    stop("`thresholds` must lie strictly between 0 and 1.", call. = FALSE)
  invisible(TRUE)
}

check_rates <- function(se, sp) {
  if (!is.numeric(se) || !is.numeric(sp) || anyNA(c(se, sp)))
    stop("Sensitivity and specificity must be numeric without missing values.", call. = FALSE)
  if (any(se < 0 | se > 1 | sp < 0 | sp > 1))
    stop("Sensitivity and specificity must lie in [0, 1].", call. = FALSE)
  if (any(se + sp <= 1))
    stop("Requires se + sp > 1.", call. = FALSE)
  invisible(TRUE)
}
