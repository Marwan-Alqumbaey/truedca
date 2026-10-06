"""Decision curve analysis for outcomes recorded with error (Python port of the R package truedca)."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import optimize, stats

__all__ = [
    "nb_curve", "dca_window", "warp_threshold", "nb_adjust", "nb_bounds",
    "make_strata", "twophase_dca", "twophase_design", "design_gain",
    "nb_bias_analysis", "pilot_v", "plot_dca",
]


# ---------------------------------------------------------------------------
# checks
def _as_arrays(risk, y, thresholds):
    risk = np.asarray(risk, dtype=float)
    y = np.asarray(y, dtype=float)
    thresholds = np.atleast_1d(np.asarray(thresholds, dtype=float))
    if risk.size == 0 or np.isnan(risk).any() or (risk < 0).any() or (risk > 1).any():
        raise ValueError("`risk` must be numeric in [0, 1] without missing values.")
    if y.shape != risk.shape or np.isnan(y).any() or not np.isin(y, (0, 1)).all():
        raise ValueError("`y` must be 0/1, the same length as `risk`, without missing values.")
    if thresholds.size == 0 or np.isnan(thresholds).any() or (thresholds <= 0).any() or (thresholds >= 1).any():
        raise ValueError("`thresholds` must lie strictly between 0 and 1.")
    return risk, y, thresholds


def _check_rates(se, sp):
    se, sp = np.asarray(se, dtype=float), np.asarray(sp, dtype=float)
    if np.isnan(se).any() or np.isnan(sp).any():
        raise ValueError("Sensitivity and specificity must be numeric without missing values.")
    if (se < 0).any() or (se > 1).any() or (sp < 0).any() or (sp > 1).any():
        raise ValueError("Sensitivity and specificity must lie in [0, 1].")
    if (se + sp <= 1).any():
        raise ValueError("Requires se + sp > 1.")


def _default_thresholds():
    return np.round(np.arange(1, 51) / 100, 10)


# ---------------------------------------------------------------------------
# net benefit on recorded outcomes
def nb_curve(risk, y, thresholds=None):
    """Decision curve on the recorded outcome.

    Returns a DataFrame with threshold, nb, nb_all, delta_all (model minus
    treat all), q (share flagged), y1 and y0 (recorded event rates among
    flagged and unflagged patients). Patients with risk >= threshold are flagged.
    """
    risk, y, thresholds = _as_arrays(risk, y, _default_thresholds() if thresholds is None else thresholds)
    ybar = y.mean()
    rows = []
    for p in thresholds:
        t = risk >= p
        q = t.mean()
        rows.append(dict(
            threshold=p,
            nb=np.mean(t * (y - p)) / (1 - p),
            nb_all=(ybar - p) / (1 - p),
            delta_all=np.mean((1 - t) * (p - y)) / (1 - p),
            q=q,
            y1=y[t].mean() if q > 0 else np.nan,
            y0=y[~t].mean() if q < 1 else np.nan,
        ))
    return pd.DataFrame(rows)


def dca_window(se, sp):
    """Floor (1 - sp), ceiling (se) and fixed point of the informative threshold window."""
    _check_rates(se, sp)
    fp = (1 - sp) / ((1 - se) + (1 - sp)) if se + sp < 2 else np.nan
    return {"floor": 1 - sp, "ceiling": se, "fixed_point": fp}


def warp_threshold(p, se, sp):
    """Threshold at which recorded outcomes judge a decision: (p - (1 - sp)) / (se + sp - 1)."""
    _check_rates(se, sp)
    return (np.asarray(p, dtype=float) - (1 - sp)) / (se + sp - 1)


# ---------------------------------------------------------------------------
# correction and bounds
def nb_adjust(risk, y, thresholds=None, se=1.0, sp=1.0, se0=None, sp0=None, se1=None, sp1=None):
    """Net benefit corrected for known outcome error, possibly different in flagged (1) and unflagged (0) patients."""
    se0 = se if se0 is None else se0
    sp0 = sp if sp0 is None else sp0
    se1 = se if se1 is None else se1
    sp1 = sp if sp1 is None else sp1
    _check_rates([se0, se1], [sp0, sp1])
    obs = nb_curve(risk, y, thresholds)
    p, q = obs.threshold.to_numpy(), obs.q.to_numpy()
    pi1 = (obs.y1.to_numpy() - 1 + sp1) / (se1 + sp1 - 1)
    pi0 = (obs.y0.to_numpy() - 1 + sp0) / (se0 + sp0 - 1)
    both = np.concatenate([pi1, pi0])
    both = both[~np.isnan(both)]
    if ((both < 0) | (both > 1)).any():
        warnings.warn("Corrected event rates fall outside [0, 1]; the assumed error rates conflict with the data.")
    with np.errstate(invalid="ignore"):
        nb = np.where(q > 0, q * (pi1 - p) / (1 - p), 0.0)
        delta = np.where(q < 1, (1 - q) * (p - pi0) / (1 - p), 0.0)
    return pd.DataFrame(dict(threshold=p, nb=nb, nb_all=nb - delta, delta_all=delta, q=q, pi1=pi1, pi0=pi0))


def _wilson(m, s, nrev, z, cap=False):
    m = np.clip(np.asarray(m, dtype=float), 0, 1)
    s = np.asarray(s, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ok = (m > 0) & (m < 1) & (s > 0)
        ne = np.where(ok, m * (1 - m) / np.where(ok, s, 1) ** 2, np.maximum(nrev, 1))
    if cap:  # never more information than the charts actually reviewed
        ne = np.maximum(np.minimum(ne, nrev), 1)
    den = 1 + z ** 2 / ne
    ctr = (m + z ** 2 / (2 * ne)) / den
    half = z * np.sqrt(m * (1 - m) / ne + z ** 2 / (4 * ne ** 2)) / den
    return np.column_stack([np.maximum(ctr - half, 0), np.minimum(ctr + half, 1)])


def nb_bounds(risk, y, thresholds=None, se0=(0.5, 1), sp0=(0.95, 1), se1=None, sp1=None,
              tip_sp0=1.0, tip_se1=1.0, tip_sp1=1.0, level=0.95):
    """Sharp bounds on true net benefit, tipping points, and (if level is set) confidence limits."""
    se1 = se0 if se1 is None else se1
    sp1 = sp0 if sp1 is None else sp1
    for b in (se0, sp0, se1, sp1):
        if len(b) != 2 or b[0] > b[1] or min(b) < 0 or max(b) > 1:
            raise ValueError("Each interval must be (lower, upper) within [0, 1].")
    if not (0 <= tip_sp0 <= 1 and 0 <= tip_se1 <= 1 and 0 <= tip_sp1 <= 1):
        raise ValueError("`tip_sp0`, `tip_se1` and `tip_sp1` must lie in [0, 1].")
    if se0[0] + sp0[0] <= 1 or se1[0] + sp1[0] <= 1:
        raise ValueError("Requires lower(se) + lower(sp) > 1.")
    risk_a = np.asarray(risk, dtype=float)
    obs = nb_curve(risk, y, thresholds)
    p, q = obs.threshold.to_numpy(), obs.q.to_numpy()
    y1, y0 = obs.y1.to_numpy(), obs.y0.to_numpy()

    def lo(yk, se, sp):
        return np.maximum(0, (yk - 1 + sp[0]) / (se[1] + sp[0] - 1))

    def hi(yk, se, sp):
        return np.minimum(1, (yk - 1 + sp[1]) / (se[0] + sp[1] - 1))

    def feasible(yk, se, sp):
        return (yk >= 1 - sp[1]) & (yk <= se[1])

    with np.errstate(invalid="ignore"):
        bad1 = ~(feasible(y1, se1, sp1) | np.isnan(y1))
        bad0 = ~(feasible(y0, se0, sp0) | np.isnan(y0))
        k1, k0 = q / (1 - p), (1 - q) / (1 - p)
        out = pd.DataFrame(dict(
            threshold=p,
            nb_lower=np.where(bad1, np.nan, np.where(q > 0, k1 * (lo(y1, se1, sp1) - p), 0.0)),
            nb_upper=np.where(bad1, np.nan, np.where(q > 0, k1 * (hi(y1, se1, sp1) - p), 0.0)),
            delta_lower=np.where(bad0, np.nan, np.where(q < 1, k0 * (p - hi(y0, se0, sp0)), 0.0)),
            delta_upper=np.where(bad0, np.nan, np.where(q < 1, k0 * (p - lo(y0, se0, sp0)), 0.0)),
            se_tip=1 - tip_sp0 + (y0 - 1 + tip_sp0) / p,
            sp_tip=(1 - y1 - p * (1 - tip_se1)) / (1 - p),
            se1_tip=1 - tip_sp1 + (y1 - 1 + tip_sp1) / p,
            conflict=bad1 | bad0,
        ))
        out["beats_none"] = out.nb_lower.notna() & (out.nb_lower > 0)
        out["beats_all"] = out.delta_lower.notna() & (out.delta_lower > 0)
        if level is not None:
            z = stats.norm.ppf(1 - (1 - level) / 2)
            n1 = np.array([np.sum(risk_a >= pp) for pp in p], dtype=float)
            n0 = risk_a.size - n1

            def ci(m, nn):
                m0 = np.nan_to_num(m, nan=0.0)
                return _wilson(m0, np.sqrt(np.maximum(m0 * (1 - m0), 0) / np.maximum(nn, 1)), nn, z)

            c1, c0 = ci(y1, n1), ci(y0, n0)
            out["nb_lower_ci"] = np.where(q > 0, k1 * (lo(c1[:, 0], se1, sp1) - p), 0.0)
            out["nb_upper_ci"] = np.where(q > 0, k1 * (hi(c1[:, 1], se1, sp1) - p), 0.0)
            out["delta_lower_ci"] = np.where(q < 1, k0 * (p - hi(c0[:, 1], se0, sp0)), 0.0)
            out["delta_upper_ci"] = np.where(q < 1, k0 * (p - lo(c0[:, 0], se0, sp0)), 0.0)
            out["se_tip_upper"] = 1 - tip_sp0 + (c0[:, 1] - 1 + tip_sp0) / p
            out["sp_tip_upper"] = (1 - c1[:, 0] - p * (1 - tip_se1)) / (1 - p)
            out["se1_tip_lower"] = 1 - tip_sp1 + (c1[:, 0] - 1 + tip_sp1) / p
            out["se1_tip_upper"] = 1 - tip_sp1 + (c1[:, 1] - 1 + tip_sp1) / p
    return out


# ---------------------------------------------------------------------------
# two-phase chart review
def make_strata(risk, y, breaks=5):
    """Strata crossing risk groups with the recorded outcome (labels 'group:y')."""
    risk = np.asarray(risk, dtype=float)
    y = np.asarray(y, dtype=float).astype(int)
    if np.ndim(breaks) == 0:
        cuts = np.unique(np.quantile(risk, np.linspace(0, 1, int(breaks) + 1)))
        if cuts.size < 2:
            g = np.zeros(risk.size, dtype=int)
        else:  # right-closed intervals, lowest included (as R's cut(include.lowest = TRUE))
            g = np.clip(np.searchsorted(cuts, risk, side="left"), 1, cuts.size - 1) - 1
    else:  # left-closed intervals [a, b) over (-inf, breaks, inf)
        g = np.searchsorted(np.sort(np.asarray(breaks, dtype=float)), risk, side="right")
    return pd.Categorical([f"{a}:{b}" for a, b in zip(g, y)])


def _weighted_logit(X, t, w, maxit=25, eps=1e-8):
    """IRLS for a weighted logistic regression, matching R's glm defaults. Returns (beta, converged)."""
    beta = np.zeros(X.shape[1])
    mu = (w * t + 0.5) / (w + 1)          # R's binomial initialisation
    eta = np.log(mu / (1 - mu))
    dev_old = np.inf
    for _ in range(maxit):
        mu = 1 / (1 + np.exp(-eta))
        var = np.maximum(mu * (1 - mu), 1e-12)
        z = eta + (t - mu) / var
        W = w * var
        beta = np.linalg.lstsq(X * np.sqrt(W)[:, None], z * np.sqrt(W), rcond=None)[0]
        eta = X @ beta
        mu = 1 / (1 + np.exp(-eta))
        with np.errstate(divide="ignore", invalid="ignore"):
            dev = 2 * np.sum(w * (np.where(t > 0, t * np.log(t / mu), 0) +
                                  np.where(t < 1, (1 - t) * np.log((1 - t) / (1 - mu)), 0)))
        if np.abs(dev - dev_old) / (np.abs(dev) + 0.1) < eps:
            return beta, True
        dev_old = dev
    return beta, False


def _design(lr, y, both):
    cols = [np.ones_like(lr), lr]
    if both:
        cols += [y, lr * y]
    return np.column_stack(cols)


def _logistic_mu(risk, y, r, tt, prob, folds, rng):
    if r.sum() < 40 or tt[r].sum() < 10 or (1 - tt[r]).sum() < 10:
        return None
    lr = stats.logistic.ppf(np.clip(risk, 1e-6, 1 - 1e-6))
    n = risk.size
    fold = rng.permutation(np.resize(np.arange(folds), n))
    mu = np.zeros(n)
    for f in range(folds):
        tr = r & (fold != f)
        both = np.unique(y[tr]).size > 1
        beta, ok = _weighted_logit(_design(lr[tr], y[tr], both), tt[tr], 1 / prob[tr])
        if not ok:
            return None
        te = fold == f
        mu[te] = 1 / (1 + np.exp(-_design(lr[te], y[te], both) @ beta))
    return mu


def _sup_crit(phi, se, level, B, rng):
    n = phi.shape[0]
    keep = se > 0
    if not keep.any():
        return stats.norm.ppf(1 - (1 - level) / 2)
    phi = phi[:, keep]
    sc = se[keep] * n
    chunk = max(1, min(100, int(1e7 // n)))
    draws = []
    while len(draws) < B:
        m = min(chunk, B - len(draws))
        xi = rng.standard_normal((m, n))
        draws.extend(np.abs((xi @ phi) / sc).max(axis=1))
    return float(np.quantile(draws, level))


def twophase_dca(risk, y, truth, prob, thresholds=None, strata=None, level=0.95, B=1000,
                 interval="split", working="logistic", folds=5, rng=None):
    """Decision curve from a two-phase chart review (AIPW estimator).

    ``truth`` holds the reviewed true outcome and NaN for patients not reviewed;
    ``prob`` is each patient's known review probability. Returns estimates,
    standard errors, pointwise (``*_lo``, ``*_hi``) and simultaneous
    (``*_slo``, ``*_shi``) limits for net benefit and the gain over treat all.
    """
    if interval not in ("split", "wilson", "wald") or working not in ("logistic", "strata"):
        raise ValueError("Unknown `interval` or `working`.")
    rng = np.random.default_rng(rng)
    risk, y, thresholds = _as_arrays(risk, y, _default_thresholds() if thresholds is None else thresholds)
    n = risk.size
    truth = np.asarray(truth, dtype=float)
    prob = np.asarray(prob, dtype=float)
    strata = make_strata(risk, y) if strata is None else strata
    if truth.size != n or prob.size != n or len(strata) != n:
        raise ValueError("`truth`, `prob` and `strata` must have the same length as `risk`.")
    if np.isnan(prob).any() or (prob <= 0).any() or (prob > 1).any():
        raise ValueError("`prob` must lie in (0, 1].")
    r = ~np.isnan(truth)
    if not np.isin(truth[r], (0, 1)).all():
        raise ValueError("`truth` must be 0/1, with NaN for patients not reviewed.")
    if (r[prob == 1] == False).any():  # noqa: E712
        raise ValueError("Patients selected with `prob` = 1 must have `truth` recorded.")
    if not 0 < level < 1:
        raise ValueError("`level` must be a single number in (0, 1).")
    codes, _ = pd.factorize(pd.Series(strata).astype(str))
    if (codes < 0).any():
        raise ValueError("`strata` must not contain missing values.")
    tt = np.where(r, truth, 0.0)
    H = codes.max() + 1
    wsum = np.bincount(codes, r / prob, H)
    num = np.bincount(codes, r * tt / prob, H)
    empty = (wsum == 0).any()
    with np.errstate(invalid="ignore", divide="ignore"):
        mu_h = np.where(wsum > 0, num / wsum, np.sum(r * tt / prob) / np.sum(r / prob))
    mu = mu_h[codes]
    nh = np.bincount(codes, r.astype(float), H)
    cf = np.where(nh > 1, nh / np.maximum(nh - 1, 1), 1.0)[codes]
    fit_mu = _logistic_mu(risk, y, r, tt, prob, folds, rng) if working == "logistic" else None
    if fit_mu is not None:
        mu = fit_mu
        cf = np.full(n, r.sum() / max(r.sum() - 4, 1))
    if empty and fit_mu is None:
        warnings.warn("Some strata have no reviewed patients; their working mean uses all reviewed patients.")
    tstar = mu + r * (tt - mu) / prob
    tvar = mu + r * (tt - mu) / prob * np.sqrt(cf)
    K = thresholds.size
    phi_nb, phi_d = np.zeros((n, K)), np.zeros((n, K))
    est_nb, est_d, q, n1, n0 = (np.zeros(K) for _ in range(5))
    sub = np.zeros((K, 2, 2, 4))  # threshold, group (flagged, unflagged), recorded y: share, rate, se, reviewed
    for k, p in enumerate(thresholds):
        flag = risk >= p
        for g, grp in enumerate((flag, ~flag)):
            for j in (0, 1):
                idx = grp & (y == j)
                if idx.any():
                    tv = tvar[idx]
                    sub[k, g, j] = (idx.sum() / grp.sum(), tstar[idx].mean(),
                                    np.sqrt(np.sum((tv - tv.mean()) ** 2)) / idx.sum(), r[idx].sum())
        est_nb[k] = np.mean(flag * (tstar - p)) / (1 - p)
        est_d[k] = np.mean(~flag * (p - tstar)) / (1 - p)
        av = flag * (tvar - p) / (1 - p)
        bv = ~flag * (p - tvar) / (1 - p)
        phi_nb[:, k], phi_d[:, k] = av - av.mean(), bv - bv.mean()
        q[k], n1[k], n0[k] = flag.mean(), np.sum(r & flag), np.sum(r & ~flag)
    se_nb = np.sqrt((phi_nb ** 2).sum(0)) / n
    se_d = np.sqrt((phi_d ** 2).sum(0)) / n
    z = stats.norm.ppf(1 - (1 - level) / 2)
    crit = {"nb": _sup_crit(phi_nb, se_nb, level, B, rng), "delta_all": _sup_crit(phi_d, se_d, level, B, rng)}
    p = thresholds
    if interval == "split":
        # MOVER over the two parts (recorded outcome 0 or 1) and the phase-one share of each part
        def lims(g, zz):
            cc0, cc1 = sub[:, g, 0, 0], sub[:, g, 1, 0]
            m = cc0 * sub[:, g, 0, 1] + cc1 * sub[:, g, 1, 1]
            ng = q * n if g == 0 else (1 - q) * n
            share = zz ** 2 * cc1 * (1 - cc1) * (sub[:, g, 1, 1] - sub[:, g, 0, 1]) ** 2 / np.maximum(ng, 1)
            lo, hi = share.copy(), share.copy()
            for j in (0, 1):
                cc, pp, ss, nn = (sub[:, g, j, i] for i in range(4))
                w = _wilson(pp, ss, nn, zz, cap=True)
                lo += (cc * (pp - w[:, 0])) ** 2
                hi += (cc * (w[:, 1] - pp)) ** 2
            return m, np.maximum(m - np.sqrt(lo), 0), np.minimum(m + np.sqrt(hi), 1)

        # map to net benefit, adding the phase-one error in the share flagged q
        def to_scale(lim, zz, flagged):
            m, lo, hi = lim
            k = (q if flagged else 1 - q) / (1 - p)
            sgn = 1 if flagged else -1
            est, a_, b_ = (k * sgn * (x - p) for x in (m, lo, hi))
            c0, c1 = np.minimum(a_, b_), np.maximum(a_, b_)
            vq = zz ** 2 * q * (1 - q) / n * ((m - p) / (1 - p)) ** 2
            return np.column_stack([est - np.sqrt((est - c0) ** 2 + vq), est + np.sqrt((c1 - est) ** 2 + vq)])

        with np.errstate(divide="ignore", invalid="ignore"):
            nb_pw, nb_sb = to_scale(lims(0, z), z, True), to_scale(lims(0, crit["nb"]), crit["nb"], True)
            d_pw = to_scale(lims(1, z), z, False)
            d_sb = to_scale(lims(1, crit["delta_all"]), crit["delta_all"], False)
        nb_pw[q == 0], nb_sb[q == 0] = 0, 0
        d_pw[q == 1], d_sb[q == 1] = 0, 0
    elif interval == "wald":
        def lim(e, s, zz):
            return np.column_stack([e - zz * s, e + zz * s])
        nb_pw, nb_sb = lim(est_nb, se_nb, z), lim(est_nb, se_nb, crit["nb"])
        d_pw, d_sb = lim(est_d, se_d, z), lim(est_d, se_d, crit["delta_all"])
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            m1, s1 = p + est_nb * (1 - p) / q, se_nb * (1 - p) / q
            m0, s0 = p - est_d * (1 - p) / (1 - q), se_d * (1 - p) / (1 - q)

            def to_nb(m):
                return q[:, None] * (m - p[:, None]) / (1 - p[:, None])

            def to_d(m):
                return ((1 - q[:, None]) * (p[:, None] - m) / (1 - p[:, None]))[:, ::-1]

            nb_pw, nb_sb = to_nb(_wilson(m1, s1, n1, z)), to_nb(_wilson(m1, s1, n1, crit["nb"]))
            d_pw, d_sb = to_d(_wilson(m0, s0, n0, z)), to_d(_wilson(m0, s0, n0, crit["delta_all"]))
        nb_pw[q == 0], nb_sb[q == 0] = 0, 0
        d_pw[q == 1], d_sb[q == 1] = 0, 0
    out = pd.DataFrame(dict(
        threshold=thresholds, nb=est_nb, nb_se=se_nb, nb_lo=nb_pw[:, 0], nb_hi=nb_pw[:, 1],
        nb_slo=nb_sb[:, 0], nb_shi=nb_sb[:, 1], delta_all=est_d, delta_se=se_d,
        delta_lo=d_pw[:, 0], delta_hi=d_pw[:, 1], delta_slo=d_sb[:, 0], delta_shi=d_sb[:, 1]))
    out.attrs["crit"] = crit
    return out


def _threshold_weight(risk, rng_, comparator):
    pl, pu = rng_
    h_all = np.where(risk < pu, (1 / (1 - pu) - 1 / (1 - np.maximum(risk, pl))) / (pu - pl), 0.0)
    h_none = np.where(risk >= pl, (1 / (1 - np.minimum(risk, pu)) - 1 / (1 - pl)) / (pu - pl), 0.0)
    return {"all": h_all, "none": h_none, "both": h_all + h_none}[comparator]


def twophase_design(risk, y, n_review, range=(0.05, 0.3), comparator="all", v=None,
                    floor=0.01, defensive=0.3):
    """Review probabilities minimising the threshold-averaged chart review variance (proportional to sqrt(v h))."""
    risk = np.asarray(risk, dtype=float)
    n = risk.size
    if comparator not in ("all", "none", "both"):
        raise ValueError("`comparator` must be 'all', 'none' or 'both'.")
    if np.isnan(risk).any() or (risk < 0).any() or (risk > 1).any() or len(y) != n:
        raise ValueError("`risk` must be in [0, 1] and `y` the same length.")
    if not 0 <= floor < 1:
        raise ValueError("`floor` must lie in [0, 1).")
    if not 0 <= defensive <= 1:
        raise ValueError("`defensive` must lie in [0, 1].")
    v = np.asarray(v, dtype=float)
    if v.size != n or np.isnan(v).any() or (v < 0).any():
        raise ValueError("`v` must be a non-negative vector, one value per patient.")
    if n_review <= floor * n or n_review > n:
        raise ValueError("`n_review` must exceed floor * n and be at most n.")
    if not 0 < range[0] < range[1] < 1:
        raise ValueError("`range` must satisfy 0 < lower < upper < 1.")
    a = np.sqrt(v * _threshold_weight(risk, range, comparator))
    if (a == 0).all():
        return np.full(n, n_review / n)
    pos = a > 0
    if n_review >= pos.sum() + floor * (~pos).sum():
        opt = np.where(pos, 1.0, (n_review - pos.sum()) / max((~pos).sum(), 1))
    else:
        def f(lam):
            return np.minimum(1, np.maximum(floor, lam * a)).sum() - n_review
        lam = optimize.brentq(f, 0, 1 / a[pos].min(), xtol=1e-14, rtol=1e-14, maxiter=500)
        opt = np.minimum(1, np.maximum(floor, lam * a))
    return (1 - defensive) * opt + defensive * n_review / n


def design_gain(risk, y, n_review, range=(0.05, 0.3), comparator="all", v=None, floor=0.001):
    """Variance ratio (SRS / optimal) for the chart review part of the variance, its small-fraction bound, and the fraction."""
    risk = np.asarray(risk, dtype=float)
    f = n_review / risk.size
    a = np.asarray(v, dtype=float) * _threshold_weight(risk, range, comparator)
    srs = a.mean() * (1 / f - 1)
    opt = twophase_design(risk, y, n_review, range, comparator, v=v, floor=floor, defensive=0)
    return {"gain": srs / np.mean(a * (1 / opt - 1)), "bound": a.mean() / np.mean(np.sqrt(a)) ** 2, "fraction": f}


# ---------------------------------------------------------------------------
# probabilistic bias analysis
def nb_bias_analysis(risk, y, thresholds=None, se0=(80, 20), sp0=(990, 10), se1=None, sp1=None,
                     draws=2000, level=0.95, rng=None):
    """Probabilistic bias analysis with Beta priors on sensitivity and specificity."""
    rng = np.random.default_rng(rng)
    risk, y, thresholds = _as_arrays(risk, y, _default_thresholds() if thresholds is None else thresholds)
    n = y.size
    shared = se1 is None and sp1 is None
    se1 = se0 if se1 is None else se1
    sp1 = sp0 if sp1 is None else sp1
    s0, c0 = rng.beta(*se0, draws), rng.beta(*sp0, draws)
    s1, c1 = (s0, c0) if shared else (rng.beta(*se1, draws), rng.beta(*sp1, draws))
    al = (1 - level) / 2
    rows = []
    for p in thresholds:
        flag = risk >= p
        n1 = flag.sum()
        x1, x0 = y[flag].sum(), y[~flag].sum()
        y1 = rng.beta(x1 + 0.5, n1 - x1 + 0.5, draws)
        y0 = rng.beta(x0 + 0.5, n - n1 - x0 + 0.5, draws)
        p1 = (y1 - 1 + c1) / (s1 + c1 - 1)
        p0 = (y0 - 1 + c0) / (s0 + c0 - 1)
        ok = (p1 >= 0) & (p1 <= 1) & (p0 >= 0) & (p0 <= 1) & (s1 + c1 > 1) & (s0 + c0 > 1)
        q = n1 / n
        nb = (q * (p1 - p) / (1 - p))[ok]
        d = ((1 - q) * (p - p0) / (1 - p))[ok]

        def qq(x):
            return np.quantile(x, [0.5, al, 1 - al]) if x.size else np.full(3, np.nan)

        a, b = qq(nb), qq(d)
        rows.append(dict(threshold=p, nb=a[0], nb_lo=a[1], nb_hi=a[2], delta_all=b[0], delta_lo=b[1],
                         delta_hi=b[2], p_beats_none=np.mean(nb > 0) if nb.size else np.nan,
                         p_beats_all=np.mean(d > 0) if d.size else np.nan, kept=ok.mean()))
    return pd.DataFrame(rows)


def pilot_v(risk, y, truth, prob=None):
    """Guess Var(T | risk, y) from a pilot review with the logistic working model (input ``v`` of the design)."""
    risk = np.asarray(risk, dtype=float)
    y = np.asarray(y, dtype=float)
    truth = np.asarray(truth, dtype=float)
    prob = np.ones_like(risk) if prob is None else np.asarray(prob, dtype=float)
    ok = ~np.isnan(truth)
    if ok.sum() < 10:
        raise ValueError("Need at least 10 reviewed pilot charts.")
    lr = stats.logistic.ppf(np.clip(risk, 1e-6, 1 - 1e-6))
    t, l, yy, w = [truth[ok]], [lr[ok]], [y[ok]], [1 / prob[ok]]
    for j in np.unique(y):
        t.append(np.array([0.0, 1.0])); l.append(np.repeat(np.median(lr[y == j]), 2))
        yy.append(np.array([j, j])); w.append(np.array([0.5, 0.5]))
    t, l, yy, w = (np.concatenate(a) for a in (t, l, yy, w))
    both = np.unique(y).size > 1
    X = _design(l, yy, both)
    if both and np.linalg.matrix_rank(X) < X.shape[1]:   # interaction not estimable: drop it (as R does)
        X, cols = X[:, :3], 3
    else:
        cols = X.shape[1]
    beta, _ = _weighted_logit(X, t, w)
    m = 1 / (1 + np.exp(-_design(lr, y, both)[:, :cols] @ beta))
    return m * (1 - m)


def plot_dca(x, ax=None):
    """Plot a decision curve (model, treat all, treat none, and a simultaneous band if present)."""
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    x = x.sort_values("threshold")
    nb_all = x["nb_all"] if "nb_all" in x else x["nb"] - x["delta_all"]
    if "nb_slo" in x:
        ax.fill_between(x.threshold, x.nb_slo, x.nb_shi, color="grey", alpha=0.4, lw=0)
    ax.plot(x.threshold, x.nb, lw=2, label="Model")
    ax.plot(x.threshold, nb_all, ls="--", label="Treat all")
    ax.axhline(0, ls=":", color="k", label="Treat none")
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.legend(frameon=False)
    return ax
