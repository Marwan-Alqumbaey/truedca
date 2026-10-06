"""Properties from the theory, checked by simulation."""
import numpy as np
import pytest

import truedca as td


def test_window_and_warp():
    w = td.dca_window(0.9, 0.97)
    assert w["floor"] == pytest.approx(0.03) and w["ceiling"] == 0.9
    assert td.warp_threshold(w["fixed_point"], 0.9, 0.97) == pytest.approx(w["fixed_point"])
    with pytest.raises(ValueError):
        td.dca_window(0.4, 0.5)


def test_perfect_recording_gives_identity():
    rng = np.random.default_rng(1)
    risk = rng.uniform(size=2000)
    y = rng.binomial(1, risk)
    a, b = td.nb_curve(risk, y, [0.1, 0.3]), td.nb_adjust(risk, y, [0.1, 0.3])
    np.testing.assert_allclose(a.nb, b.nb)


def test_uninformative_score_looks_useful_when_cases_are_missed():
    # Theorem 2: with Sp0 = 1, missed cases among unflagged patients inflate the gain over treat all
    rng = np.random.default_rng(2)
    n = 200000
    truth = rng.binomial(1, 0.2, n)
    risk = rng.uniform(size=n)
    y = truth * rng.binomial(1, np.where(risk >= 0.15, 0.95, 0.4))
    obs = td.nb_curve(risk, y, [0.15]).delta_all.iloc[0]
    tru = td.nb_curve(risk, truth, [0.15]).delta_all.iloc[0]
    assert tru < 0 < obs


def test_twophase_unbiased_and_covers():
    rng = np.random.default_rng(3)
    n, hits, ests = 5000, 0, []
    risk0 = rng.beta(1, 6, n)
    for i in range(60):
        truth = rng.binomial(1, risk0)
        y = truth * rng.binomial(1, 0.7, n)
        prob = np.where(y == 1, 0.3, 0.1)
        tt = np.where(rng.binomial(1, prob) == 1, truth, np.nan)
        f = td.twophase_dca(risk0, y, tt, prob, [0.1], B=50, rng=i)
        tr = td.nb_curve(risk0, truth, [0.1]).nb.iloc[0]
        ests.append(f.nb.iloc[0] - tr)
        hits += f.nb_lo.iloc[0] <= tr <= f.nb_hi.iloc[0]
    assert abs(np.mean(ests)) < 3 * np.std(ests) / np.sqrt(60) + 1e-3
    assert hits >= 50


def test_bias_analysis_shapes():
    rng = np.random.default_rng(4)
    risk = rng.uniform(0, 0.5, 3000)
    y = rng.binomial(1, risk * 0.7)
    out = td.nb_bias_analysis(risk, y, [0.1, 0.2], se0=(70, 30), sp0=(995, 5), draws=500, rng=1)
    assert out.shape[0] == 2 and (out.kept > 0.5).all()
    assert (out.nb_lo <= out.nb).all() and (out.nb <= out.nb_hi).all()


def test_bounds_ci_contain_bounds():
    rng = np.random.default_rng(5)
    risk = rng.uniform(size=5000)
    y = rng.binomial(1, risk * 0.7)
    b = td.nb_bounds(risk, y, [0.1, 0.2], se0=(0.6, 1), sp0=(0.98, 1))
    assert (b.nb_lower_ci <= b.nb_lower).all() and (b.nb_upper <= b.nb_upper_ci).all()
    assert (b.delta_lower_ci <= b.delta_lower).all() and (b.se_tip <= b.se_tip_upper).all()
