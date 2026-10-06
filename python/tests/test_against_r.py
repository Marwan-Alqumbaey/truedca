"""The Python port must reproduce the R package (fixtures from tests/fixtures/make_fixtures.R)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import truedca as td
from truedca._core import _design, _weighted_logit

FX = Path(__file__).parent / "fixtures"
D = pd.read_csv(FX / "data.csv")
TH = [0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
RISK, Y, TRUTH, PROB, ST = D.risk.values, D.y.values, D.truth.values, D.prob.values, D.strata.values


def same(py, r, cols=None, tol=1e-8):
    cols = cols or [c for c in r.columns if c in py.columns]
    for c in cols:
        a = py[c].astype(float).to_numpy()
        b = r[c].astype(float).to_numpy()
        np.testing.assert_allclose(a, b, rtol=tol, atol=tol, equal_nan=True, err_msg=c)


def test_nb_curve():
    same(td.nb_curve(RISK, Y, TH), pd.read_csv(FX / "nb_curve.csv"))


def test_nb_adjust():
    py = td.nb_adjust(RISK, Y, TH, se0=0.6, sp0=0.99, se1=0.85, sp1=0.98)
    same(py, pd.read_csv(FX / "nb_adjust.csv"))


def test_nb_bounds():
    py = td.nb_bounds(RISK, Y, TH, se0=(0.4, 1), sp0=(0.97, 1), se1=(0.7, 1), sp1=(0.98, 1), tip_sp0=0.99, tip_sp1=0.995)
    r = pd.read_csv(FX / "nb_bounds.csv")
    same(py, r, [c for c in r.columns if c not in ("conflict", "beats_none", "beats_all")])
    for c in ("conflict", "beats_none", "beats_all"):
        assert (py[c].to_numpy() == r[c].to_numpy()).all(), c


def test_strata_match():
    for col, brk in (("strata", [0.05, 0.1, 0.2]), ("strata5", 5)):
        py = pd.factorize(pd.Series(td.make_strata(RISK, Y, brk)).astype(str))[0]
        r = D[col].values
        assert pd.crosstab(py, r).gt(0).sum(axis=1).eq(1).all()
        assert pd.crosstab(py, r).gt(0).sum(axis=0).eq(1).all()


def test_twophase_intervals():
    py = td.twophase_dca(RISK, Y, TRUTH, PROB, TH, strata=ST, B=10, working="strata", rng=1)
    same(py, pd.read_csv(FX / "twophase_strata.csv"))
    for iv in ("wald", "wilson"):
        pw = td.twophase_dca(RISK, Y, TRUTH, PROB, TH, strata=ST, B=10, working="strata", interval=iv, rng=1)
        same(pw, pd.read_csv(FX / f"twophase_{iv}.csv"))


def test_design_and_gain():
    r = pd.read_csv(FX / "design.csv")
    np.testing.assert_allclose(td.twophase_design(RISK, Y, 300, (0.02, 0.3), "all", v=r.v, floor=0.01), r["all"], rtol=1e-6)
    p = td.twophase_design(RISK, Y, 300, (0.02, 0.3), "both", v=r.v, floor=0.01, defensive=0)
    np.testing.assert_allclose(p, r["both"], rtol=1e-6)
    assert abs(p.sum() - 300) < 1e-6
    g = td.design_gain(RISK, Y, 300, (0.02, 0.3), "both", v=r.v)
    rg = pd.read_csv(FX / "design_gain.csv").iloc[0]
    for k in ("gain", "bound", "fraction"):
        assert g[k] == pytest.approx(rg[k], rel=1e-6)


def test_weighted_logit_matches_glm():
    seen = ~np.isnan(TRUTH)
    lr = np.log(RISK / (1 - RISK))
    beta, ok = _weighted_logit(_design(lr[seen], Y[seen], True), TRUTH[seen], 1 / PROB[seen])
    assert ok
    np.testing.assert_allclose(beta, pd.read_csv(FX / "glm_coef.csv").to_numpy().ravel(), rtol=1e-6)


def test_pilot_v():
    r = pd.read_csv(FX / "pilot_v.csv").v.to_numpy()
    tr = pd.read_csv(FX / "pilot_truth.csv", skip_blank_lines=False).truth.to_numpy()
    np.testing.assert_allclose(td.pilot_v(RISK, Y, tr), r, rtol=1e-6)


def test_pilot_v_rank_deficient():
    r = pd.read_csv(FX / "pilot_v_edge.csv", skip_blank_lines=False)
    np.testing.assert_allclose(td.pilot_v(RISK, Y, r.truth.to_numpy()), r.v.to_numpy(), rtol=1e-6, atol=1e-9)
