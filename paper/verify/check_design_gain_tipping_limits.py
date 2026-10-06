"""Theorem 4(iv) design gain, the flagged-group tipping point, and coverage of the confidence limits
for the bounds (Section 4). Needs the Python package. Run from the paper folder:
python verify/check_design_gain_tipping_limits.py"""
import numpy as np
import truedca as td

rng = np.random.default_rng(5)

# Theorem 4(iv): closed form vs direct computation of the chart review variance ratio
for trial in range(5):
    n = 20000
    a = rng.gamma(0.5 + trial, 1.0, n) * (rng.uniform(size=n) < 0.7)
    f = 0.03
    lam = f * n / np.sqrt(a).sum()
    rho = lam * np.sqrt(a)
    assert rho.max() < 1, "pick f so that no cap binds"
    pos = a > 0
    direct = (a.mean() * (1 / f - 1)) / np.mean(np.where(pos, a * (1 / np.where(pos, rho, 1) - 1), 0))
    closed = (1 - f) * a.mean() / (np.sqrt(a).mean() ** 2 - f * a.mean())
    assert abs(direct / closed - 1) < 1e-10, (direct, closed)
print("Theorem 4(iv) closed form: OK")

# Flagged-group tipping point: NB > 0 iff Se1 < tau1, over random populations
for _ in range(2000):
    pi1, sp1, p = rng.uniform(0.01, 0.99), rng.uniform(0.8, 1), rng.uniform(0.02, 0.6)
    se1 = rng.uniform(1.0001 - sp1, 1)
    y1 = se1 * pi1 + (1 - sp1) * (1 - pi1)
    tau1 = 1 - sp1 + (y1 - 1 + sp1) / p
    assert (pi1 > p) == (se1 < tau1)
print("Flagged-group tipping point: OK")

# Confidence limits of the bounds cover the identified set (and the truth) in >= 95% of samples
n, reps, p = 3000, 1000, 0.15
hit = 0
for _ in range(reps):
    risk = rng.uniform(0, 0.5, n)
    truth = rng.binomial(1, risk)
    se = np.where(risk < p, 0.7, 0.9)
    y = np.where(truth == 1, rng.binomial(1, se), rng.binomial(1, 0.01, n))
    b = td.nb_bounds(risk, y, [p], se0=(0.6, 1), sp0=(0.98, 1), se1=(0.6, 1), sp1=(0.98, 1))
    # population quantities: pi0 = E(risk | risk < p) = p/2; true D
    pi0 = p / 2
    q = 1 - p / 0.5
    D = (1 - q) * (p - pi0) / (1 - p)
    hit += b.delta_lower_ci.iloc[0] <= D <= b.delta_upper_ci.iloc[0]
print(f"Coverage of the confidence limits for the gain over treat-all: {hit / reps:.3f} (target >= 0.95)")
assert hit / reps >= 0.94
