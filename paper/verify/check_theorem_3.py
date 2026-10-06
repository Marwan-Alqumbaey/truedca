"""Theorem 3: closed-form bounds versus a brute-force search over the (Se, Sp) box, and the tipping points.
Run from the paper folder: python verify/check_theorem_3.py"""
import numpy as np
rng = np.random.default_rng(2)

def closed(y, a, b, c, d):
    if not (1-d <= y <= b): return None
    return max(0, (y-1+c)/(b+c-1)), min(1, (y-1+d)/(a+d-1))

def brute(y, a, b, c, d, m=801):
    Se, Sp = np.meshgrid(np.linspace(a, b, m), np.linspace(c, d, m))
    pi = (y-1+Sp)/(Se+Sp-1)
    ok = (pi >= -1e-12) & (pi <= 1+1e-12)
    if not ok.any(): return None
    v = np.sort(pi[ok]); gap = np.max(np.diff(v)) if v.size > 1 else 0
    return v[0], v[-1], gap

worst, mism, n_empty, n_clipL, n_clipU, tested = 0, 0, 0, 0, 0, 0
for _ in range(3000):
    a, b = np.sort(rng.uniform(0.4, 1, 2)); c, d = np.sort(rng.uniform(0.4, 1, 2))
    if a+c <= 1: continue
    y = rng.uniform(0, 1); tested += 1
    cf, bf = closed(y, a, b, c, d), brute(y, a, b, c, d)
    if (cf is None) != (bf is None): mism += 1; continue
    if cf is None: n_empty += 1; continue
    n_clipL += (y-1+c) < 0; n_clipU += (y-1+d)/(a+d-1) > 1
    worst = max(worst, abs(cf[0]-bf[0]), abs(cf[1]-bf[1]))
print(f"tested {tested}: emptiness mismatches={mism}, empty={n_empty}, clipped L={n_clipL}, clipped U={n_clipU}, "
      f"max |closed-brute| = {worst:.2e} (grid resolution)")

# tipping points: D>0 iff pi0<p iff Se0>Se_tip (Sp0 fixed); NB>0 iff Sp1>Sp_tip (Se1=1)
bad = 0
for _ in range(20000):
    p, y0, Sp0, Se0 = rng.uniform(0.02, 0.98, 4)
    if Se0+Sp0 <= 1 or not (1-Sp0 <= y0 <= Se0): continue
    pi0 = (y0-1+Sp0)/(Se0+Sp0-1); tip = 1-Sp0+(y0-1+Sp0)/p
    bad += (pi0 < p) != (Se0 > tip)
    y1, Sp1 = rng.uniform(0, 1, 2)
    if 1-Sp1 <= y1:
        pi1 = (y1-1+Sp1)/Sp1; bad += (pi1 > p) != (Sp1 > (1-y1)/(1-p))
print("tipping-point equivalence failures:", bad)
assert mism == 0 and worst < 1e-2 and bad == 0
print("all checks passed")
