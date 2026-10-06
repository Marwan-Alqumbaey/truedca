"""Numerical checks of equation (1), Theorem 1, Corollary 1 (identity, floor, ceiling, direction) and Theorem 2.
Run from the paper folder: python verify/check_theorems_1_2_corollary_1.py"""
import numpy as np
rng = np.random.default_rng(1)
n = 10**6

def nb(t, out, p):            # empirical net benefit of policy t against outcome 'out'
    w = p/(1-p); return np.mean(t*out) - w*np.mean(t*(1-out))

def regions(r, T, Y, p):
    t = (r >= p).astype(float); q = t.mean()
    k1, k0 = t == 1, t == 0
    pi1, pi0 = T[k1].mean(), T[k0].mean()
    y1, y0 = Y[k1].mean(), Y[k0].mean()
    Se1 = Y[k1 & (T == 1)].mean(); Sp1 = 1-Y[k1 & (T == 0)].mean()
    Se0 = Y[k0 & (T == 1)].mean(); Sp0 = 1-Y[k0 & (T == 0)].mean()
    return t, q, pi1, pi0, y1, y0, Se1, Sp1, Se0, Sp0

# --- differential, region-specific misclassification (also depends on covariate x)
x = rng.normal(size=n)
r = 1/(1+np.exp(-(-1.5 + 1.2*x + 0.5*rng.normal(size=n))))
T = rng.binomial(1, 1/(1+np.exp(-(-1.5+1.4*x)))).astype(float)
Se_i = np.clip(0.5 + 0.4*r + 0.05*x, 0, 1)          # depends on r and x
Sp_i = np.clip(0.98 - 0.1*r, 0, 1)
Y = np.where(T == 1, rng.binomial(1, Se_i), rng.binomial(1, 1-Sp_i)).astype(float)
maxerr = 0
for p in [0.05, 0.1, 0.2, 0.3, 0.5]:
    t, q, pi1, pi0, y1, y0, Se1, Sp1, Se0, Sp0 = regions(r, T, Y, p)
    NB, NBs = nb(t, T, p), nb(t, Y, p)
    NBall, NBalls = nb(np.ones(n), T, p), nb(np.ones(n), Y, p)
    D, Ds = NB-NBall, NBs-NBalls
    errs = [NB - q*(pi1-p)/(1-p), D - (1-q)*(p-pi0)/(1-p),
            y1 - (Se1*pi1+(1-Sp1)*(1-pi1)), y0 - (Se0*pi0+(1-Sp0)*(1-pi0)),
            (NBs-NB) - q/(1-p)*((1-Sp1)*(1-pi1)-(1-Se1)*pi1),
            (Ds-D) - (1-q)/(1-p)*((1-Se0)*pi0-(1-Sp0)*(1-pi0))]
    maxerr = max(maxerr, max(abs(e) for e in errs))
print("Equation (1) and Theorem 1 (differential errors), max abs identity error:", maxerr)
assert maxerr < 1e-10

# --- Corollary 1: non-differential Se, Sp
Se, Sp = 0.8, 0.9; J = Se+Sp-1
Yn = np.where(T == 1, rng.binomial(1, Se, n), rng.binomial(1, 1-Sp, n)).astype(float)
print("Corollary 1 identity (population Se, Sp; differences are Monte Carlo error):")
for p in [0.12, 0.2, 0.3, 0.5, 0.7]:
    t, q, pi1, pi0, y1, y0, *_ = regions(r, T, Yn, p)
    d = (p-(1-Sp))/J; c = (Se-p)/(1-p)
    NBs = nb(t, Yn, p); Ds = NBs - nb(np.ones(n), Yn, p)
    NBt = q*(pi1-d)/(1-d); Dt = (1-q)*(d-pi0)/(1-d)
    print(f"  p={p}: NB*={NBs:.5f} cNB_t={c*NBt:.5f} | D*={Ds:.5f} cD_t={c*Dt:.5f} delta={d:.3f}")
# exact algebra on a grid (no MC), incl. p outside (1-Sp,Se)
g = np.linspace(0.01, 0.99, 99); bad = 0; outside = 0
for Se_ in g:
    for Sp_ in g:
        J_ = Se_+Sp_-1
        if J_ <= 0: continue
        for p in g:
            for pik in [0.0, 0.1, 0.5, 0.9, 1.0]:
                yk = J_*pik + 1-Sp_; d = (p-1+Sp_)/J_
                if abs(1-d) < 1e-12: continue
                lhs = (yk-p)/(1-p); rhs = (Se_-p)/(1-p)*(pik-d)/(1-d)
                bad += abs(lhs-rhs) > 1e-9
                outside += not (0 < d < 1)
print("Corollary 1 exact identity failures:", bad, "| cases with delta outside (0,1):", outside)
assert bad == 0

# floor/ceiling with perfect model r=T
for p in [0.05, 0.1, 0.8, 0.85]:
    t = (T >= p).astype(float)
    NBs = nb(t, Yn, p); Ds = NBs - nb(np.ones(n), Yn, p)
    print(f"  perfect model p={p}: D*={Ds:.4f} (floor 1-Sp=0.1)  NB*={NBs:.4f} (ceiling Se=0.8)")

# direction claims (c): sign comparisons over a grid of (pi_k, p)
pstar = (1-Sp)/((1-Se)+(1-Sp)); print("p* =", pstar)
viol = {"lo_all": 0, "lo_none": 0, "hi_all": 0, "hi_none": 0}; cnt = dict.fromkeys(viol, 0)
for p in np.linspace(0.11, 0.79, 69):
    d = (p-1+Sp)/J
    for pik in np.linspace(0, 1, 1001):
        yk = J*pik+1-Sp
        Dpos, Dspos = pik < p, yk < p          # sign of D and D*
        NBpos, NBspos = pik > p, yk > p
        if p < pstar:
            viol["lo_all"] += Dspos and not Dpos      # conservative: D*>0 => D>0
            viol["lo_none"] += NBpos and not NBspos   # anti-conservative: NB>0 => NB*>0
            cnt["lo_all"] += Dpos and not Dspos; cnt["lo_none"] += NBspos and not NBpos
        elif p > pstar:
            viol["hi_all"] += Dpos and not Dspos
            viol["hi_none"] += NBspos and not NBpos
            cnt["hi_all"] += Dspos and not Dpos; cnt["hi_none"] += NBpos and not NBspos
print("direction violations:", viol, " sign errors in the stated direction:", cnt)
assert sum(viol.values()) == 0

# --- Theorem 2: Sp0 = 1, uninformative score
ru = rng.uniform(size=n); prev = T.mean(); Se0 = 0.6
Yu = np.where(T == 1, rng.binomial(1, Se0, n), 0).astype(float)
print(f"Theorem 2: pi={prev:.4f}, Se0*pi={Se0*prev:.4f}")
worst2 = 0
for p in [0.08, 0.12, 0.15, 0.18, 0.21, 0.25]:
    t, q, pi1, pi0, y1, y0, Se1, Sp1, se0_hat, Sp0 = regions(ru, T, Yu, p)
    D = nb(t, T, p)-nb(np.ones(n), T, p); Ds = nb(t, Yu, p)-nb(np.ones(n), Yu, p)
    f = (1-q)*(1-se0_hat)*pi0/(1-p); worst2 = max(worst2, abs(Ds-D-f))
    print(f"  p={p}: D={D:+.4f} D*={Ds:+.4f} D*-D={Ds-D:.5f} formula={f:.5f}")
assert Sp0 == 1 and worst2 < 1e-10
print("all checks passed")
