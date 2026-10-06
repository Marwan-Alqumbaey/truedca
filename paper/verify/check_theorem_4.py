"""Theorem 4: closed forms of h, optimality of the design rho* (iii), and unbiasedness and variance of the
chart review estimator (i)-(ii) against exact targets computed by quadrature.
Run from the paper folder: python verify/check_theorem_4.py"""
import numpy as np
from scipy.integrate import quad
from scipy.stats import beta
rng = np.random.default_rng(3)
pL, pU = 0.05, 0.30

# --- closed forms of h vs quadrature (g uniform on [pL,pU])
def h_all(r): return 0.0 if r >= pU else (1/(1-pU)-1/(1-max(r, pL)))/(pU-pL)
def h_none(r): return 0.0 if r < pL else (1/(1-min(r, pU))-1/(1-pL))/(pU-pL)
e = 0
for r in np.linspace(0, 1, 401):
    qa = quad(lambda p: (r < p)/(1-p)**2/(pU-pL), pL, pU, points=[r] if pL < r < pU else None)[0]
    qn = quad(lambda p: (r >= p)/(1-p)**2/(pU-pL), pL, pU, points=[r] if pL < r < pU else None)[0]
    e = max(e, abs(qa-h_all(r)), abs(qn-h_none(r)))
print("h_all/h_none closed form vs quadrature, max abs err:", e)
assert e < 1e-10

# --- data-generating model with known mu(r,Y)=E[T|r,Y]
def gen(n):
    r = rng.beta(2, 6, n); T = rng.binomial(1, r)
    Se = 0.5+0.4*r; Sp = 0.97
    Y = np.where(T == 1, rng.binomial(1, Se), rng.binomial(1, 1-Sp))
    m1 = r*Se/(r*Se+(1-r)*(1-Sp)); m0 = r*(1-Se)/(r*(1-Se)+(1-r)*Sp)
    mu = np.where(Y == 1, m1, m0); return r, Y, T, mu

grid = np.linspace(pL, pU, 26); gw = np.full(grid.size, 1/grid.size)  # discretised uniform g
def h_all_disc(r): return ((r[:, None] < grid[None, :])/(1-grid[None, :])**2 * gw).sum(1)

N = 2*10**6; rB, YB, TB, muB = gen(N); vB = muB*(1-muB); hB = h_all_disc(rB)
kappa = 0.10
def lam_design(score):               # pi = min(1, lam*score) with E[pi]=kappa
    lo, hi = 0, 1e6
    for _ in range(200):
        m = (lo+hi)/2; (lo, hi) = (m, hi) if np.minimum(1, m*score).mean() < kappa else (lo, m)
    return lambda s: np.minimum(1, (lo+hi)/2*s)
def phase2(pi_vals): # E[v h (1/pi - 1)], pi=0 allowed where v*h=0
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(vB*hB > 0, vB*hB*(1/pi_vals-1), 0.0)
    return term.mean()
f_star = lam_design(np.sqrt(vB*hB)); pis = f_star(np.sqrt(vB*hB))
designs = {"SRS": np.full(N, kappa),
           "prop. to h": lam_design(hB)(hB),
           "prop. to sqrt(v)": lam_design(np.sqrt(vB))(np.sqrt(vB)),
           "model-neg only (r<pU)": lam_design((rB < pU)*1.0)((rB < pU)*1.0),
           "pi*": pis}
for s in [0.1, 0.3]:   # perturbed pi*, renormalised to same E[pi]
    z = np.sqrt(vB*hB)*np.exp(s*rng.normal(size=N)); designs[f"pi* perturbed s={s}"] = lam_design(z)(z)
print(f"min pi* = {pis.min():.3g}, share pi*=0: {(pis == 0).mean():.3f}, share capped: {(pis == 1).mean():.3f}")
for k, pv in designs.items(): print(f"  {k:24s} E[pi]={pv.mean():.4f}  phase-2 term={phase2(pv):.5f}")
# small deterministic perturbations around pi* (same mean): directional derivative check
worse = 0
for _ in range(50):
    u = rng.normal(size=N)*(pis > 0)*(pis < 1); u -= u[(pis > 0) & (pis < 1)].mean()*((pis > 0) & (pis < 1))
    pv = np.clip(pis*(1+0.05*u), 1e-6, 1); pv *= kappa/pv.mean()
    worse += phase2(pv) >= phase2(pis)
print("random same-mean perturbations not improving on pi*:", worse, "/ 50")
assert worse == 50 and all(phase2(pv) >= phase2(pis) - 1e-12 for pv in designs.values())


# --- exact targets: D(p) and V(p) by quadrature; Monte Carlo under SRS and a stratified design
rng = np.random.default_rng(4)
zs = []
p, n, reps = 0.15, 20000, 6000
Sp = 0.97
f = beta(2, 6).pdf
def mus(r):
    Se = 0.5+0.4*r
    py1 = r*Se+(1-r)*(1-Sp)
    return py1, r*Se/py1, r*(1-Se)/(1-py1)
def design(r, Y):  # stratified on (r<p, Y): oversample model-negatives with Y=0
    return np.where(r < p, np.where(Y == 0, 0.25, 0.6), 0.02)
# exact D and V
D = quad(lambda r: (p-r)*f(r), 0, p)[0]/(1-p)
E2 = quad(lambda r: (p**2-2*p*r+r)*f(r), 0, p)[0]/(1-p)**2   # E[psi^2], E[T^2|r]=r
Vpsi = E2-D**2
def ph2(pi1, pi0):
    def g(r):
        py1, m1, m0 = mus(r)
        return f(r)*(py1*m1*(1-m1)*(1-pi1)/pi1+(1-py1)*m0*(1-m0)*(1-pi0)/pi0)
    return quad(g, 0, p)[0]/(1-p)**2
designs = {"SRS 0.1": (lambda r, Y: np.full(r.size, 0.1), ph2(0.1, 0.1)),
           "stratified": (design, ph2(0.6, 0.25))}
for name, (dfn, P2) in designs.items():
    V = Vpsi+P2; est = np.empty(reps)
    for b in range(reps):
        r = rng.beta(2, 6, n); T = rng.binomial(1, r)
        Y = np.where(T == 1, rng.binomial(1, 0.5+0.4*r), rng.binomial(1, 1-Sp, n))
        _, m1, m0 = mus(r); mu = np.where(Y == 1, m1, m0)
        pi = dfn(r, Y); R = rng.uniform(size=n) < pi
        Tt = mu+R*(T-mu)/pi
        est[b] = np.mean((r < p)*(p-Tt)/(1-p))
    se_var = est.var()*np.sqrt(2/(reps-1))
    print(f"{name}: mean={est.mean():.6f} D={D:.6f} z={(est.mean()-D)/(est.std()/np.sqrt(reps)):+.2f} | "
          f"emp var={est.var():.4e} V/n={V/n:.4e} z={(est.var()-V/n)/se_var:+.2f}")
    zs += [(est.mean()-D)/(est.std()/np.sqrt(reps)), (est.var()-V/n)/se_var]
assert all(abs(z) < 4 for z in zs)
print("all checks passed")
