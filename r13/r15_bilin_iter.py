# Round 15 candidate B (FLOAT): bilinear moment kernels F(x,y) >= f(x)^T M f(y), f = (1, x, h, phi, R). Integrated against mu x mu:
# EF >= mbar^T M mbar, mbar = (1, s, e, t, Rbar): N3 stays 4-D. Loop: sup of the relaxation (round-14 cuts, SU kernels + K2 + the
# kernels so far) -> point P -> a law with moments close to P (grid LP) -> its exact moments P' -> grid LP for the best M at P'
# (bounded, <= EF of that law) -> add M -> repeat. Grid-LP kernels are only grid-valid (float estimate, not a certificate).
# Usage: r15_bilin_iter.py eps a0 theta iters
import sys, json, numpy as np
from scipy.optimize import linprog, minimize
args = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import *
from fractions import Fraction as Fr
PAIRS = [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]; RMAX = 0.2150
cuts = [np.array([float(Fr(v)) for v in row]) for row in json.load(open('r14_cuts_rat.json'))]
idx = [(i, j) for i in range(5) for j in range(i, 5)]
def fb(xs): return np.stack([np.ones_like(xs), xs, h(xs), phi(xs), R(xs)], axis=1)
xl = np.unique(np.concatenate([np.linspace(0, 1, 401), np.linspace(0.55, 0.75, 201), np.logspace(-12, -1, 60), 1 - np.logspace(-12, -1, 60), [r]]))
fl = fb(xl); XI, YI = np.triu_indices(len(xl))
Abil = np.stack([fl[XI, i] * fl[YI, j] + (fl[XI, j] * fl[YI, i] if i != j else 0 * fl[XI, i]) for i, j in idx], axis=1); bbil = F(xl[XI], xl[YI])
xg = np.unique(np.concatenate([np.linspace(0, 1, 4001), np.logspace(-10, -2, 200), 1 - np.logspace(-10, -2, 200), [r]])); Gg = fb(xg).T
def realise(P):
    n = len(xg); res = linprog(np.concatenate([np.zeros(n), np.ones(10)]), A_eq=np.hstack([Gg, np.eye(5), -np.eye(5)]), b_eq=P,
                               bounds=[(0, None)] * (n + 10), method='highs')
    w = res.x[:n]; return Gg @ w, res.fun
def bestM(P):
    obj = -np.array([P[i] * P[j] * (2 if i != j else 1) for i, j in idx])
    res = linprog(obj, A_ub=Abil, b_ub=bbil, bounds=[(-50, 50)] * len(idx), method='highs')
    M = np.zeros((5, 5))
    for (i, j), v in zip(idx, res.x): M[i, j] = M[j, i] = v
    return M, -res.fun
KER = []
def lam(S, E, BB, RB, cfg, theta):
    T = (S - E / p - BB) / c
    ok = (T > 0) & (T <= 1) & (E <= 0.8 * T + 27 / 400) & (E <= h(S))
    for a in cuts: ok &= (a[0] + a[1] * S + a[2] * E + a[3] * T + a[4] * RB >= 0)
    vb = E - hr * T; x = E - es
    EF = np.max([2 * T * (a_ * vb * vb + k_ * RB) for a_, k_ in PAIRS], axis=0)
    mv = [np.ones_like(S), S, E, T, RB]
    for M in KER: EF = np.maximum(EF, sum(M[i, j] * mv[i] * mv[j] for i in range(5) for j in range(5)))
    Jp = np.maximum(0.0, q * BB * (2 * T / 5 - 1 / (25 * L)) - 11 * q / 25 * RB)
    U = EF - q * L * T * T + 2 * S * E
    lx = (cfg.A0 + cfg.Bp) * np.maximum(x, 0) + (cfg.A0 - cfg.Bp) * np.maximum(-x, 0)
    lb = np.maximum(0.0, bm * BB + (1 + cfg.C) * RB + bm * (m - S))
    beta = EF + Jp + 2 * (S - m + cfg.eps) * E
    return np.where(ok, (lx + lb + cfg.mu * (astar * E + d0 * x - U) - theta * beta) / E, -np.inf)
def f1(w, cfg, th):
    s, e, bb, Rb = w
    if not (m - cfg.eps <= s <= 1 and cfg.z <= e <= L and 0 <= bb <= 1 and 0 <= Rb <= RMAX): return -np.inf
    return float(lam(np.array([s]), np.array([e]), np.array([bb]), np.array([Rb]), cfg, th)[0])
def sup(cfg, th):
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.01, 14), np.linspace(m + 0.01, 0.95, 14)]); E = np.linspace(cfg.z, L, 36)
    B = np.concatenate([[0.0], np.logspace(-7, 0, 22)]); Rg = np.concatenate([[0.0], np.logspace(-8, np.log10(RMAX), 26)])
    SS, EE, BB, RR = np.meshgrid(S, E, B, Rg, indexing='ij'); V = lam(SS, EE, BB, RR, cfg, th).ravel(); best = (-np.inf, None)
    for i in np.argsort(-V)[:6]:
        if not np.isfinite(V[i]): break
        v0 = np.array([SS.ravel()[i], EE.ravel()[i], BB.ravel()[i], RR.ravel()[i]])
        res = minimize(lambda w: -f1(w, cfg, th) if np.isfinite(f1(w, cfg, th)) else 1e9, v0, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=5000))
        for cand in (res.x, v0):
            vv = f1(cand, cfg, th)
            if vv > best[0]: best = (vv, tuple(cand))
    return best
eps, a0, th, iters = float(args[0]), float(args[1]), float(args[2]), int(args[3])
cfg = Cfg(eps=eps, mu=2.8, a0=a0); probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1)
for it in range(iters):
    v, (s, e, bb, Rb) = sup(cfg, th); t = (s - e / p - bb) / c
    P = np.array([1, s, e, t, Rb]); P2, resid = realise(P); M, val = bestM(P2)
    print(f"it={it} nker={len(KER)}: psiG={v:.6f} G0={cfg.gain - max(v, probe) - 2*eps*th:+.6f} at s-m={s-m:+.2e} e={e:.4f} b={bb:.2e} R={Rb:.2e}; "
          f"realised (resid {resid:.1e}), new kernel value there {val:.4e}", flush=True)
    KER.append(M); np.save('r15_bilin_kernels.npy', np.array(KER))
