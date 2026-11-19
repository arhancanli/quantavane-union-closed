# Round 14 step 2 (FLOAT): moment-hull cuts for the G-side relaxation (N3). Each cut alpha . (1, x, h(x), phi(x), R(x)) >= 0 on [0,1]
# is a 1-D inequality; integrated it gives alpha . (1, s, e, t, Rbar) >= 0, which removes unrealisable points of D_G.
# Cutting-plane loop: sup of Lambda_G/e over D_G and the cuts (vectorised grid + Nelder-Mead), separate the maximiser, add the cut.
# Library form of r14_cutplane.py (functions only).   (kernels: SU's two + K2 = (3/20, 51/2000), as certified at 2.7e-4)
import sys, json, numpy as np
from scipy.optimize import linprog, minimize
sys.argv = sys.argv[:1]
from gateA import *
PAIRS = [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]; RMAX = 0.2150
xs = np.unique(np.concatenate([np.linspace(0, 1, 20001), np.logspace(-15, -2, 300), 1 - np.logspace(-15, -2, 300)])); xs = xs[(xs >= 0) & (xs <= 1)]
Gc = np.stack([np.ones_like(xs), xs, h(xs), phi(xs), R(xs)], axis=1)
def sep(P):
    res = linprog(P, A_ub=-Gc, b_ub=np.zeros(len(xs)), bounds=[(-1, 1)] * 5, method='highs'); return res.fun, res.x
def lam(S, E, BB, RB, cfg, theta, cuts):
    T = (S - E / p - BB) / c
    ok = (T > 0) & (T <= 1) & (E <= 0.8 * T + 27 / 400) & (E <= h(S))
    for a in cuts: ok &= (a[0] + a[1] * S + a[2] * E + a[3] * T + a[4] * RB >= 0)
    vb = E - hr * T; x = E - es
    EF = np.max([2 * T * (a_ * vb * vb + k_ * RB) for a_, k_ in PAIRS], axis=0)
    Jp = np.maximum(0.0, q * BB * (2 * T / 5 - 1 / (25 * L)) - 11 * q / 25 * RB)
    U = EF - q * L * T * T + 2 * S * E
    lx = (cfg.A0 + cfg.Bp) * np.maximum(x, 0) + (cfg.A0 - cfg.Bp) * np.maximum(-x, 0)
    lb = np.maximum(0.0, bm * BB + (1 + cfg.C) * RB + bm * (m - S))
    beta = EF + Jp + 2 * (S - m + cfg.eps) * E
    return np.where(ok, (lx + lb + cfg.mu * (astar * E + d0 * x - U) - theta * beta) / E, -np.inf)
def f1(w, cfg, theta, cuts):
    s, e, bb, Rb = w
    if not (m - cfg.eps <= s <= 1 and cfg.z <= e <= L and 0 <= bb <= 1 and 0 <= Rb <= RMAX): return -np.inf
    return float(lam(np.array([s]), np.array([e]), np.array([bb]), np.array([Rb]), cfg, theta, cuts)[0])
def sup(cfg, theta, cuts):
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.01, 16), np.linspace(m + 0.01, 0.95, 16)]); E = np.linspace(cfg.z, L, 40)
    B = np.concatenate([[0.0], np.logspace(-7, 0, 24)]); Rg = np.concatenate([[0.0], np.logspace(-8, np.log10(RMAX), 30)])
    SS, EE, BB, RR = np.meshgrid(S, E, B, Rg, indexing='ij'); V = lam(SS, EE, BB, RR, cfg, theta, cuts).ravel()
    idx = np.argsort(-V)[:8]; best = (-np.inf, None)
    for i in idx:
        if not np.isfinite(V[i]): break
        v0 = np.array([SS.ravel()[i], EE.ravel()[i], BB.ravel()[i], RR.ravel()[i]])
        res = minimize(lambda w: -f1(w, cfg, theta, cuts) if np.isfinite(f1(w, cfg, theta, cuts)) else 1e9, v0, method='Nelder-Mead',
                       options=dict(xatol=1e-13, fatol=1e-15, maxiter=6000))
        for cand in (res.x, v0):
            vv = f1(cand, cfg, theta, cuts)
            if vv > best[0]: best = (vv, tuple(cand))
    return best
