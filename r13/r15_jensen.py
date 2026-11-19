# Round 15 candidate A (FLOAT): drop the Jensen step E v^2 >= vbar^2. Add the moment w = E v(X)^2 and cuts on the curve
# g6(x) = (1, x, h, phi, R, v^2). Lambda_G decreases in w (EF_lo increases in w, coefficient -(mu + theta)), so N3 needs only
# w_lo(s, e, t, Rbar) = max(vbar^2, max over cuts with alpha_w > 0 of -(alpha_0..4 . (1, s, e, t, Rbar)) / alpha_w): N3 stays 4-D.
# Cutting plane on the 6-D hull. Usage: r15_jensen.py eps a0 theta [theta ...]
import sys, json, numpy as np
from scipy.optimize import linprog, minimize
args = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import *
from fractions import Fraction as Fr
PAIRS = [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]; RMAX = 0.2150
xs = np.unique(np.concatenate([np.linspace(0, 1, 20001), np.logspace(-15, -2, 300), 1 - np.logspace(-15, -2, 300)])); xs = xs[(xs >= 0) & (xs <= 1)]
G6 = np.stack([np.ones_like(xs), xs, h(xs), phi(xs), R(xs), vfun(xs) ** 2], axis=1)
def sep(P):
    res = linprog(P, A_ub=-G6, b_ub=np.zeros(len(xs)), bounds=[(-1, 1)] * 6, method='highs'); return res.fun, res.x
cuts = [np.append(np.array([float(Fr(v)) for v in row]), 0.0) for row in json.load(open('r14_cuts_rat.json'))]   # round-14 cuts (alpha_w = 0)
import os
if os.environ.get('R15_CONTINUE'): cuts = [np.array(c_) for c_ in json.load(open('r15_cuts6.json'))]   # continue from saved cuts
def lam(S, E, BB, RB, cfg, theta, cuts):
    T = (S - E / p - BB) / c
    ok = (T > 0) & (T <= 1) & (E <= 0.8 * T + 27 / 400) & (E <= h(S))
    vb = E - hr * T; Wlo = vb * vb; Whi = np.full_like(S, np.inf)
    for a in cuts:
        lin = a[0] + a[1] * S + a[2] * E + a[3] * T + a[4] * RB
        if a[5] > 1e-12: Wlo = np.maximum(Wlo, -lin / a[5])
        elif a[5] < -1e-12: Whi = np.minimum(Whi, lin / (-a[5]))
        else: ok &= (lin >= 0)
    ok &= (Wlo <= Whi)
    x = E - es
    EF = np.max([2 * T * (a_ * Wlo + k_ * RB) for a_, k_ in PAIRS], axis=0)
    Jp = np.maximum(0.0, q * BB * (2 * T / 5 - 1 / (25 * L)) - 11 * q / 25 * RB)
    U = EF - q * L * T * T + 2 * S * E
    lx = (cfg.A0 + cfg.Bp) * np.maximum(x, 0) + (cfg.A0 - cfg.Bp) * np.maximum(-x, 0)
    lb = np.maximum(0.0, bm * BB + (1 + cfg.C) * RB + bm * (m - S))
    beta = EF + Jp + 2 * (S - m + cfg.eps) * E
    return np.where(ok, (lx + lb + cfg.mu * (astar * E + d0 * x - U) - theta * beta) / E, -np.inf), Wlo, T
def f1(w, cfg, theta, cuts):
    s, e, bb, Rb = w
    if not (m - cfg.eps <= s <= 1 and cfg.z <= e <= L and 0 <= bb <= 1 and 0 <= Rb <= RMAX): return -np.inf, 0, 0
    v, W, T = lam(np.array([s]), np.array([e]), np.array([bb]), np.array([Rb]), cfg, theta, cuts); return float(v[0]), float(W[0]), float(T[0])
def sup(cfg, theta, cuts):
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.01, 12), np.linspace(m + 0.01, 0.95, 12)]); E = np.linspace(cfg.z, L, 36)
    B = np.concatenate([[0.0], np.logspace(-7, 0, 22)]); Rg = np.concatenate([[0.0], np.logspace(-8, np.log10(RMAX), 26)])
    SS, EE, BB, RR = np.meshgrid(S, E, B, Rg, indexing='ij'); V = lam(SS, EE, BB, RR, cfg, theta, cuts)[0].ravel()
    best = (-np.inf, None)
    for i in np.argsort(-V)[:6]:
        if not np.isfinite(V[i]): break
        v0 = np.array([SS.ravel()[i], EE.ravel()[i], BB.ravel()[i], RR.ravel()[i]])
        res = minimize(lambda w: -f1(w, cfg, theta, cuts)[0] if np.isfinite(f1(w, cfg, theta, cuts)[0]) else 1e9, v0, method='Nelder-Mead',
                       options=dict(xatol=1e-13, fatol=1e-15, maxiter=5000))
        for cand in (res.x, v0):
            vv = f1(cand, cfg, theta, cuts)[0]
            if vv > best[0]: best = (vv, tuple(cand))
    return best
eps, a0 = float(args[0]), float(args[1])
cfg = Cfg(eps=eps, mu=2.8, a0=a0); probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1)
for th in [float(v) for v in args[2:]]:
    for it in range(30):
        v, arg = sup(cfg, th, cuts); s, e, bb, Rb = arg; _, W, T = f1(arg, cfg, th, cuts)
        sv, a = sep(np.array([1, s, e, T, Rb, W]))
        print(f"eps={eps:.2e} a0={a0} theta={th} it={it} ncuts={len(cuts)}: psiG={v:.6f} G0={cfg.gain - max(v, probe) - 2*eps*th:+.6f} at s-m={s-m:+.2e} e={e:.4f} b={bb:.2e} R={Rb:.2e} w={W:.3e} sep={sv:+.2e}", flush=True)
        if sv > -1e-7: break
        cuts.append(a)
json.dump([list(map(float, c_)) for c_ in cuts], open('r15_cuts6.json', 'w'))
