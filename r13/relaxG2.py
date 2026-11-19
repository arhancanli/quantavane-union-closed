# Refined G-side relaxation sup (grid + Nelder-Mead from the 8 best grid points), with Rbar <= 0.2150 and e <= h(s). FLOAT.
import sys, numpy as np
from scipy.optimize import minimize
from gateA import *
RMAX = 0.2150
def Fval(v, cfg, theta):
    s, e, bb, Rb = v
    if not (m - cfg.eps <= s <= 1 and cfg.z <= e <= L and 0 <= bb <= 1 and 0 <= Rb <= RMAX): return -np.inf
    t = (s - e / p - bb) / c
    if not (0 < t <= 1) or e > 0.8 * t + 27 / 400 or e > h(s): return -np.inf
    vb = e - hr * t; x = e - es
    EF = max(2 * 0.22 * t * vb * vb, 2 * 0.05 * t * Rb)
    Jp = max(0.0, q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    U = EF - q * L * t * t + 2 * s * e
    lx = (cfg.A0 + cfg.Bp) * max(x, 0) + (cfg.A0 - cfg.Bp) * max(-x, 0)
    lb = max(0.0, bm * bb + (1 + cfg.C) * Rb + bm * (m - s))
    beta = EF + Jp + 2 * (s - m + cfg.eps) * e
    return (lx + lb + cfg.mu * (astar * e + d0 * x - U) - theta * beta) / e
def sup(cfg, theta):
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.01, 20), np.linspace(m + 0.01, 0.95, 20)]); E = np.linspace(cfg.z, L, 50)
    B = np.concatenate([[0.0], np.logspace(-7, 0, 30)]); Rg = np.concatenate([[0.0], np.logspace(-8, np.log10(RMAX), 30)])
    pts = []
    for s in S:
        for e in E:
            for bb in B:
                for Rb in Rg:
                    pts.append((Fval((s, e, bb, Rb), cfg, theta), (s, e, bb, Rb)))
    pts.sort(key=lambda a: -a[0])
    best = pts[0]
    for val, v0 in pts[:8]:
        res = minimize(lambda w: -Fval(w, cfg, theta) if np.isfinite(Fval(w, cfg, theta)) else 1e9, np.array(v0), method='Nelder-Mead',
                       options=dict(xatol=1e-13, fatol=1e-15, maxiter=8000))
        vv = Fval(res.x, cfg, theta)
        if vv > best[0]: best = (vv, tuple(res.x))
    return best
for eps in [float(a) for a in sys.argv[1:]]:
    cfg = Cfg(eps=eps)
    for th in [50.0, 60.0, 70.0]:
        v, arg = sup(cfg, th)
        print(f"eps={eps:.2e} theta={th}: psiG={v:+.6f} G0={cfg.gain - max(v, -0.0052) - 2*eps*th:+.6f} at s-m={arg[0]-m:+.2e} e={arg[1]:.5f} b={arg[2]:.2e} R={arg[3]:.2e}", flush=True)
