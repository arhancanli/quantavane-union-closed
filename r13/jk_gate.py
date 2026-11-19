# Step 3 FLOAT gate: G-side relaxation with EF_lo = max over grid-valid joint kernels (a, kappa') of 2 t (a vbar^2 + kappa' Rbar),
# including SU's (0.22, 0) and (0, 0.05). Same everything else as relaxG2. Reports the best G0 over (a0, theta) at several eps.
import sys, itertools, numpy as np
from scipy.optimize import minimize
from gateA import *
PAIRS = [(0.22, 0.0), (0.0, 0.05), (0.22, 0.0005), (0.21, 0.0040), (0.20, 0.0078), (0.18, 0.0150), (0.15, 0.0255), (0.10, 0.0418),
         (0.05, 0.0535), (0.0, 0.0583)]
RMAX = 0.2150
def Fval(v, cfg, theta, pairs):
    s, e, bb, Rb = v
    if not (m - cfg.eps <= s <= 1 and cfg.z <= e <= L and 0 <= bb <= 1 and 0 <= Rb <= RMAX): return -np.inf
    t = (s - e / p - bb) / c
    if not (0 < t <= 1) or e > 0.8 * t + 27 / 400 or e > h(s): return -np.inf
    vb = e - hr * t; x = e - es
    EF = max(2 * t * (a * vb * vb + k * Rb) for a, k in pairs)
    Jp = max(0.0, q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    U = EF - q * L * t * t + 2 * s * e
    lx = (cfg.A0 + cfg.Bp) * max(x, 0) + (cfg.A0 - cfg.Bp) * max(-x, 0)
    lb = max(0.0, bm * bb + (1 + cfg.C) * Rb + bm * (m - s))
    beta = EF + Jp + 2 * (s - m + cfg.eps) * e
    return (lx + lb + cfg.mu * (astar * e + d0 * x - U) - theta * beta) / e
def sup(cfg, theta, pairs):
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.01, 16), np.linspace(m + 0.01, 0.95, 16)]); E = np.linspace(cfg.z, L, 40)
    B = np.concatenate([[0.0], np.logspace(-7, 0, 24)]); Rg = np.concatenate([[0.0], np.logspace(-8, np.log10(RMAX), 30)])
    pts = sorted(((Fval((s, e, bb, Rb), cfg, theta, pairs), (s, e, bb, Rb)) for s in S for e in E for bb in B for Rb in Rg), key=lambda a_: -a_[0])
    best = pts[0]
    for val, v0 in pts[:6]:
        res = minimize(lambda w: -Fval(w, cfg, theta, pairs) if np.isfinite(Fval(w, cfg, theta, pairs)) else 1e9, np.array(v0),
                       method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=6000))
        vv = Fval(res.x, cfg, theta, pairs)
        if vv > best[0]: best = (vv, tuple(res.x))
    return best
for eps in ([] if __name__ != "__main__" else [float(a_) for a_ in sys.argv[1:]]):
    for label, pairs in [('SU kernels', PAIRS[:2]), ('joint frontier', PAIRS)]:
        best = None
        for a0, th in itertools.product([0.18, 0.21, 0.24], [36.0, 40.0, 44.0]):
            cfg = Cfg(eps=eps, mu=2.8, a0=a0); v, arg = sup(cfg, th, pairs)
            probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1); G0 = cfg.gain - max(v, probe) - 2 * eps * th
            if best is None or G0 > best[0]: best = (G0, a0, th, v, arg)
        G0, a0, th, v, arg = best
        print(f"eps={eps:.2e} [{label}]: best G0={G0:+.6f} at a0={a0} theta={th}; psiG={v:.6f} at s-m={arg[0]-m:+.1e} e={arg[1]:.4f} b={arg[2]:.1e} R={arg[3]:.2e}", flush=True)
