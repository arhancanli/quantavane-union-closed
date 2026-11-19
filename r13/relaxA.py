# Round 13, G-side moment relaxation built ONLY from SU's certified local inputs (FLOAT scan; NOT a certificate).
# Per G-law moments: s in [m-eps, 1], e in [z, L], bbar in [0,1], Rbar >= 0. Identities (SU (2.12)):
#   t = (s - e/p - bbar)/c,  vbar = e - h_r t,  x = e - e* = gam (s-m) + alpha vbar - gam bbar.
# Certified inputs used: EF >= 2 kE t vbar^2 (F >= kE[v^2 phi + phi v^2] + Jensen); EF >= 2 kF t Rbar;
#   Jpart := q E_K J + qL Var(t_Z) >= q bbar (2t/5 - 1/(25L)) - (11q/25) Rbar  and Jpart >= 0  (localJ + b in [0,1]);
#   -sigma' <= b_m bbar + Rbar + b_m (m - s)  (point1 + (2.12));  e <= (4/5) t + 27/400  (hregular);  t <= 1.
# Exact: g = EF + Jpart,  U = EF - qL t^2 + 2 s e,  beta = g + 2 (s - m + eps) e.
# Loss  l = (A0+B') x_+ + (A0-B') x_- + max(0, b_m bbar + (1+C) Rbar + b_m (m-s)) + mu (a* e + d0 x - U).
# The adversary takes EF and Jpart minimal (both enter l - theta beta with negative coefficients -(mu+theta), -theta).
import numpy as np
from gateA import *

def relaxG_value(s, e, bb, Rb, cfg, theta):
    t = (s - e / p - bb) / c
    if t <= 0 or t > 1 or e > 0.8 * t + 27 / 400 + 1e-15: return -np.inf
    vb = e - hr * t
    x = e - es
    kE, kF = 11 / 50, 1 / 20
    EF = max(2 * kE * t * vb * vb, 2 * kF * t * Rb)
    Jp = max(0.0, q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    U = EF - q * L * t * t + 2 * s * e
    beta = EF + Jp + 2 * (s - m + cfg.eps) * e
    lx = (cfg.A0 + cfg.Bp) * max(x, 0) + (cfg.A0 - cfg.Bp) * max(-x, 0)
    lb = max(0.0, bm * bb + (1 + cfg.C) * Rb + bm * (m - s))
    lU = astar * e + d0 * x - U
    return (lx + lb + cfg.mu * lU - theta * beta) / e

def psiG_relax(theta, cfg, n=(60, 60, 40, 40), refine=3):
    best = (-np.inf, None)
    S = np.linspace(m - cfg.eps, 0.9, n[0]); E = np.linspace(cfg.z, L, n[1])
    B = np.concatenate([[0.0], np.logspace(-6, -0.5, n[2] - 1)]); Rg = np.concatenate([[0.0], np.logspace(-7, -1, n[3] - 1)])
    for s in S:
        for e in E:
            for bb in B:
                for Rb in Rg:
                    v = relaxG_value(s, e, bb, Rb, cfg, theta)
                    if v > best[0]: best = (v, (s, e, bb, Rb))
    # local refinement (coordinate search)
    from scipy.optimize import minimize
    x0 = np.array(best[1])
    lo_ = np.array([m - cfg.eps, cfg.z, 0.0, 0.0]); hi_ = np.array([1.0, L, 1.0, 1.0])
    f = lambda v: -relaxG_value(*np.clip(v, lo_, hi_), cfg, theta) if np.isfinite(relaxG_value(*np.clip(v, lo_, hi_), cfg, theta)) else 1e9
    for _ in range(refine):
        res = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
        x0 = np.clip(res.x, lo_, hi_)
    val = relaxG_value(*x0, cfg, theta)
    return (val, tuple(x0)) if val > best[0] else best
