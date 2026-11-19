import sys, numpy as np
from gateA import *
from relaxA import relaxG_value
from scipy.optimize import minimize
def vec_vals(cfg, theta, S, E, B, Rg):
    s, e, bb, Rb = np.meshgrid(S, E, B, Rg, indexing='ij')
    t = (s - e / p - bb) / c
    ok = (t > 0) & (t <= 1) & (e <= 0.8 * t + 27 / 400)
    vb = e - hr * t; x = e - es
    EF = np.maximum(2 * 0.22 * t * vb * vb, 2 * 0.05 * t * Rb)
    Jp = np.maximum(0.0, q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    U = EF - q * L * t * t + 2 * s * e
    beta = EF + Jp + 2 * (s - m + cfg.eps) * e
    lx = (cfg.A0 + cfg.Bp) * np.maximum(x, 0) + (cfg.A0 - cfg.Bp) * np.maximum(-x, 0)
    lb = np.maximum(0.0, bm * bb + (1 + cfg.C) * Rb + bm * (m - s))
    val = (lx + lb + cfg.mu * (astar * e + d0 * x - U) - theta * beta) / e
    val[~ok] = -np.inf
    i = np.unravel_index(np.argmax(val), val.shape)
    return val[i], (S[i[0]], E[i[1]], B[i[2]], Rg[i[3]])
eps = float(sys.argv[1]) if len(sys.argv) > 1 else 2e-4
cfg = Cfg(eps=eps)
S = np.concatenate([np.linspace(m - eps, m + 0.02, 30), np.linspace(m + 0.02, 0.95, 30)]); E = np.linspace(cfg.z, L, 70)
B = np.concatenate([[0.0], np.logspace(-7, -0.3, 45)]); Rg = np.concatenate([[0.0], np.logspace(-8, -0.5, 45)])
lo_ = np.array([m - eps, cfg.z, 0.0, 0.0]); hi_ = np.array([1.0, L, 1.0, 1.0])
for theta in [20.0, 30.0, 40.0, 50.0, 70.0, 100.0]:
    v, arg = vec_vals(cfg, theta, S, E, B, Rg)
    f = lambda w: -relaxG_value(*np.clip(w, lo_, hi_), cfg, theta)
    x0 = np.array(arg)
    for _ in range(4):
        res = minimize(lambda w: f(w) if np.isfinite(f(w)) else 1e9, x0, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=6000))
        if -res.fun > v: v, x0 = -res.fun, np.clip(res.x, lo_, hi_)
    s_, e_, b_, R_ = x0
    print(f"theta={theta:6.1f} psiG_relax={v:+.6f} at s-m={s_-m:+.2e} e={e_:.5f} bbar={b_:.2e} Rbar={R_:.2e} | G0(if T ok)={cfg.gain - v - 2*eps*theta:+.6f}", flush=True)
