# T-side relaxation with sigma >= max_k (a_k t + b_k e + c_k s + d_k) (LP minorants, rounded DOWN by 1e-4 for slack). FLOAT.
import sys, numpy as np
from scipy.optimize import linprog
from gateA import *
xs = np.unique(np.concatenate([np.linspace(0, 1, 20001), np.logspace(-14, -2, 400), 1 - np.logspace(-14, -2, 400)]))
f = np.sqrt(np.maximum(xs * h(xs), 0)); A = np.stack([phi(xs), h(xs), xs, np.ones_like(xs)], axis=1)
pts = [(0.7296, 0.385, m), (0.70, 0.37, m), (0.65, 0.33, m), (0.55, 0.27, m), (0.4, 0.2, m), (0.25, 0.12, m), (0.75, 0.385, 0.63), (0.6, 0.385, 0.66)]
MIN = []
for (t0, e0, s0) in pts:
    res = linprog(-np.array([t0, e0, s0, 1.0]), A_ub=A, b_ub=f - 1e-4, bounds=[(None, None)] * 3 + [(None, 0)], method='highs')
    MIN.append(res.x)
MIN = np.array(MIN)
np.save('sigma_minorants.npy', MIN)
print("minorants (a,b,c,d):\n", np.round(MIN, 5))
Db = 11 / (1 - 1 / (44 * q * L)); kE = 0.22
def psiT_relax(cfg, theta, n=(140, 200, 120)):
    E = np.logspace(-6, np.log10(cfg.z), n[0]); T = np.linspace(0, 1, n[1])
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.02, n[2] // 2), np.linspace(m + 0.02, 1, n[2] // 2)])
    e, t, s = np.meshgrid(E, T, S, indexing='ij')
    ok = (e <= 0.8 * t + 27 / 400) & (e <= 1.2 * t + 21 / 10000) & (e <= h(s) + 1e-15)
    U0 = np.maximum(809 / 500 * s * e, 2 * s * e - q * L * t * t)
    bb = s - e / p - c * t
    gLB = np.maximum.reduce([U0 + q * L * t * t - 2 * s * e, 2 * kE * t * (e - hr * t) ** 2,
                             t * t / Db * np.maximum(bb, 0), 4 / 25 * np.maximum(np.sqrt(L) * t - 1.4 * e, 0) ** 2])
    sig = np.zeros_like(e)
    for (a, b, c_, d) in MIN: sig = np.maximum(sig, a * t + b * e + c_ * s + d)
    beta = gLB + 2 * (s - m + cfg.eps) * e
    val = ((m + cfg.gain + cfg.mu) * e - sig ** 2 - cfg.mu * U0 - theta * beta) / e
    val[~ok] = -np.inf
    i = np.unravel_index(np.argmax(val), val.shape)
    return val[i], (E[i[0]], T[i[1]], S[i[2]] - m)
for eps in [1.5e-4, 2e-4]:
    cfg = Cfg(eps=eps)
    for th in [40, 50, 60, 70]:
        v, arg = psiT_relax(cfg, th)
        print(f"eps={eps:.1e} theta={th}: psiT_relax={v:+.6f} at (e,t,s-m)=({arg[0]:.4g},{arg[1]:.4f},{arg[2]:+.2e})", flush=True)
