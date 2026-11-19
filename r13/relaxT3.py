# T-side relaxation (same as N4's Lambda_T with the eight minorants of params.json) at lower theta (FLOAT grid).
import json, numpy as np
from fractions import Fraction as Fr
from gateA import *
MIN = [[float(Fr(c_)) for c_ in row] for row in json.load(open('params.json'))['minorants']]
Db = 11 / (1 - 1 / (44 * q * L)); kE = 0.22
def psiT(cfg, theta):
    E = np.logspace(-6, np.log10(cfg.z), 140); T = np.linspace(0, 1, 220)
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.02, 60), np.linspace(m + 0.02, 1, 60)])
    e, t, s = np.meshgrid(E, T, S, indexing='ij')
    ok = (e <= 0.8 * t + 27 / 400) & (e <= 1.2 * t + 21 / 10000) & (e <= h(s) + 1e-15)
    U0 = np.maximum(809 / 500 * s * e, 2 * s * e - q * L * t * t); bb = s - e / p - c * t
    gLB = np.maximum.reduce([U0 + q * L * t * t - 2 * s * e, 2 * kE * t * (e - hr * t) ** 2, t * t / Db * np.maximum(bb, 0),
                             4 / 25 * np.maximum(np.sqrt(L) * t - 1.4 * e, 0) ** 2])
    sig = np.zeros_like(e)
    for (a, b, c_, d) in MIN: sig = np.maximum(sig, a * t + b * e + c_ * s + d)
    val = ((m + cfg.gain + cfg.mu) * e - sig ** 2 - cfg.mu * U0 - theta * (gLB + 2 * (s - m + cfg.eps) * e)) / e
    val[~ok] = -np.inf; i = np.unravel_index(np.argmax(val), val.shape)
    return val[i], (E[i[0]], T[i[1]], S[i[2]] - m)
for th in [26.0, 30.0, 34.0]:
    cfg = Cfg(eps=2.8e-4, mu=2.8, a0=0.18); v, arg = psiT(cfg, th)
    print(f"theta={th}: psiT_relax={v:+.6f} at (e,t,s-m)=({arg[0]:.4g},{arg[1]:.4f},{arg[2]:+.2e})  [psiG ~ 0.0153 at theta 30]", flush=True)
