# Round 14 (FLOAT): T-side sup of Lambda_T/e on D_T (as relaxT2.psiT_relax) with the CERTIFIED minorants of params_e270.json.
# Usage: r14_Tside.py eps a0 theta [theta ...]
import sys, json, numpy as np
from fractions import Fraction as Fr
args = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import *
MIN = np.array([[float(Fr(v)) for v in row] for row in json.load(open('params_e270.json'))['minorants']])
Db = 11 / (1 - 1 / (44 * q * L)); kE = 0.22
def psiT(cfg, theta, n=(160, 240, 140)):
    E = np.logspace(-6, np.log10(cfg.z), n[0]); T = np.linspace(0, 1, n[1])
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.02, n[2] // 2), np.linspace(m + 0.02, 1, n[2] // 2)])
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
eps, a0 = float(args[0]), float(args[1]); cfg = Cfg(eps=eps, mu=2.8, a0=a0)
for th in [float(a) for a in args[2:]]:
    v, arg = psiT(cfg, th)
    print(f"eps={eps:.2e} a0={a0} theta={th}: psiT={v:+.6f} at (e,t,s-m)=({arg[0]:.4g},{arg[1]:.4f},{arg[2]:+.2e})", flush=True)
