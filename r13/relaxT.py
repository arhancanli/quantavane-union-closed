# T-side relaxation in (e, t, s) from SU's certified tail inequalities (5.x) plus Phi >= (c1 t)^2 (FLOAT scan).
import sys, numpy as np
from gateA import *
xs = np.linspace(1e-9, 1 - 1e-9, 400001)
rat = np.sqrt(xs * h(xs)) / np.maximum(phi(xs), 1e-300)
print("min sqrt(xh)/phi on (0,1):", rat.min(), "at", xs[np.argmin(rat)])
c1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.47
Db = 11 / (1 - 1 / (44 * q * L)); kE = 0.22
def psiT_relax(cfg, theta, n=(120, 160, 120)):
    E = np.concatenate([np.logspace(-6, np.log10(cfg.z), n[0])]); T = np.linspace(0, 1, n[1])
    S = np.concatenate([np.linspace(m - cfg.eps, m + 0.02, n[2] // 2), np.linspace(m + 0.02, 1, n[2] // 2)])
    e, t, s = np.meshgrid(E, T, S, indexing='ij')
    ok = (e <= 0.8 * t + 27 / 400) & (e <= 1.2 * t + 21 / 10000) & (e <= h(s) + 1e-15)
    U0 = np.maximum(809 / 500 * s * e, 2 * s * e - q * L * t * t)
    bb = s - e / p - c * t
    gLB = np.maximum.reduce([U0 + q * L * t * t - 2 * s * e, 2 * kE * t * (e - hr * t) ** 2,
                             t * t / Db * np.maximum(bb, 0), 4 / 25 * np.maximum(np.sqrt(L) * t - 1.4 * e, 0) ** 2])
    beta = gLB + 2 * (s - m + cfg.eps) * e
    val = ((m + cfg.gain + cfg.mu) * e - (c1 * t) ** 2 - cfg.mu * U0 - theta * beta) / e
    val[~ok] = -np.inf
    i = np.unravel_index(np.argmax(val), val.shape)
    return val[i], (E[i[0]], T[i[1]], S[i[2]] - m)
for eps in [1.5e-4, 1.7e-4, 2e-4]:
    cfg = Cfg(eps=eps)
    for th in [50, 60, 70]:
        v, arg = psiT_relax(cfg, th)
        print(f"eps={eps:.1e} theta={th}: psiT_relax={v:+.6f} at (e,t,s-m)=({arg[0]:.4g},{arg[1]:.4f},{arg[2]:+.2e})", flush=True)
