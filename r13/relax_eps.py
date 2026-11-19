# eps* of the G-side relaxation (SU certified inputs; Rbar <= Rmax; bbar <= 1), assuming the T side holds (checked separately).
import numpy as np
from gateA import *
xs = np.linspace(0, 1, 200001); Rmax = float(R(xs).max()); print("Rmax =", Rmax, "at x =", xs[np.argmax(R(xs))])
def G0_of(eps, thetas):
    cfg = Cfg(eps=eps)
    S = np.concatenate([np.linspace(m - eps, m + 0.01, 25), np.linspace(m + 0.01, 0.95, 25)]); E = np.linspace(cfg.z, L, 60)
    B = np.concatenate([[0.0], np.logspace(-7, 0, 40)]); Rg = np.concatenate([[0.0], np.logspace(-8, np.log10(Rmax), 40)])
    s, e, bb, Rb = np.meshgrid(S, E, B, Rg, indexing='ij')
    t = (s - e / p - bb) / c
    ok = (t > 0) & (t <= 1) & (e <= 0.8 * t + 27 / 400)
    vb = e - hr * t; x = e - es
    EF = np.maximum(2 * 0.22 * t * vb * vb, 2 * 0.05 * t * Rb)
    Jp = np.maximum(0.0, q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    U = EF - q * L * t * t + 2 * s * e
    lx = (cfg.A0 + cfg.Bp) * np.maximum(x, 0) + (cfg.A0 - cfg.Bp) * np.maximum(-x, 0)
    lb = np.maximum(0.0, bm * bb + (1 + cfg.C) * Rb + bm * (m - s))
    base = lx + lb + cfg.mu * (astar * e + d0 * x - U)
    beta = EF + Jp + 2 * (s - m + cfg.eps) * e
    out = []
    for th in thetas:
        val = (base - th * beta) / e; val[~ok] = -np.inf
        out.append((cfg.gain - val.max() - 2 * eps * th, th))
    return max(out)
for eps in [1.0e-4, 1.4e-4, 1.7e-4, 1.9e-4, 2.0e-4, 2.2e-4]:
    G0, th = G0_of(eps, [30, 40, 50, 60, 70, 80, 100, 130])
    print(f"eps={eps:.1e}: best G0 = {G0:+.6f} at theta={th}", flush=True)
