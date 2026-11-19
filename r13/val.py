import numpy as np
from gateA import *
cfg = Cfg(eps=2e-4)
print(f"D={D:.6f} t0={t0:.6f} e*={es:.6f} a*={astar:.6f} d0={d0:.6f} eta={cfg.eta:.6f} W={cfg.W:.6f} B'={cfg.Bp:.6f} mu={cfg.mu:.6f} gain={cfg.gain:.6f}")
def st_of(xs, ws, Zc=1):
    x = np.array(xs, float); M = np.array([ws], float)
    return law_stats(x, M, None, cfg)
phs = 1 - m / r
s = st_of([0, r], [phs, 1 - phs]); print("contact:", {k_: round(v, 10) for k_, v in s.items() if k_ != 'w'}, "lG", lossG(s, cfg))
for dl in [-1e-2, -1e-3, 1e-3]:
    y = r + dl; pi = m / y
    s = st_of([0, y], [1 - pi, pi]); print(f"perturbed d={dl}: x/sqrt(g)={s['x']/np.sqrt(s['g']):.4f} lG/sqrt(beta)={lossG(s,cfg)/np.sqrt(s['beta']):.4f} x={s['x']:.3e} g={s['g']:.3e}")
s = st_of([0, 0.801], [1 - 0.7715, 0.7715]); print("two-point 0.801:", round(s['U'], 5), round(s['g'], 6), round(s['e'], 5), "lT", lossT(s, cfg))
s = st_of([m], [1.0]); print("point mass m: g =", round(s['g'], 6), "x =", round(s['x'], 5))
