# Step 3 refinement (FLOAT): ideal per-law G0 near its zero; plus an exact-Phi T-side check at the best point.
import sys, numpy as np
from gateA import *
from ideal_eps import psiG_ideal
mu = 2.855121
best = {}
for eps in [3.6e-4, 3.7e-4, 3.75e-4]:
    for a0 in [0.12, 0.15, 0.18]:
        cfg = Cfg(eps=eps, mu=mu, a0=a0); probe = m + cfg.gain - mu * (2 * (m - eps) - 1)
        for th in [21.0, 23.0, 25.0]:
            v = psiG_ideal(cfg, th); G0 = cfg.gain - max(v, probe) - 2 * eps * th
            print(f"eps={eps:.2e} a0={a0} theta={th}: gain={cfg.gain:.6f} psiG={v:.6f} probe={probe:+.5f} G0={G0:+.6f}", flush=True)
            if eps not in best or G0 > best[eps][0]: best[eps] = (G0, a0, th)
    print(f"BEST eps={eps:.2e}: G0={best[eps][0]:+.6f} (a0={best[eps][1]}, theta={best[eps][2]})", flush=True)
# exact-Phi T side at the best point of the smallest eps
eps = 3.6e-4; G0, a0, th = best[eps]; cfg = Cfg(eps=eps, mu=mu, a0=a0, Tmode='Phi')
phs = 1 - m / r
vT = -np.inf
for (k, Zc, ns) in [(3, 2, 60), (4, 2, 40)]:
    seeds = [seed_vec(a, b, k, Zc, k, cfg, 'T') for a, b in [([0, 0.801], [1 - 0.7715, 0.7715]), ([0, 1, 0.7], [0.38, 0.55, 0.07]), ([0, 0.5, 0.9], [0.2, 0.3, 0.5])] if len(a) <= k]
    val, v = psi(th, cfg, 'T', k=k, Zc=Zc, Zd=k, nstart=ns, seed=5 + k, seeds=seeds)
    vT = max(vT, val)
print(f"T side (exact Phi) at eps={eps:.2e} a0={a0} theta={th}: psiT={vT:+.6f}", flush=True)
