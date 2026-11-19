# Step 3 (FLOAT): size of route JB. IDEAL per-law version: exact U, V (best CIID coupling, up to 3 latent states), g, exact
# per-law losses; G0 = gain - max(psiG, probeT) - 2 eps theta. psiT (exact Phi) was far below psiG in every run; probes bound it.
import sys, numpy as np
from gateA import *
phs = 1 - m / r
def psiG_ideal(cfg, theta):
    best = -np.inf
    for (k, Zc, ns) in [(3, 2, 60), (4, 3, 40)]:
        seeds = [seed_vec(a, b, k, Zc, 0, cfg, 'G') for a, b in
                 [([0, r], [phs, 1 - phs]), ([0, 0.66], [0.065, 0.935]), ([0, 0.65], [0.06, 0.94]), ([0, 0.64], [0.055, 0.945]),
                  ([m], [1.0]), ([0, 0.63, 0.66], [0.05, 0.5, 0.45])] if len(a) <= k]
        val, v = psi(theta, cfg, 'G', k=k, Zc=Zc, Zd=0, nstart=ns, seed=int(1000 * cfg.eps * 1e4) + k, seeds=seeds)
        best = max(best, val)
    return best
mu = 2.855121
for eps in [float(a) for a in sys.argv[1:]]:
    rows = []
    for a0 in [0.24, 0.30, 0.375]:
        cfg = Cfg(eps=eps, mu=mu, a0=a0)
        probe = m + cfg.gain - mu * (2 * (m - eps) - 1)
        for th in [28.0, 33.0, 38.0, 44.0]:
            v = psiG_ideal(cfg, th)
            G0 = cfg.gain - max(v, probe) - 2 * eps * th
            rows.append((G0, a0, th, v))
            print(f"eps={eps:.2e} a0={a0} theta={th}: gain={cfg.gain:.6f} psiG={v:.6f} probe={probe:+.5f} G0={G0:+.6f}", flush=True)
    b = max(rows)
    print(f"BEST eps={eps:.2e}: G0={b[0]:+.6f} (a0={b[1]}, theta={b[2]})", flush=True)
