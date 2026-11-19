# Robust gate scan at a given eps: psiG, psiT over a theta grid, several (k, Zc), seeds; prints G0(theta).
import sys, time, json, numpy as np
from gateA import *
eps = float(sys.argv[1]); mu = float(sys.argv[2]) if len(sys.argv) > 2 else None; Tmode = sys.argv[3] if len(sys.argv) > 3 else 'Phi'
cfg = Cfg(eps=eps, mu=mu, Tmode=Tmode)
print(f"eps={eps} mu={cfg.mu:.4f} gain={cfg.gain:.6f} B'={cfg.Bp:.5f} Tmode={Tmode}", flush=True)
phs = 1 - m / r
def seeds_for(k, Zc, Zd, side):
    S = [([0, r], [phs, 1 - phs]), ([0, 0.66], [0.065, 0.935]), ([0, r - 0.02], [1 - m/(r-0.02), m/(r-0.02)]),
         ([m], [1.0]), ([0, 0.801], [1 - 0.7715, 0.7715]), ([0, 1, 0.7], [0.38, 0.55, 0.07]), ([0, 1, 0.7], [0.3, 0.3, 0.4]),
         ([0, 0.5, 0.9], [0.1, 0.3, 0.6]), ([0.2, 0.7, 0.95], [0.2, 0.5, 0.3])]
    return [seed_vec(a, b, k, Zc, Zd, cfg, side) for a, b in S if len(a) <= k]
out = {}
for theta in [25.0, 35.0, 45.0, 60.0, 80.0]:
    res = {}
    for side in ['G', 'T']:
        best = (-np.inf, None, None)
        for (k, Zc, ns) in [(3, 2, 80), (4, 2, 80), (4, 3, 60), (5, 3, 50)]:
            Zd = k if (side == 'T' and Tmode == 'Phi') else 0
            val, v = psi(theta, cfg, side, k=k, Zc=Zc, Zd=Zd, nstart=ns, seed=int(theta) + k + 7 * Zc, seeds=seeds_for(k, Zc, Zd, side))
            if val > best[0]: best = (val, describe(v, k, Zc, Zd, cfg), (k, Zc))
        res[side] = best
        print(f"  theta={theta:5.1f} {side}: psi={best[0]:+.6f} via k,Zc={best[2]} {best[1]}", flush=True)
    G0 = cfg.gain - max(res['G'][0], res['T'][0]) - 2 * eps * theta
    print(f"theta={theta:5.1f}  G0={G0:+.6f}", flush=True)
