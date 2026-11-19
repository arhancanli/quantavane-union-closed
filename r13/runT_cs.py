# Semi-ideal T side: exact per-law U, V (best CIID coupling), g, but Phi replaced by the provable Phi_CS = (E sqrt(X h(X)))^2.
import sys, numpy as np
import gateA
from gateA import *
def lossT_cs(st, cfg, x=None, w=None):
    return m * st['e'] - st['PhiCS'] - cfg.mu * (st['U'] - st['e']) + cfg.gain * st['e']
_orig = gateA.law_stats
def law_stats_cs(x, M, P, cfg):
    st = _orig(x, M, None, cfg)
    st['PhiCS'] = float(st['w'] @ np.sqrt(np.maximum(x * h(x), 0))) ** 2
    return st
gateA.law_stats = law_stats_cs
gateA.lossT = lossT_cs
eps = float(sys.argv[1]); cfg = Cfg(eps=eps, Tmode='CS')
phs = 1 - m / r
for theta in [float(a) for a in sys.argv[2:]]:
    best = (-np.inf, None)
    for (k, Zc, ns) in [(3, 2, 80), (4, 2, 60), (4, 3, 50), (5, 3, 40)]:
        seeds = [seed_vec(a, b, k, Zc, 0, cfg, 'G') for a, b in [([0, 0.801], [1-0.7715, 0.7715]), ([0, 1, 0.7], [0.38, 0.55, 0.07]), ([0, 1, 0.7], [0.3, 0.3, 0.4]), ([0, 0.5, 0.9], [0.2, 0.3, 0.5]), ([0.05, 0.9, 1.0], [0.3, 0.3, 0.4])] if len(a) <= k]
        val, v = psi(theta, cfg, 'T', k=k, Zc=Zc, Zd=0, nstart=ns, seed=11 + k + Zc, seeds=seeds)
        if val > best[0]: best = (val, describe(v, k, Zc, 0, cfg))
    print(f"theta={theta:6.1f} psiT_CS={best[0]:+.6f} {best[1]}", flush=True)
