import time, numpy as np
from gateA import *
cfg = Cfg(eps=2e-4)
phs = 1 - m / r
for theta in [30.0, 44.0, 60.0]:
    for side in ['G', 'T']:
        k, Zc = 4, 2
        Zd = k if side == 'T' else 0
        seeds = [seed_vec([0, r], [phs, 1 - phs], k, Zc, Zd, cfg, side),
                 seed_vec([0, r - 0.01], [1 - m/(r-0.01), m/(r-0.01)], k, Zc, Zd, cfg, side),
                 seed_vec([m], [1.0], k, Zc, Zd, cfg, side),
                 seed_vec([0, 0.801], [1 - 0.7715, 0.7715], k, Zc, Zd, cfg, side),
                 seed_vec([0, 1, 0.7], [0.38, 0.55, 0.07], k, Zc, Zd, cfg, side)]
        t1 = time.time()
        val, v = psi(theta, cfg, side, k=k, Zc=Zc, Zd=Zd, nstart=60, seed=1, seeds=seeds)
        print(f"theta={theta} side={side}: psi={val:.6f} ({time.time()-t1:.1f}s)", describe(v, k, Zc, Zd, cfg) if v is not None else None, flush=True)
