# Relaxation (provable ingredients) margin at larger eps: G0 = gain - max(psiG, probeT) - 2 eps theta (FLOAT).
import sys, itertools
from gateA import *
from relaxG2 import sup
for eps in [float(a) for a in sys.argv[1:]]:
    best = None
    for a0, mu, th in itertools.product([0.21, 0.24, 0.27], [2.8], [44.0, 47.0, 50.0]):
        cfg = Cfg(eps=eps, mu=mu, a0=a0)
        v, arg = sup(cfg, th)
        probe = m + cfg.gain - mu * (2 * (m - eps) - 1)
        G0 = cfg.gain - max(v, probe) - 2 * eps * th
        print(f"eps={eps:.3e} a0={a0} mu={mu} theta={th}: psiG={v:.6f} probe={probe:+.5f} G0={G0:+.6f}", flush=True)
        if best is None or G0 > best[0]: best = (G0, a0, mu, th, v)
    print(f"BEST eps={eps:.3e}: G0={best[0]:+.6f} a0={best[1]} mu={best[2]} theta={best[3]} psiG={best[4]:.6f}", flush=True)
