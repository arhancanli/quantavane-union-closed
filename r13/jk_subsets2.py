import sys, itertools
argv = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import *
from jk_gate import sup
SU = [(0.22, 0.0), (0.0, 0.05)]
subsets = {'SU+K1': SU + [(0.20, 0.0078)], 'SU+K2': SU + [(0.15, 0.0255)], 'SU+K1+K2': SU + [(0.20, 0.0078), (0.15, 0.0255)],
           'SU+K1+K2+K3': SU + [(0.20, 0.0078), (0.15, 0.0255), (0.05, 0.0535)]}
for eps in [float(a) for a in argv]:
    for name, pairs in subsets.items():
        best = None
        for a0, th in itertools.product([0.18, 0.21], [32.0, 36.0, 40.0, 44.0]):
            cfg = Cfg(eps=eps, mu=2.8, a0=a0); v, arg = sup(cfg, th, pairs)
            probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1); G0 = cfg.gain - max(v, probe) - 2 * eps * th
            if best is None or G0 > best[0]: best = (G0, a0, th, v)
        print(f"eps={eps:.2e} {name:14s}: G0={best[0]:+.6f} (a0={best[1]}, theta={best[2]}, psiG={best[3]:.6f})", flush=True)
