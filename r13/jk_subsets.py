import sys, itertools
sys_argv = sys.argv[:]; sys.argv = sys.argv[:1]
from gateA import *
from jk_gate import sup, PAIRS
SU = [(0.22, 0.0), (0.0, 0.05)]
subsets = {'SU': SU, 'SU+(0.20,.0078)': SU + [(0.20, 0.0078)], 'SU+(0.18,.0150)': SU + [(0.18, 0.0150)],
           'SU+(0.15,.0255)': SU + [(0.15, 0.0255)], 'SU+(.18,.015)+(.10,.0418)': SU + [(0.18, 0.0150), (0.10, 0.0418)],
           'SU+(.20,.0078)+(.15,.0255)+(.05,.0535)': SU + [(0.20, 0.0078), (0.15, 0.0255), (0.05, 0.0535)], 'full': PAIRS}
eps = float(sys_argv[1])
for name, pairs in subsets.items():
    best = None
    for th in [28.0, 30.0, 32.0]:
        cfg = Cfg(eps=eps, mu=2.8, a0=0.18); v, arg = sup(cfg, th, pairs)
        probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1); G0 = cfg.gain - max(v, probe) - 2 * eps * th
        if best is None or G0 > best[0]: best = (G0, th, v)
    print(f"eps={eps:.2e} {name:42s}: G0={best[0]:+.6f} (theta={best[1]}, psiG={best[2]:.6f})", flush=True)
