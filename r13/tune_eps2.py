import itertools
from gateA import *
import sys as _s; _argv = _s.argv[1:]; _s.argv = _s.argv[:1]
from relaxG2 import sup
for eps in [float(a) for a in _argv]:
    best = None
    for a0, th in itertools.product([0.18, 0.21], [40.0, 42.0, 44.0]):
        cfg = Cfg(eps=eps, mu=2.8, a0=a0); v, arg = sup(cfg, th)
        probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1); G0 = cfg.gain - max(v, probe) - 2 * eps * th
        print(f"eps={eps:.3e} a0={a0} theta={th}: gain={cfg.gain:.6f} psiG={v:.6f} probe={probe:+.5f} G0={G0:+.6f}", flush=True)
        if best is None or G0 > best[0]: best = (G0, a0, th, v, cfg.gain)
    print(f"BEST eps={eps:.3e}: G0={best[0]:+.6f} a0={best[1]} theta={best[2]} psiG={best[3]:.6f} gain={best[4]:.6f}", flush=True)
