import numpy as np, itertools
from gateA import *
from relaxG2 import sup
eps = 2e-4
for a0, mu, th in itertools.product([0.24, 0.27, 0.30, 0.33], [2.80], [45.0, 50.0, 55.0, 60.0]):
    cfg = Cfg(eps=eps, mu=mu, a0=a0)
    v, arg = sup(cfg, th)
    probe = m + cfg.gain - mu * (2 * (m - eps) - 1)
    print(f"a0={a0} mu={mu} theta={th}: gain={cfg.gain:.6f} psiG={v:+.6f} probeT={probe:+.6f} G0={cfg.gain - max(v, probe) - 2*eps*th:+.6f} arg s-m={arg[0]-m:+.2e} e={arg[1]:.4f} b={arg[2]:.2e} R={arg[3]:.2e}", flush=True)
