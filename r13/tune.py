import numpy as np, itertools
from gateA import *
import relaxG2
from relaxG2 import sup
eps = 2e-4
rows = []
for a0, mu, th in itertools.product([0.33, 0.375, 0.42], [2.78, 2.855, 3.0], [55.0, 60.0, 65.0]):
    cfg = Cfg(eps=eps, mu=mu, a0=a0)
    v, arg = sup(cfg, th)
    probe = m + cfg.gain - mu * (2 * (m - eps) - 1)   # T-side probe value (per unit entropy); must be <= Psi
    G0 = cfg.gain - max(v, probe) - 2 * eps * th
    rows.append((G0, a0, mu, th, v, probe))
    print(f"a0={a0} mu={mu} theta={th}: gain={cfg.gain:.6f} psiG={v:+.6f} probeT={probe:+.6f} G0={G0:+.6f}", flush=True)
rows.sort(reverse=True); print("BEST", rows[:3])
