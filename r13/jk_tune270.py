import sys, itertools
argv = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import *
from jk_gate import sup
PAIRS = [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]
for eps in [float(a) for a in argv]:
    rows = []
    for a0, th in itertools.product([0.19, 0.21, 0.23], [32.0, 34.0, 36.0, 38.0]):
        cfg = Cfg(eps=eps, mu=2.8, a0=a0); v, arg = sup(cfg, th, PAIRS)
        probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1); G0 = cfg.gain - max(v, probe) - 2 * eps * th
        rows.append((G0, a0, th, v, cfg.gain)); print(f"eps={eps:.2e} a0={a0} theta={th}: psiG={v:.6f} gain={cfg.gain:.6f} G0={G0:+.6f}", flush=True)
    b = max(rows); print(f"BEST eps={eps:.2e}: G0={b[0]:+.6f} a0={b[1]} theta={b[2]} psiG={b[3]:.6f} gain={b[4]:.6f}", flush=True)
