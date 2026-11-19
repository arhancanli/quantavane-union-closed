# Round 14 (FLOAT): sup of Lambda_G/e on D_G with the RATIONAL cuts of a json file (as certified), for (eps, a0, mu, theta) lists.
# Usage: r14_psi_rat.py cuts.json eps a0 mu theta [theta ...]
import sys, json, numpy as np
from fractions import Fraction as Fr
args = sys.argv[1:]
from r14_lib import *
cuts = [np.array([float(Fr(v)) for v in row]) for row in json.load(open(args[0]))]
eps, a0, mu = float(args[1]), float(args[2]), float(args[3])
cfg = Cfg(eps=eps, mu=mu, a0=a0); probe = m + cfg.gain - mu * (2 * (m - eps) - 1)
for th in [float(v) for v in args[4:]]:
    v, arg = sup(cfg, th, cuts); s, e, bb, Rb = arg
    print(f"rational cuts: eps={eps:.3e} a0={a0} mu={mu} theta={th}: psiG={v:.7f} G0={cfg.gain - max(v, probe) - 2*eps*th:+.6f} gain={cfg.gain:.7f} "
          f"at s-m={s-m:+.2e} e={e:.5f} b={bb:.3e} R={Rb:.3e}", flush=True)
