# Round 14 step 1 (FLOAT): what binds at eps = 2.8e-4 with kernel K2. For each theta: sup of Lambda_G/e over D_G (grid + Nelder-Mead,
# jk_gate.sup), its location, and G0 = gain - max(psiG, probeT) - 2 eps theta. Also the N5 margin of red270_eps280.json.
import sys; sys.argv = sys.argv[:1]
from fractions import Fraction as Fr
from gateA import *
from jk_gate import sup
PAIRS = [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]
cfg0 = Cfg(eps=2.8e-4, mu=2.8, a0=0.21)
print(f"red270_eps280 (theta 35, Psi 0.01537): N5 = Gamma0 - Psi - 2 eps theta = {cfg0.gain - 0.01537 - 2*2.8e-4*35:+.6e}")
for a0 in [0.21]:
    cfg = Cfg(eps=2.8e-4, mu=2.8, a0=a0)
    probe = m + cfg.gain - 2.8 * (2 * (m - 2.8e-4) - 1)
    for th in [30.0, 32.0, 33.0, 34.0, 35.0, 36.0, 38.0]:
        v, arg = sup(cfg, th, PAIRS)
        reg = 'CORNER (large Rbar, bbar)' if arg[3] > 0.05 else 'contact region'
        print(f"a0={a0} theta={th}: psiG={v:.6f} at s-m={arg[0]-m:+.2e} e={arg[1]:.4f} bbar={arg[2]:.3e} Rbar={arg[3]:.3e} [{reg}]  G0={cfg.gain - max(v, probe) - 2*2.8e-4*th:+.6f}", flush=True)
