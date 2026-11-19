# Round 15 prep (FLOAT): ideal per-law sup of l_G - theta beta over laws (k atoms, CIID couplings with Zc latent values) at the
# round-14 parameters (eps 2.9e-4, a0 0.19, mu 2.8, theta 31), against the certified relaxation's 0.014856. The difference is the
# slack that stronger local lemmas could recover at these parameters. Also the ideal sup at theta 23-27 (where the ideal optimum sits).
import sys; sys.argv = sys.argv[:1]
from gateA import *
from ideal_eps import psiG_ideal
for th in [31.0, 27.0, 23.0]:
    cfg = Cfg(eps=2.9e-4, mu=2.8, a0=0.19); probe = m + cfg.gain - 2.8 * (2 * (m - 2.9e-4) - 1)
    v = psiG_ideal(cfg, th)
    print(f"eps=2.9e-4 a0=0.19 mu=2.8 theta={th}: ideal psiG={v:.6f}  G0_ideal={cfg.gain - max(v, probe) - 2*2.9e-4*th:+.6f}", flush=True)
