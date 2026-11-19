# Round 15 (FLOAT): sup of Lambda_G/e with the RATIONAL cuts (round 14, r14_cuts_rat.json) and 6-cuts (r15_cuts6_rat.json) as certified
# by certA_v2/certB_v2: w_lo rule in EF_lo. Usage: r15_psi6.py eps a0 theta [theta ...]
import sys, json, numpy as np
from fractions import Fraction as Fr
args = sys.argv[1:]
src = open('r15_jensen.py').read().split("eps, a0 = float(args[0])")[0].replace("args = sys.argv[1:]; sys.argv = sys.argv[:1]", "sys.argv = sys.argv[:1]")
ns = {'__name__': 'r15lib'}; exec(src, ns)
cuts = ns['cuts'] + [np.array([float(Fr(v)) for v in row]) for row in json.load(open('r15_cuts6_rat.json'))]
Cfg, m, sup = ns['Cfg'], ns['m'], ns['sup']
eps, a0 = float(args[0]), float(args[1]); cfg = Cfg(eps=eps, mu=2.8, a0=a0); probe = m + cfg.gain - 2.8 * (2 * (m - eps) - 1)
for th in [float(v) for v in args[2:]]:
    v, arg = sup(cfg, th, cuts)
    print(f"rational cuts+cuts6: eps={eps:.3e} a0={a0} theta={th}: psiG={v:.7f} G0={cfg.gain - max(v, probe) - 2*eps*th:+.6f} at s-m={arg[0]-m:+.2e} e={arg[1]:.5f} b={arg[2]:.3e} R={arg[3]:.3e}", flush=True)
