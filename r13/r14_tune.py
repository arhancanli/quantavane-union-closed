# Round 14 (FLOAT): parameter scan with the moment-hull cuts (cutting plane per point; cuts shared through r14_cuts.json).
# Usage: r14_tune.py eps a0,a0.. mu,mu.. theta,theta..
import sys, json, itertools, numpy as np
args = sys.argv[1:]
from r14_lib import *
eps = float(args[0]); A0s, MUs, THs = ([float(v) for v in a.split(',')] for a in args[1:4])
cuts = [np.array(a) for a in json.load(open('r14_cuts.json'))]
for a0, mu, th in itertools.product(A0s, MUs, THs):
    cfg = Cfg(eps=eps, mu=mu, a0=a0); probe = m + cfg.gain - mu * (2 * (m - eps) - 1)
    for it in range(20):
        v, arg = sup(cfg, th, cuts); s, e, bb, Rb = arg; t = (s - e / p - bb) / c; sv, a = sep(np.array([1, s, e, t, Rb]))
        if sv > -1e-7: break
        cuts.append(a); json.dump([list(map(float, c_)) for c_ in cuts], open('r14_cuts.json', 'w'))
    print(f"eps={eps:.3e} a0={a0} mu={mu} theta={th}: psiG={v:.6f} G0={cfg.gain - max(v, probe) - 2*eps*th:+.6f} gain={cfg.gain:.6f} at s-m={s-m:+.2e} e={e:.4f} b={bb:.2e} R={Rb:.2e} ncuts={len(cuts)}", flush=True)
