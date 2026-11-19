# Round 15: rational 6-component cuts alpha . (1, x, h, phi, R, v^2) from r15_cuts6.json (only those not already in r14_cuts_rat.json,
# i.e. the new ones; the first 5 entries of r15_cuts6.json are the round-14 cuts). alpha1..5 to denominators <= 10^7, alpha0 raised to a
# multiple of 1e-7 leaving margin >= DELTA on a 2,008,000-point grid. Validity on [0,1] is N7. Usage: r15_ratcuts6.py DELTA out.json
import sys, json, numpy as np
from fractions import Fraction as Fr
args = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import h, phi, R, vfun
DELTA = float(args[0])
xs = np.unique(np.concatenate([np.linspace(0, 1, 2000001), np.logspace(-16, -1, 4000), 1 - np.logspace(-16, -1, 4000), [0.0, 1.0]]))
xs = xs[(xs >= 0) & (xs <= 1)]
G = np.stack([xs, h(xs), phi(xs), R(xs), vfun(xs) ** 2], axis=1)
out = []
for a in json.load(open('r15_cuts6.json'))[5:]:
    rest = [Fr(v).limit_denominator(10**7) for v in a[1:]]
    mn = float(np.min(G @ np.array([float(v) for v in rest])))
    a0 = Fr(int(np.ceil((DELTA - mn) * 1e7)), 10**7)
    out.append([str(a0)] + [str(v) for v in rest])
    print(f"cut6: alpha0 {a[0]:+.3e} -> {float(a0):+.7f}; alpha5 = {float(rest[-1]):+.4f}")
json.dump(out, open(args[1], 'w'), indent=0)
