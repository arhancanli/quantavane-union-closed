# Round 14: turn the float cuts of r14_cuts.json into exact rational cuts for params['cuts'].
# alpha1..alpha4 rounded to denominators <= 10^7; alpha0 = smallest multiple of 1e-7 with alpha0 + min_x (rest) >= DELTA on a fine grid
# (2e6 uniform points plus 4000 log-spaced points at each end). Validity on all of [0,1] is then N6 (both programs).
# Usage: r14_ratcuts.py DELTA out.json
import sys, json, numpy as np
from fractions import Fraction as Fr
args = sys.argv[1:]; sys.argv = sys.argv[:1]
from gateA import h, phi, R
DELTA = float(args[0])
xs = np.unique(np.concatenate([np.linspace(0, 1, 2000001), np.logspace(-16, -1, 4000), 1 - np.logspace(-16, -1, 4000), [0.0, 1.0]]))
xs = xs[(xs >= 0) & (xs <= 1)]
G = np.stack([xs, h(xs), phi(xs), R(xs)], axis=1)
out = []
for a in json.load(open('r14_cuts.json')):
    rest = [Fr(v).limit_denominator(10**7) for v in a[1:]]
    mn = float(np.min(G @ np.array([float(v) for v in rest])))
    a0 = Fr(int(np.ceil((DELTA - mn) * 1e7)), 10**7)
    gap = float(a0) - (a[0])
    out.append([str(a0)] + [str(v) for v in rest])
    print(f"cut: float alpha0 {a[0]:+.3e} -> {float(a0):+.7f} (weakened by {gap:+.2e}); rest {[float(v) for v in rest]}")
json.dump(out, open(args[1], 'w'), indent=0)
