# Program B point control for red270_theta32.json (N3): search a small grid near the large-Rbar corner for a point that is feasible in
# certB_jk.py's own arithmetic (its search domain starts at s = SMIN, 3e-16 below m - eps) with a strictly positive lower bound.
import sys, importlib.util, itertools
from fractions import Fraction
sys.argv = ['certB_jk.py', 'red270_theta32.json']
spec = importlib.util.spec_from_file_location('cbj', 'certB_jk.py'); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
X = B.X; best = None
for s, e, bb, Rb in itertools.product([B.SMIN, B.SMIN + 1e-6, 0.6175, 0.6177], [0.385 + 1e-9, 0.3851, 0.3855, 0.387, 0.39],
                                      [0.40, 0.45, 0.48, 0.50, 0.51, 0.515], [0.2, 0.21, 0.2149]):
    val, t = B.lamG(X(s), X(e), X(bb), X(Rb))
    feas = t.a > 0 and t.b <= 1 and (X(e) - B.F45 * t - B.cv(Fraction(27, 400))).b < 0 and (X(e) - B.ent_box(s, s)).b < 0 and e >= B.ZG and Rb <= B.RMX
    if feas and (best is None or val.a > best[0]): best = (val.a, (s, e, bb, Rb))
v, pt = best
print(f"red270_theta32.json: best feasible point (s,e,b,R) = {pt} (s = SMIN is in Program B's search domain): Program B lower value {v:.4e} -> {'cannot PASS (red OK)' if v > 0 else 'no conclusion'}")
