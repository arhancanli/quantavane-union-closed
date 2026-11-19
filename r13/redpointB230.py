# Program B red controls at eps = 2.3e-4: point controls (B's own lamG at the float maximiser) for N3.
import sys, importlib.util
from fractions import Fraction
sys.argv = ['x']
from gateA import Cfg, m
from relaxG2 import sup
for params, th in [('red230_Psi_low.json', 42.0), ('red230_theta38.json', 38.0)]:
    cfg = Cfg(eps=2.3e-4, mu=2.8, a0=0.21)
    v, (s, e, bb, Rb) = sup(cfg, th)
    sys.argv = ['certB.py', params]
    spec = importlib.util.spec_from_file_location('cb', 'certB.py'); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
    X = B.X; s = max(s, B.SMIN)
    val, t = B.lamG(X(s), X(e), X(max(bb, 0.0)), X(max(Rb, 0.0)))
    feas = t.a > 0 and t.b <= 1 and (X(e) - B.F45 * t - B.cv(Fraction(27, 400))).b < 0 and (X(e) - B.ent_box(s, s)).b < 0
    print(f"{params}: float sup Lambda_G/e = {v:.6f} at (s-m, e, b, R) = ({s - m:+.2e}, {e:.5f}, {bb:.2e}, {Rb:.2e}); "
          f"Program B value in [{val.a:.4e}, {val.b:.4e}], feasible={feas} -> {'cannot PASS (red OK)' if val.a > 0 and feas else 'no conclusion'}", flush=True)
