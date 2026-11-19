# Program B red controls at eps = 2.7e-4 (with kernel K2): point controls with certB_jk.py's own lamG at the float maximiser.
import sys, importlib.util
from fractions import Fraction
sys.argv = ['x']
from gateA import Cfg, m
from jk_gate import sup
cases = [('red270_Psi_low.json', 35.0, [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]),
         ('red270_noK2.json', 35.0, [(0.22, 0.0), (0.0, 0.05)]),
         ('red270_theta32.json', 32.0, [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)])]
for params, th, pairs in cases:
    cfg = Cfg(eps=2.7e-4, mu=2.8, a0=0.21)
    v, (s, e, bb, Rb) = sup(cfg, th, pairs)
    sys.argv = ['certB_jk.py', params]
    spec = importlib.util.spec_from_file_location('cbj', 'certB_jk.py'); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
    X = B.X; s = max(s, B.SMIN); Rb = min(max(Rb, 0.0), B.RMX)
    val, t = B.lamG(X(s), X(e), X(max(bb, 0.0)), X(Rb))
    feas = t.a > 0 and t.b <= 1 and (X(e) - B.F45 * t - B.cv(Fraction(27, 400))).b < 0 and (X(e) - B.ent_box(s, min(s, 1.0))).b < 0
    print(f"{params}: float sup Lambda_G/e = {v:.6f} at (s-m, e, b, R) = ({s - m:+.2e}, {e:.5f}, {bb:.2e}, {Rb:.2e}); "
          f"Program B value in [{val.a:.4e}, {val.b:.4e}], feasible={feas} -> {'cannot PASS (red OK)' if val.a > 0 and feas else 'no conclusion'}", flush=True)
