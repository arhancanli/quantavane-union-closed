# Point red controls for Program B: if B's own evaluator gives a strictly positive LOWER bound at one feasible point, no box
# containing that point is ever accepted, so the bisection cannot PASS. Uses certB.py's functions unchanged.
import sys, json, importlib.util
from fractions import Fraction
def load(params):
    sys.argv = ['certB.py', params]
    spec = importlib.util.spec_from_file_location('cb_' + params.replace('.', '_'), 'certB.py'); mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod); return mod
# N3 point: the float maximiser region (s = s_min, e ~ 0.6094, bbar = 0, Rbar ~ 0.00323), Program A's failing box for Psi = 0.0140
for params, pt in [('red_Psi_low.json', (None, 0.6093838610045, 0.0, 0.0032293947226)),
                   ('red_Psi_0130.json', (None, 0.6093838610045, 0.0, 0.0032293947226)),
                   ('red_theta40.json', (None, 0.6317235248655, 0.0, 0.0062422642041)),
                   ('params.json', (None, 0.6093838610045, 0.0, 0.0032293947226))]:
    B = load(params); X = B.X
    s = B.SMIN; e, bb, Rb = pt[1:]
    s_, e_, b_, R_ = X(s), X(e), X(bb), X(Rb)
    v, t = B.lamG(s_, e_, b_, R_)
    feas = t.a > 0 and t.b <= 1 and (e_ - B.F45 * t - B.cv(Fraction(27, 400))).b < 0 and (e_ - B.ent_box(s, s)).b < 0
    print(f"{params}: N3 point value in [{v.a:.6e}, {v.b:.6e}], feasible={feas} -> {'cannot PASS (red OK)' if v.a > 0 and feas else 'no conclusion'}")
# N1 point: R(0.984) > 0.2148 ;  N2 point: minorant 0 with d = 0.01 at x = 0.9147
B = load('red_Rmax_low.json'); xx = 0.984
Rv = (B.X(B.ent_at(B.r.a * xx).a) + B.q * B.L * B.X(B.phi_at(xx).a) - B.cv(xx) * B.hr - B.r * B.X(B.ent_at(xx).b)) / (B.q * B.L)
print(f"red_Rmax_low.json: R(0.984) >= {Rv.a:.6f} vs Rmax 0.2148 -> {'cannot PASS (red OK)' if Rv.a > 0.2148 else 'no conclusion'}")
B = load('red_minor_up.json'); xx = 0.9147; a_, b_, c_, d_ = B.MINS[0]
hb = B.ent_at(xx); pb = B.phi_at(xx)
rhs = a_ * pb + b_ * hb + c_ * B.cv(xx) + d_
import math
lhs_hi = math.nextafter(math.sqrt(math.nextafter(xx * hb.b, math.inf)), math.inf)
print(f"red_minor_up.json: at x = 0.9147 rhs >= {rhs.a:.6f} > lhs <= {lhs_hi:.6f} -> {'cannot PASS (red OK)' if rhs.a > lhs_hi else 'no conclusion'}")
