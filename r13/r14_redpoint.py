# Round 14 point controls for N3 red runs, in BOTH programs (certA_cut.py's LamG in Arb, certB_cut.py's lamG in binary64 outward).
# At an exact point P = (s, e, bbar, Rbar) of D_G: if every domain test and every cut of the params file holds at P (so no box
# containing P is dropped) and Lambda_G - Psi e has a strictly positive lower bound at P (so no box containing P is accepted),
# the N3 run on that params file cannot PASS. (Domain-box tests z <= e <= L, 0 <= bbar <= 1, Rbar >= 0, s <= 1 added after the
# round-14 referee, NOTE 6.)  Usage: r14_redpoint.py params.json s e bbar Rbar   (decimal strings, exact)
import sys, importlib.util
from fractions import Fraction as Fr
params, pt = sys.argv[1], [Fr(v) for v in sys.argv[2:6]]
sys.argv = ['certA_cut.py', params]
src = open('certA_cut.py').read().split("if __name__ == '__main__':")[0]
ns = {'__name__': 'x'}; exec(src, ns)
A, LamG, lo, hi, CUTS, SMIN = ns['A'], ns['LamG'], ns['lo'], ns['hi'], ns['CUTS'], ns['SMIN']
s, e, bb, Rb = pt; s = max(s, SMIN)
val, t = LamG(A(s), A(e), A(bb), A(Rb))
feasA = (lo(t) > 0 and hi(t) <= 1 and hi(A(e) - A(Fr(4, 5)) * t - A(Fr(27, 400))) <= 0 and hi(A(e) - ns['h_pt'](A(s))) <= 0
         and A(Rb) <= A(ns['RMAX']) and Rb >= 0 and 0 <= bb <= 1 and s <= 1 and ns['Z'] <= e <= ns['LHI'] and all(lo(A(a0) + A(a1) * A(s) + A(a2) * A(e) + A(a3) * t + A(a4) * A(Rb)) >= 0 for (a0, a1, a2, a3, a4) in CUTS))
okA = feasA and lo(val) > 0
print(f"Program A {params} at (s_min={float(s):.12f}, {float(e)}, {float(bb)}, {float(Rb)}): Lambda_G - Psi e in {val.str(8)}, "
      f"in domain and all {len(CUTS)} cuts hold: {feasA} -> {'cannot PASS (red OK)' if okA else 'no conclusion'}")
sys.argv = ['certB_cut.py', params]
spec = importlib.util.spec_from_file_location('cbc', 'certB_cut.py'); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
X, cv = B.X, B.cv
sB = max(float(s), B.SMIN)   # SMIN (float) <= m - eps; s_min of Program A is also <= m - eps
xs, xe, xb, xR = X(sB), cv(e), cv(bb), cv(Rb)
v, tB = B.lamG(xs, xe, xb, xR)
feasB = (B.ZG <= float(e) <= B.L.b and 0 <= bb <= 1 and Rb >= 0 and sB <= 1 and tB.a > 0 and tB.b <= 1 and (xe - B.F45 * tB - cv(Fr(27, 400))).b <= 0 and (xe - B.ent_box(sB, sB)).b <= 0 and xR.b <= B.RMX
         and all((a0 + a1 * xs + a2 * xe + a3 * tB + a4 * xR).a >= 0 for (a0, a1, a2, a3, a4) in B.CUTS))
okB = feasB and v.a > 0
print(f"Program B {params} at (s={sB!r}, same e, bbar, Rbar): Lambda_G - Psi e in [{v.a:.6e}, {v.b:.6e}], in domain and all "
      f"{len(B.CUTS)} cuts hold: {feasB} -> {'cannot PASS (red OK)' if okB else 'no conclusion'}")
