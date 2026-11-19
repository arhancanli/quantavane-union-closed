# Program A point control for red270_Psi_low.json (N3): certA_jk.py's own LamG at the float maximiser, a feasible point
# (s = s_min, e = 0.60635, bbar = 0, Rbar = 0.00181: t = 0.939, e <= (4/5)t + 27/400, e <= h(s), Rbar <= R_max). A strictly positive
# lower bound there means no box containing the point can be accepted, so the N3 run cannot PASS.
import sys
sys.argv = ['certA_jk.py', 'red270_Psi_low.json']
src = open('certA_jk.py').read().split("if __name__ == '__main__':")[0]
ns = {'__name__': 'x'}; exec(src, ns)
A, LamG, Fr = ns['A'], ns['LamG'], ns['Fr']
pt = (Fr(ns['SMIN']), Fr('0.60635'), Fr(0), Fr('0.00181'))
val, t = LamG(*(A(v) for v in pt))
print(f"Program A point control red270_Psi_low at (s_min, 0.60635, 0, 0.00181): Lambda_G - Psi e in {val.str(8)}, t in {t.str(8)} -> {'cannot PASS (red OK)' if val.lower() > 0 else 'no conclusion'}")
