# Round 14: Program B point control for the N6 red file red285_cutbad.json: B's own enclosure of alpha . g at the single point
# x0 (e.g. x0 = 0, where alpha . g = alpha0 = -1e-6; Program A located its failure near 0.6538, where the cut crosses 0). If its UPPER end is < 0, no box containing x0 can pass B's test
# (the box enclosures contain the point value), so B's N6 run cannot PASS.  Usage: r14_redpointN6.py red.json x0
import sys, importlib.util
params, x0 = sys.argv[1], float(sys.argv[2])
sys.argv = ['certB_cut.py', params]
spec = importlib.util.spec_from_file_location('cbc', 'certB_cut.py'); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
X = B.X
a0_, a1_, a2_, a3_, a4_ = B.CUTS[0]
v = a0_ + a1_ * X(x0) + a2_ * B.ent_box(x0, x0) + a3_ * B.phi_box(x0, x0) + a4_ * B.R_box(x0, x0)
print(f"Program B {params} cut 0 at x0 = {x0!r}: alpha . g in [{v.a:.6e}, {v.b:.6e}] -> {'cannot PASS (red OK)' if v.b < 0 else 'no conclusion'}")
