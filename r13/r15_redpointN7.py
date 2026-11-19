# Round 15: Program B point control for the N7 red file red300_cut6bad.json: B's enclosure of alpha . (1, x, h, phi, R, v^2) for 6-cut 0
# at a point x0 (x0 = 0 gives alpha0 = -1e-6). An upper end < 0 means no box containing x0 passes B's N7 test. Usage: r15_redpointN7.py red.json x0
import sys, importlib.util
params, x0 = sys.argv[1], float(sys.argv[2])
sys.argv = ['certB_v2.py', params]
spec = importlib.util.spec_from_file_location('cbv', 'certB_v2.py'); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
X = B.X
a0_, a1_, a2_, a3_, a4_, a5_ = B.CUTS6[0]
hb, pb = B.ent_box(x0, x0), B.phi_box(x0, x0)
v = a0_ + a1_ * X(x0) + a2_ * hb + a3_ * pb + a4_ * B.R_box(x0, x0) + a5_ * B.sqr(hb - B.hr * pb)
print(f"Program B {params} 6-cut 0 at x0 = {x0!r}: alpha . g6 in [{v.a:.6e}, {v.b:.6e}] -> {'cannot PASS (red OK)' if v.b < 0 else 'no conclusion'}")
