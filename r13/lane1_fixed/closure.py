# (1) Re-evaluate every edge [SC] condition with its threshold enlarged so the edge region swallows the
#     float sliver (factor 1+1e-9 on small thresholds). (2) Check numerically that each contact function and its
#     gradient contain 0 at the r-ball (the analytic identities the Hessian argument relies on).
import sys; sys.argv = ['x', 'SC']
from progA_core import *
import progA_checks as P
E = Q(1, 10**9) + 1
ell, Phiu, hq, ap, kappa4 = P.ell, P.Phiu, P.hq, P.ap, P.kappa4
def hp(x): return h_point(x)
res = []
def ok(n, c): res.append((n, bool(c))); print(('PASS ' if c else 'FAIL ') + n)
ok("K1 small at 3/20*(1+e)", (2 * ell(Q(3, 20) * E) / L).sqrt() > rho0)
ok("K1/K3b/K4a near 1 at a=(1/100)(1+e)", hp(Q(1, 100) * E) < qL ** 2 / (r * r * L))
l = ell(Q(1, 200) * E)
ok("K2 small at (1/200)(1+e)", (hr / ap - r) * l + (1 - hr * hr / ap) * (2 * l / L).sqrt() - rho0 > 0)
a2 = Q(1, 200) * E; kp = r / qL - (1 - r + hr / ap - a2)
ok("K2 near 1 at (1/200)(1+e)", hp(a2).sqrt() < (1 - hr * hr / ap) / (kp * sqrtL))
l = ell(Q(1, 500) * E); Pu = Phiu(l)
ok("K3 small at (1/500)(1+e)", 1 + l / ap - c * Pu > 0 and 1 / ap - c * (1 + tau) / (L * Pu) > 0)
ok("K3 near 1 at (1/100)(1+e)", 1 - Q(1, 100) * E - c * (hp(Q(2, 100) * E) / L).sqrt() > 0)
ok("K3 1-b small at (1/100)(1+e)", 1 - Q(1, 100) * E - hp(Q(1, 100) * E) / ap > 0)
ok("K4a at (1/100)(1+e)", hp(Q(1, 100) * E) < Q(27, 400))
ok("K4b at (1/5000)(1+e)", hp(Q(1, 5000) * E) < Q(21, 10000))
ok("W1 small at (1/10)(1+e)", Q(2, 5) * (2 * ell(Q(1, 10) * E) / L).sqrt() > 1)
ok("W1 near 1 at (1/25)(1+e)", hp(Q(1, 25) * E) <= Q(4, 25) / L)
xx = Q(1, 1000) * E; l = ell(xx); phk = phi(Jet(xx, n=1, order=0)).v
ok("W3 small at (1/1000)(1+e)", xx * (l + 1) ** 2 / (2 * l / L).sqrt() <= Q(9, 100) and hr * hr * phk <= Q(9, 100))
BETA = Q(4, 10**5) * E
ok("K11 at beta(1+e) (closes corner and W2/X1 slivers)", kappa4 * (2 - (1 + tau) * BETA) * ell(2 * BETA) / L >= 2 * L)
D1 = Q(2, 10**7) * E
ok("X1 at delta1(1+e)", hp(D1) <= (P.MU * kappa4) ** 2 / L)
ok("X0 at delta0(1+e)", kappa4 * (2 * ell(Q(1, 10**10) * E) / L).sqrt() >= Q(2, 5))
Phi1 = (hp(2 * D1) / L).sqrt(); la = ell(P.YA); hb = la / Phiu(la); yb = (L / (2 * la)).sqrt()
val = (1 - D1 - Phi1 / (5 * ap)) * hb - L * Phi1 - Q(1, 5) * (1 + hp(D1) / ap) - Q(1, 5) * Phi1 * yb
ok("Z3 at delta1(1+e)", val > 0)
l1 = ell(Q(1, 1000) * E); rho1 = 1 - rho0 * (L / (2 * l1)).sqrt(); beta1 = Q(1, 1000) * E + hp(Q(1, 1000) * E) / ap
ok("Z4a at y1(1+e)", Q(11, 50) * rho1 - beta1 / 5 > 0)
print("Z2 B&B with delta1(1+e): rerun below; Z4b with delta0'(1+e): rerun below")
print("---- contact identities at the r-ball (value and gradient must contain 0)")
X = Jet.var(r, 0, 2, 1); Y = Jet.var(r, 1, 2, 1); x1 = Jet.var(r, 0, 1, 1)
for nm, G in [("G1", P.G1), ("G2", P.G2), ("G3", P.G3), ("G4", P.G4)]:
    g = G(X, Y); ok(f"contact {nm}: value and gradient contain 0", g.v.contains(0) and all(d.contains(0) for d in g.d))
S1b = lambda x: (x - r) * jh(x) - vfun(x) * (hr / p) + R(x)
for nm, f in [("R", R), ("S1b", S1b), ("b", bfun)]:
    g = f(x1); ok(f"contact {nm}: value and derivative contain 0", g.v.contains(0) and g.d[0].contains(0))
print("ALL EDGE CLOSURES PASS" if all(c for _, c in res) else "SOME CLOSURE FAILED")
