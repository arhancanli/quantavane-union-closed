# Quantavane Program A: all checks of EDGE_LEMMAS.md with Arb balls. Usage:
#   python progA_checks.py <group> [red]      groups: SC S6 S5 K W Z C H G1 G2 G3 G4 all
# 'red' runs the deliberately broken variant of the group; it must FAIL.
import sys, time, math
from progA_core import *

RED = 'red' in sys.argv[2:]
GROUP = sys.argv[1] if len(sys.argv) > 1 else 'SC'
LOG = []
def out(s):
    print(s, flush=True); LOG.append(s)

def ok(name, cond, val=None):
    out(f"{'PASS' if cond else 'FAIL'} {name}" + ("" if val is None else f"  [{val}]"))
    if not cond:
        global FAILED; FAILED = True
FAILED = False

def ell(x): return (1 / x).log()
def Phiu(l): return ((1 + tau) * (2 * l + 1 - (1 + tau).log()) / L).sqrt()
def hq(n, d=1): return h_point(Q(n, d))

kappa4 = L * (q - Q(4, 25))
ap = -p                      # |p|
E_ = 1 + Q(1, 10**9)
D0P, D1P, BETA, Z1, YA, Y1 = Q(1, 10**10) * E_, Q(2, 10**7) * E_, Q(4, 10**5), Q(5, 10**5), Q(1, 10**6), Q(1, 1000)
MU = Q(32, 1000)

# ---------------------------------------------------------------- paper constants (Sections 2-7)
eps = Q(92, 10**6) if not (RED and GROUP == 'S6') else Q(93, 10**6)
z = Q(77, 200); lam = Q(1232, 1000); Gam = Q(63); kF = Q(1, 20); kE = Q(11, 50)
Db = 11 / (1 - 1 / (44 * qL))
gam = hr / r; alpha = 1 - gam / p; t0 = r - m; bm = -hr * gam / (p * alpha); Dc = -hr / (p * alpha) - t0
es = gam * m; tmin = (z - Q(27, 400)) / Q(4, 5); Cb = Db / tmin ** 2; a_s = m / r
kap = Q(39, 100); Ccov = Q(13, 25); a0 = Q(3, 8); bY = Q(2, 5)

def scalar_constants(eps):
    k = {}
    ellc = eps / m
    k0 = (1 / a_s + eps / (c * tmin * a_s)) / es
    k1 = k0 + 1 / (ap * c * tmin * a_s)
    k2 = Cb / (c * tmin * a_s)
    BE = (alpha / (1 - ellc)) ** 2 / (2 * kE)
    Ds = arb.union(L - es, es - z) if False else None
    Ds = (L - es) if (L - es) > (es - z) else (es - z)
    assert ((L - es) > (es - z)) or ((es - z) > (L - es))
    astar = m * h_point(r * r) / (r * hr)
    d0 = 2 * m - astar + 2 * qL * a_s / (p * c)
    D0 = lam - astar; rD = -d0 / D0
    beta = 1 - Dc / kap; b0 = (beta * a0.log()).exp(); W = (kap - Dc) / (1 - b0)
    eta = (kap - Dc) / Dc * (b0 - a0) / (1 - b0)
    kc = t0 - bm * eps / es; kh = eta * kc; A0 = W / 2; Bc = W * (a0 - HALF)
    Cm = bm * Cb + (1 + Ccov) / (2 * kF * tmin)
    Btot = eta * (kc + Dc) + Bc + (kh + m) * rD
    k.update(locals())
    return k

def M_of_v(K, v):
    # u = 2 eps v^2, sqrt(u) = sqrt(2 eps) v ; v is an arb ball inside [0,1]
    eps_ = K['eps']; u = 2 * eps_ * v * v; su = (2 * eps_).sqrt() * v
    C0 = (gam * (eps_ - u / 2) / z + gam * eps_ / es + gam * Cb * u) / (1 - K['ellc'])
    BE, k0, k1, k2 = K['BE'], K['k0'], K['k1'], K['k2']
    Ebar = C0 + BE * k1 * u / 2 + BE.sqrt() * su * (k0 + k2 * u + k1 * C0 + BE * k1 * k1 * u / 4).sqrt()
    n = (1 + Ebar) / es
    Kf = k0 + k1 * Ebar + k2 * u
    A2 = gam * K['Ds'] / (1 - bY) * ((eps_ - u / 2) / z + eps_ * n) + u / (1 - bY) * (alpha ** 2 / (8 * bY * kE * tmin) + gam * (es - z) * Cb)
    B0 = qL / (p * p * c * c) * A2 + 2 * eps_ + qL * eps_ ** 2 / (c * c * z)
    rho0_ = (1 - K['astar'] + B0 + Gam * (2 * eps_ - u)) / K['D0']
    A0, Bt = K['A0'], K['Btot']
    Qf = gam * (A0 + Bt) * (eps_ - u / 2) / z + gam * (A0 - Bt) * eps_ * n + (A0 + Bt) * alpha * su * (Kf / (2 * kE)).sqrt() + (A0 - Bt) * gam * Cb * u
    return K['kh'] - (K['kh'] + m) * rho0_ - Qf - K['Cm'] * u


UNDER = [0]
def exrad(x0, x1, xm):
    d1 = arb(xm) - arb(x0); d2 = arb(x1) - arb(xm)   # exact at prec 128 for doubles in [1e-10, 1]
    rad = d1 if d1 >= d2 else d2
    if rad > arb((x1 - x0) / 2): UNDER[0] += 1
    return rad
# ---------------------------------------------------------------- generic branch and bound
def iv(a, b): return arb.union(arb(a), arb(b))

def split_point(a, b):
    if a > 0 and b / a > 4: return math.sqrt(a * b)
    if b < 1 and (1 - a) / (1 - b) > 4: return 1 - math.sqrt((1 - a) * (1 - b))
    return (a + b) / 2

def bb1(name, f, a, b, contact=None, maxn=2_000_000, minw=1e-15):
    # f(Jet) -> Jet ; certify f > 0 on [a,b] (contact interval accepted separately)
    st = [(a, b)]; n = 0; mn = None; t = time.time()
    while st:
        x0, x1 = st.pop(); n += 1
        if contact and contact[0] <= x0 and x1 <= contact[1]:
            continue
        X = Jet(iv(x0, x1), n=1, order=0)
        try: v = f(X).v
        except (ValueError, ZeroDivisionError): v = NAN()
        if v > 0:
            lbv = v.lower(); mn = lbv if mn is None or lbv < mn else mn; continue
        try:
            xm = (x0 + x1) / 2
            vc = f(Jet(arb(xm), n=1, order=0)).v
            if vc < 0:
                ok(f"{name} on [{a},{b}]", False, f"counterexample: value at x={xm!r} is {vc}"); return False
            g = f(Jet.var(iv(x0, x1), 0, 1, 1)).d[0]
            lb = vc if g is None else vc - abs(g).upper() * exrad(x0, x1, xm)
            if lb > 0:
                lbv = lb.lower(); mn = lbv if mn is None or lbv < mn else mn; continue
        except (ValueError, ZeroDivisionError):
            pass
        if x1 - x0 < minw or n > maxn:
            ok(f"{name} on [{a},{b}]", False, f"stuck at [{x0!r},{x1!r}] value {v}")
            return False
        s = split_point(x0, x1); st.append((x0, s)); st.append((s, x1))
    ok(f"{name} on [{a},{b}]", True, f"boxes={n} min lower={float(mn) if mn is not None else 'contact only':.3e} t={time.time()-t:.1f}s")
    return True

def bb2(name, G, boxes, contact=None, sym=True, maxn=20_000_000, minw=1e-14, report=200000):
    st = list(boxes); n = 0; acc = 0; mn = None; t = time.time()
    while st:
        x0, x1, y0, y1 = st.pop(); n += 1
        if n % report == 0: out(f"   .. {name}: {n} boxes, stack {len(st)}, {time.time()-t:.0f}s")
        if sym and x0 > y1: continue
        if contact and contact[0] <= x0 and x1 <= contact[1] and contact[0] <= y0 and y1 <= contact[1]:
            continue
        X = Jet(iv(x0, x1), n=2, order=0); Y = Jet(iv(y0, y1), n=2, order=0)
        try: v = G(X, Y).v
        except (ValueError, ZeroDivisionError): v = NAN()
        if v > 0:
            lbv = v.lower(); mn = lbv if mn is None or lbv < mn else mn; acc += 1; continue
        good = False
        try:
            xm = (x0 + x1) / 2; ym = (y0 + y1) / 2
            vc = G(Jet(arb(xm), n=2, order=0), Jet(arb(ym), n=2, order=0)).v
            if vc < 0:
                ok(f"{name}", False, f"counterexample: value at (x,y)=({xm!r},{ym!r}) is {vc}"); return False
            g = G(Jet.var(iv(x0, x1), 0, 2, 1), Jet.var(iv(y0, y1), 1, 2, 1))
            lb = vc
            if g.d[0] is not None: lb = lb - abs(g.d[0]).upper() * exrad(x0, x1, xm)
            if g.d[1] is not None: lb = lb - abs(g.d[1]).upper() * exrad(y0, y1, ym)
            if lb > 0:
                good = True; lbv = lb.lower(); mn = lbv if mn is None or lbv < mn else mn; acc += 1
        except (ValueError, ZeroDivisionError):
            pass
        if good: continue
        wx, wy = x1 - x0, y1 - y0
        if max(wx, wy) < minw or n > maxn:
            ok(f"{name}", False, f"stuck at x=[{x0!r},{x1!r}] y=[{y0!r},{y1!r}] value {v}")
            return False
        # split the relatively wider side (relative to distance from the nearer edge)
        rx = wx / max(min(x0, 1 - x1), wx); ry = wy / max(min(y0, 1 - y1), wy)
        if rx >= ry:
            s = split_point(x0, x1); st.append((x0, s, y0, y1)); st.append((s, x1, y0, y1))
        else:
            s = split_point(y0, y1); st.append((x0, x1, y0, s)); st.append((x0, x1, s, y1))
    ok(f"{name}", True, f"boxes={n} accepted={acc} min lower={float(mn) if mn is not None else float('nan'):.3e} t={time.time()-t:.1f}s")
    return True

def hess_pd(name, G, k0, k1, nvar, ngrid=None):
    # G is convex on K = [k0,k1]^nvar if its interval Hessian is positive definite on every cell of a grid
    # covering K. Together with G(z*) = 0 and grad G(z*) = 0 at z* = (r,..) in K (analytic) this gives G >= 0 on K.
    assert k0 < float(r.lower()) and float(r.upper()) < k1
    ng = ngrid or (64 if nvar == 1 else 24)
    pts = [k0 + (k1 - k0) * i / ng for i in range(ng + 1)]
    worst = None
    for i in range(ng):
        if nvar == 1:
            Hh = G(Jet.var(iv(pts[i], pts[i + 1]), 0, 1, 2)).H.get((0, 0))
            if Hh is None or not Hh > 0:
                ok(f"{name} second derivative > 0 on [{k0},{k1}]", False, f"cell {i}: [{Hh.lower()},{Hh.upper()}]"); return False
            w = Hh.lower(); worst = w if worst is None or w < worst else worst
            continue
        for j in range(ng):
            st = [(pts[i], pts[i + 1], pts[j], pts[j + 1], 0)]
            while st:
                a0_, a1_, b0_, b1_, dep = st.pop()
                X = Jet.var(iv(a0_, a1_), 0, 2, 2); Y = Jet.var(iv(b0_, b1_), 1, 2, 2)
                Hm = G(X, Y).H
                h11, h12, h22 = Hm.get((0, 0)), Hm.get((0, 1), ZERO), Hm.get((1, 1))
                det = h11.lower() * h22.lower() - abs(h12).upper() ** 2
                if h11 > 0 and h22 > 0 and det > 0:
                    worst = det if worst is None or det < worst else worst; continue
                if dep >= 6:
                    ok(f"{name} Hessian PD on [{k0},{k1}]^2", False, f"cell x[{a0_},{a1_}] y[{b0_},{b1_}]: H11 [{h11.lower()},{h11.upper()}] |H12|<={abs(h12).upper()} H22 [{h22.lower()},{h22.upper()}]"); return False
                am, bm_ = (a0_ + a1_) / 2, (b0_ + b1_) / 2
                st += [(a0_, am, b0_, bm_, dep + 1), (am, a1_, b0_, bm_, dep + 1), (a0_, am, bm_, b1_, dep + 1), (am, a1_, bm_, b1_, dep + 1)]
    ok(f"{name} {'second derivative' if nvar == 1 else 'Hessian'} positive definite on [{k0},{k1}]^{nvar} ({ng}^{nvar} cells)", True, f"min lower {'f2' if nvar == 1 else 'det'} = {float(worst):.4e}")
    return True

# ---------------------------------------------------------------- the kernel two-variable functions
def G1(x, y):
    co = Q(1, 20) if not RED else Q(1, 5)
    return F(x, y) - (phi(x) * R(y) + phi(y) * R(x)) * co
def G2(x, y):
    co = Q(11, 50) if not RED else Q(1, 10)
    return Jk(x, y) - (phi(x) * bfun(y) + phi(y) * bfun(x)) * Q(1, 5) + (R(x) + R(y)) * co
def G3(x, y):
    co = Q(11, 50) if not RED else Q(1, 1)
    vx, vy = vfun(x), vfun(y)
    return F(x, y) - (vx * vx * phi(y) + vy * vy * phi(x)) * co
def G4(x, y):
    co = Q(4, 25) if not RED else Q(1, 1)
    dx = phi(x) * sqrtL - jh(x) * Q(7, 5); dy = phi(y) * sqrtL - jh(y) * Q(7, 5)
    return F(x, y) - jpos(dx) * jpos(dy) * co

def compact_boxes(d0, d1, corner):
    # [d0, 1-d1]^2 minus [1-corner, 1]^2, as a list of rectangles; coarse initial grid
    cuts = [d0, 1e-6, 1e-3, 0.05, 0.2, 0.4, 0.6, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.8, 0.9, 0.99, 1 - 1e-3, 1 - corner, 1 - d1]
    bx = []
    for i in range(len(cuts) - 1):
        for j in range(len(cuts) - 1):
            x0, x1, y0, y1 = cuts[i], cuts[i + 1], cuts[j], cuts[j + 1]
            if x0 >= 1 - corner - 1e-18 and y0 >= 1 - corner - 1e-18: continue
            bx.append((x0, x1, y0, y1))
    return bx

# ================================================================ groups
def group_SC():
    out("== SC0: constants and sign conditions (paper Sections 2-7)")
    ok("r bracket width < 1e-30", (r.upper() - r.lower()) < Q(1, 10**30), str(r))
    ok("0 < tau < 1/2", tau > 0 and tau < HALF, str(tau))
    ok("q > 5/22", q > Q(5, 22), str(q))
    ok("44 qL > 1", 44 * qL > 1)
    ok("c > 0, p < 0", c > 0 and p < 0)
    K = scalar_constants(eps)
    for nm, cond in [("alpha>0", alpha > 0), ("b_m>0", bm > 0), ("D>0", Dc > 0), ("D<kappa", Dc < kap),
                     ("e*>z", es > z), ("k_c>0", K['kc'] > 0), ("D0>0", K['D0'] > 0), ("k_h+m>0", K['kh'] + m > 0),
                     ("Btot>0", K['Btot'] > 0), ("Btot<A0", K['Btot'] < K['A0']), ("2(m-eps)>lambda", 2 * (m - eps) > lam),
                     ("Gamma>25q/2", Gam > 25 * q / 2), ("df/ds>0", 2 * z - 2 * qL / c ** 2 * (1 + L / ap) > 0),
                     ("tmin>0", tmin > 0), ("W>0", K['W'] > 0), ("eta>0", K['eta'] > 0),
                     ("G empty case: lambda-2 Gamma eps>1", lam - 2 * Gam * eps > 1), ("e1>0", 2 * (m - eps) - lam > 0)]:
        ok("constant " + nm, cond)
    out("== SC1: edge scalar conditions (EDGE_LEMMAS.md)")
    ok("K1 small: sqrt(2 log(20/3)/L) > rho0", (2 * ell(Q(3, 20)) / L).sqrt() > rho0)
    ok("K1 near 1: h(1/100) < (qL)^2/(r^2 L)", hq(1, 100) < qL ** 2 / (r * r * L))
    l = ell(Q(1, 200))
    ok("K2 small: (hr/|p|-r) ell + (1-hr^2/|p|) sqrt(2ell/L) - rho0 > 0, coefficient hr/|p|-r>0",
       (hr / ap - r) * l + (1 - hr * hr / ap) * (2 * l / L).sqrt() - rho0 > 0 and hr / ap - r > 0 and 1 - hr * hr / ap > 0)
    kp = r / qL - (1 - r + hr / ap - Q(1, 200))
    ok("K2 near 1: sqrt(h(1/200)) < (1-hr^2/|p|)/(kp sqrt L)", hq(1, 200).sqrt() < (1 - hr * hr / ap) / (kp * sqrtL))
    l = ell(Q(1, 500)); Pu = Phiu(l)
    ok("K3 small: 1 + ell/|p| - c Phiu(ell) > 0 and derivative > 0 at ell(1/500)",
       1 + l / ap - c * Pu > 0 and 1 / ap - c * (1 + tau) / (L * Pu) > 0 and (1 + tau) * Q(1, 500) ** 2 <= HALF)
    ok("K3 near 1: 1 - 1/100 - c sqrt(h(1/50)/L) > 0", 1 - Q(1, 100) - c * (hq(1, 50) / L).sqrt() > 0)
    ok("K3 1-b small: 1 - 1/100 - h(1/100)/|p| > 0", 1 - Q(1, 100) - hq(1, 100) / ap > 0)
    ok("K3 1-b near 1: sqrt(L) < c|p|/sqrt(L)", sqrtL < c * ap / sqrtL)
    ok("K4 h(1/100) < 27/400", hq(1, 100) < Q(27, 400))
    ok("K4 h(1/5000) < 21/10000", hq(1, 5000) < Q(21, 10000))
    ok("W1 small: (2/5) sqrt(2 log 10 / L) > 1", Q(2, 5) * (2 * ell(Q(1, 10)) / L).sqrt() > 1)
    ok("W1 near 1: h(1/25) <= 4/(25L)", hq(1, 25) <= Q(4, 25) / L)
    ok("W2 small: sqrt(2 log(100/69)/L) > 0.032 and 69/100 < r", (2 * ell(Q(69, 100)) / L).sqrt() > MU and Q(69, 100) < r)
    l = ell(Q(1, 1000)); xx = Q(1, 1000)
    phx = Jet(xx, n=1, order=0); phi_k = phi(phx).v
    ok("W3 small: x(ell+1)^2/sqrt(2ell/L) <= 9/100 and hr^2 phi(1/1000) <= 9/100 (x(ell+1)^2/sqrt(2ell/L) increasing)",
       xx * (l + 1) ** 2 / (2 * l / L).sqrt() <= Q(9, 100) and hr * hr * phi_k <= Q(9, 100))
    ok("W3 near 1: sqrt(L) h(1/400)^1.5 <= 9/100 and hr^2 sqrt(h(1/200)/L) <= 9/100",
       sqrtL * hq(1, 400) * hq(1, 400).sqrt() <= Q(9, 100) and hr * hr * (hq(1, 200) / L).sqrt() <= Q(9, 100))
    ok("reduction: kappa4 < qL - 1/10 < qL - 99/2500, kappa4 > 0", kappa4 > 0 and kappa4 < qL - Q(1, 10) and qL - Q(1, 10) < qL - Q(99, 2500))
    ok("X0: kappa4 sqrt(2 log(1e10)/L) >= 2/5", kappa4 * (2 * ell(D0P) / L).sqrt() >= Q(2, 5))
    ok("X1: h(2e-7) <= (0.032 kappa4)^2 / L", h_point(D1P) <= (MU * kappa4) ** 2 / L)
    ok("K11: kappa4 (2-(1+tau)beta) log(1/(2beta))/L >= 2L", kappa4 * (2 - (1 + tau) * BETA) * ell(2 * BETA) / L >= 2 * L)
    ok("Z1: h(5e-5) <= (qL/(11 r))^2 / L", h_point(Z1) <= (qL / (11 * r)) ** 2 / L)
    Phi1 = (h_point(2 * D1P) / L).sqrt(); la = ell(YA)
    hb = la / Phiu(la); yb = (L / (2 * la)).sqrt()
    val = (1 - D1P - Phi1 / (5 * ap)) * hb - L * Phi1 - Q(1, 5) * (1 + h_point(D1P) / ap) - Q(1, 5) * Phi1 * yb
    ok("Z3: bound at ell(y_a) > 0 and 1-delta1-Phi1/(5|p|) > 0", val > 0 and 1 - D1P - Phi1 / (5 * ap) > 0, str(val))
    l1 = ell(Y1); rho1 = 1 - rho0 * (L / (2 * l1)).sqrt(); beta1 = Y1 + h_point(Y1) / ap
    ok("Z4a: (11/50) rho1 - beta1/5 > 0", Q(11, 50) * rho1 - beta1 / 5 > 0, str(Q(11, 50) * rho1 - beta1 / 5))
    ok("Z4: (1+tau) x^2 <= 1/2 and 3x/2 <= 1/e at x = 1e-10", (1 + tau) * D0P ** 2 <= HALF and Q(3, 2) * D0P < Q(1, 3))
    d = Q(1, 2 ** 16); dl = Q(8, 10000)
    ok("C1 small x: (6/100) log(100/3) < kappa - 1/12", Q(6, 100) * ell(Q(3, 100)) < kap - Q(1, 12))
    ok("C2 x<=d: delta(1/12 - d) - (kappa - 1/12 + d) h(d) > 0", dl * (Q(1, 12) - d) - (kap - Q(1, 12) + d) * h_point(d) > 0)
    ok("C2 x>=63/64: kappa (3/2) h(2/3)(63/64) - delta 11/12 - (1+kappa) h(1/64) > 0",
       kap * Q(3, 2) * hq(2, 3) * Q(63, 64) - dl * Q(11, 12) - (1 + kap) * hq(1, 64) > 0)
    ok("C4: h(9/10) <= kappa (4/5) log 9", hq(9, 10) <= kap * Q(4, 5) * arb(9).log())
    ok("K5 Boppana: 809/500 < (1+sqrt5)/2", Q(809, 500) < (1 + arb(5).sqrt()) / 2)

def group_S6():
    out("== S6: M(u) > 0 for u in [0, 2 eps] (paper Prop 7.1/S6)" + (" [RED: eps = 93e-6]" if RED else ""))
    K = scalar_constants(eps)
    st = [(0.0, 1.0)]; n = 0; mn = None
    while st:
        a, b = st.pop(); n += 1
        val = M_of_v(K, iv(a, b))
        if val > 0:
            mn = val.lower() if mn is None or val.lower() < mn else mn; continue
        if b - a < 1e-9 or n > 10**6:
            ok("S6 M>0", False, f"at v in [{a},{b}] M in {val}"); return
        s = (a + b) / 2; st += [(a, s), (s, b)]
    ok("S6 M(u)>0 on [0,2eps]", True, f"intervals={n} min lower={float(mn):.4e}")

def group_S5():
    out("== S5: Psi(e,t,s) >= 0 on the constrained box (paper Lemma 5.2/Prop S5)" + (" [RED: lambda=1.240]" if RED else ""))
    lam_ = lam if not RED else Q(1240, 1000)
    sL = sqrtL; e1 = (2 * (m - eps) - lam) / (2 * q * Q(49, 25))
    def Psi(e, t, s):
        U0 = (Q(809, 500) * s * e).max(2 * s * e - qL * t * t)
        a1 = U0 + qL * t * t - 2 * s * e
        a2 = 2 * kE * t * (e - hr * t) ** 2
        a3 = t * t / Db * (s - e / p - c * t)
        w = sL * t - Q(7, 5) * e
        a4 = Q(4, 25) * w.max(ZERO) ** 2 if w.upper() > 0 else ZERO
        mx = a1.max(a2).max(a3).max(a4)
        return U0 + Gam * mx - lam_ * e
    e1l = float(e1.lower()); zf = 77 / 200; sl = float((m - eps).lower())
    st = [(e1l, zf, 0.0, 1.0, sl, 1.0)]; n = 0; disc = 0; mn = None; t = time.time()
    while st:
        e0, e1_, t0_, t1, s0, s1 = st.pop(); n += 1
        E = iv(e0, e1_); T = iv(t0_, t1); S = iv(s0, s1)
        if (arb(e0) - Q(4, 5) * arb(t1) - Q(27, 400)) > 0 or (arb(e0) - Q(6, 5) * arb(t1) - Q(21, 10000)) > 0:
            disc += 1; continue
        v = Psi(E, T, S)
        if v.upper() < 0 and (arb(e0) - Q(4, 5) * arb(t0_) - Q(27, 400)) < 0 and (arb(e0) - Q(6, 5) * arb(t0_) - Q(21, 10000)) < 0:
            # corner (e0, t0_, s0) region: confirm with a point evaluation inside the domain
            pv = Psi(arb(e0), arb(t1), arb(s0))
            if pv < 0 and (arb(e0) - Q(4, 5) * arb(t1) - Q(27, 400)) <= 0 and (arb(e0) - Q(6, 5) * arb(t1) - Q(21, 10000)) <= 0:
                ok("S5 Psi>=0", False, f"counterexample at e={e0}, t={t1}, s={s0}: Psi = {pv}"); return
        if v > 0:
            mn = v.lower() if mn is None or v.lower() < mn else mn; continue
        if n > 30_000_000 or max(e1_ - e0, t1 - t0_, s1 - s0) < 1e-13:
            ok("S5 Psi>=0", False, f"stuck e=[{e0},{e1_}] t=[{t0_},{t1}] s=[{s0},{s1}] val {v}"); return
        # split the dimension with largest scaled width
        we, wt, ws = (e1_ - e0) / 0.4, (t1 - t0_) / 1.0, (s1 - s0) / 0.4
        if we >= wt and we >= ws:
            md = (e0 + e1_) / 2; st += [(e0, md, t0_, t1, s0, s1), (md, e1_, t0_, t1, s0, s1)]
        elif wt >= ws:
            md = (t0_ + t1) / 2; st += [(e0, e1_, t0_, md, s0, s1), (e0, e1_, md, t1, s0, s1)]
        else:
            md = (s0 + s1) / 2; st += [(e0, e1_, t0_, t1, s0, md), (e0, e1_, t0_, t1, md, s1)]
    ok("S5 Psi>0 on domain (e in [e1_lower, z])", True, f"boxes={n} discarded={disc} min lower={float(mn):.3e} t={time.time()-t:.0f}s")
    ok("S5 e1 lower endpoint used is <= e1", arb(e1l) <= e1)

def group_K():
    out("== K: one-variable kernel inequalities (paper S1, S3)")
    rl, ru = float(r.lower()), float(r.upper())
    if RED:
        bb1("RED K1 R(x) - 0.001 >= 0", lambda x: R(x) - Q(1, 1000), 0.15, 0.99)
        return
    hess_pd("H1 R", R, 0.66, 0.72, 1)
    bb1("K1 R>0", R, 0.15, 0.99, contact=(0.66, 0.72))
    S1b = lambda x: (x - r) * jh(x) - vfun(x) * (hr / p) + R(x)
    hess_pd("H2 S1b", S1b, 0.66, 0.72, 1)
    bb1("K2 S1b>0", S1b, 0.005, 0.995, contact=(0.66, 0.72))
    hess_pd("H3 b", bfun, 0.66, 0.72, 1)
    bb1("K3 b>0", bfun, 0.002, 0.99, contact=(0.66, 0.72))
    bb1("K3 1-b>0", lambda x: 1 - bfun(x), 0.01, 0.75)
    bb1("K4 (4/5)phi+27/400-h>0", lambda x: phi(x) * Q(4, 5) + Q(27, 400) - jh(x), 0.01, 0.99)
    bb1("K4 (6/5)phi+21/10000-h>0", lambda x: phi(x) * Q(6, 5) + Q(21, 10000) - jh(x), 1 / 5000, 1 - 1 / 5000)
    bb1("K5 h(t^2)-(809/500) t h(t) > 0", lambda x: jh(x * x) - x * jh(x) * Q(809, 500), 0.001, 0.999)
    bb1("K6 7h/5 - sqrt(L) phi > 0 on [0.62,0.719]", lambda x: jh(x) * Q(7, 5) - phi(x) * sqrtL, 0.62, 0.719)

def group_W():
    out("== W: helper one-variable inequalities")
    if RED:
        bb1("RED W1 (1/5)phi - omega", lambda y: phi(y) * Q(1, 5) + (1 - y) * jlog(1 - y), 0.1, 0.96); return
    bb1("W1 (2/5)phi - omega > 0", lambda y: phi(y) * Q(2, 5) + (1 - y) * jlog(1 - y), 0.1, 0.96)
    bb1("W2 phi - 0.032 y > 0", lambda y: phi(y) - y * MU, 0.69, 1 - 4e-5)
    bb1("W3 (9/100)phi - v^2 > 0", lambda x: phi(x) * Q(9, 100) - vfun(x) * vfun(x), 0.001, 1 - 1 / 400)

def group_Z():
    out("== Z: edge checks for G2 (J inequality)")
    Phi1 = (h_point(2 * D1P) / L).sqrt(); hd1 = h_point(D1P)
    def Z2(y):
        val = jh(y) - hd1 - phi(y) * (L * Phi1) - bfun(y) * (Phi1 / 5) - phi(y) * ((1 + hd1 / ap) / 5) + R(y) * Q(11, 50)
        return val if not RED else val - Q(1, 1000)
    bb1("Z2 on [y_a, 1-z1]", Z2, 1e-6, 1 - 5e-5)
    if RED: return
    l0 = ell(D0P); Pu0 = Phiu(l0)
    def Z4b(y0, y1):
        Y = iv(y0, y1); X = Jet(Y, n=1, order=0)
        ph = phi(X).v; bb = bfun(X).v; RR = R(X).v
        Klo = 1 + tau * (1 - D0P) * (1 - Y)
        yK = Y * Klo
        A = yK - ph / (5 * ap)
        P = Q(11, 50) - L * ph - bb / 5
        lyK = -(yK.log())
        base = yK * (l0 + lyK) - ph / 5 * (1 + (l0 + 1) / ap) - Q(11, 50) * rho0
        pm = P * (2 * l0 / L).sqrt(); pu = P * Pu0
        Lam = base + arb(min(pm.lower(), pu.lower()))
        Pneg = arb(max(0, (-P).upper()))
        beta_ok = Lam > 0 and (A - Pneg * (1 + tau) / (L * Pu0)) >= 0
        C0 = abs(yK * yK.log()) + ph / 5 * (1 + 1 / ap) + Q(11, 50) * rho0
        alpha_ok = (Q(11, 50) * RR - D0P * (abs(A) * l0 + C0 + abs(P) * Pu0)) > 0
        return beta_ok or alpha_ok
    st = [(1e-3, 1.0)]; n = 0
    while st:
        a, b = st.pop(); n += 1
        if Z4b(a, b): continue
        if b - a < 1e-13: ok("Z4b on [y1,1]", False, f"stuck at [{a},{b}]"); return
        s = (a + b) / 2; st += [(a, s), (s, b)]
    ok("Z4b on [y1,1] (alpha or beta on every piece)", True, f"intervals={n}")

def group_C():
    out("== C: covariance inequality (paper Prop S4)")
    kp = kap if not RED else Q(37, 100)
    bb1("C1 kappa - 1/12 - k(x) > 0", lambda x: kp - Q(1, 12) - (2 * x * (1 - x) * jlog((1 - x) / x) - x), 0.03, 0.5)
    dl = Q(8, 10000)
    def Gcov(a, x):
        return kp * jh(a) * x / a + dl * (a - x) - (x + kp - a) * jh(x) + R(x) * Ccov
    bb2("C2 Gcov>0 on [1/12,2/3]x[2^-16,63/64]", Gcov, [(1 / 12, 2 / 3, 2.0 ** -16, 63 / 64)], sym=False)
    def Pkl(a, x):
        KL = -jh(x) - x * jlog(a) - (1 - x) * jlog(1 - a)
        return KL * kp - (a - x) * (jh(a) - jh(x))
    # domain x <= 1-a: boxes with x0 >= 1 - a0 lie outside it (and there the left side is <= 0 anyway)
    st = [(2 / 3, 0.9, 0.0, 0.34)]; n = 0; mn = None
    while st:
        a0_, a1_, x0, x1 = st.pop(); n += 1
        if arb(x0) >= 1 - arb(a0_): continue
        A = Jet(iv(a0_, a1_), n=2, order=0); X = Jet(iv(x0, x1), n=2, order=0)
        try: v = Pkl(A, X).v
        except (ValueError, ZeroDivisionError): v = NAN()
        if v > 0: mn = v.lower() if mn is None or v.lower() < mn else mn; continue
        try:
            vc = Pkl(Jet(arb((a0_ + a1_) / 2), n=2, order=0), Jet(arb((x0 + x1) / 2), n=2, order=0)).v
            g = Pkl(Jet.var(iv(a0_, a1_), 0, 2, 1), Jet.var(iv(x0, x1), 1, 2, 1))
            lb = vc - abs(g.d[0]).upper() * exrad(a0_, a1_, (a0_ + a1_) / 2) - abs(g.d[1]).upper() * exrad(x0, x1, (x0 + x1) / 2)
            if lb > 0: mn = lb.lower() if mn is None or lb.lower() < mn else mn; continue
        except (ValueError, ZeroDivisionError, TypeError):
            pass
        if max(a1_ - a0_, x1 - x0) < 1e-13 or n > 5_000_000:
            ok("C3 KL", False, f"stuck a=[{a0_},{a1_}] x=[{x0},{x1}] v={v}"); return
        if a1_ - a0_ > x1 - x0:
            s = (a0_ + a1_) / 2; st += [(a0_, s, x0, x1), (s, a1_, x0, x1)]
        else:
            s = (x0 + x1) / 2; st += [(a0_, a1_, x0, s), (a0_, a1_, s, x1)]
    ok("C3 kappa KL - (a-x)(h(a)-h(x)) > 0 on {2/3<=a<=9/10, x<=1-a}", True, f"boxes={n} min lower={float(mn):.3e}")

def group_H():
    out("== RB: contact identities at the r-ball (value and gradient enclosures must contain 0)")
    X = Jet.var(r, 0, 2, 1); Y = Jet.var(r, 1, 2, 1); x1 = Jet.var(r, 0, 1, 1)
    for nm, G in [("G1", G1), ("G2", G2), ("G3", G3), ("G4", G4)]:
        g = G(X, Y)
        ok(f"RB {nm}(r,r) = 0 and grad = 0 (enclosures contain 0)", g.v.contains(0) and all(d is None or d.contains(0) for d in g.d))
    S1b = lambda x: (x - r) * jh(x) - vfun(x) * (hr / p) + R(x)
    for nm, f in [("R", R), ("S1b", S1b), ("b", bfun)]:
        g = f(x1)
        ok(f"RB {nm}(r) = 0 and {nm}'(r) = 0 (enclosures contain 0)", g.v.contains(0) and g.d[0].contains(0))
    out("== H: Hessians at the contact for the two-variable inequalities")
    for nm, G, k0, k1 in [("H4 G1", G1, 0.675, 0.705), ("H5 G2", G2, 0.675, 0.705), ("H6 G3", G3, 0.675, 0.705), ("H7 G4", G4, 0.675, 0.705)]:
        hess_pd(nm, G, k0, k1, 2)

def group_G(which):
    G = {'G1': G1, 'G2': G2, 'G3': G3, 'G4': G4}[which]
    corner = 5e-5 if which == 'G2' else 4e-5
    out(f"== {which}: compact domain [1e-10, 1-2e-7]^2 minus corner [1-{corner}, 1]^2, contact box [0.675,0.705]^2" + (" [RED]" if RED else ""))
    bb2(f"{which} > 0 off contact", G, compact_boxes(1e-10, 2e-7, corner), contact=(0.675, 0.705))

if __name__ == '__main__':
    out(f"Program A (Arb, python-flint 0.7.1, prec {ctx.prec}) group {GROUP} red={RED} start {time.strftime('%Y-%m-%d %H:%M:%S %z')}")
    groups = {'SC': group_SC, 'S6': group_S6, 'S5': group_S5, 'K': group_K, 'W': group_W, 'Z': group_Z, 'C': group_C, 'H': group_H,
              'G1': lambda: group_G('G1'), 'G2': lambda: group_G('G2'), 'G3': lambda: group_G('G3'), 'G4': lambda: group_G('G4')}
    if GROUP == 'all':
        for g in ['SC', 'S6', 'K', 'W', 'Z', 'C', 'H', 'S5', 'G1', 'G2', 'G3', 'G4']: groups[g]()
    else:
        groups[GROUP]()
    out(f"radius underestimates by the original code (boxes): {UNDER[0]}")
    out(f"END group {GROUP} red={RED} {'FAILED' if FAILED else 'ALL PASS'} {time.strftime('%Y-%m-%d %H:%M:%S %z')}")
