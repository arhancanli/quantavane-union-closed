# Round 13, Program B (written separately from certA.py). Binary64 intervals with one-ulp outward widening after every
# operation (math.nextafter); constants and transcendental point values from mpmath.iv (60 digits), rounded outward.
# Certifies N1-N5 of Theorem JB for a params JSON file. Quantavane union-closed lane, 18 Nov 2026.
import sys, json, math, time
from collections import deque
from fractions import Fraction
import mpmath
from mpmath import iv
iv.dps = 60
INF = math.inf
dn = lambda v: math.nextafter(v, -INF)
up = lambda v: math.nextafter(v, INF)

class X:
    __slots__ = ('a', 'b')
    def __init__(s, a, b=None):
        s.a = a; s.b = a if b is None else b
    def __add__(s, o):
        o = cv(o); return X(dn(s.a + o.a), up(s.b + o.b))
    __radd__ = __add__
    def __neg__(s): return X(-s.b, -s.a)
    def __sub__(s, o):
        o = cv(o); return X(dn(s.a - o.b), up(s.b - o.a))
    def __rsub__(s, o): return cv(o) - s
    def __mul__(s, o):
        o = cv(o); ps = (s.a * o.a, s.a * o.b, s.b * o.a, s.b * o.b)
        return X(dn(min(ps)), up(max(ps)))
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = cv(o)
        assert o.a > 0 or o.b < 0, "division by interval containing 0"
        qs = (s.a / o.a, s.a / o.b, s.b / o.a, s.b / o.b)
        return X(dn(min(qs)), up(max(qs)))
    def __rtruediv__(s, o): return cv(o) / s
    def __repr__(s): return f"[{s.a!r}, {s.b!r}]"
def cv(v):
    if isinstance(v, X): return v
    if isinstance(v, int) and abs(v) < 2**52: return X(float(v))
    if isinstance(v, Fraction):
        f = float(v)
        return X(f) if Fraction(f) == v else X(dn(f), up(f))
    if isinstance(v, float): return X(v)
    raise TypeError(type(v))
def fromiv(w):       # mpmath interval -> outward binary64 interval
    return X(dn(float(w.a)), up(float(w.b)))
def hull(*xs): return X(min(x.a for x in xs), max(x.b for x in xs))
def vmax(x, y): return X(max(x.a, y.a), max(x.b, y.b))
def ppos(x): return X(max(x.a, 0.0), max(x.b, 0.0))
def sqr(x):
    lo_ = 0.0 if x.a <= 0 <= x.b else min(abs(x.a), abs(x.b))
    hi_ = max(abs(x.a), abs(x.b))
    return X(dn(lo_ * lo_) if lo_ > 0 else 0.0, up(hi_ * hi_))
def xsqrt(x):
    assert x.a >= 0 or x.a > -1e-300
    return X(dn(math.sqrt(max(x.a, 0.0))) if x.a > 0 else 0.0, up(math.sqrt(x.b)))

# ---- constants via mpmath.iv -------------------------------------------------------------------------------------
def ent_iv(x):   # x: mpmath iv point/interval strictly inside (0,1)
    return -x * iv.log(x) - (1 - x) * iv.log(1 - x)
def rootfun(x): return -(1 - x) * iv.log(1 - x) - x * iv.log(1 + x)
mpmath.mp.dps = 60
r_f = mpmath.findroot(lambda x: -(1 - x) * mpmath.log(1 - x) - x * mpmath.log(1 + x), 0.69)
rlo, rhi = r_f - mpmath.mpf(10) ** -50, r_f + mpmath.mpf(10) ** -50
assert rootfun(iv.mpf(rlo)).a > 0 and rootfun(iv.mpf(rhi)).b < 0, "root bracket"
R_ = iv.mpf([rlo, rhi]); L_ = iv.log(2)
M_ = iv.mpf(1) / 2 + (1 - R_) * iv.log(1 + R_) / (2 * L_); Q_ = 2 * M_ - 1
TAU_ = (iv.mpf(1) / 2 - R_ ** 2) / (R_ ** 2 * (1 - R_) ** 2)
HR_ = ent_iv(R_); P_ = iv.log((1 - R_) / R_); C_ = R_ - HR_ / P_
GAM_ = HR_ / R_; ALP_ = 1 - GAM_ / P_; T0_ = R_ - M_; BM_ = -HR_ * GAM_ / (P_ * ALP_)
D_ = HR_ / (abs(P_) * ALP_) - T0_          # = -h_r/(p alpha) - t0, written with |p|
ES_ = GAM_ * M_; AST_ = M_ * ent_iv(R_ ** 2) / (R_ * HR_)
D0_ = 2 * M_ - AST_ - 2 * Q_ * L_ * (M_ / R_) / (abs(P_) * C_)
DB_ = iv.mpf(11) / (1 - 1 / (44 * Q_ * L_))
KAP_ = iv.mpf(39) / 100; CC_ = iv.mpf(13) / 25; KE_ = iv.mpf(11) / 50; KF_ = iv.mpf(1) / 20

Pj = json.load(open(sys.argv[1]))
fr = lambda k: Fraction(Pj[k])
EPSf, Zf, A0f, MUf, THf, PSIf, RMAXf, E1f = (fr(k) for k in ['eps', 'z', 'a0', 'mu', 'theta', 'Psi', 'Rmax', 'e1'])
mpq = lambda f: iv.mpf(f.numerator) / f.denominator
EPS_, A0_, MU_, TH_, PSI_ = mpq(EPSf), mpq(A0f), mpq(MUf), mpq(THf), mpq(PSIf)
BETA_ = 1 - D_ / KAP_
B0_ = iv.exp(BETA_ * iv.log(A0_))
W_ = (KAP_ - D_) / (1 - B0_); ETA_ = (KAP_ - D_) / D_ * (B0_ - A0_) / (1 - B0_)
AA_ = W_ / 2; BP_ = ETA_ * (T0_ + D_) + W_ * (A0_ - iv.mpf(1) / 2) - MU_ * D0_
G0_ = ETA_ * T0_ - MU_ * (1 - AST_)
assert BETA_.a > 0 and (AA_ - BP_).a > 0

r, L, m, q, hr, p, c = map(fromiv, (R_, L_, M_, Q_, HR_, P_, C_))
tau, gam, t0, bm, D, es, astar, d0, Db = map(fromiv, (TAU_, GAM_, T0_, BM_, D_, ES_, AST_, D0_, DB_))
kE, kF, Ccov = fromiv(KE_), fromiv(KF_), fromiv(CC_)
eps, mu, th, Psi = map(fromiv, (EPS_, MU_, TH_, PSI_))
Apl, Amn, Gamma0 = fromiv(AA_ + BP_), fromiv(AA_ - BP_), fromiv(G0_)
sqrtL = fromiv(iv.sqrt(L_))
MINS = [tuple(cv(Fraction(v)) for v in row) for row in Pj['minorants']]
F45, F65, F75, F9825 = cv(Fraction(4, 5)), cv(Fraction(6, 5)), cv(Fraction(7, 5)), cv(Fraction(98, 25))

# ---- point values and monotone enclosures ------------------------------------------------------------------------
def ent_at(x):          # x a float in [0,1]; returns outward interval of h(x)
    if x <= 0.0 or x >= 1.0: return X(0.0)
    return fromiv(ent_iv(iv.mpf(x)))
_cache_phi = {}
def phi_at(x):
    if x <= 0.0 or x >= 1.0: return X(0.0)
    if x in _cache_phi: return _cache_phi[x]
    xx = iv.mpf(x); Ixx = xx * xx * (1 + TAU_ * (1 - xx) ** 2)
    v = fromiv(iv.sqrt(ent_iv(Ixx) / L_)); _cache_phi[x] = v; return v
def ent_box(xa, xb):    # h increasing on [0,1/2], decreasing on [1/2,1]
    ha, hb = ent_at(xa), ent_at(xb)
    if xb <= 0.5: return X(ha.a, hb.b)
    if xa >= 0.5: return X(hb.a, ha.b)
    return X(min(ha.a, hb.a), L.b)
def phi_box(xa, xb):    # phi increasing on [0,r], decreasing on [r,1]  (I(x,x) increasing; h(I) peaks at I = 1/2 <=> x = r)
    pa, pb = phi_at(xa), phi_at(xb)
    if xb < r.a: return X(pa.a, pb.b)
    if xa > r.b: return X(pb.a, pa.b)
    return X(min(pa.a, pb.a), 1.0)

def run1d(name, ok, minw=1e-13, cap=3_000_000):
    q_ = deque([(0.0, 1.0)]); n = 0
    while q_:
        a_, b_ = q_.popleft(); n += 1
        if n > cap: raise RuntimeError(f"{name}: cap")
        if ok(a_, b_): continue
        if b_ - a_ < minw: raise RuntimeError(f"{name}: FAIL near {a_!r}")
        md = 0.5 * (a_ + b_); q_.append((a_, md)); q_.append((md, b_))
    return n

def N1():
    Rm = cv(RMAXf); qL = q * L
    def ok(a_, b_):
        rx_a, rx_b = (r * cv(a_)).a, (r * cv(b_)).b
        hrx_hi = L.b if rx_a <= 0.5 <= rx_b else max(ent_at(max(rx_a, 0.0)).b, ent_at(min(rx_b, 1.0)).b)
        num = X(hrx_hi) + qL * X(phi_box(a_, b_).b) - cv(a_) * hr - r * X(ent_box(a_, b_).a)
        return (num / qL).b < Rm.a
    return run1d('N1', ok)

def N2():
    tot = 0
    for k, (a_, b_, c_, d_) in enumerate(MINS):
        def ok(xa, xb, a_=a_, b_=b_, c_=c_, d_=d_):
            hb = ent_box(xa, xb); pb = phi_box(xa, xb); xx = X(xa, xb)
            rhs = a_ * pb + b_ * hb + c_ * xx + d_
            if rhs.b < 0: return True
            pl = dn(xa * hb.a) if (xa > 0 and hb.a > 0) else 0.0
            lhs_lo = dn(math.sqrt(pl)) if pl > 0 else 0.0
            return lhs_lo > rhs.b
        tot += run1d(f'N2[{k}]', ok)
    return tot

SMIN = (m - eps).a          # float <= m - eps
ZG = dn(float(Zf)) if Fraction(float(Zf)) > Zf else float(Zf)      # <= z
ZT = up(float(Zf)) if Fraction(float(Zf)) < Zf else float(Zf)      # >= z
E1 = dn(float(E1f)) if Fraction(float(E1f)) > E1f else float(E1f)  # <= e1
RMX = up(float(RMAXf))

def lamG(s, e, bb, Rb):
    t = (s - e / p - bb) / c
    vb = e - hr * t; x = e - es
    EF = vmax(2 * kE * t * sqr(vb), 2 * kF * t * Rb)
    Jp = ppos(q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    lb = ppos(bm * bb + (1 + Ccov) * Rb + bm * (m - s))
    v = (Apl * ppos(x) + Amn * ppos(-x) + lb + mu * (astar * e + d0 * x + q * L * sqr(t) - 2 * s * e)
         - (mu + th) * EF - th * Jp - 2 * th * (s - m + eps) * e - Psi * e)
    return v, t

def N3(cap=60_000_000):
    dom = ((SMIN, 1.0), (ZG, L.b), (0.0, 1.0), (0.0, RMX))
    span = [b_ - a_ for a_, b_ in dom]
    q_ = deque([dom]); n = 0
    while q_:
        B = q_.pop(); n += 1
        if n > cap: raise RuntimeError("N3 cap")
        s, e, bb, Rb = (X(a_, b_) for a_, b_ in B)
        v, t = lamG(s, e, bb, Rb)
        if t.b <= 0 or t.a > 1: continue
        if (e - F45 * t - cv(Fraction(27, 400))).a > 0: continue
        if (e - ent_box(B[0][0], min(B[0][1], 1.0))).a > 0: continue
        if v.b <= 0: continue
        # split: largest (width / span) scaled by a fixed sensitivity table (different from Program A's rule)
        sc = [60.0, 1.0, 6.0, 6.0]
        k = max(range(4), key=lambda i: (B[i][1] - B[i][0]) * sc[i])
        if B[k][1] - B[k][0] < 1e-13: raise RuntimeError(f"N3 FAIL box {B} upper {v.b}")
        md = 0.5 * (B[k][0] + B[k][1])
        B1 = list(B); B2 = list(B); B1[k] = (B[k][0], md); B2[k] = (md, B[k][1])
        q_.append(tuple(B1)); q_.append(tuple(B2))
    return n

def N4(cap=60_000_000):
    Pm = mu * (2 * (m - eps) - 1) - m - Gamma0 + Psi
    assert (th - 12.5 * mu * q).a > 0, "theta >= 12.5 mu q"
    assert Pm.a > 0 and (cv(E1f) * F9825 * mu * q - Pm).b < 0, "e1 small enough"
    dom = ((E1, ZT), (0.0, 1.0), (SMIN, 1.0))
    q_ = deque([dom]); n = 0
    while q_:
        B = q_.pop(); n += 1
        if n > cap: raise RuntimeError("N4 cap")
        e, t, s = (X(a_, b_) for a_, b_ in B)
        if (e - F45 * t - cv(Fraction(27, 400))).a > 0: continue
        if (e - F65 * t - cv(Fraction(21, 10000))).a > 0: continue
        if (e - ent_box(B[2][0], min(B[2][1], 1.0))).a > 0: continue
        tt = sqr(t)
        U0 = vmax(cv(Fraction(809, 500)) * s * e, 2 * s * e - q * L * tt)
        g1 = U0 + q * L * tt - 2 * s * e
        g2 = 2 * kE * t * sqr(e - hr * t)
        g3 = tt / Db * ppos(s - e / p - c * t)
        g4 = cv(Fraction(4, 25)) * sqr(ppos(sqrtL * t - F75 * e))
        glo = vmax(vmax(g1, g2), vmax(g3, g4))
        sig = X(0.0)
        for (a_, b_, c_, d_) in MINS: sig = vmax(sig, a_ * t + b_ * e + c_ * s + d_)
        v = (Gamma0 + m + mu - Psi) * e - sqr(sig) - mu * U0 - th * (glo + 2 * (s - m + eps) * e)
        if v.b <= 0: continue
        sc = [3.0, 1.5, 60.0]
        k = max(range(3), key=lambda i: (B[i][1] - B[i][0]) * sc[i])
        if B[k][1] - B[k][0] < 1e-13: raise RuntimeError(f"N4 FAIL box {B} upper {v.b}")
        md = 0.5 * (B[k][0] + B[k][1])
        B1 = list(B); B2 = list(B); B1[k] = (B[k][0], md); B2[k] = (md, B[k][1])
        q_.append(tuple(B1)); q_.append(tuple(B2))
    return n

def N5():
    v = G0_ - PSI_ - 2 * EPS_ * TH_
    assert v.a > 0, f"N5 FAIL {v}"
    return float(v.a)

if __name__ == '__main__':
    which = sys.argv[2:] or ['N5', 'N1', 'N2', 'N4', 'N3']
    print(f"Program B (binary64 outward + mpmath.iv); params {sys.argv[1]}")
    print(f"m in {m}, Gamma0 in {Gamma0}, A0+B' in {Apl}, A0-B' in {Amn}, SMIN = {SMIN!r}")
    for w in which:
        t1 = time.time()
        try:
            print(f"{w}: PASS ({globals()[w]()}) [{time.time()-t1:.1f}s]", flush=True)
        except (AssertionError, RuntimeError) as ex:
            print(f"{w}: FAIL: {ex} [{time.time()-t1:.1f}s]", flush=True)
