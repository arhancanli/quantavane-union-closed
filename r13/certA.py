# Round 13, Program A (python-flint / Arb ball arithmetic). Certifies the finite inequalities N1-N5 of Theorem JB
# (NOTEBOOK 13.4-13.5) for the parameters in params.json. Quantavane union-closed lane, 18 Nov 2026.
#   N1  R(x) <= RMAX on [0,1]
#   N2  sqrt(x h(x)) >= a_k phi(x) + b_k h(x) + c_k x + d_k on [0,1] for every minorant k
#   N3  Lambda_G(s,e,bb,Rb) - Psi e <= 0 on D_G                       (4-D bisection)
#   N4  Lambda_T(e,t,s) - Psi e <= 0 on D_T  (analytic for e <= e1, bisection on [e1, z])
#   N5  Gamma0 - Psi - 2 eps theta > 0
# Usage: python certA.py params.json [N1 N2 N3 N4 N5] [--red]   (--red: no special meaning; red controls use other params files)
import sys, json, time
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 160

def A(v):
    if isinstance(v, arb): return v
    if isinstance(v, Fr): return arb(v.numerator) / v.denominator
    if isinstance(v, int): return arb(v)
    if isinstance(v, float): return arb(Fr(v).numerator) / Fr(v).denominator   # exact binary value
    return arb(v)
def iv(lo, hi):
    lo, hi = A(lo), A(hi)
    return lo.union(hi)
def lo(a): return a.lower()
def hi(a): return a.upper()
def amax(a, b):
    l = lo(a) if lo(a) >= lo(b) else lo(b)
    h_ = hi(a) if hi(a) >= hi(b) else hi(b)
    return l.union(h_)
def pos(a): return amax(a, arb(0))
def sq(a):
    al, ah = a.abs_lower(), a.abs_upper()
    return (al * al).union(ah * ah)

# ---- constants -------------------------------------------------------------------------------------------------
L = arb(2).log()
def h_pt(x):            # binary entropy at an exact point (arb ball allowed; 0 and 1 exactly handled)
    x = A(x)
    if x == 0 or x == 1: return arb(0)
    return -x * x.log() - (1 - x) * (1 - x).log()
def froot(x): x = A(x); return -(1 - x) * (1 - x).log() - x * (1 + x).log()
rl, rh = Fr('0.69090773872540775'), Fr('0.69090773872540776')
assert froot(rl) > 0 and froot(rh) < 0, "r enclosure"
for _ in range(80):     # refine the enclosure; each step certified by a sign
    mid = (rl + rh) / 2
    v = froot(mid)
    if v > 0: rl = mid
    elif v < 0: rh = mid
    else: break
assert froot(rl) > 0 and froot(rh) < 0
r = iv(rl, rh)
m = arb(1) / 2 + (1 - r) * (1 + r).log() / (2 * L); q = 2 * m - 1
tau = (arb(1) / 2 - r * r) / (r * r * (1 - r) ** 2)
hr = h_pt(r); p = ((1 - r) / r).log(); c = r - hr / p
gam = hr / r; alpha = 1 - gam / p; t0 = r - m; bm = -hr * gam / (p * alpha); D = -hr / (p * alpha) - t0
es = gam * m; astar = m * h_pt(r * r) / (r * hr); a_s = m / r; d0 = 2 * m - astar + 2 * q * L * a_s / (p * c)
Db = arb(11) / (1 - 1 / (44 * q * L)); kE = arb(11) / 50; kF = arb(1) / 20; kap = arb(39) / 100; Ccov = arb(13) / 25
assert (b_m_check := (bm - (t0 + D) * gam)).abs_upper() < arb('1e-40'), "identity b_m = (t0+D) gamma"

P = json.load(open(sys.argv[1]))
EPS = Fr(P['eps']); Z = Fr(P['z']); A0f = Fr(P['a0']); MU = Fr(P['mu']); TH = Fr(P['theta']); PSI = Fr(P['Psi'])
RMAX = Fr(P['Rmax']); E1 = Fr(P['e1']); MINS = [tuple(Fr(c_) for c_ in row) for row in P['minorants']]
eps, z, a0, mu, theta, Psi = A(EPS), A(Z), A(A0f), A(MU), A(TH), A(PSI)
beta_ = 1 - D / kap; b0 = (beta_ * a0.log()).exp()
W = (kap - D) / (1 - b0); eta = (kap - D) / D * (b0 - a0) / (1 - b0); A0 = W / 2
Bp = eta * (t0 + D) + W * (a0 - arb(1) / 2) - mu * d0; Gamma0 = eta * t0 - mu * (1 - astar)
assert beta_ > 0 and A0 + Bp > 0 and A0 - Bp > 0

# ---- monotone enclosures ---------------------------------------------------------------------------------------
def h_iv(xl, xh):       # xl <= xh exact (Fractions); h increasing on [0,1/2], decreasing on [1/2,1]
    vl, vh = h_pt(xl), h_pt(xh)
    if xh <= Fr(1, 2): return lo(vl).union(hi(vh))
    if xl >= Fr(1, 2): return lo(vh).union(hi(vl))
    lw = lo(vl) if lo(vl) <= lo(vh) else lo(vh)
    return lw.union(hi(L))
def Ixx_pt(x): x = A(x); return x * x * (1 + tau * (1 - x) ** 2)
def phi_pt(x):
    x = A(x)
    if x == 0 or x == 1: return arb(0)
    hv = h_pt(Ixx_pt(x))
    return (hv / L).sqrt()
# I(x,x) is increasing on [0,1] (derivative 2x[1 + tau(1-x)(1-2x)] > 0 as tau < 1/2), so phi increases on [0,r], decreases on [r,1].
def phi_iv(xl, xh):
    vl, vh = phi_pt(xl), phi_pt(xh)
    if A(xh) < r: return lo(vl).union(hi(vh))
    if A(xl) > r: return lo(vh).union(hi(vl))
    lw = lo(vl) if lo(vl) <= lo(vh) else lo(vh)
    return lw.union(arb(1))

def bisect1d(name, test, a=Fr(0), b=Fr(1), minw=Fr(1, 10**12), cap=2_000_000):
    stack = [(a, b)]; n = 0
    while stack:
        xl, xh = stack.pop(); n += 1
        if n > cap: raise RuntimeError(f"{name}: box cap exceeded")
        if test(xl, xh): continue
        if xh - xl < minw: raise RuntimeError(f"{name}: FAIL near x = {float(xl):.15g}")
        mid = (xl + xh) / 2; stack += [(xl, mid), (mid, xh)]
    return n

def N1():
    # R(x) = [h(r x) + qL phi(x) - x h_r - r h(x)] / (qL)   (phi(r) = 1 exactly since I(r,r) = 1/2)
    rmax = A(RMAX)
    def test(xl, xh):
        rxl, rxh = lo(r * A(xl)), hi(r * A(xh))
        # h(rx) upper bound: h is increasing on [0,1/2], decreasing after; r x in [0, r]
        vals = [h_pt(rxl), h_pt(rxh)]
        hrx_hi = hi(L) if (rxl <= arb(1) / 2 <= rxh) else (hi(vals[0]) if hi(vals[0]) >= hi(vals[1]) else hi(vals[1]))
        ph = hi(phi_iv(xl, xh)); hx_lo = lo(h_iv(xl, xh))
        up = (hrx_hi + q * L * ph - A(xl) * hr - r * hx_lo) / (q * L)
        return up < rmax
    return bisect1d('N1', test)

def N2():
    tot = 0
    for k, (a_, b_, c_, d_) in enumerate(MINS):
        def test(xl, xh, a_=a_, b_=b_, c_=c_, d_=d_):
            hiv = h_iv(xl, xh); ph = phi_iv(xl, xh)
            lhs_lo = (A(xl) * lo(hiv)).sqrt() if lo(hiv) > 0 else arb(0)
            rhs = A(a_) * ph + A(b_) * hiv + A(c_) * iv(xl, xh) + A(d_)
            return lo(lhs_lo) > hi(rhs) or hi(rhs) < 0
        tot += bisect1d(f'N2[{k}]', test)
    return tot

SMIN = Fr(lo(m - eps).str(40, radius=False)) - Fr(1, 10**35)   # rational <= m - eps
assert A(SMIN) <= m - eps
LHI = Fr(hi(L).str(40, radius=False)) + Fr(1, 10**35)

def LamG(s, e, bb, Rb):
    t = (s - e / p - bb) / c
    vb = e - hr * t; x = e - es
    EF = amax(2 * kE * t * sq(vb), 2 * kF * t * Rb)
    Jp = pos(q * bb * (2 * t / 5 - 1 / (25 * L)) - 11 * q / 25 * Rb)
    lb = pos(bm * bb + (1 + Ccov) * Rb + bm * (m - s))
    val = ((A0 + Bp) * pos(x) + (A0 - Bp) * pos(-x) + lb + mu * (astar * e + d0 * x + q * L * sq(t) - 2 * s * e)
           - (mu + theta) * EF - theta * Jp - 2 * theta * (s - m + eps) * e - Psi * e)
    return val, t

def N3(cap=30_000_000):
    dom = [(SMIN, Fr(1)), (Z, LHI), (Fr(0), Fr(1)), (Fr(0), RMAX)]
    wts = [Fr(80), Fr(2), Fr(8), Fr(8)]
    stack = [dom]; n = 0; worst = None
    while stack:
        B = stack.pop(); n += 1
        if n > cap: raise RuntimeError("N3 box cap exceeded")
        s, e, bb, Rb = (iv(a_, b_) for a_, b_ in B)
        val, t = LamG(s, e, bb, Rb)
        # infeasible boxes: t <= 0, t > 1, e > (4/5)t + 27/400, e > h(s)
        if hi(t) <= 0 or lo(t) > 1: continue
        if lo(e - arb(4) / 5 * t - arb(27) / 400) > 0: continue
        if lo(e - h_iv(B[0][0], B[0][1])) > 0: continue
        if hi(val) <= 0: continue
        k = max(range(4), key=lambda i: (B[i][1] - B[i][0]) * wts[i])
        if (B[k][1] - B[k][0]) < Fr(1, 10**13):
            raise RuntimeError(f"N3 FAIL at box {[(float(a_), float(b_)) for a_, b_ in B]} upper {hi(val)}")
        mid = (B[k][0] + B[k][1]) / 2
        B1 = list(B); B2 = list(B); B1[k] = (B[k][0], mid); B2[k] = (mid, B[k][1])
        stack += [B1, B2]
    return n

def N4(cap=30_000_000):
    # analytic part, e <= e1: needs theta >= (25/2) mu q and e1 <= Pm / ((98/25) mu q), Pm = mu (2 (m - eps) - 1) - m - Gamma0 + Psi
    Pm = mu * (2 * (m - eps) - 1) - m - Gamma0 + Psi
    assert theta - arb(25) / 2 * mu * q > 0, "N4 analytic: theta >= 12.5 mu q"
    assert Pm > 0 and A(E1) * arb(98) / 25 * mu * q - Pm < 0, "N4 analytic: e1 <= Pm/(3.92 mu q)"
    dom = [(E1, Z), (Fr(0), Fr(1)), (SMIN, Fr(1))]
    wts = [Fr(4), Fr(2), Fr(80)]
    stack = [dom]; n = 0
    rL = L.sqrt()
    while stack:
        B = stack.pop(); n += 1
        if n > cap: raise RuntimeError("N4 box cap exceeded")
        e, t, s = (iv(a_, b_) for a_, b_ in B)
        if lo(e - arb(4) / 5 * t - arb(27) / 400) > 0: continue
        if lo(e - arb(6) / 5 * t - arb(21) / 10000) > 0: continue
        if lo(e - h_iv(B[2][0], B[2][1])) > 0: continue
        U0 = amax(arb(809) / 500 * s * e, 2 * s * e - q * L * sq(t))
        g_lo = amax(amax(U0 + q * L * sq(t) - 2 * s * e, 2 * kE * t * sq(e - hr * t)),
                    amax(sq(t) / Db * pos(s - e / p - c * t), arb(4) / 25 * sq(pos(rL * t - arb(7) / 5 * e))))
        sig = arb(0)
        for (a_, b_, c_, d_) in MINS:
            sig = amax(sig, A(a_) * t + A(b_) * e + A(c_) * s + A(d_))
        val = (Gamma0 + m + mu - Psi) * e - sq(sig) - mu * U0 - theta * (g_lo + 2 * (s - m + eps) * e)
        if hi(val) <= 0: continue
        k = max(range(3), key=lambda i: (B[i][1] - B[i][0]) * wts[i])
        if (B[k][1] - B[k][0]) < Fr(1, 10**13):
            raise RuntimeError(f"N4 FAIL at box {[(float(a_), float(b_)) for a_, b_ in B]} upper {hi(val)}")
        mid = (B[k][0] + B[k][1]) / 2
        B1 = list(B); B2 = list(B); B1[k] = (B[k][0], mid); B2[k] = (mid, B[k][1])
        stack += [B1, B2]
    return n

def N5():
    v = Gamma0 - Psi - 2 * eps * theta
    assert v > 0, f"N5 FAIL {v}"
    return v

if __name__ == '__main__':
    which = [a for a in sys.argv[2:] if not a.startswith('--')] or ['N5', 'N1', 'N2', 'N4', 'N3']
    print(f"Program A (Arb, prec {ctx.prec}); params {sys.argv[1]}: {P}")
    print(f"r in [{float(rl):.17f}, {float(rh):.17f}], m = {m.str(20)}, c* + eps = {(1 - m + eps).str(20)}")
    print(f"eta = {eta.str(12)} W = {W.str(12)} A0 = {A0.str(12)} B' = {Bp.str(12)} Gamma0 = {Gamma0.str(12)}")
    for w in which:
        t1 = time.time()
        try:
            out = globals()[w]()
            print(f"{w}: PASS ({out if not isinstance(out, arb) else out.str(12)}) [{time.time()-t1:.1f}s]", flush=True)
        except (AssertionError, RuntimeError) as ex:
            print(f"{w}: FAIL: {ex} [{time.time()-t1:.1f}s]", flush=True)
