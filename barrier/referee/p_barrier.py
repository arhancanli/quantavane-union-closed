# Referee's independent interval check of Proposition P (product-measure barrier).
# Arb balls (python-flint), own certified bracket for r. Nothing from the lab's scripts is used.
import sys
from flint import arb, ctx
ctx.prec = 256

def h(x):
    # binary entropy in nats on an arb ball strictly inside (0,1)
    return -x * x.log() - (1 - x) * (1 - x).log()

def rootfun(r):
    return -(1 - r) * (1 - r).log() - r * (1 + r).log()

# --- certified bracket for r: bisection on exact rationals, signs certified by arb
lo, hi = arb(6) / 10, arb(8) / 10
assert rootfun(lo) > 0 and rootfun(hi) < 0
lo_q, hi_q = 0.6, 0.8
from fractions import Fraction
lo_f, hi_f = Fraction(6, 10), Fraction(8, 10)
def A(fr):
    return arb(fr.numerator) / arb(fr.denominator)
for _ in range(200):
    mid = (lo_f + hi_f) / 2
    v = rootfun(A(mid))
    if v > 0:
        lo_f = mid
    elif v < 0:
        hi_f = mid
    else:
        break
assert rootfun(A(lo_f)) > 0 and rootfun(A(hi_f)) < 0
r = A(lo_f).union(A(hi_f))
print("r ball:", r.str(25), " width ~", float(hi_f - lo_f))
L = arb(2).log()
m = arb(1) / 2 + (1 - r) * (1 + r).log() / (2 * L)
q = 2 * m - 1
tau = (arb(1) / 2 - r * r) / (r * r * (1 - r) ** 2)
cstar = 1 - m
print("m =", m.str(22), " c* =", cstar.str(22), " q =", q.str(15), " tau =", tau.str(15))

def I(x, y):
    return x * y * (1 + tau * (1 - x) * (1 - y))

def quantities(eps):
    s = m - eps
    e = h(s)
    U = h(s * s) - e                 # per-coordinate U-budget term U - e
    V = h(I(s, s)) - e               # per-coordinate V-budget term V - e
    g = h(s * s) + q * h(I(s, s)) - 2 * s * e
    C = g + 2 * (s - m) * e          # combined budget per coordinate
    C2 = h(s * s) + q * h(I(s, s)) - (1 + q) * e   # same thing written as U + qV - (1+q)e
    F = (s - m) * e                  # fibre inequality (4.2) term
    return dict(s=s, e=e, g=g, C=C, C2=C2, U=U, V=V, F=F)

def show(name, epsfr):
    eps = A(epsfr)
    Q = quantities(eps)
    print(f"--- eps = {epsfr} ({float(epsfr):.7e}) ---")
    for k in ("C", "C2", "U", "F", "V", "g"):
        x = Q[k]
        sgn = "NEG" if x < 0 else ("POS" if x > 0 else "UNDECIDED")
        print(f"  {k:3s} = {x.str(14)}  {sgn}")
    print("  abundance 1-s =", (1 - Q["s"]).str(14), " c*+eps =", (cstar + eps).str(14))
    return Q

show("main", Fraction(1205, 10**6))
show("red1", Fraction(119, 10**5))
show("red2", Fraction(12044, 10**7))
show("lab 121e-5", Fraction(121, 10**5))
show("eps 12045e-7", Fraction(12045, 10**7))

# --- eps_P by bisection on the sign of C
a_f, b_f = Fraction(12044, 10**7), Fraction(1205, 10**6)
assert quantities(A(a_f))["C"] > 0 and quantities(A(b_f))["C"] < 0
for _ in range(80):
    mid = (a_f + b_f) / 2
    c = quantities(A(mid))["C"]
    if c > 0:
        a_f = mid
    elif c < 0:
        b_f = mid
    else:
        break
print("eps_P in [", float(a_f), ",", float(b_f), "]  ->  c* + eps_P =", (cstar + A(a_f)).str(16))

# --- monotonicity: dC/ds on the bracket [12044e-7, 1205e-6] (s = m - eps): C decreasing in eps <=> dC/ds > 0
def hp(x):
    return ((1 - x) / x).log()
def dC_ds(s):
    Ip = 2 * s * (1 + tau * (1 - s) * (1 - 2 * s))   # d/ds I(s,s)
    return 2 * s * hp(s * s) + q * hp(I(s, s)) * Ip - (1 + q) * hp(s)
epsbox = A(Fraction(12044, 10**7)).union(A(Fraction(1205, 10**6)))
sbox = m - epsbox
d = dC_ds(sbox)
print("dC/ds on s-box of [12044e-7,1205e-6]:", d.str(10), " (positive => C strictly decreasing in eps there, root unique)")

# --- C < 0 on [1205e-6, m - 1/2] (s in [1/2, m - 1205e-6]): interval bisection in s
def C_of_s(s):
    return h(s * s) + q * h(I(s, s)) - (1 + q) * h(s)
stack = [(Fraction(1, 2), None)]
# work with arb endpoints: s from 1/2 to m - 1205e-6 (upper end is an arb; use rational just below)
s_hi_ball = m - A(Fraction(1205, 10**6))
# choose rational s_hi with s_hi >= upper end of the ball, cover [1/2, s_hi] but treat last box specially
s_hi = Fraction(int(float(s_hi_ball.mid()) * 10**12) + 1, 10**12)
boxes = [(Fraction(1, 2), s_hi)]
ok, nbox = True, 0
while boxes:
    a, b = boxes.pop()
    nbox += 1
    sb = A(a).union(A(b))
    # only the part with s <= m - 1205e-6 matters; points above it have eps < 1205e-6
    val = C_of_s(sb)
    if val < 0:
        continue
    if b - a < Fraction(1, 10**11):
        # box straddles s = m - 1205e-6 only if it contains it; check its lower point
        if A(a) > s_hi_ball:  # entirely beyond the range of interest
            continue
        ok = False
        print("  undecided box", float(a), float(b), val.str(8))
        break
    mid = (a + b) / 2
    boxes.append((a, mid)); boxes.append((mid, b))
print(f"C < 0 for every s in [1/2, m - 1205e-6] (eps in [1205e-6, m-1/2]): {'PASS' if ok else 'FAIL'} ({nbox} boxes)")

# --- U-budget and V-budget sign over the same range (for the 'for eps >= eps_P' statement)
def U_of_s(s):
    return h(s * s) - h(s)
boxes = [(Fraction(1, 2), s_hi)]; okU = True; nb = 0
while boxes:
    a, b = boxes.pop(); nb += 1
    v = U_of_s(A(a).union(A(b)))
    if v < 0: continue
    if b - a < Fraction(1, 10**11):
        if A(a) > s_hi_ball: continue
        okU = False; print("  U undecided", float(a), float(b)); break
    mid = (a + b) / 2; boxes += [(a, mid), (mid, b)]
print(f"h(s^2) < h(s) for s in [1/2, m - 1205e-6]: {'PASS' if okU else 'FAIL'} ({nb} boxes)")
print("golden-ratio point where h(s^2) = h(s): s = (sqrt5-1)/2 =", ((arb(5).sqrt() - 1) / 2).str(15), "; m - that =", (m - (arb(5).sqrt() - 1) / 2).str(10))
