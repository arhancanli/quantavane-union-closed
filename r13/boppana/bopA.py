# Program A (Arb): G(t) = h(t^2) - c t h(t) >= 0 on [0,1], c = 809/500 (default), i.e. the input of SU Lemma lem:prodred.
# Edges (analytic, constants checked here): t <= 1/100: G >= t^2 [(2-c) log(1/t) - c];
# u = 1 - t <= 1/10000: G >= u [(2-u-c) log(1/u) - (2-u) log 2 - c] >= u [(2 - 1e-4 - c) log(1e4) - 2 log 2 - c].
# Middle [1/100, 1 - 1/10000]: bisection with monotone enclosures of h (increasing on [0,1/2], decreasing on [1/2,1]).
import sys
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 128
c = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(809, 500)
A = lambda q: arb(q.numerator) / q.denominator
def hp(x):
    x = A(x) if isinstance(x, Fr) else x
    if x == 0 or x == 1: return arb(0)
    return -x * x.log() - (1 - x) * (1 - x).log()
def hbox(a, b):   # enclosure of h on [a, b], a <= b exact Fractions
    va, vb = hp(a), hp(b)
    if b <= Fr(1, 2): return va.lower(), vb.upper()
    if a >= Fr(1, 2): return vb.lower(), va.upper()
    return (va.lower() if va.lower() <= vb.lower() else vb.lower()), arb(2).log().upper()
L2 = arb(2).log(); C = A(c)
e1 = (2 - C) * arb(100).log() - C
e2 = (2 - A(Fr(1, 10000)) - C) * arb(10000).log() - 2 * L2 - C
print(f"c = {c} = {float(c)}; edge constants: left {e1.str(10)} (must be > 0), right {e2.str(10)} (must be > 0)")
ok_edges = e1 > 0 and e2 > 0
lo_, hi_ = Fr(1, 100), 1 - Fr(1, 10000)
stack = [(lo_, hi_)]; n = 0; fail = None
while stack:
    a, b = stack.pop(); n += 1
    h2lo, _ = hbox(a * a, b * b)
    _, hhi = hbox(a, b)
    G_lo = h2lo - C * A(b) * hhi
    if G_lo > 0: continue
    if b - a < Fr(1, 10**12): fail = (a, b); break
    mid = (a + b) / 2; stack += [(a, mid), (mid, b)]
print(f"middle: {'PASS' if fail is None else 'FAIL near ' + str(float(fail[0]))} ({n} intervals); edges {'PASS' if ok_edges else 'FAIL'}")
print("RESULT:", "PASS" if (fail is None and ok_edges) else "FAIL")
