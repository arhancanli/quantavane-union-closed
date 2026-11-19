# Proposition P certificate, Program A (python-flint / Arb balls). 19 Nov 2026.
# Claim: at eps = EPS (rational), s = m - eps, the product measure Bern(s)^n satisfies, per coordinate,
#   (C) g(s) + 2 (s - m) h(s) <= 0,  g(s) = h(s^2) + q h(I(s,s)) - 2 s h(s)   [SU's combined budget (4.1)]
#   (U) h(s^2) - h(s) <= 0                                                      [sum U <= H]
#   (F) (s - m) h(s) <= 0                                                       [fibre inequality (4.2), every fibre]
# and violates (V) h(I(s,s)) - h(s) <= 0 (printed, must be > 0). r is enclosed by a certified sign change of
# phi(x) = -(1-x) log(1-x) - x log(1+x) (concave, phi(0) = 0, phi'(0) = 1, phi(1-) = -log 2).
# Usage: d2_certP_A.py <eps numerator> <eps denominator>
import sys
from flint import arb, ctx
ctx.prec = 200
num, den = int(sys.argv[1]), int(sys.argv[2])
def phi(x): return -(1 - x) * (1 - x).log() - x * (1 + x).log()
lo, hi = arb('0.6909077387'), arb('0.6909077388')
for _ in range(120):
    mid = (lo + hi) / 2
    if phi(mid) > 0: lo = mid
    elif phi(mid) < 0: hi = mid
    else: raise SystemExit('undecided sign at midpoint')
assert phi(lo) > 0 and phi(hi) < 0, 'r bracket not certified'
r = lo.union(hi)
L = arb(2).log()
def h(x): return -x * x.log() - (1 - x) * (1 - x).log()
m = arb(1) / 2 + (1 - r) * (1 + r).log() / (2 * L); q = 2 * m - 1
tau = (arb(1) / 2 - r * r) / (r * r * (1 - r) ** 2)
eps = arb(num) / den; s = m - eps
I = s * s * (1 + tau * (1 - s) ** 2)
g = h(s * s) + q * h(I) - 2 * s * h(s)
C = g + 2 * (s - m) * h(s); U = h(s * s) - h(s); Fb = (s - m) * h(s); V = h(I) - h(s)
print(f"Program A (Arb, 200 bits): eps = {num}/{den}; r in {r.str(20)}; m = {m.str(20)}")
for name, val in [('C', C), ('U', U), ('F', Fb), ('V', V)]:
    print(f"  {name} = {val.str(15)}  upper bound {'< 0' if val < 0 else ('> 0 (lower bound positive)' if val > 0 else 'UNDECIDED')}")
ok = bool(C < 0) and bool(U < 0) and bool(Fb < 0) and bool(V > 0)
print('RESULT:', 'PASS' if ok else 'FAIL')
