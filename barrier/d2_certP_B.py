# Proposition P certificate, Program B (mpmath.iv interval arithmetic, independent of Program A). 19 Nov 2026.
# Same claims as d2_certP_A.py, written differently: r by interval bisection on [1/2, 9/10] from scratch (sign of the concave
# root function), g expanded as h(s^2) + q h(I) - 2 s h(s) with h(x) = -x log x - (1 - x) log(1 - x) evaluated as
# -(x log x + log(1-x) - x log(1-x)).
# Usage: d2_certP_B.py <eps numerator> <eps denominator>
import sys
from mpmath import iv, mpf
iv.dps = 50
num, den = int(sys.argv[1]), int(sys.argv[2])
def rootfun(x): return -(1 - x) * iv.log(1 - x) - x * iv.log(1 + x)
a, b = mpf('0.5'), mpf('0.9')
for _ in range(160):
    c = (a + b) / 2
    v = rootfun(iv.mpf(c))
    if v.a > 0: a = c
    elif v.b < 0: b = c
    else: break
assert rootfun(iv.mpf(a)).a > 0 and rootfun(iv.mpf(b)).b < 0, 'r bracket not certified'
r = iv.mpf([a, b])
def h(x): return -(x * iv.log(x) + iv.log(1 - x) - x * iv.log(1 - x))
L = iv.log(2)
m = iv.mpf(1) / 2 + (1 - r) * iv.log(1 + r) / (2 * L); q = 2 * m - 1
tau = (iv.mpf(1) / 2 - r ** 2) / (r ** 2 * (1 - r) ** 2)
eps = iv.mpf(num) / den; s = m - eps
I = s ** 2 + tau * s ** 2 * (1 - s) ** 2
C = h(s ** 2) + q * h(I) - 2 * s * h(s) + 2 * (s - m) * h(s)
U = h(s ** 2) - h(s); Fb = (s - m) * h(s); V = h(I) - h(s)
print(f"Program B (mpmath.iv, 50 digits): eps = {num}/{den}; r in [{iv.nstr(r, 20)}]; width {float(r.delta):.1e}")
res = {}
for name, val in [('C', C), ('U', U), ('F', Fb), ('V', V)]:
    sign = 'neg' if val.b < 0 else ('pos' if val.a > 0 else 'undecided')
    res[name] = sign
    print(f"  {name} in [{float(val.a):+.12e}, {float(val.b):+.12e}]  {sign}")
ok = res['C'] == 'neg' and res['U'] == 'neg' and res['F'] == 'neg' and res['V'] == 'pos'
print('RESULT:', 'PASS' if ok else 'FAIL')
