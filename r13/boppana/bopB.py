# Program B (written separately): binary64 intervals with one-ulp outward widening, h values from mpmath.iv (60 digits).
# Same statement as bopA.py; breadth-first queue; own edge-constant evaluation.
import sys, math
from collections import deque
from fractions import Fraction
from mpmath import iv
iv.dps = 60
dn = lambda v: math.nextafter(v, -math.inf); up = lambda v: math.nextafter(v, math.inf)
cq = Fraction(sys.argv[1]) if len(sys.argv) > 1 else Fraction(809, 500)
civ = iv.mpf(cq.numerator) / cq.denominator
def hval(x):   # x float in [0,1] -> (lo, hi)
    if x <= 0.0 or x >= 1.0: return 0.0, 0.0
    X = iv.mpf(x); v = -X * iv.log(X) - (1 - X) * iv.log(1 - X)
    return dn(float(v.a)), up(float(v.b))
def henc(a, b):
    (la, ua), (lb, ub) = hval(a), hval(b)
    if b <= 0.5: return la, ub
    if a >= 0.5: return lb, ua
    return min(la, lb), up(math.log(2))
edgeL = (2 - civ) * iv.log(100) - civ
edgeR = (2 - iv.mpf(1) / 10000 - civ) * iv.log(10000) - 2 * iv.log(2) - civ
print(f"c = {cq}; edge constants: left [{float(edgeL.a):.10f}, ...], right [{float(edgeR.a):.10f}, ...]")
ok_edges = edgeL.a > 0 and edgeR.a > 0
cl, ch = dn(float(civ.a)), up(float(civ.b))
q = deque([(0.01, 0.9999)]); n = 0; fail = None   # 0.01 and 0.9999 are not exact: widen the middle outward
q = deque([(dn(0.01), up(0.9999))])
while q:
    a, b = q.popleft(); n += 1
    a2, b2 = dn(a * a), up(b * b)
    h2lo = henc(max(a2, 0.0), min(b2, 1.0))[0]
    hhi = henc(a, b)[1]
    Ghi_sub = up(up(ch * b) * hhi)          # upper bound of c t h(t)
    if dn(h2lo - Ghi_sub) > 0: continue
    if b - a < 1e-12: fail = a; break
    md = 0.5 * (a + b); q.append((a, md)); q.append((md, b))
    if n > 5_000_000: fail = -1.0; break
print(f"middle: {'PASS' if fail is None else 'FAIL near ' + repr(fail)} ({n} intervals); edges {'PASS' if ok_edges else 'FAIL'}")
print("RESULT:", "PASS" if (fail is None and ok_edges) else "FAIL")
