# Lemma F, numerical check B (written separately from check A, referee's lemma_f.py). 21 Nov 2026.
# For random small union-closed families (all ordered pairs (A, B) enumerated exactly), conditional absence probabilities
# Q_i(x) are exact rationals (fractions.Fraction); logarithms are evaluated with mpmath at 40 significant digits.
# Checks, for a random coordinate order:
#   (I)  log|F| - sum_i U_i  ==  sum_i E d(Q_i(a) Q_i(b) || Q_i(a u b))         (identity; |difference| < 1e-30)
#   (F1) sum_i E[Q_i(c) 1{Q_i(a) Q_i(b) = 0}] <= sum_i E[-log(1 - Q_i(c)) 1{Q_i(a) Q_i(b) = 0}] <= log|F| - sum_i U_i
# with U_i = E h(Q_i(a) Q_i(b)), c = a u b (presence-union of the prefix patterns), d = binary relative entropy.
# Red control (--red): the same computation on families that are NOT union-closed; there A u B can leave F, the identity has
# no reason to hold, and the program must report the identity or (F1) failing on some family.
import sys, random
from fractions import Fraction
from mpmath import mp, mpf, log
mp.dps = 40

def M(x):
    # exact rational -> mpf (40 digits)
    if isinstance(x, Fraction): return mpf(x.numerator) / x.denominator
    return mpf(x)

def h(p):
    p = M(p)
    if p == 0 or p == 1: return mpf(0)
    return -p * log(p) - (1 - p) * log(1 - p)

def d(p, q):
    p, q = M(p), M(q)
    t = mpf(0)
    if p > 0:
        if q == 0: return mpf('inf')            # impossible for union-closed families (A u B lies in F)
        t += p * log(p / q)
    if p < 1:
        if q == 1: return mpf('inf')
        t += (1 - p) * log((1 - p) / (1 - q))
    return t

def closure(gens):
    F = set(gens) | {0}
    changed = True
    while changed:
        changed = False
        for a in list(F):
            for b in list(F):
                if a | b not in F: F.add(a | b); changed = True
    return sorted(F)

PERTURB = None   # red control 2: multiply the union-prefix conditional Q(c) by this factor (identity must then fail)
def check(F, order):
    N = len(F); Fset = set(F)
    pos = {e: k for k, e in enumerate(order)}
    def prefix(X, k): return frozenset(e for e in order[:k] if (X >> e) & 1)
    H = log(mpf(N))
    sumU = mpf(0); sumD = mpf(0); lhs1 = mpf(0); lhs2 = mpf(0); escaped = 0
    for k, i in enumerate(order):
        cnt = {}; absn = {}
        for X in F:
            p = prefix(X, k); cnt[p] = cnt.get(p, 0) + 1
            if not (X >> i) & 1: absn[p] = absn.get(p, 0) + 1
        Q = {p: Fraction(absn.get(p, 0), cnt[p]) for p in cnt}
        for A in F:
            a = prefix(A, k)
            for B in F:
                b = prefix(B, k); c = a | b
                if c not in Q: escaped += 1; continue          # union prefix not realised (only possible without closure)
                P = Q[a] * Q[b]
                Qc = Q[c] if PERTURB is None else min(Q[c] * PERTURB, Fraction(1))
                sumU += h(P) / (N * N)
                sumD += d(P, Qc) / (N * N)
                if P == 0:
                    lhs1 += M(Q[c]) / (N * N)
                    lhs2 += (-log(1 - M(Q[c])) if Q[c] < 1 else mpf('inf')) / (N * N)
    return H - sumU, sumD, lhs1, lhs2, escaped

if __name__ == '__main__':
    red = '--red' in sys.argv
    if '--red2' in sys.argv: red = False; PERTURB = Fraction(999, 1000)
    random.seed(20261121 + (7 if red else 0))
    worst_id = mpf(0); bad = 0; tested = 0
    for t in range(300 if not red else 200):
        n = random.randint(3, 6)
        gens = [random.getrandbits(n) for _ in range(random.randint(2, 6))]
        F = closure(gens) if not red else sorted(set(gens) | {0})
        if red and closure(F) == F: continue                   # red control uses only non-union-closed families
        if len(F) < 2: continue
        order = list(range(n)); random.shuffle(order)
        gap, sumD, l1, l2, esc = check(F, order)
        tested += 1
        idf = abs(gap - sumD) if esc == 0 else mpf('inf')
        ok = (idf < mpf('1e-30')) and (l1 <= l2 + mpf('1e-30')) and (l2 <= gap + mpf('1e-30'))
        if not ok: bad += 1
        if esc == 0: worst_id = max(worst_id, idf)
    tag = 'RED CONTROL 2 (Q at union prefixes x 0.999)' if PERTURB is not None else ('RED CONTROL (non-union-closed families)' if red else 'union-closed families')
    print(f"{tag}: tested {tested}, "
          f"failures {bad}, max |identity residual| {mp.nstr(worst_id, 5)}")
    if red or PERTURB is not None:
        print('RESULT:', 'FAIL (as a red control must)' if bad > 0 else 'PASS (red control did NOT fail: investigate)')
    else:
        print('RESULT:', 'PASS' if bad == 0 else 'FAIL')
