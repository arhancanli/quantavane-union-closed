# Exhaustive |R| in {3, 4} (R containing the identity) for the non-abelian groups of c_groups.py: min over R of avg|S|/|G| for
# F_L(R), excluding the exact-coset cases (avg = |G|/2 exactly), and the min of avg_L(R) + avg_L(R^{-1}) - |G| (the provable part).
import sys, itertools
from fractions import Fraction
exec(open("c_groups.py").read().split("\nif __name__")[0])
MAXF = int(sys.argv[1]) if len(sys.argv) > 1 else 600_000
sizes = [int(a) for a in sys.argv[2].split(',')] if len(sys.argv) > 2 else [3, 4]
import os
GF = os.environ.get("GROUPS", "").split(",")
for name, (n, mul) in groups().items():
    if GF != [""] and name not in GF: continue
    inv = [next(b for b in range(n) if mul[a][b] == 0) for a in range(n)]
    for k in sizes:
        best = (Fraction(10), None); bestsum = (Fraction(10**6), None); tested = skipped = 0; nontriv_min = (Fraction(10), None)
        for rest in itertools.combinations(range(1, n), k - 1):
            Rm = 1 | sum(1 << x for x in rest)
            FL = family(mul, n, Rm, 'L')
            if FL is None: skipped += 1; continue
            Rinv = sum(1 << inv[x] for x in range(n) if (Rm >> x) & 1)
            FLi = family(mul, n, Rinv, 'L')
            if FLi is None: skipped += 1; continue
            tested += 1
            a, ai = avg_size(FL), avg_size(FLi)
            if a / n < best[0]: best = (a / n, (rest, len(FL)))
            if a / n > Fraction(1, 2) and a / n < nontriv_min[0]: nontriv_min = (a / n, (rest, len(FL)))
            if a + ai - n < bestsum[0]: bestsum = (a + ai - n, rest)
        print(f"{name} |G|={n} |R|={k}: tested {tested}, skipped {skipped}; min avg/|G| = {float(best[0]):.6f} at {best[1]}; "
              f"min above 1/2 = {float(nontriv_min[0]):.6f} at {nontriv_min[1]}; min avg_L(R)+avg_L(R^-1)-|G| = {float(bestsum[0]):.4f}",
              flush=True)
