# Tier 1 lane, union-closed, direction (c) (exact integer computation). 19 Nov 2026, 22:2x +04.
# For a finite group G (given by its multiplication table) and R subset of G, the families
#   F_L(R) = { A R : A subset G }  (unions of left translates gR),   F_R(R) = { R B } (unions of right translates Rg).
# G acts transitively on F_L(R) by left multiplication, so every element has the same abundance avg|S| / |G|, and Frankl for
# F_L(R) is equivalent to avg_{S in F_L(R)} |S| >= |G|/2. Aaronson-Ellis-Leader (EJC 28(2) 2021, arXiv 2010.08795) prove this for
# ABELIAN G. Here: exhaustive or sampled check for non-abelian G; record min of avg|S| / |G| and whether avg_L = avg_R.
import sys, itertools, random
from fractions import Fraction

def perm_group_table(gens, deg):
    # closure of permutation generators -> elements as tuples; return list and multiplication (composition) table
    idp = tuple(range(deg)); elems = [idp]; idx = {idp: 0}; frontier = [idp]
    def comp(p, q): return tuple(p[q[i]] for i in range(deg))   # (p*q)(i) = p(q(i))
    while frontier:
        nf = []
        for x in frontier:
            for g in gens:
                y = comp(x, g)
                if y not in idx: idx[y] = len(elems); elems.append(y); nf.append(y)
        frontier = nf
    n = len(elems)
    mul = [[idx[comp(elems[a], elems[b])] for b in range(n)] for a in range(n)]
    return n, mul

def groups():
    G = {}
    def cyc(n): return tuple((i + 1) % n for i in range(n))
    def refl(n): return tuple((-i) % n for i in range(n))
    for n in [7, 8, 9, 10, 11, 12]:
        G[f'D{n}'] = perm_group_table([cyc(n), refl(n)], n)            # dihedral of order 2n
    # S4 (24), A4 (12)
    G['S4'] = perm_group_table([(1, 2, 3, 0), (1, 0, 2, 3)], 4)
    # Z7 x| Z3 (order 21): x -> 2x mod 7 has order 3, acting on Z7 as affine maps
    G['F21'] = perm_group_table([tuple((i + 1) % 7 for i in range(7)), tuple((2 * i) % 7 for i in range(7))], 7)
    # Z5 x| Z4 (order 20): affine group AGL(1,5)
    G['F20'] = perm_group_table([tuple((i + 1) % 5 for i in range(5)), tuple((2 * i) % 5 for i in range(5))], 5)
    # quaternion-like order 16 groups via permutations: generalized quaternion Q16 as regular rep is heavy; use SL(2,3) (24)
    # SL(2,3) acting on the 8 nonzero vectors of F3^2
    vecs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
    def mat(M): return tuple(vecs.index(((M[0] * v[0] + M[1] * v[1]) % 3, (M[2] * v[0] + M[3] * v[1]) % 3)) for v in vecs)
    G['SL23'] = perm_group_table([mat((1, 1, 0, 1)), mat((1, 0, 1, 1))], 8)
    # Heisenberg group mod 3 (order 27) acting on F3^3 affine? use upper unitriangular 3x3 over F3 acting on 9 vectors (x, y) -> regular rep
    els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    def hmul(p, q): return ((p[0] + q[0]) % 3, (p[1] + q[1]) % 3, (p[2] + q[2] + p[0] * q[1]) % 3)
    idx = {e: i for i, e in enumerate(els)}
    G['Heis27'] = (27, [[idx[hmul(els[a], els[b])] for b in range(27)] for a in range(27)])
    G['D13'] = perm_group_table([tuple((i + 1) % 13 for i in range(13)), tuple((-i) % 13 for i in range(13))], 13)
    # Q16 (generalized quaternion) as pairs (k, e) with x^8 = 1, y^2 = x^4, y x y^-1 = x^-1: elements x^k y^e
    def qmul(p, q):
        k1, e1 = p; k2, e2 = q
        if e1 == 0: return ((k1 + k2) % 8, e2)
        # y x^k2 = x^-k2 y ; y y = x^4
        k = (k1 - k2) % 8
        return ((k + 4) % 8, 0) if e2 == 1 else (k, 1)
    qe = [(k, e) for e in range(2) for k in range(8)]; qi = {v: i for i, v in enumerate(qe)}
    G['Q16'] = (16, [[qi[qmul(qe[a], qe[b])] for b in range(16)] for a in range(16)])
    return G

def family(mul, n, Rmask, side):
    # translates: left gR = {g r}, right Rg = {r g}
    trans = set()
    for g in range(n):
        m = 0
        for r in range(n):
            if (Rmask >> r) & 1:
                x = mul[g][r] if side == 'L' else mul[r][g]
                m |= 1 << x
        trans.add(m)
    F = {0}
    for t in trans:
        F |= {S | t for S in F}
        if len(F) > MAXF: return None
    return F

def avg_size(F):
    return Fraction(sum(bin(S).count('1') for S in F), len(F))

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    maxF = int(sys.argv[2]) if len(sys.argv) > 2 else 300_000
    MAXF = maxF; NSAMP = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    for name, (n, mul) in groups().items():
        ab = all(mul[a][b] == mul[b][a] for a in range(n) for b in range(n))
        worst = (Fraction(10), None); neqLR = 0; tested = 0; skipped = 0
        def rnd():
            d = random.uniform(0.15, 0.75)
            return 1 | sum(1 << i for i in range(1, n) if random.random() < d)
        Rs = (rnd() for _ in range(NSAMP))
        for Rm in Rs:
            FL = family(mul, n, Rm, 'L')
            if FL is None: skipped += 1; continue
            FR = family(mul, n, Rm, 'R')
            if FR is None: skipped += 1; continue
            aL, aR = avg_size(FL), avg_size(FR)
            tested += 1
            if aL != aR: neqLR += 1
            for a, F in [(aL, FL), (aR, FR)]:
                if a / n < worst[0]: worst = (a / n, (bin(Rm), len(F)))
        print(f"{name}: |G|={n} abelian={ab} tested R={tested} skipped(>{maxF})={skipped}; min avg|S|/|G| = "
              f"{float(worst[0]):.6f} at R={worst[1]}; #R with avg_L != avg_R: {neqLR}", flush=True)
