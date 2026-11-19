# Referee check of the non-abelian translate-family claims (NOTES (c), "Proved here").
# F_L(R) = {A R : A subset of G}. Exact integer arithmetic on bitmasks.
import sys, itertools, random
import numpy as np
from fractions import Fraction

# ---------------- groups as multiplication tables ----------------
def from_mult(elements, mul):
    idx = {g: i for i, g in enumerate(elements)}
    n = len(elements)
    T = np.zeros((n, n), dtype=np.int64)
    for a in elements:
        for b in elements:
            T[idx[a], idx[b]] = idx[mul(a, b)]
    return T

def closure_perm(gens):
    gens = [tuple(g) for g in gens]
    ident = tuple(range(len(gens[0])))
    S = {ident}; frontier = [ident]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = tuple(x[g[i]] for i in range(len(g)))   # x o g
                if y not in S:
                    S.add(y); new.append(y)
        frontier = new
    els = sorted(S)
    return els, (lambda a, b: tuple(a[b[i]] for i in range(len(b))))

def dihedral(k):
    rot = [(i + 1) % k for i in range(k)]; ref = [(-i) % k for i in range(k)]
    return closure_perm([rot, ref])

def dicyclic(mm):
    els = [(a, x) for a in range(2 * mm) for x in range(2)]
    def mul(u, v):
        (k1, j1), (k2, j2) = u, v
        return ((k1 + (k2 if j1 == 0 else -k2) + (mm if (j1 and j2) else 0)) % (2 * mm), j1 ^ j2)
    return els, mul

def alt4():
    els = [p for p in itertools.permutations(range(4)) if sum(1 for i in range(4) for j in range(i) if p[j] > p[i]) % 2 == 0]
    return els, (lambda a, b: tuple(a[b[i]] for i in range(4)))

def sym(k):
    els = list(itertools.permutations(range(k)))
    return els, (lambda a, b: tuple(a[b[i]] for i in range(k)))

def z7z3():
    els = [(a, b) for a in range(7) for b in range(3)]
    return els, (lambda u, v: ((u[0] + pow(2, u[1], 7) * v[0]) % 7, (u[1] + v[1]) % 3))

def cyclic(k):
    return list(range(k)), (lambda a, b: (a + b) % k)

# ---------------- group utilities ----------------
def basics(T):
    n = T.shape[0]
    e = int([i for i in range(n) if all(T[i, j] == j for j in range(n))][0])
    inv = np.array([int(np.where(T[i] == e)[0][0]) for i in range(n)])
    return n, e, inv

def automorphisms(T):
    n, e, inv = basics(T)
    order = []
    for g in range(n):
        x, k = g, 1
        while x != e:
            x = T[x, g]; k += 1
        order.append(k)
    # find a small generating set greedily
    def gen_closure(gs):
        S = {e}; fr = [e]
        while fr:
            new = []
            for x in fr:
                for g in gs:
                    y = int(T[x, g])
                    if y not in S: S.add(y); new.append(y)
            fr = new
        return S
    gens = []
    while len(gen_closure(gens)) < n:
        best = max(range(n), key=lambda g: len(gen_closure(gens + [g])))
        gens.append(best)
    auts = []
    for imgs in itertools.product(range(n), repeat=len(gens)):
        if any(order[imgs[i]] != order[gens[i]] for i in range(len(gens))):
            continue
        mp_ = {e: e}; fr = [e]; ok = True
        while fr and ok:
            new = []
            for x in fr:
                for g, y in zip(gens, imgs):
                    xg = int(T[x, g]); im = int(T[mp_[x], y])
                    if xg in mp_:
                        if mp_[xg] != im: ok = False; break
                    else:
                        mp_[xg] = im; new.append(xg)
                if not ok: break
            fr = new
        if not ok or len(mp_) != n or len(set(mp_.values())) != n:
            continue
        a = np.array([mp_[i] for i in range(n)])
        if np.all(a[T] == T[np.ix_(a, a)]):
            auts.append(a)
    return auts

def twisted_translation_perms(T, auts):
    # all permutations x -> u * beta(x) * v  (beta in Aut, u, v in G); deduplicated
    n = T.shape[0]
    P = set()
    for a in auts:
        for u in range(n):
            ua = T[u][a]                   # u*beta(x)
            for v in range(n):
                P.add(tuple(T[ua, v]))
    return np.array(sorted(P), dtype=np.int64)

# ---------------- family computations ----------------
POP = np.bitwise_count

def family(trmask, n):
    # all unions of the sets trmask[x] (x in G): array over all A of OR_{x in A} trmask[x]; return unique
    arr = np.zeros(1, dtype=np.int64)
    for x in range(n):
        arr = np.concatenate([arr, arr | np.int64(trmask[x])])
    return np.unique(arr)

def check_R(T, n, e, inv, Rset, perms):
    Rl = sorted(Rset); Rinv = sorted({int(inv[x]) for x in Rl})
    trR = [sum(1 << int(T[x, rr]) for rr in Rl) for x in range(n)]        # mask of xR
    trRi = [sum(1 << int(T[x, rr]) for rr in Rinv) for x in range(n)]     # mask of xR^{-1}
    FL = family(trR, n); FLi = family(trRi, n)
    full = (1 << n) - 1
    # Int(S) = {x : xR subset S}
    Int = np.zeros_like(FL)
    for x in range(n):
        Int |= np.where((FL & trR[x]) == trR[x], np.int64(1 << x), np.int64(0))
    phi = full & ~Int
    # (G \ S) R^{-1}
    comp = full & ~FL
    alt = np.zeros_like(FL)
    for x in range(n):
        alt |= np.where((comp >> x) & 1 == 1, np.int64(trRi[x]), np.int64(0))
    res = {}
    res['formula'] = bool(np.all(alt == phi))
    res['into'] = bool(np.all(np.isin(phi, FLi)))
    res['inj'] = len(np.unique(phi)) == len(FL)
    res['bij'] = res['into'] and res['inj'] and len(FL) == len(FLi)
    res['sizeineq'] = bool(np.all(POP(FL) + POP(phi) >= n))
    sumL = int(POP(FL).sum()); sumLi = int(POP(FLi).sum()); N = len(FL)
    res['avgsum'] = (sumL + sumLi) >= n * N if len(FL) == len(FLi) else None
    # transitivity: all element frequencies equal
    freqs = [int(((FL >> x) & 1).sum()) for x in range(n)]
    res['transitive'] = len(set(freqs)) == 1
    res['frankl'] = 2 * freqs[0] >= N
    # equivalence: #{T in F_L(R^{-1}) : R^{-1} subset T} == #{S in F_L(R): e not in S}
    Rimask = sum(1 << x for x in Rinv)
    cntT = int(((FLi & Rimask) == Rimask).sum()); cntS = int(N - freqs[e])
    res['equiv'] = (cntT == cntS)
    # criterion R^{-1} = alpha(g R h)
    imgs = (np.int64(1) << perms[:, Rl]).sum(axis=1)
    res['crit'] = bool(np.any(imgs == Rimask))
    res['avg_eq'] = (sumL * len(FLi) == sumLi * N)
    res['ratio'] = Fraction(sumL, N * n)
    res['coset'] = None
    res['e_in_R'] = e in Rset
    return res

def is_left_coset_of_subgroup(T, e, inv, Rset):
    # R = g H for a subgroup H  <=>  g^{-1}R is a subgroup for some (any) g in R
    g = next(iter(Rset))
    H = {int(T[inv[g], x]) for x in Rset}
    return all(int(T[a, b]) in H for a in H for b in H)

def run(name, els_mul, mode, nsample=None, seed=0):
    els, mul = els_mul
    T = from_mult(els, mul)
    n, e, inv = basics(T)
    auts = automorphisms(T)
    perms = twisted_translation_perms(T, auts)
    abelian = bool(np.all(T == T.T))
    rng = random.Random(seed)
    if mode == 'all':
        Rs = [frozenset(x for x in range(n) if (mask >> x) & 1) for mask in range(1, 1 << n)]
    else:
        Rs = []
        while len(Rs) < nsample:
            k = rng.randint(1, n - 1)
            Rs.append(frozenset(rng.sample(range(n), k)))
    agg = dict(total=0, fail=[], noE=0, crit=0, avg_neq=0, crit_avg_neq=0, frankl_fail=0, nontrans=0)
    minr_all = None; minr_noncrit = None; minr_noncoset = None
    for Rset in Rs:
        res = check_R(T, n, e, inv, Rset, perms)
        agg['total'] += 1
        for key in ('formula', 'into', 'inj', 'bij', 'sizeineq', 'avgsum', 'equiv'):
            if not res[key]:
                agg['fail'].append((key, sorted(Rset)))
        if not res['transitive']: agg['nontrans'] += 1
        if not res['e_in_R']: agg['noE'] += 1
        if res['crit']: agg['crit'] += 1
        if not res['avg_eq']: agg['avg_neq'] += 1
        if res['crit'] and not res['avg_eq']: agg['crit_avg_neq'] += 1
        if not res['frankl']: agg['frankl_fail'] += 1
        r_ = res['ratio']
        minr_all = r_ if minr_all is None else min(minr_all, r_)
        if not res['crit']:
            minr_noncrit = r_ if minr_noncrit is None else min(minr_noncrit, r_)
        if not is_left_coset_of_subgroup(T, e, inv, Rset):
            minr_noncoset = r_ if minr_noncoset is None else min(minr_noncoset, r_)
    print(f"{name}: |G|={n} abelian={abelian} |Aut|={len(auts)} twisted-translation perms={len(perms)} mode={mode}")
    print(f"   R tested {agg['total']} (e not in R: {agg['noE']}); step failures: {len(agg['fail'])} {agg['fail'][:3]}")
    print(f"   non-transitive frequency vectors: {agg['nontrans']}; Frankl failures: {agg['frankl_fail']}")
    print(f"   R in class R^-1 = alpha(gRh): {agg['crit']}; avg_L(R) != avg_L(R^-1): {agg['avg_neq']}; in class but avg differ: {agg['crit_avg_neq']}")
    print(f"   min avg/|G|: all R {minr_all} = {float(minr_all):.4f}; R outside class: "
          f"{(minr_noncrit, float(minr_noncrit)) if minr_noncrit is not None else None}; non-coset R: "
          f"{(minr_noncoset, float(minr_noncoset)) if minr_noncoset is not None else None}")
    sys.stdout.flush()

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'small':
        run('Z6 (abelian control)', cyclic(6), 'all')
        run('S3', sym(3), 'all')
        run('D4 (order 8)', dihedral(4), 'all')
        run('Q8', dicyclic(2), 'all')
        run('D5 (order 10)', dihedral(5), 'all')
    elif which == 'twelve':
        run('A4', alt4(), 'all')
        run('Dic3 (order 12)', dicyclic(3), 'all')
        run('D6 (order 12)', dihedral(6), 'all')
    elif which == 'big':
        run('Q16 = Dic4', dicyclic(4), 'sample', nsample=300, seed=1)
        run('D8 (order 16)', dihedral(8), 'sample', nsample=300, seed=2)
        run('Z7:Z3', z7z3(), 'sample', nsample=40, seed=3)
