# Referee check of Lemma F and of the per-coordinate identity
#   log|F| - sum_i U_i = sum_i E d(Q_i(a)Q_i(b) || Q_i(a u b))
# on random union-closed families, exact enumeration over all ordered pairs (A, B), float64.
import sys, math, random
import numpy as np

def xlogy(a, b):
    # a*log(b) with 0*log(0) = 0
    out = np.zeros_like(a)
    nz = a > 0
    out[nz] = a[nz] * np.log(b[nz])
    return out

def hbin(x):
    return -(xlogy(x, x) + xlogy(1 - x, 1 - x))

def dkl(a, b):
    # binary KL d(a||b); inf where a>0,b=0 or a<1,b=1
    out = np.zeros_like(a)
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = np.where(a > 0, a * np.log(np.where(a > 0, a, 1) / np.where(b > 0, b, 0)), 0.0)
        t2 = np.where(a < 1, (1 - a) * np.log(np.where(a < 1, 1 - a, 1) / np.where(b < 1, 1 - b, 0)), 0.0)
    return t1 + t2

def closure(gens):
    F = set(gens)
    frontier = list(F)
    while frontier:
        new = []
        cur = list(F)
        for x in frontier:
            for y in cur:
                z = x | y
                if z not in F:
                    F.add(z); new.append(z)
        frontier = new
    return sorted(F)

def analyse(F, n, perm):
    Fa = np.array(F, dtype=np.int64)
    N = len(Fa)
    A = np.repeat(Fa, N); B = np.tile(Fa, N); S = A | B
    H = math.log(N)
    sumU = 0.0; sumD = 0.0; sumJ0 = 0.0; sumD0 = 0.0
    lhs1 = 0.0; lhs2 = 0.0; lhs2_two = 0.0
    pm = 0
    for k in range(n):
        i = perm[k]
        absent = ((Fa >> i) & 1) == 0
        pre = Fa & pm
        keys, inv = np.unique(pre, return_inverse=True)
        cnt = np.bincount(inv).astype(float)
        cab = np.bincount(inv, weights=absent.astype(float))
        Qk = cab / cnt                       # Q_i for each prefix key
        def Q(x):
            idx = np.searchsorted(keys, x & pm)
            assert np.all(keys[idx] == (x & pm))
            return Qk[idx]
        qa, qb, qc = Q(A), Q(B), Q(S)
        P = qa * qb
        U = hbin(P).mean(); sumU += U
        dd = dkl(P, qc)
        assert np.all(np.isfinite(dd)), "infinite KL: union prefix inconsistent"
        sumD += dd.mean()
        # J0 / D0 split: pi(c) = E[P | c]
        ckeys, cinv = np.unique(S & pm, return_inverse=True)
        pic = np.bincount(cinv, weights=P) / np.bincount(cinv)
        pi_at = pic[cinv]
        sumJ0 += dkl(P, pi_at).mean(); sumD0 += dkl(pi_at, qc).mean()
        forcedA = qa == 0
        lhs1 += np.where(forcedA, qc, 0.0).mean()
        lhs2 += np.where(forcedA, -np.log1p(-qc), 0.0).mean()
        lhs2_two += np.where(P == 0, -np.log1p(-qc), 0.0).mean()
        pm |= (1 << i)
    return dict(H=H, sumU=sumU, sumD=sumD, J0=sumJ0, D0=sumD0, lhs1=lhs1, lhs2=lhs2, lhs2_two=lhs2_two, N=len(F))

def main(seed, trials):
    rng = random.Random(seed)
    worst_id = 0.0; viol = 0; viol_two = 0; maxratio = 0.0; nfam = 0; worst_split = 0.0
    for t in range(trials):
        n = rng.randint(3, 8)
        k = rng.randint(2, 9)
        dens = rng.choice([0.2, 0.35, 0.5, 0.65, 0.8])
        gens = [sum(1 << j for j in range(n) if rng.random() < dens) for _ in range(k)]
        if rng.random() < 0.3:
            gens.append(0)                    # sometimes include the empty set
        F = closure(gens)
        if len(F) < 2 or len(F) > 300:
            continue
        perm = list(range(n)); rng.shuffle(perm)
        R = analyse(F, n, perm)
        nfam += 1
        resid = abs((R['H'] - R['sumU']) - R['sumD'])
        worst_id = max(worst_id, resid)
        worst_split = max(worst_split, abs(R['J0'] + R['D0'] - R['sumD']))
        tol = 1e-10
        if not (R['lhs1'] <= R['lhs2'] + tol and R['lhs2'] <= R['H'] - R['sumU'] + tol):
            viol += 1
        if not (R['lhs2_two'] <= R['H'] - R['sumU'] + tol):
            viol_two += 1
        if R['H'] - R['sumU'] > 1e-12:
            maxratio = max(maxratio, R['lhs2'] / (R['H'] - R['sumU']))
    print(f"seed {seed}: families {nfam}; max |(log|F| - sum U) - sum E d| = {worst_id:.2e}; "
          f"max |J0+D0 - sum E d| = {worst_split:.2e}; Lemma F violations {viol}; two-sided version violations {viol_two}; "
          f"max lhs2/(H - sum U) = {maxratio:.4f}")

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
