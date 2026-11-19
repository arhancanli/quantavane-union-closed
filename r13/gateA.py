# Round 13 gate, Program A (FLOAT; NOT a certificate). Quantavane union-closed lane, 18 Nov 2026.
#
# Route JB: SU's scalar chain with every loss written PER LAW, all in one uniformly random compatible order,
# against ONE budget beta = g + 2 (s - m + eps) e (sum over all coordinates <= 2 eps H in every order).
# Dual certificate: find theta >= 0, mu >= 0 with
#     G0 = eta t0 - mu (1 - a*) - max(psiG(theta), psiT(theta)) - 2 eps theta > 0,
#     psiG(theta) = sup over laws with e >= z  of (lG - theta beta) / e,
#     psiT(theta) = sup over laws with e <= z  of (lT - theta beta) / e,
# lG = (A0+B') x_+ + (A0-B') x_- + max(0, -sigma' + C Rbar) + mu (a* e + d0 x - U),   x = e - e*,
#      sigma' = E[(X-m)h(X)] + D e - b_m m,
# lT = m e - Phi - mu (U - e) + (eta t0 - mu (1 - a*)) e   [Phi = inf over decompositions of E p e; 'SU' mode: Phi = 0].
# Laws: k atoms; CIID couplings with Zc latent components (matrix M[z, j]); decompositions with Zd components.
# The sup is a nonconvex maximisation: multi-start SLSQP plus seeded starts. A sup that is missed makes the gate too
# OPTIMISTIC, so the seeds include every binding law known from rounds 1-12 and lane 1.
import sys, json, time
import numpy as np
from scipy.optimize import minimize
from scipy.special import xlogy

L = np.log(2.0)
def h(x):
    x = np.clip(x, 0.0, 1.0)
    return -xlogy(x, x) - xlogy(1 - x, 1 - x)
# r: root of -(1-r)log(1-r) - r log(1+r) = 0 in (0,1)
lo, hi = 0.5, 0.9
for _ in range(200):
    mid = (lo + hi) / 2
    val = -(1 - mid) * np.log1p(-mid) - mid * np.log1p(mid)
    lo, hi = (mid, hi) if val > 0 else (lo, mid)
r = (lo + hi) / 2
m = 0.5 + (1 - r) * np.log1p(r) / (2 * L); q = 2 * m - 1; cstar = 1 - m
tau = (0.5 - r * r) / (r * r * (1 - r) ** 2)
hr = h(r); p = np.log((1 - r) / r); c = r - hr / p
gam = hr / r; alpha = 1 - gam / p; t0 = r - m; bm = -hr * gam / (p * alpha); D = -hr / (p * alpha) - t0
es = gam * m; a_s = m / r; astar = m * h(r * r) / (r * hr); d0 = 2 * m - astar + 2 * q * L * a_s / (p * c)

def I(x, y): return x * y * (1 + tau * (1 - x) * (1 - y))
def phi(x): return np.sqrt(np.maximum(h(I(x, x)), 0.0) / L)
def F(x, y): return h(x * y) + q * L * phi(x) * phi(y) - y * h(x) - x * h(y)
def R(x): return F(r, x) / (q * L)
def vfun(x): return h(x) - hr * phi(x)
def bfun(x): return x - h(x) / p - c * phi(x)

class Cfg:
    def __init__(self, eps=2e-4, z=77/200, kap=0.39, C=13/25, a0=3/8, mu=None, Tmode='Phi'):
        self.eps, self.z, self.kap, self.C, self.a0, self.Tmode = eps, z, kap, C, a0, Tmode
        beta_ = 1 - D / kap; b0 = a0 ** beta_
        self.W = (kap - D) / (1 - b0); self.eta = (kap - D) / D * (b0 - a0) / (1 - b0)
        self.A0 = self.W / 2
        self.mu = mu if mu is not None else (self.eta * t0 + m) / (1.232 - astar)
        self.Bp = self.eta * (t0 + D) + self.W * (a0 - 0.5) - self.mu * d0
        self.gain = self.eta * t0 - self.mu * (1 - astar)

def unpack(v, k, Zc, Zd):
    x = v[:k]
    M = v[k:k + Zc * k].reshape(Zc, k)
    P = v[k + Zc * k:k + Zc * k + Zd * k].reshape(Zd, k) if Zd else None
    return x, M, P

def law_stats(x, M, P, cfg):
    w = M.sum(axis=0)                       # atom weights
    om = M.sum(axis=1)                      # latent weights
    s = float(w @ x); hx = h(x); e = float(w @ hx); t = float(w @ phi(x))
    X, Y = np.meshgrid(x, x, indexing='ij')
    U = float(w @ h(X * Y) @ w)
    KI = h(I(X, Y))
    V = 0.0
    for zz in range(M.shape[0]):
        if om[zz] > 1e-15:
            V += float(M[zz] @ KI @ M[zz]) / om[zz]
    g = U + q * V - 2 * s * e
    beta = g + 2 * (s - m + cfg.eps) * e
    xdev = e - es
    sig = float(w @ ((x - m) * hx)) + D * e - bm * m
    Rb = float(w @ R(x))
    Phi = 0.0
    if P is not None:
        for a in range(P.shape[0]):
            pa = P[a].sum()
            if pa > 1e-15:
                Phi += float(P[a] @ x) * float(P[a] @ hx) / pa
    return dict(s=s, e=e, t=t, U=U, V=V, g=g, beta=beta, x=xdev, sig=sig, Rb=Rb, Phi=Phi, w=w)

def lossG(st, cfg):
    x = st['x']
    lx = (cfg.A0 + cfg.Bp) * max(x, 0) + (cfg.A0 - cfg.Bp) * max(-x, 0)
    lb = max(0.0, -st['sig'] + cfg.C * st['Rb'])
    lU = astar * st['e'] + d0 * x - st['U']
    return lx + lb + cfg.mu * lU

def lossT(st, cfg):
    Phi = st['Phi'] if cfg.Tmode == 'Phi' else 0.0
    return m * st['e'] - Phi - cfg.mu * (st['U'] - st['e']) + cfg.gain * st['e']

def psi(theta, cfg, side, k=4, Zc=2, Zd=0, nstart=200, seed=0, seeds=()):
    rng = np.random.default_rng(seed)
    if side == 'T' and cfg.Tmode == 'Phi' and Zd == 0: Zd = k
    nv = k + Zc * k + Zd * k
    def obj(v):
        x, M, P = unpack(v, k, Zc, Zd)
        st = law_stats(x, M, P, cfg)
        if st['e'] < 1e-9: return 0.0
        l = lossG(st, cfg) if side == 'G' else lossT(st, cfg)
        return -(l - theta * st['beta']) / st['e']
    cons = [{'type': 'eq', 'fun': lambda v: v[k:k + Zc * k].sum() - 1.0},
            {'type': 'ineq', 'fun': lambda v: law_stats(*unpack(v, k, Zc, Zd), cfg)['s'] - (m - cfg.eps)}]
    if side == 'G':
        cons.append({'type': 'ineq', 'fun': lambda v: law_stats(*unpack(v, k, Zc, Zd), cfg)['e'] - cfg.z})
    else:
        cons.append({'type': 'ineq', 'fun': lambda v: cfg.z - law_stats(*unpack(v, k, Zc, Zd), cfg)['e']})
        cons.append({'type': 'ineq', 'fun': lambda v: law_stats(*unpack(v, k, Zc, Zd), cfg)['e'] - 1e-4})
    if Zd:
        def dec_eq(v):
            x, M, P = unpack(v, k, Zc, Zd)
            return P.sum(axis=0) - M.sum(axis=0)
        cons.append({'type': 'eq', 'fun': dec_eq})
    bounds = [(0.0, 1.0)] * nv
    best = (-np.inf, None)
    starts = []
    for sd in seeds:
        starts.append(sd)
    for _ in range(nstart):
        x0 = rng.random(k)
        if rng.random() < 0.5: x0[0] = 0.0
        if rng.random() < 0.3: x0[-1] = 1.0
        M0 = rng.random((Zc, k)); M0 /= M0.sum()
        v0 = np.concatenate([x0, M0.ravel()])
        if Zd:
            wv = M0.sum(axis=0); Pf = rng.random((Zd, k)); Pf = Pf / Pf.sum(axis=0) * wv
            v0 = np.concatenate([v0, Pf.ravel()])
        starts.append(v0)
    for v0 in starts:
        v0 = np.asarray(v0, float)
        if len(v0) != nv: continue
        try:
            res = minimize(obj, v0, method='SLSQP', bounds=bounds, constraints=cons,
                           options=dict(maxiter=400, ftol=1e-12))
        except Exception:
            continue
        v = np.clip(res.x, 0, 1)
        ok = all((c_['fun'](v) >= -1e-7) if c_['type'] == 'ineq' else np.all(np.abs(c_['fun'](v)) < 1e-7) for c_ in cons)
        if not ok: continue
        val = -obj(v)
        if val > best[0]: best = (val, v)
    return best

def seed_vec(xs, ws, k, Zc, Zd, cfg, side):
    xs = list(xs) + [0.5] * (k - len(xs)); ws = list(ws) + [0.0] * (k - len(ws))
    M = np.zeros((Zc, k)); M[0] = ws
    v = np.concatenate([np.array(xs, float), M.ravel()])
    if side == 'T' and cfg.Tmode == 'Phi':
        Zd = Zd or k
        P = np.zeros((Zd, k)); P[0] = ws
        v = np.concatenate([v, P.ravel()])
    return v

def describe(v, k, Zc, Zd, cfg):
    x, M, P = unpack(np.clip(v, 0, 1), k, Zc, Zd)
    st = law_stats(x, M, P, cfg)
    return dict(atoms=[round(float(a), 5) for a in x], weights=[round(float(a), 5) for a in st['w']],
                s=round(st['s'], 6), e=round(st['e'], 6), g=float(st['g']), beta=float(st['beta']), x=round(st['x'], 6))

if __name__ == '__main__':
    cfg = Cfg(eps=float(sys.argv[1]) if len(sys.argv) > 1 else 2e-4)
    print(f"r={r:.15f} m={m:.15f} c*={cstar:.15f} tau={tau:.6f} D={D:.6f} t0={t0:.6f} e*={es:.6f} a*={astar:.6f} d0={d0:.6f}")
    print(f"eta={cfg.eta:.6f} W={cfg.W:.6f} A0={cfg.A0:.6f} B'={cfg.Bp:.6f} mu={cfg.mu:.6f} gain={cfg.gain:.6f}")
