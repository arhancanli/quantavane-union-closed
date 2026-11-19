# Round 15 prep (FLOAT): best bilinear moment kernel F(x,y) >= f(x)^T M f(y), f = (1, x, h, phi, R), at the relaxation maximiser P
# (eps 2.9e-4), on successively finer grids (the coarse-grid LP of r15_slack.py over-fits between grid points). Each M is then
# checked on a much finer test grid (min of F - f^T M f) and against the realising law of r15_slack.py.
import sys; sys.argv = sys.argv[:1]
import numpy as np
from scipy.optimize import linprog
from gateA import *
eps = 2.9e-4; s0, e0, bb0, R0 = m - eps, 0.61414, 8.597e-4, 0.002428
t0_ = (s0 - e0 / p - bb0) / c; P = np.array([1.0, s0, e0, t0_, R0])
idx = [(i, j) for i in range(5) for j in range(i, 5)]
def fb(xs): return np.stack([np.ones_like(xs), xs, h(xs), phi(xs), R(xs)], axis=1)
def lp(xs):
    f = fb(xs); X, Y = np.triu_indices(len(xs))
    A = np.stack([f[X, i] * f[Y, j] + (f[X, j] * f[Y, i] if i != j else 0 * f[X, i]) for i, j in idx], axis=1)
    obj = -np.array([P[i] * P[j] * (2 if i != j else 1) for i, j in idx])
    res = linprog(obj, A_ub=A, b_ub=F(xs[X], xs[Y]), bounds=[(-50, 50)] * len(idx), method='highs')
    M = np.zeros((5, 5))
    for (i, j), v in zip(idx, res.x): M[i, j] = M[j, i] = v
    return -res.fun, M
def worst(M, xs):
    f = fb(xs); K = f @ M @ f.T; D = F(xs[:, None], xs[None, :]) - K; i = np.unravel_index(np.argmin(D), D.shape)
    return D[i], (xs[i[0]], xs[i[1]])
test = np.unique(np.concatenate([np.linspace(0, 1, 3001), np.logspace(-12, -1, 120), 1 - np.logspace(-12, -1, 120)]))
lawx = np.array([0.0, 0.306, 0.65, 0.6505]); laww = np.array([5.1170e-02, 5.2000e-04, 4.2599e-01, 5.2232e-01])
for n, dens in [(201, 101), (401, 201), (801, 301)]:
    xs = np.unique(np.concatenate([np.linspace(0, 1, n), np.linspace(0.55, 0.75, dens), np.logspace(-12, -1, 60), 1 - np.logspace(-12, -1, 60), [r]]))
    val, M = lp(xs); wv, wp = worst(M, test)
    fl = fb(lawx); EFl = float(laww @ F(lawx[:, None], lawx[None, :]) @ laww); Kl = float(laww @ (fl @ M @ fl.T) @ laww)
    print(f"grid {len(xs)}: LP bound at P = {val:.4e}; min(F - K) on test grid = {wv:+.3e} at {np.round(wp, 5)}; law: EF {EFl:.4e} vs kernel {Kl:.4e}", flush=True)
