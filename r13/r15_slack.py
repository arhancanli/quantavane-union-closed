# Round 15 prep (FLOAT): where is the slack at the certified relaxation's maximiser (eps 2.9e-4, a0 0.19, theta 31)?
# P = (s, e, t, Rbar) = (m - eps, 0.61414, t, 0.002428), bbar = 8.597e-4.
# (a) EF_lo there (SU's kernels + K2) and with every grid-valid joint kernel of the frontier (kappa_K > kappa4 or not);
# (b) min of EF = E_{mu x mu} F over laws with exactly these four moments (SLSQP, k atoms, many starts) and that law's Var v(X);
# (c) best BILINEAR moment kernel: max sum M_ij P_i P_j s.t. f(x)^T M f(y) <= F(x,y) on a grid, f = (1, x, h, phi, R) (LP).
import sys; sys.argv = sys.argv[:1]
import numpy as np
from scipy.optimize import minimize, linprog
from gateA import *
eps = 2.9e-4; s0, e0, bb0, R0 = m - eps, 0.61414, 8.597e-4, 0.002428
t0_ = (s0 - e0 / p - bb0) / c; vb0 = e0 - hr * t0_
P = np.array([1.0, s0, e0, t0_, R0])
K2 = [(0.22, 0.0), (0.0, 0.05), (0.15, 0.0255)]
FR = [(0.22, 0.0005), (0.21, 0.0040), (0.20, 0.0078), (0.18, 0.0150), (0.15, 0.0255), (0.10, 0.0418), (0.05, 0.0535), (0.0, 0.0583)]
EF = lambda pairs, w2: max(2 * t0_ * (a * w2 + k * R0) for a, k in pairs)
print(f"P: s={s0:.6f} e={e0} t={t0_:.6f} R={R0} vbar={vb0:+.3e}")
print(f"(a) EF_lo certified (K2 set) = {EF(K2, vb0**2):.4e}; with whole frontier = {EF(K2 + FR, vb0**2):.4e}")
# (b) laws realising P (up to the 2e-6 slack of the rational cuts): LP on a grid, minimise the L1 moment residual, then report the law
xg = np.unique(np.concatenate([np.linspace(0, 1, 4001), np.logspace(-10, -2, 200), 1 - np.logspace(-10, -2, 200), [r]]))
G = np.stack([np.ones_like(xg), xg, h(xg), phi(xg), R(xg)])
n = len(xg)
# variables: w (n) >= 0, slack+ (5), slack- (5); G w + sp - sm = P; minimise sum(sp + sm)
Aeq = np.hstack([G, np.eye(5), -np.eye(5)])
res = linprog(np.concatenate([np.zeros(n), np.ones(10)]), A_eq=Aeq, b_eq=P, bounds=[(0, None)] * (n + 10), method='highs')
w = res.x[:n]; sel = w > 1e-10; x, w = xg[sel], w[sel]
vx = vfun(x); EFlaw = float(w @ F(x[:, None], x[None, :]) @ w); varv = float(w @ vx**2 - (w @ vx) ** 2)
print(f"(b) a law with these moments (L1 residual {res.fun:.2e}): atoms {np.round(x, 5)} weights {np.round(w, 5)}")
print(f"    its EF = {EFlaw:.4e}; Var v = {varv:.3e}; EF_lo with E v^2 in place of vbar^2 = {EF(K2, float(w @ vx**2)):.4e}")
xs = np.unique(np.concatenate([np.linspace(0, 1, 121), np.logspace(-8, -1.3, 25), 1 - np.logspace(-8, -1.3, 25), [r]]))
fb = np.stack([np.ones_like(xs), xs, h(xs), phi(xs), R(xs)], axis=1)
X, Y = np.meshgrid(np.arange(len(xs)), np.arange(len(xs)), indexing='ij'); X, Y = X.ravel(), Y.ravel(); keep = X <= Y; X, Y = X[keep], Y[keep]
idx = [(i, j) for i in range(5) for j in range(i, 5)]
A = np.stack([(fb[X, i] * fb[Y, j] + (fb[X, j] * fb[Y, i] if i != j else 0) * (1 if i != j else 0)) if i != j else fb[X, i] * fb[Y, i] for i, j in idx], axis=1)
b_ = F(xs[X], xs[Y])
obj = -np.array([P[i] * P[j] * (2 if i != j else 1) for i, j in idx])
res = linprog(obj, A_ub=A, b_ub=b_, bounds=[(-50, 50)] * len(idx), method='highs')
print(f"(c) best bilinear moment bound on EF at P (grid LP, |M_ij| <= 50): {(-res.fun if res.status == 0 else float('nan')):.4e} (status {res.status})")
