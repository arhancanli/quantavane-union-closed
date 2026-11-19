# Round 14 (FLOAT): is the binding point of the N3 relaxation realisable? Moments (1, s, e, t, Rbar) of a law X on [0,1] lie in the
# convex hull of g(x) = (1, x, h(x), phi(x), R(x)). b-bar is not a free coordinate: b-bar = s - e/p - c t.
# Separation LP: min alpha . P  s.t. alpha . g(x) >= 0 on a fine grid, |alpha_j| <= 1. A negative value = P outside the hull (a cut).
import sys; sys.argv = sys.argv[:1]
import numpy as np
from scipy.optimize import linprog
from gateA import h, phi, R, bfun, m, p, c, L
xs = np.unique(np.concatenate([np.linspace(0, 1, 20001), np.logspace(-15, -2, 300), 1 - np.logspace(-15, -2, 300)]))
xs = xs[(xs >= 0) & (xs <= 1)]
G = np.stack([np.ones_like(xs), xs, h(xs), phi(xs), R(xs)], axis=1)
def sep(P):
    res = linprog(P, A_ub=-G, b_ub=np.zeros(len(xs)), bounds=[(-1, 1)] * 5, method='highs')
    return res.fun, res.x
bz = bfun(xs); print("b(x) min", bz.min(), "zeros of b near:", xs[bz < 1e-6][:3], xs[bz < 1e-6][-3:])
print("R max", R(xs).max(), "at x =", xs[np.argmax(R(xs))])
for (s, e, bb, Rb, lab) in [(m - 2.8e-4, 0.385, 0.4962, 0.2150, 'theta30 corner'), (m - 2.8e-4, 0.6063, 0.0, 1.81e-3, 'theta35 contact'),
                            (m - 2.8e-4, 0.6063, 0.01, 1.81e-3, 'contact b=0.01'), (m - 2.8e-4, 0.6063, 0.03, 0.0, 'contact b=0.03')]:
    t = (s - e / p - bb) / c; P = np.array([1, s, e, t, Rb]); v, a = sep(P)
    print(f"{lab:18s} t={t:.4f}: separation value {v:+.4e}  ({'OUTSIDE hull' if v < -1e-9 else 'inside'})  alpha={np.round(a, 5)}")
