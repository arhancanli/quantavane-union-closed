# Round 14 (FLOAT LP): best linear majorant R(x) <= a0 + a1 x + a2 h(x) + a3 phi(x) + a4 b(x) on [0,1] at a given moment point
# (1, s, e, t, bbar). Its value is the largest Rbar any law with those moments can have (LP duality, up to grid resolution).
import numpy as np
from scipy.optimize import linprog
from gateA import h, phi, R, bfun, m, p, c
xs = np.unique(np.concatenate([np.linspace(0, 1, 40001), np.logspace(-14, -2, 400), 1 - np.logspace(-14, -2, 400)]))
xs = xs[(xs >= 0) & (xs <= 1)]
Fm = np.stack([np.ones_like(xs), xs, h(xs), phi(xs), bfun(xs)], axis=1); Rv = R(xs)
def best(s, e, bb):
    t = (s - e / p - bb) / c; mom = np.array([1, s, e, t, bb])
    res = linprog(mom, A_ub=-Fm, b_ub=-(Rv + 1e-5), bounds=[(None, None)] * 5, method='highs')   # minimise <alpha, mom> s.t. Fm alpha >= R + 1e-5
    if res.status != 0: return float("nan"), [res.status, res.message[:60]], t
    return res.fun, res.x, t
for (s, e, bb, lab) in [(m - 2.8e-4, 0.385, 0.496, 'CORNER theta 30'), (m - 2.8e-4, 0.385, 0.40, 'corner b=0.40'), (m - 2.8e-4, 0.385, 0.30, 'corner b=0.30'),
                        (m - 2.8e-4, 0.606, 0.0, 'contact region'), (m - 2.8e-4, 0.45, 0.2, 'middle')]:
    v, a, t = best(s, e, bb)
    print(f"{lab:16s} (s-m={s-m:+.1e}, e={e}, t={t:.4f}, bbar={bb}): max Rbar <= {v:.5f}  alpha = {a if isinstance(a, list) else np.round(a, 5)}")
