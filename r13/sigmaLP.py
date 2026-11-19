# Best linear minorant sqrt(x h(x)) >= a phi(x) + b h(x) + c x + d on [0,1] for the moment point (t0, e0, s0) (LP, FLOAT).
import numpy as np
from scipy.optimize import linprog
from gateA import *
xs = np.unique(np.concatenate([np.linspace(0, 1, 20001), np.logspace(-14, -2, 400), 1 - np.logspace(-14, -2, 400)]))
xs = xs[(xs >= 0) & (xs <= 1)]
f = np.sqrt(np.maximum(xs * h(xs), 0))
A = np.stack([phi(xs), h(xs), xs, np.ones_like(xs)], axis=1)
for (t0, e0, s0) in [(0.7296, 0.385, m - 1.5e-4), (0.744, 0.385, 0.618), (0.60, 0.30, 0.62), (0.40, 0.20, 0.62)]:
    res = linprog(-np.array([t0, e0, s0, 1.0]), A_ub=A, b_ub=f, bounds=[(None, None)] * 3 + [(None, 0)], method='highs')
    a, b, c_, d = res.x
    print(f"point (t,e,s)=({t0},{e0},{s0:.4f}): sigma >= {-res.fun:.5f}  with a={a:.5f} b={b:.5f} c={c_:.5f} d={d:.5f}; min slack on grid {np.min(f - A@res.x):.2e}")
