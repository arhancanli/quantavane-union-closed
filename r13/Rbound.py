# FLOAT: smallest c1 with R(x) <= c1 v(x)^2 + c2 b(x) on [0,1], for several c2 (grid incl. log-refined edges and the contact r).
import numpy as np
from gateA import R, vfun, bfun, r
xs = np.unique(np.clip(np.concatenate([np.logspace(-14, -1, 400), np.linspace(0.1, 0.9, 40001), 1 - np.logspace(-14, -1, 400),
                                       r + np.linspace(-0.01, 0.01, 4001)]), 0, 1))
Rv, v2, bv = R(xs), vfun(xs)**2, bfun(xs)
for c2 in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
    rem = Rv - c2 * bv
    need = np.where(v2 > 1e-300, rem / np.maximum(v2, 1e-300), np.where(rem > 1e-15, np.inf, 0.0))
    i = np.argmax(need)
    print(f"c2={c2:.2f}: min c1 = {need[i]:.4f} (binding x = {xs[i]:.6g}); R_max={Rv.max():.4f}")
