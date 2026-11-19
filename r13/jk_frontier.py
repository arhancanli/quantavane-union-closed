# Step 3 FLOAT: valid joint kernels F(x,y) >= a[v^2(x)phi(y)+v^2(y)phi(x)] + aw[w(x)phi(y)+w(y)phi(x)] + k[phi(x)R(y)+phi(y)R(x)],
# w = v^2/phi (0 where phi = 0). Grid check on [0,1]^2 (log-refined at 0, 1 and near r). For each (a, aw) report the max k.
import numpy as np
from gateA import h, phi, F, R, vfun, r
g1 = np.unique(np.clip(np.concatenate([np.logspace(-12, -1, 90), np.linspace(0.1, 0.9, 1201), 1 - np.logspace(-12, -1, 90)[::-1],
                                       r + np.linspace(-0.03, 0.03, 301)]), 0, 1))
X, Y = np.meshgrid(g1, g1, indexing='ij')
Fv = F(X, Y); px = phi(X); py = phi(Y); vx = vfun(X); vy = vfun(Y); Rx = R(X); Ry = R(Y)
wx = np.where(px > 0, vx**2 / np.maximum(px, 1e-300), 0.0); wy = np.where(py > 0, vy**2 / np.maximum(py, 1e-300), 0.0)
Av = vx**2 * py + vy**2 * px; Aw = wx * py + wy * px; AR = px * Ry + py * Rx
rows = []
for a in [0.0, 0.05, 0.10, 0.15, 0.18, 0.20, 0.22]:
    for aw in [0.0, 0.05, 0.10, 0.15, 0.18, 0.20, 0.22]:
        slack = Fv - a * Av - aw * Aw
        bad = slack[AR <= 1e-14]
        if bad.size and bad.min() < -1e-12: continue
        sel = AR > 1e-14
        k = float(np.min(slack[sel] / AR[sel]))
        if k < 0: continue
        rows.append((a, aw, round(k, 5)))
        print(f"a={a:.2f} aw={aw:.2f}: max kappa' = {k:.5f}", flush=True)
np.save('jk_frontier.npy', np.array(rows))
