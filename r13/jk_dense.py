# Denser validation of the three key joint kernels: max kappa' for a = 0.20, 0.15, 0.05 on a grid refined near 0, 1, r and (1,1).
import numpy as np
from gateA import h, phi, F, R, vfun, r
edge = np.logspace(-15, -1, 300)
g1 = np.unique(np.clip(np.concatenate([edge, np.linspace(0.1, 0.9, 2401), 1 - edge[::-1], r + np.linspace(-0.02, 0.02, 801),
                                       1 - np.logspace(-6, -2, 400)]), 0, 1))
X, Y = np.meshgrid(g1, g1, indexing='ij')
Fv = F(X, Y); px = phi(X); py = phi(Y); vx = vfun(X); vy = vfun(Y); Rx = R(X); Ry = R(Y)
Av = vx**2 * py + vy**2 * px; AR = px * Ry + py * Rx
print("grid", len(g1))
for a in [0.20, 0.15, 0.05, 0.0]:
    slack = Fv - a * Av; sel = AR > 1e-14
    ratio = slack[sel] / AR[sel]; i = np.argmin(ratio)
    print(f"a={a:.2f}: max kappa' = {ratio[i]:.6f} at (x,y)=({X[sel][i]:.8f},{Y[sel][i]:.8f})", flush=True)
