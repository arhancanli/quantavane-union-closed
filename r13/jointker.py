# FLOAT scan: for which (kE', kap') does F(x,y) >= kE'[v(x)^2 phi(y) + v(y)^2 phi(x)] + kap'[phi(x)R(y) + phi(y)R(x)]
# hold on a grid of [0,1]^2 (log-refined near 0, 1 and r)? Reports the max kap' for each kE'.
import numpy as np
from gateA import h, phi, F, R, vfun, r
g1 = np.concatenate([np.logspace(-12, -1, 120), np.linspace(0.1, 0.9, 1601), 1 - np.logspace(-12, -1, 120)[::-1], r + np.linspace(-0.02, 0.02, 401)])
g1 = np.unique(np.clip(g1, 0, 1))
X, Y = np.meshgrid(g1, g1, indexing='ij')
Fv = F(X, Y); px = phi(X); py = phi(Y); vx = vfun(X); vy = vfun(Y); Rx = R(X); Ry = R(Y)
Av = vx**2 * py + vy**2 * px; AR = px * Ry + py * Rx
mask = (Av + AR) > 1e-14
print("grid", len(g1), "min F", Fv.min(), "min R", min(Rx.min(), 0))
for kE in [0.0, 0.10, 0.15, 0.18, 0.20, 0.21, 0.22, 0.23, 0.24]:
    # need F - kE Av - k AR >= 0 everywhere; for AR > 0: k <= (F - kE Av)/AR
    slack = Fv - kE * Av
    if np.any(slack[(AR <= 1e-14) & mask] < -1e-12):
        print(f"kE'={kE}: FAILS even with kap'=0 (min slack {slack[(AR<=1e-14)&mask].min():.3e})"); continue
    sel = AR > 1e-14
    ratio = slack[sel] / AR[sel]
    i = np.argmin(ratio)
    print(f"kE'={kE:.2f}: max kap' = {ratio[i]:.5f} at (x,y)=({X[sel][i]:.5f},{Y[sel][i]:.5f}); min slack at kap'=0: {slack.min():.3e}")
