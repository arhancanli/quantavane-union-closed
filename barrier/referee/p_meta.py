# Referee check of the logical part of Prop P against Theorem JB (SELF_CONTAINED.md Sec. 6), and of an input
# that the product measure violates besides the V-budget. mpmath 40 digits (FLOAT-level, sanity only).
import mpmath as mp
mp.mp.dps = 40
L = mp.log(2)
h = lambda x: -x*mp.log(x) - (1-x)*mp.log(1-x)
r = mp.findroot(lambda r: -(1-r)*mp.log(1-r) - r*mp.log(1+r), 0.69)
m = mp.mpf(1)/2 + (1-r)*mp.log(1+r)/(2*L); q = 2*m-1
tau = (mp.mpf(1)/2 - r**2)/(r**2*(1-r)**2)
I = lambda x, y: x*y*(1+tau*(1-x)*(1-y))
phi = lambda x: mp.sqrt(h(I(x, x))/L)
Fxy = lambda x, y: h(x*y) + q*L*phi(x)*phi(y) - y*h(x) - x*h(y)
R = lambda x: Fxy(r, x)/(q*L)
hr = h(r); p = mp.log((1-r)/r); c = r - hr/p
gamma = hr/r; alpha = 1 - gamma/p; t0 = r - m
b_m = -hr*gamma/(p*alpha); D = -hr/(p*alpha) - t0
estar = gamma*m; astar = m*h(r**2)/(r*hr); a_s = m/r
d0 = 2*m - astar + 2*q*L*a_s/(p*c)
kappa = mp.mpf(39)/100; Cc = mp.mpf(13)/25
print("D", mp.nstr(D, 8), "t0", mp.nstr(t0, 8), "e*", mp.nstr(estar, 8), "a*", mp.nstr(astar, 8), "d0", mp.nstr(d0, 8))

def consts(a0, mu):
    be = 1 - D/kappa; b0 = a0**be
    W = (kappa - D)/(1 - b0); eta = ((kappa - D)/D)*(b0 - a0)/(1 - b0)
    A0 = W/2; Bp = eta*(t0 + D) + W*(a0 - mp.mpf(1)/2) - mu*d0; G0 = eta*t0 - mu*(1 - astar)
    return A0, Bp, G0, eta, W

A0, Bp, G0, eta, W = consts(mp.mpf(19)/100, mp.mpf(14)/5)
print("JB params a0=19/100 mu=14/5: eta", mp.nstr(eta, 12), "W", mp.nstr(W, 12), "Gamma0", mp.nstr(G0, 12), "B'", mp.nstr(Bp, 12))

def lG_point(s, a0, mu):
    A0, Bp, G0, eta, W = consts(a0, mu)
    e = h(s); x = e - estar
    sig = (s - m)*e + D*e - b_m*m
    lg = (A0+Bp)*max(x, 0) + (A0-Bp)*max(-x, 0) + max(Cc*R(s) - sig, 0) + mu*(astar*e + d0*x - h(s*s))
    return lg, G0*e

# (1) the chain of JB applied to the product measure must give l_G(delta_s) >= Gamma0 h(s) for EVERY (a0, mu)
worst = None
for eps in (mp.mpf(29)/100000, mp.mpf(1205)/10**6, mp.mpf('0.01')):
    s = m - eps
    for a0 in [mp.mpf(k)/100 for k in range(1, 100)]:
        for mu in [mp.mpf(k)/10 for k in range(0, 101)]:
            lg, rhs = lG_point(s, a0, mu)
            d = lg - rhs
            if worst is None or d < worst[0]:
                worst = (d, eps, a0, mu)
print("min over grid of l_G(delta_s) - Gamma0 h(s) =", mp.nstr(worst[0], 8), "at eps,a0,mu =", [mp.nstr(v, 5) for v in worst[1:]])

# (2) direct JB contradiction at eps: (P_G) at delta_s + (M) need l_G - Gamma0 e < theta*C(s)
for eps in (mp.mpf(29)/100000, mp.mpf('0.0012'), mp.mpf(1205)/10**6):
    s = m - eps
    g = h(s*s) + q*h(I(s, s)) - 2*s*h(s)
    Cs = g - 2*eps*h(s)
    lg, rhs = lG_point(s, mp.mpf(19)/100, mp.mpf(14)/5)
    print(f"eps={mp.nstr(eps,6)}: C(delta_s)={mp.nstr(Cs,6)}  g={mp.nstr(g,8)}  l_G-Gamma0*e (JB params)={mp.nstr(lg-rhs,8)}")

# (3) inputs that hold for EVERY uniform measure on a union-closed family but are not in the list (i)-(iv):
#     Gibbs-gap identity H - sum U_i = sum_i E d(Q_i(a)Q_i(b) || Q_i(a u b))   (X0 = 0)
d = lambda a, b: a*mp.log(a/b) + (1-a)*mp.log((1-a)/(1-b))
for eps in (mp.mpf(0), mp.mpf(1205)/10**6):
    s = m - eps
    print(f"eps={mp.nstr(eps,6)} per coordinate: h(s)-h(s^2) = {mp.nstr(h(s)-h(s*s),6)};  d(s^2||s) = {mp.nstr(d(s*s, s),6)};"
          f"  2 eps h(s) = {mp.nstr(2*eps*h(s),6)};  X0/n = {mp.nstr(h(s) - h(s*s) - d(s*s, s),6)}")
    # cross-entropy -E log mu(A u B) vs H(mu) per coordinate
    ce = -(s*s*mp.log(s) + (1-s*s)*mp.log(1-s))
    print(f"   cross-entropy -E log mu(AuB)/n = {mp.nstr(ce,8)} vs H/n = {mp.nstr(h(s),8)} (uniform measures: equal)")
# (4) the older point-mass estimate in NOTES (b): g(delta_m)/(2 h(m))
gm = h(m*m) + q*h(I(m, m)) - 2*m*h(m)
print("g(delta_m) =", mp.nstr(gm, 8), " g/(2h(m)) =", mp.nstr(gm/(2*h(m)), 8), " C at eps=1.2e-3 is listed above")
