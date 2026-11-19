# Tier 1 lane, union-closed (mpmath, 40 digits; FLOAT-level, not an interval certificate). 19 Nov 2026, 22:4x +04.
# Product-measure barrier: mu = Bern(s)^n on absence vectors, s = m - eps. Every coordinate law is the point mass at s in every
# order and every fibre (fibres of products are products). Check, as functions of eps:
#   combined budget (SU (4.1)):  g(s) + 2 (s - m) h(s) <= 0,  g(s) = h(s^2) + q h(I(s,s)) - 2 s h(s)
#   U-budget:                     h(s^2) <= h(s)
#   fibre inequality (4.2):       (s - m) h(s) <= 0
#   V-budget (NOT used by SU/JB): h(I(s,s)) <= h(s)
# eps_P := root of g(m - eps) = 2 eps h(m - eps).
from mpmath import mp, mpf, log, sqrt, findroot
mp.dps = 40
def h(x): return -x*log(x) - (1-x)*log(1-x)
r = findroot(lambda x: -(1-x)*log(1-x) - x*log(1+x), mpf('0.69'))
L = log(2); m = mpf(1)/2 + (1-r)*log(1+r)/(2*L); q = 2*m - 1; tau = (mpf(1)/2 - r*r)/(r*r*(1-r)**2)
I = lambda x, y: x*y*(1 + tau*(1-x)*(1-y))
g = lambda s: h(s*s) + q*h(I(s, s)) - 2*s*h(s)
epsP = findroot(lambda e: g(m - e) - 2*e*h(m - e), mpf('0.00117'))
print('r =', r); print('m =', m, ' c* =', 1 - m)
print('g(delta_m) =', g(m), ' g/(2h) at m =', g(m)/(2*h(m)))
print('eps_P =', epsP, ' c* + eps_P =', 1 - m + epsP)
s = m - epsP
print('at eps_P: U-budget h(s^2) - h(s) =', h(s*s) - h(s), '; V-budget h(I(s,s)) - h(s) =', h(I(s, s)) - h(s),
      '; fibre term (s-m)h(s) =', (s - m)*h(s))
for e in [mpf('1e-3'), mpf('1.2e-3'), mpf('1.5e-3'), mpf('3e-3')]:
    s = m - e
    print(f"eps={float(e):.1e}: combined {float(g(s) + 2*(s-m)*h(s)):+.3e}  U {float(h(s*s)-h(s)):+.3e}  V {float(h(I(s,s))-h(s)):+.3e}")
