# Point red controls for bopB.py: with B's own enclosures, an upper bound < 0 of G at one point t means bopB cannot PASS.
import sys, math, importlib.util
from fractions import Fraction
for cq in ['1618035/1000000', '7/4']:
    sys.argv = ['bopB.py', cq]
    src = open('bopB.py').read().split("cl, ch = dn")[0]      # constants, hval, henc, edge constants only (no bisection)
    ns = {}; exec(src, ns)
    t = 0.6180339887498949  # ~ 1/phi
    a2 = ns['dn'](t * t); b2 = ns['up'](t * t)
    h2hi = ns['henc'](a2, b2)[1]; hlo = ns['hval'](t)[0]
    civ = ns['civ']; clo = ns['dn'](float(civ.a))
    Gup = ns['up'](h2hi - ns['dn'](ns['dn'](clo * t) * hlo))
    print(f"c = {cq}: edges ok = {ns['ok_edges']}; G(1/phi) <= {Gup:.3e} -> {'cannot PASS (red OK)' if (Gup < 0 or not ns['ok_edges']) else 'no conclusion'}")
