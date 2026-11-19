# Quantavane Program A, core: Arb ball arithmetic (python-flint 0.7.1) with forward-mode jets
# (value, gradient, Hessian) and monotone endpoint enclosures for h, h', h'', sqrt, log.
# Written 16 Nov 2026 for the replay of Sadhu-Ursin (Zenodo 23322347). Independent of the authors' code
# and of Program B (C, binary64 intervals + MPFR).
from flint import arb, fmpq, ctx
ctx.prec = 128

def Q(n, d=1):
    return arb(fmpq(n, d))

ZERO = arb(0); ONE = arb(1); HALF = Q(1, 2)
L = arb.const_log2()
INF = arb(0, float('inf')) if False else None

def NAN():
    return arb.nan()

def hull(*xs):
    out = xs[0]
    for x in xs[1:]:
        out = arb.union(out, x)
    return out

def lo(x): return x.lower()
def hi(x): return x.upper()

def h_point(u):
    # u: arb (thin or not), assumed inside [0,1]; exact 0 and 1 handled
    if u.is_zero():
        return ZERO
    if (u - 1).is_zero():
        return ZERO
    return -u * u.log() - (1 - u) * (1 - u).log()

def clamp01(u):
    a = lo(u); b = hi(u)
    if a < 0: a = ZERO
    if b > 1: b = ONE
    if a > b:
        raise ValueError("empty after clamp")
    return arb(a), arb(b)

def h_range(u):
    a, b = clamp01(u)
    ha, hb = h_point(a), h_point(b)
    if b <= HALF:
        return hull(ha, hb)
    if a >= HALF:
        return hull(ha, hb)
    return hull(ha, hb, L)

def hp_range(u):
    # h'(u) = log((1-u)/u), strictly decreasing on (0,1)
    a, b = clamp01(u)
    if a <= 0 or b >= 1:
        return NAN()
    return hull(((1 - b) / b).log(), ((1 - a) / a).log())

def hpp_range(u):
    # h''(u) = -1/(u(1-u))
    a, b = clamp01(u)
    if a <= 0 or b >= 1:
        return NAN()
    ga, gb = a * (1 - a), b * (1 - b)
    g = hull(ga, gb)
    if a <= HALF and b >= HALF:
        g = hull(g, Q(1, 4))
    return -1 / g

def sqrt_range(w):
    a = lo(w); b = hi(w)
    if a < 0: a = ZERO
    if b < 0:
        raise ValueError("sqrt of negative")
    return hull(arb(a).sqrt(), arb(b).sqrt())

def log_range(w):
    a = lo(w); b = hi(w)
    if a <= 0:
        return NAN()
    return hull(arb(a).log(), arb(b).log())

class Jet:
    """Forward-mode jet: value v, gradient d (list, None = structural zero), Hessian H (dict (i,j)->arb, i<=j)."""
    __slots__ = ('v', 'd', 'H', 'n', 'order')
    def __init__(self, v, d=None, H=None, n=0, order=0):
        self.v = v; self.n = n; self.order = order
        self.d = d if d is not None else [None] * n
        self.H = H if H is not None else {}
    @staticmethod
    def var(val, i, n, order):
        d = [None] * n
        if order >= 1:
            d[i] = ONE
        return Jet(val, d, {}, n, order)
    def const(self, c):
        return Jet(c if isinstance(c, arb) else arb(c), [None] * self.n, {}, self.n, self.order)
    def _lift(self, o):
        if isinstance(o, Jet): return o
        return self.const(o)
    def __add__(self, o):
        o = self._lift(o)
        d = [_add(a, b) for a, b in zip(self.d, o.d)]
        H = dict(self.H)
        for k, val in o.H.items():
            H[k] = H[k] + val if k in H else val
        return Jet(self.v + o.v, d, H, self.n, self.order)
    __radd__ = __add__
    def __neg__(self):
        return Jet(-self.v, [None if a is None else -a for a in self.d], {k: -v for k, v in self.H.items()}, self.n, self.order)
    def __sub__(self, o):
        return self + (-self._lift(o))
    def __rsub__(self, o):
        return self._lift(o) + (-self)
    def __mul__(self, o):
        if not isinstance(o, Jet):
            c = o if isinstance(o, arb) else arb(o)
            return Jet(self.v * c, [None if a is None else a * c for a in self.d], {k: v * c for k, v in self.H.items()}, self.n, self.order)
        a, b = self, o
        d = []
        for i in range(self.n):
            t = None
            if b.d[i] is not None: t = a.v * b.d[i]
            if a.d[i] is not None: t = _add(t, b.v * a.d[i])
            d.append(t)
        H = {}
        if self.order >= 2:
            for k, val in a.H.items(): H[k] = H.get(k, None); H[k] = _add(H[k], b.v * val)
            for k, val in b.H.items(): H[k] = _add(H.get(k, None), a.v * val)
            for i in range(self.n):
                for j in range(i, self.n):
                    t = None
                    if a.d[i] is not None and b.d[j] is not None: t = a.d[i] * b.d[j]
                    if a.d[j] is not None and b.d[i] is not None: t = _add(t, a.d[j] * b.d[i])
                    if t is not None: H[(i, j)] = _add(H.get((i, j), None), t)
            H = {k: v for k, v in H.items() if v is not None}
        return Jet(a.v * b.v, d, H, self.n, self.order)
    __rmul__ = __mul__
    def compose(self, f0, f1, f2):
        # f(u) with f0 = f(u.v), f1 = f'(u.v), f2 = f''(u.v) enclosures
        d = [None if a is None else f1 * a for a in self.d]
        H = {}
        if self.order >= 2:
            for k, val in self.H.items(): H[k] = f1 * val
            for i in range(self.n):
                for j in range(i, self.n):
                    if self.d[i] is not None and self.d[j] is not None:
                        H[(i, j)] = _add(H.get((i, j), None), f2 * self.d[i] * self.d[j])
        return Jet(f0, d, H, self.n, self.order)
    def has_deriv(self):
        return any(a is not None for a in self.d)
    def inv(self):
        f0 = 1 / self.v
        if self.order == 0 or not (self.has_deriv() or self.H):
            return Jet(f0, [None] * self.n, {}, self.n, self.order)
        return self.compose(f0, -f0 * f0, 2 * f0 * f0 * f0)
    def __truediv__(self, o):
        if not isinstance(o, Jet):
            c = o if isinstance(o, arb) else arb(o)
            return self * (1 / c)
        return self * o.inv()
    def __rtruediv__(self, o):
        return self._lift(o) * self.inv()

def _add(a, b):
    if a is None: return b
    if b is None: return a
    return a + b

def need(u):
    return u.order >= 1 and (u.has_deriv() or u.H)

def jh(u):
    f0 = h_range(u.v)
    if not need(u): return Jet(f0, [None] * u.n, {}, u.n, u.order)
    return u.compose(f0, hp_range(u.v), hpp_range(u.v) if u.order >= 2 else ZERO)

def jsqrt(u):
    f0 = sqrt_range(u.v)
    if not need(u): return Jet(f0, [None] * u.n, {}, u.n, u.order)
    a = lo(u.v)
    if a <= 0:
        return u.compose(f0, NAN(), NAN())
    s = f0
    return u.compose(f0, 1 / (2 * s), -1 / (4 * s * s * s) if u.order >= 2 else ZERO)

def jlog(u):
    f0 = log_range(u.v)
    if not need(u): return Jet(f0, [None] * u.n, {}, u.n, u.order)
    return u.compose(f0, 1 / u.v, -1 / (u.v * u.v) if u.order >= 2 else ZERO)

def jpos(u):
    # positive part; Lipschitz, derivative in hull{0, u'} (enough for mean-value bounds); no Hessian
    if hi(u.v) <= 0:
        return u.const(ZERO)
    if lo(u.v) >= 0:
        return u
    v = hull(ZERO, arb(hi(u.v)))
    d = [None if a is None else hull(ZERO, a) for a in u.d]
    H = {k: NAN() for k in u.H} if u.order >= 2 else {}
    if u.order >= 2:
        for i in range(u.n):
            for j in range(i, u.n): H[(i, j)] = NAN()
    return Jet(v, d, H, u.n, u.order)

# ---------------- constants ----------------
def f_root(x):
    return -(1 - x) * (1 - x).log() - x * (1 + x).log()

def root_bracket(bits=110):
    # f concave on (0,1), f(0)=0, f'(0)=1, f(1)=-log2 < 0: unique root r; f>0 on (0,r), f<0 on (r,1).
    a, b = fmpq(69, 100), fmpq(70, 100)
    assert f_root(arb(a)) > 0 and f_root(arb(b)) < 0
    for _ in range(bits):
        mid = (a + b) / 2
        fm = f_root(arb(mid))
        if fm > 0: a = mid
        elif fm < 0: b = mid
        else: raise RuntimeError("undecided sign")
    return a, b

RA, RB = root_bracket()
r = arb.union(arb(RA), arb(RB))
m = HALF + (1 - r) * (1 + r).log() / (2 * L)
q = 2 * m - 1
cstar = 1 - m
hr = h_point(r)
p = ((1 - r) / r).log()
c = r - hr / p
tau = (HALF - r * r) / (r * r * (1 - r) ** 2)
qL = q * L
sqrtL = L.sqrt()
rho0 = -(1 - r) * (1 - r).log() / qL    # R(x) >= phi(x) - rho0 x

# ---------------- kernel functions on jets ----------------
def I2(x, y):
    return x * y * (1 + tau * (1 - x) * (1 - y))

def phi(x):
    return jsqrt(jh(I2(x, x)) / L)

def R(x):
    # R(x) = F(r,x)/(qL) with phi(r) = 1 exactly (I(r,r) = 1/2)
    return phi(x) + (jh(x * r) - x * (hr + Q(1, 10**9)) - jh(x) * r) / qL   # RED: h_r + 1e-9 inside R only

def vfun(x):
    return jh(x) - phi(x) * hr

def bfun(x):
    return x - jh(x) / p - phi(x) * c

def F(x, y):
    return jh(x * y) + phi(x) * phi(y) * qL - y * jh(x) - x * jh(y)

def Jk(x, y):
    return jh(I2(x, y)) - phi(x) * phi(y) * L
