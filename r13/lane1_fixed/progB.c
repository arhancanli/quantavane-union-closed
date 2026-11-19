/* Quantavane Program B: replay of the finite inequalities of Sadhu-Ursin (Zenodo 23322347) as reduced in
 * EDGE_LEMMAS.md. Written 16 Nov 2026 separately from Program A (python-flint Arb): binary64 intervals with
 * outward nextafter after every correctly rounded IEEE operation (+ - * / sqrt), MPFR (53-bit, directed rounding)
 * for log, log1p, exp; univariate second-order Taylor jets along a direction (gradient from two passes, Hessian
 * from three passes); midpoint bisection; constants derived again from the paper's definitions.
 * Build: cc -O2 -ffp-contract=off -fno-fast-math progB.c -I$MP/include -L$MP/lib -lmpfr -lgmp -lm -Wl,-rpath,$MP/lib
 * Usage: ./progB <group> [red]   groups: SC S6 S5 K W Z C H G1 G2 G3 G4
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpfr.h>

typedef struct { double lo, hi; } I;
static int RED = 0, FAILED = 0;

static inline double dn(double x) { return nextafter(x, -INFINITY); }
static inline double up(double x) { return nextafter(x, INFINITY); }
static inline I mk(double a, double b) { I r = {a, b}; return r; }
static inline I pt(double a) { return mk(a, a); }
static inline int bad(I a) { return isnan(a.lo) || isnan(a.hi); }
static inline I add(I a, I b) { return mk(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static inline I sub(I a, I b) { return mk(dn(a.lo - b.hi), up(a.hi - b.lo)); }
static inline I neg(I a) { return mk(-a.hi, -a.lo); }
static inline double mn4(double a, double b, double c, double d) { double m = a; if (b < m) m = b; if (c < m) m = c; if (d < m) m = d; return m; }
static inline double mx4(double a, double b, double c, double d) { double m = a; if (b > m) m = b; if (c > m) m = c; if (d > m) m = d; return m; }
static inline I mul(I a, I b) {
    if (bad(a) || bad(b)) return mk(NAN, NAN);
    if ((a.lo == 0 && a.hi == 0) || (b.lo == 0 && b.hi == 0)) return pt(0);
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    return mk(dn(mn4(p1, p2, p3, p4)), up(mx4(p1, p2, p3, p4)));
}
static inline I inv(I a) {
    if (bad(a)) return mk(NAN, NAN);
    if (a.lo <= 0 && a.hi >= 0) return mk(NAN, NAN);
    return mk(dn(1.0 / a.hi), up(1.0 / a.lo));
}
static inline I dvd(I a, I b) { return mul(a, inv(b)); }
static inline I sc(double c, I a) { return mul(pt(c), a); }
static inline I sqr(I a) {
    if (bad(a)) return mk(NAN, NAN);
    if (a.lo >= 0) return mk(dn(a.lo * a.lo), up(a.hi * a.hi));
    if (a.hi <= 0) return mk(dn(a.hi * a.hi), up(a.lo * a.lo));
    double m = fmax(-a.lo, a.hi); return mk(0, up(m * m));
}
static inline I isq(I a) { /* sqrt, argument clipped at 0 from below */
    if (bad(a)) return mk(NAN, NAN);
    double l = a.lo < 0 ? 0 : a.lo; if (a.hi < 0) return mk(NAN, NAN);
    double sl = sqrt(l), sh = sqrt(a.hi); return mk(sl > 0 ? dn(sl) : 0, up(sh));
}
static inline I hull(I a, I b) { if (bad(a) || bad(b)) return mk(NAN, NAN); return mk(fmin(a.lo, b.lo), fmax(a.hi, b.hi)); }
static inline I iabs(I a) { if (bad(a)) return mk(NAN, NAN); if (a.lo >= 0) return a; if (a.hi <= 0) return neg(a); return mk(0, fmax(-a.lo, a.hi)); }
static inline I imax(I a, I b) { if (bad(a) || bad(b)) return mk(NAN, NAN); return mk(fmax(a.lo, b.lo), fmax(a.hi, b.hi)); }
static inline I rat(long n, long d) { return dvd(pt((double)n), pt((double)d)); }

static mpfr_t MT;
static I LN2, rr, mm, qq, hr, pp, cc, tau, qL, sL, rho0, k4;
static double mlog(double x, mpfr_rnd_t rd) { mpfr_set_d(MT, x, MPFR_RNDN); mpfr_log(MT, MT, rd); return mpfr_get_d(MT, rd); }
static double mlog1p(double x, mpfr_rnd_t rd) { mpfr_set_d(MT, x, MPFR_RNDN); mpfr_log1p(MT, MT, rd); return mpfr_get_d(MT, rd); }
static double mexp(double x, mpfr_rnd_t rd) { mpfr_set_d(MT, x, MPFR_RNDN); mpfr_exp(MT, MT, rd); return mpfr_get_d(MT, rd); }
static I ilog(I a) { if (!(a.lo > 0)) return mk(NAN, NAN); return mk(mlog(a.lo, MPFR_RNDD), mlog(a.hi, MPFR_RNDU)); }
static I ilog1p(I a) { if (!(a.lo > -1)) return mk(NAN, NAN); return mk(mlog1p(a.lo, MPFR_RNDD), mlog1p(a.hi, MPFR_RNDU)); }
static I iexp(I a) { return mk(mexp(a.lo, MPFR_RNDD), mexp(a.hi, MPFR_RNDU)); }

/* binary entropy at an exact double point, as an interval */
static I hpt(double x) {
    if (x <= 0 || x >= 1) return pt(0);
    I X = pt(x);
    I t1 = mul(X, ilog(X));                 /* x log x <= 0 */
    I t2 = mul(sub(pt(1), X), ilog1p(pt(-x))); /* (1-x) log(1-x) */
    return neg(add(t1, t2));
}
/* range of h over an interval (h increasing on [0,1/2], decreasing on [1/2,1]) */
static I hR(I u) {
    double a = u.lo < 0 ? 0 : u.lo, b = u.hi > 1 ? 1 : u.hi;
    if (a > b) return mk(NAN, NAN);
    I ha = hpt(a), hb = hpt(b);
    if (b <= 0.5 || a >= 0.5) return hull(ha, hb);
    return hull(hull(ha, hb), LN2);
}
/* h'(u) = log(1-u) - log(u), decreasing */
static I hp1(double x) { return sub(ilog1p(pt(-x)), ilog(pt(x))); }
static I hpR(I u) {
    double a = u.lo, b = u.hi; if (!(a > 0 && b < 1)) return mk(NAN, NAN);
    I A = hp1(b), B = hp1(a); return mk(A.lo, B.hi);
}
/* h''(u) = -1/(u(1-u)) */
static I hppR(I u) {
    double a = u.lo, b = u.hi; if (!(a > 0 && b < 1)) return mk(NAN, NAN);
    I g = hull(mul(pt(a), sub(pt(1), pt(a))), mul(pt(b), sub(pt(1), pt(b))));
    if (a <= 0.5 && b >= 0.5) g = hull(g, pt(0.25));
    return neg(inv(g));
}

/* ---------------- second-order directional jets ---------------- */
typedef struct { I v, d, e; } T;   /* value, first and second derivative along the seed direction */
static inline T tc(I v) { T t = {v, pt(0), pt(0)}; return t; }
static inline T tadd(T a, T b) { T t = {add(a.v, b.v), add(a.d, b.d), add(a.e, b.e)}; return t; }
static inline T tsub(T a, T b) { T t = {sub(a.v, b.v), sub(a.d, b.d), sub(a.e, b.e)}; return t; }
static inline T tsc(I c, T a) { T t = {mul(c, a.v), mul(c, a.d), mul(c, a.e)}; return t; }
static inline T tmul(T a, T b) {
    T t = {mul(a.v, b.v), add(mul(a.v, b.d), mul(a.d, b.v)), add(add(mul(a.v, b.e), mul(b.v, a.e)), sc(2.0, mul(a.d, b.d)))};
    return t;
}
static inline T tcomp(T a, I f0, I f1, I f2) {
    T t; t.v = f0;
    int zd = (a.d.lo == 0 && a.d.hi == 0), ze = (a.e.lo == 0 && a.e.hi == 0);
    t.d = zd ? pt(0) : mul(f1, a.d);
    t.e = (zd && ze) ? pt(0) : add(zd ? pt(0) : mul(f2, sqr(a.d)), ze ? pt(0) : mul(f1, a.e));
    return t;
}
static int ORD = 2; /* 0: values only, 1: first derivative, 2: second */
static T th(T a) {
    I f0 = hR(a.v);
    if (ORD == 0 || (a.d.lo == 0 && a.d.hi == 0 && a.e.lo == 0 && a.e.hi == 0)) return tc(f0);
    return tcomp(a, f0, hpR(a.v), ORD >= 2 ? hppR(a.v) : pt(0));
}
static T tsqrt(T a) {
    I f0 = isq(a.v);
    if (ORD == 0 || (a.d.lo == 0 && a.d.hi == 0 && a.e.lo == 0 && a.e.hi == 0)) return tc(f0);
    if (!(a.v.lo > 0)) { T t = {f0, mk(NAN, NAN), mk(NAN, NAN)}; return t; }
    I f1 = inv(sc(2.0, f0));
    I f2 = neg(inv(sc(4.0, mul(f0, mul(f0, f0)))));
    return tcomp(a, f0, f1, f2);
}
static T tlog(T a) {
    I f0 = ilog(a.v);
    if (ORD == 0 || (a.d.lo == 0 && a.d.hi == 0 && a.e.lo == 0 && a.e.hi == 0)) return tc(f0);
    I f1 = inv(a.v); return tcomp(a, f0, f1, neg(sqr(f1)));
}
static T tpos(T a) {
    if (a.v.hi <= 0) return tc(pt(0));
    if (a.v.lo >= 0) return a;
    T t = {mk(0, a.v.hi), hull(pt(0), a.d), mk(NAN, NAN)}; return t;
}
static T tinvT(T a) { I f0 = inv(a.v); return tcomp(a, f0, neg(sqr(f0)), sc(2.0, mul(f0, sqr(f0)))); }

/* ---------------- constants (derived from the paper's definitions) ---------------- */
static I fr(double x) { /* -(1-x)log(1-x) - x log(1+x) */
    I X = pt(x); return sub(neg(mul(sub(pt(1), X), ilog1p(pt(-x)))), mul(X, ilog1p(X)));
}
static void constants(void) {
    LN2 = ilog(pt(2.0));
    double a = 0.69, b = 0.70;
    if (!(fr(a).lo > 0 && fr(b).hi < 0)) { printf("FAIL root bracket\n"); exit(2); }
    for (int i = 0; i < 200; i++) {
        double c = 0.5 * (a + b); if (c == a || c == b) break;
        I f = fr(c);
        if (f.lo > 0) a = c; else if (f.hi < 0) b = c; else break;
    }
    rr = mk(a, b);
    mm = add(pt(0.5), dvd(mul(sub(pt(1), rr), ilog1p(rr)), sc(2.0, LN2)));
    qq = sub(sc(2.0, mm), pt(1));
    /* h(r) over the bracket: r in [a,b] subset (1/2,1): h decreasing */
    hr = mk(hpt(b).lo, hpt(a).hi);
    pp = hpR(rr);
    cc = sub(rr, dvd(hr, pp));
    tau = dvd(sub(pt(0.5), sqr(rr)), mul(sqr(rr), sqr(sub(pt(1), rr))));
    qL = mul(qq, LN2); sL = isq(LN2);
    rho0 = dvd(neg(mul(sub(pt(1), rr), ilog1p(neg(rr)))), qL);
    k4 = mul(LN2, sub(qq, rat(4, 25)));
}

/* ---------------- kernel functions on jets ---------------- */
static T Ixy(T x, T y) { /* x y (1 + tau (1-x)(1-y)) */
    T one = tc(pt(1));
    return tmul(tmul(x, y), tadd(one, tsc(tau, tmul(tsub(one, x), tsub(one, y)))));
}
static T Phi(T x) { return tsqrt(tsc(inv(LN2), th(Ixy(x, x)))); }
static T Rf(T x) { /* phi(x) + (h(rx) - r h(x) - x h_r)/(qL); phi(r) = 1 */
    T D = tsub(tsub(th(tsc(rr, x)), tsc(rr, th(x))), tsc(hr, x));
    return tadd(Phi(x), tsc(inv(qL), D));
}
static T Vf(T x) { return tsub(th(x), tsc(hr, Phi(x))); }
static T Bf(T x) { return tsub(tsub(x, tsc(inv(pp), th(x))), tsc(cc, Phi(x))); }
static T Ff(T x, T y) { return tsub(tsub(tadd(th(tmul(x, y)), tsc(qL, tmul(Phi(x), Phi(y)))), tmul(y, th(x))), tmul(x, th(y))); }
static T Jf(T x, T y) { return tsub(th(Ixy(x, y)), tsc(LN2, tmul(Phi(x), Phi(y)))); }

static T G1(T x, T y) { I co = RED ? rat(1, 5) : rat(1, 20); return tsub(Ff(x, y), tsc(co, tadd(tmul(Phi(x), Rf(y)), tmul(Phi(y), Rf(x))))); }
static T G2(T x, T y) { I co = RED ? rat(1, 10) : rat(11, 50);
    return tadd(tsub(Jf(x, y), tsc(rat(1, 5), tadd(tmul(Phi(x), Bf(y)), tmul(Phi(y), Bf(x))))), tsc(co, tadd(Rf(x), Rf(y)))); }
static T G3(T x, T y) { I co = RED ? pt(1) : rat(11, 50); T vx = Vf(x), vy = Vf(y);
    return tsub(Ff(x, y), tsc(co, tadd(tmul(tmul(vx, vx), Phi(y)), tmul(tmul(vy, vy), Phi(x))))); }
static T G4(T x, T y) { I co = RED ? pt(1) : rat(4, 25);
    T dx = tsub(tsc(sL, Phi(x)), tsc(rat(7, 5), th(x))), dy = tsub(tsc(sL, Phi(y)), tsc(rat(7, 5), th(y)));
    return tsub(Ff(x, y), tsc(co, tmul(tpos(dx), tpos(dy)))); }

typedef T (*F2)(T, T);
typedef T (*F1)(T);

static void report(const char *name, int pass, const char *extra) {
    printf("%s %s  [%s]\n", pass ? "PASS" : "FAIL", name, extra); fflush(stdout);
    if (!pass) FAILED = 1;
}

/* value of a 2-variable function on a box, and its directional derivative bounds */
static I val2(F2 G, I X, I Y) { int o = ORD; ORD = 0; T r = G(tc(X), tc(Y)); ORD = o; return r.v; }
static I der2(F2 G, I X, I Y, int dir) {
    int o = ORD; ORD = 1;
    T x = {X, pt(dir == 0 ? 1 : 0), pt(0)}, y = {Y, pt(dir == 1 ? 1 : 0), pt(0)};
    T r = G(x, y); ORD = o; return r.d;
}

static long bb2(const char *name, F2 G, double bx0, double bx1, double by0, double by1, double cx0, double cx1, int sym,
                double cornerlo /* skip boxes inside [cornerlo,1]^2 */) {
    typedef struct { double x0, x1, y0, y1; } Bx;
    long cap = 1 << 20, n = 0, sp = 1; Bx *st = malloc(sizeof(Bx) * cap);
    st[0].x0 = bx0; st[0].x1 = bx1; st[0].y0 = by0; st[0].y1 = by1;
    double minlb = INFINITY; clock_t t0 = clock();
    while (sp > 0) {
        Bx b = st[--sp]; n++;
        if (sym && b.x0 > b.y1) continue;
        if (b.x0 >= cornerlo && b.y0 >= cornerlo) continue;
        if (b.x0 >= cx0 && b.x1 <= cx1 && b.y0 >= cx0 && b.y1 <= cx1) continue;
        I X = mk(b.x0, b.x1), Y = mk(b.y0, b.y1);
        I v = val2(G, X, Y);
        if (!bad(v) && v.lo > 0) { if (v.lo < minlb) minlb = v.lo; continue; }
        double xm = 0.5 * (b.x0 + b.x1), ym = 0.5 * (b.y0 + b.y1);
        I vc = val2(G, pt(xm), pt(ym));
        if (!bad(vc) && vc.hi < 0) {
            char buf[300]; snprintf(buf, 300, "counterexample: value at (%.17g, %.17g) in [%g,%g]", xm, ym, vc.lo, vc.hi);
            report(name, 0, buf); free(st); return -1;
        }
        I gx = der2(G, X, Y, 0), gy = der2(G, X, Y, 1);
        /* mean value: G >= G(c) - |Gx| (half width, rounded up) - |Gy| (...) ; half widths bounded above */
        I rx = mk(0, up(up(xm - b.x0) > up(b.x1 - xm) ? up(xm - b.x0) : up(b.x1 - xm)));
        I ry = mk(0, up(up(ym - b.y0) > up(b.y1 - ym) ? up(ym - b.y0) : up(b.y1 - ym)));
        I lb = sub(sub(vc, mul(iabs(gx), pt(rx.hi))), mul(iabs(gy), pt(ry.hi)));
        if (!bad(lb) && lb.lo > 0) { if (lb.lo < minlb) minlb = lb.lo; continue; }
        double wx = b.x1 - b.x0, wy = b.y1 - b.y0;
        if ((wx < 1e-15 && wy < 1e-15) || n > 400000000L) {
            char buf[300]; snprintf(buf, 300, "stuck x=[%.17g,%.17g] y=[%.17g,%.17g] v=[%g,%g]", b.x0, b.x1, b.y0, b.y1, v.lo, v.hi);
            report(name, 0, buf); free(st); return -1;
        }
        if (sp + 2 >= cap) { cap *= 2; st = realloc(st, sizeof(Bx) * cap); }
        /* bisect the side that is wider relative to its distance from the nearest edge of [0,1] */
        double dx = fmin(b.x0, 1 - b.x1), dy = fmin(b.y0, 1 - b.y1);
        double rxw = wx / fmax(dx, wx), ryw = wy / fmax(dy, wy);
        if (rxw >= ryw) { Bx l = b, r = b; l.x1 = xm; r.x0 = xm; st[sp++] = l; st[sp++] = r; }
        else { Bx l = b, r = b; l.y1 = ym; r.y0 = ym; st[sp++] = l; st[sp++] = r; }
    }
    char buf[200]; snprintf(buf, 200, "boxes=%ld min lower=%.3e cpu=%.1fs", n, minlb, (double)(clock() - t0) / CLOCKS_PER_SEC);
    report(name, 1, buf); free(st); return n;
}

static I val1(F1 f, I X) { int o = ORD; ORD = 0; T r = f(tc(X)); ORD = o; return r.v; }
static I der1(F1 f, I X) { int o = ORD; ORD = 1; T x = {X, pt(1), pt(0)}; T r = f(x); ORD = o; return r.d; }
static I dd1(F1 f, I X) { int o = ORD; ORD = 2; T x = {X, pt(1), pt(0)}; T r = f(x); ORD = o; return r.e; }

static int bb1(const char *name, F1 f, double a, double b, double c0, double c1) {
    double st[4096][2]; int sp = 1; long n = 0; st[0][0] = a; st[0][1] = b; double minlb = INFINITY;
    while (sp > 0) {
        sp--; double x0 = st[sp][0], x1 = st[sp][1]; n++;
        if (x0 >= c0 && x1 <= c1) continue;
        I v = val1(f, mk(x0, x1));
        if (!bad(v) && v.lo > 0) { if (v.lo < minlb) minlb = v.lo; continue; }
        double xm = 0.5 * (x0 + x1);
        I vc = val1(f, pt(xm)), g = der1(f, mk(x0, x1));
        if (!bad(vc) && vc.hi < 0) { char buf[200]; snprintf(buf, 200, "counterexample: value at %.17g in [%g,%g]", xm, vc.lo, vc.hi); report(name, 0, buf); return 0; }
        double rad = up(fmax(up(xm - x0), up(x1 - xm)));
        I lb = sub(vc, mul(iabs(g), pt(rad)));
        if (!bad(lb) && lb.lo > 0) { if (lb.lo < minlb) minlb = lb.lo; continue; }
        if (x1 - x0 < 1e-15 || sp > 4000) { char buf[200]; snprintf(buf, 200, "stuck [%.17g,%.17g] v=[%g,%g]", x0, x1, v.lo, v.hi); report(name, 0, buf); return 0; }
        st[sp][0] = x0; st[sp][1] = xm; sp++; st[sp][0] = xm; st[sp][1] = x1; sp++;
    }
    char buf[200]; snprintf(buf, 200, "on [%g,%g] intervals=%ld min lower=%.3e", a, b, n, minlb); report(name, 1, buf); return 1;
}

/* convexity on [k0,k1] (1D) or [k0,k1]^2 (2D) by adaptive cells */
static int conv1(const char *name, F1 f, double k0, double k1) {
    if (!(k0 < rr.lo && rr.hi < k1)) { report(name, 0, "box misses r"); return 0; }
    double st[4096][2]; int sp = 1; st[0][0] = k0; st[0][1] = k1; double mnv = INFINITY; long n = 0;
    while (sp > 0) {
        sp--; double a = st[sp][0], b = st[sp][1]; n++;
        I e = dd1(f, mk(a, b));
        if (!bad(e) && e.lo > 0) { if (e.lo < mnv) mnv = e.lo; continue; }
        if (b - a < 1e-6 || sp > 4000) { report(name, 0, "second derivative not positive"); return 0; }
        double c = 0.5 * (a + b); st[sp][0] = a; st[sp][1] = c; sp++; st[sp][0] = c; st[sp][1] = b; sp++;
    }
    char buf[200]; snprintf(buf, 200, "[%g,%g] cells=%ld min f''=%.4e", k0, k1, n, mnv); report(name, 1, buf); return 1;
}
static I dirsec(F2 G, I X, I Y, double ux, double uy) {
    int o = ORD; ORD = 2; T x = {X, pt(ux), pt(0)}, y = {Y, pt(uy), pt(0)}; T r = G(x, y); ORD = o; return r.e;
}
static int conv2(const char *name, F2 G, double k0, double k1) {
    if (!(k0 < rr.lo && rr.hi < k1)) { report(name, 0, "box misses r"); return 0; }
    typedef struct { double x0, x1, y0, y1; int d; } Bx;
    Bx *st = malloc(sizeof(Bx) * 200000); int sp = 0; long n = 0; double mdet = INFINITY;
    int N = 16;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        Bx b = {k0 + (k1 - k0) * i / N, k0 + (k1 - k0) * (i + 1) / N, k0 + (k1 - k0) * j / N, k0 + (k1 - k0) * (j + 1) / N, 0};
        st[sp++] = b;
    }
    while (sp > 0) {
        Bx b = st[--sp]; n++;
        I X = mk(b.x0, b.x1), Y = mk(b.y0, b.y1);
        I a11 = dirsec(G, X, Y, 1, 0), a22 = dirsec(G, X, Y, 0, 1);
        I sp_ = dirsec(G, X, Y, 1, 1), sm_ = dirsec(G, X, Y, 1, -1);
        I a12 = sc(0.25, sub(sp_, sm_));          /* (D^2_{e1+e2} - D^2_{e1-e2})/4 */
        double det = dn(dn(a11.lo * a22.lo) - up(iabs(a12).hi * iabs(a12).hi));
        if (!bad(a11) && !bad(a22) && !bad(a12) && a11.lo > 0 && a22.lo > 0 && det > 0) { if (det < mdet) mdet = det; continue; }
        if (b.d >= 7) { char buf[300]; snprintf(buf, 300, "cell x[%g,%g] y[%g,%g] a11[%g,%g] a12[%g,%g] a22[%g,%g]", b.x0, b.x1, b.y0, b.y1, a11.lo, a11.hi, a12.lo, a12.hi, a22.lo, a22.hi); report(name, 0, buf); free(st); return 0; }
        double xm = 0.5 * (b.x0 + b.x1), ym = 0.5 * (b.y0 + b.y1);
        Bx c1 = {b.x0, xm, b.y0, ym, b.d + 1}, c2 = {xm, b.x1, b.y0, ym, b.d + 1}, c3 = {b.x0, xm, ym, b.y1, b.d + 1}, c4 = {xm, b.x1, ym, b.y1, b.d + 1};
        st[sp++] = c1; st[sp++] = c2; st[sp++] = c3; st[sp++] = c4;
    }
    char buf[200]; snprintf(buf, 200, "[%g,%g]^2 cells=%ld min det lower=%.4e", k0, k1, n, mdet); report(name, 1, buf); free(st); return 1;
}

/* ------------------------------------------------------------------ groups */
static double D0E = 1e-10, D1E = 2e-7;
static I hD(double x) { return hpt(x); }
static I ellI(double x) { return neg(ilog(pt(x))); }
static I PhiU(I l) { return isq(dvd(mul(add(pt(1), tau), add(sc(2.0, l), sub(pt(1), ilog1p(tau)))), LN2)); }
static int GT(I a, I b) { return a.lo > b.hi; }

static void group_SC(void) {
    printf("== SC: constants and edge scalars (Program B)\n");
    { char b_[120]; snprintf(b_, 120, "[%.17g, %.17g] width %.2e, f(lo)>0>f(hi) certified", rr.lo, rr.hi, rr.hi - rr.lo); report("r bracket", rr.hi - rr.lo < 1e-12, b_); }
    report("0<tau<1/2", tau.lo > 0 && tau.hi < 0.5, "");
    report("q>5/22", GT(qq, rat(5, 22)), "");
    report("44qL>1", sc(44, qL).lo > 1, "");
    report("c>0, p<0", cc.lo > 0 && pp.hi < 0, "");
    report("G empty case: lambda - 2 Gamma eps > 1", sub(rat(1232, 1000), mul(pt(126), rat(92, 1000000))).lo > 1, "");
    I ap = neg(pp);
    report("K1 small", GT(isq(dvd(sc(2, neg(ilog(rat(3, 20)))), LN2)), rho0), "");
    report("K1 near 1", GT(dvd(sqr(qL), mul(sqr(rr), LN2)), hD(0.01)), "");
    I l = ellI(0.005);
    I k2 = sub(add(mul(sub(dvd(hr, ap), rr), l), mul(sub(pt(1), dvd(sqr(hr), ap)), isq(dvd(sc(2, l), LN2)))), rho0);
    report("K2 small", k2.lo > 0 && sub(dvd(hr, ap), rr).lo > 0 && sub(pt(1), dvd(sqr(hr), ap)).lo > 0, "");
    I kp = sub(dvd(rr, qL), sub(add(sub(pt(1), rr), dvd(hr, ap)), pt(0.005)));
    report("K2 near 1", GT(dvd(sub(pt(1), dvd(sqr(hr), ap)), mul(kp, sL)), isq(hD(0.005))), "");
    l = ellI(0.002); I Pu = PhiU(l);
    report("K3 small", sub(add(pt(1), dvd(l, ap)), mul(cc, Pu)).lo > 0 && sub(inv(ap), dvd(mul(cc, add(pt(1), tau)), mul(LN2, Pu))).lo > 0, "");
    report("K3 near 1", sub(sub(pt(1), pt(0.01)), mul(cc, isq(dvd(hD(0.02), LN2)))).lo > 0, "");
    report("K3 1-b small", sub(sub(pt(1), pt(0.01)), dvd(hD(0.01), ap)).lo > 0, "");
    report("K3 1-b near 1", GT(dvd(mul(cc, ap), sL), sL), "");
    report("K4 h(1/100)<27/400", GT(rat(27, 400), hD(0.01)), "");
    report("K4 h(1/5000)<21/10000", GT(rat(21, 10000), hD(1.0 / 5000)), "");
    report("W1 small", mul(rat(2, 5), isq(dvd(sc(2, ilog(pt(10))), LN2))).lo > 1, "");
    report("W1 near 1", GT(dvd(rat(4, 25), LN2), hD(0.04)), "");
    report("W2 small", GT(isq(dvd(sc(2, neg(ilog(pt(0.69)))), LN2)), pt(0.032)) && 0.69 < rr.lo, "");
    l = ellI(0.001);
    I w3a = dvd(mul(pt(0.001), sqr(add(l, pt(1)))), isq(dvd(sc(2, l), LN2)));
    T x1 = tc(pt(0.001)); int o = ORD; ORD = 0; I ph001 = Phi(x1).v; ORD = o;
    report("W3 small", GT(rat(9, 100), w3a) && GT(rat(9, 100), mul(sqr(hr), ph001)), "");
    I h400 = hD(1.0 / 400);
    report("W3 near 1", GT(rat(9, 100), mul(sL, mul(h400, isq(h400)))) && GT(rat(9, 100), mul(sqr(hr), isq(dvd(hD(1.0 / 200), LN2)))), "");
    report("reduction kappa4 order", k4.lo > 0 && GT(sub(qL, rat(1, 10)), k4) && GT(sub(qL, rat(99, 2500)), sub(qL, rat(1, 10))), "");
    report("X0", mul(k4, isq(dvd(sc(2, ellI(D0E)), LN2))).lo >= 0.4, "");
    report("X1", GT(dvd(sqr(mul(pt(0.032), k4)), LN2), hD(D1E)), "");
    I beta = pt(4e-5);
    report("K11", GT(dvd(mul(mul(k4, sub(pt(2), mul(add(pt(1), tau), beta))), ellI(8e-5)), LN2), sc(2, LN2)), "");
    report("Z1", GT(dvd(sqr(dvd(qL, sc(11, rr))), LN2), hD(5e-5)), "");
    I Phi1 = isq(dvd(hD(2 * D1E), LN2)); I la = ellI(1e-6);
    I hb = dvd(la, PhiU(la)), yb = isq(dvd(LN2, sc(2, la)));
    I c1 = sub(sub(pt(1), pt(D1E)), dvd(Phi1, sc(5, ap)));
    I z3 = sub(sub(sub(mul(c1, hb), mul(LN2, Phi1)), mul(rat(1, 5), add(pt(1), dvd(hD(D1E), ap)))), mul(rat(1, 5), mul(Phi1, yb)));
    report("Z3", z3.lo > 0 && c1.lo > 0, "");
    I l1 = ellI(1e-3); I rho1 = sub(pt(1), mul(rho0, isq(dvd(LN2, sc(2, l1))))); I beta1 = add(pt(1e-3), dvd(hD(1e-3), ap));
    report("Z4a", sub(mul(rat(11, 50), rho1), dvd(beta1, pt(5))).lo > 0, "");
    report("Z4 range", mul(add(pt(1), tau), pt(1e-20)).hi <= 0.5, "");
    I kap = rat(39, 100), dl = rat(8, 10000), d = pt(1.0 / 65536);
    report("C1 small", GT(sub(kap, rat(1, 12)), mul(pt(0.06), neg(ilog(pt(0.03))))), "");
    report("C2 x<=d", sub(mul(dl, sub(rat(1, 12), d)), mul(add(sub(kap, rat(1, 12)), d), hD(1.0 / 65536))).lo > 0, "");
    report("C2 x>=63/64", sub(sub(mul(mul(kap, pt(1.5)), mul(hR(rat(2, 3)), rat(63, 64))), mul(dl, rat(11, 12))), mul(add(pt(1), kap), hD(1.0 / 64))).lo > 0, "");
    report("C4", GT(mul(mul(kap, rat(4, 5)), ilog(pt(9))), hR(rat(9, 10))), "");
    report("K5 Boppana 809/500 < golden ratio", GT(dvd(add(pt(1), isq(pt(5))), pt(2)), rat(809, 500)), "");
}


/* 18 Nov 2026 (code-audit fix): edge closures. Each edge set is enlarged by the factor E9 so that it overlaps the
   double-precision end of the adjacent bisection range; each condition is evaluated at the end where it is weakest. */
static const double E9 = 1.000000001;
static void group_SCX(void) {
    printf("== SCX: edge conditions with thresholds enlarged by 1+1e-9 (Program B)\n");
    I ap = neg(pp);
    report("K1 small x<=3/20*E9", GT(isq(dvd(sc(2, ellI(0.15 * E9)), LN2)), rho0), "");
    report("K1 near 1 a<=E9/100", GT(dvd(sqr(qL), mul(sqr(rr), LN2)), hD(0.01 * E9)), "");
    I l = ellI(0.005 * E9);
    I k2 = sub(add(mul(sub(dvd(hr, ap), rr), l), mul(sub(pt(1), dvd(sqr(hr), ap)), isq(dvd(sc(2, l), LN2)))), rho0);
    report("K2 small x<=E9/200", k2.lo > 0, "");
    I kp = sub(dvd(rr, qL), sub(add(sub(pt(1), rr), dvd(hr, ap)), pt(0.005 * E9)));
    report("K2 near 1 a<=E9/200", kp.lo > 0 && GT(dvd(sub(pt(1), dvd(sqr(hr), ap)), mul(kp, sL)), isq(hD(0.005 * E9))), "");
    l = ellI(0.002 * E9); I Pu = PhiU(l);
    report("K3 small x<=E9/500", sub(add(pt(1), dvd(l, ap)), mul(cc, Pu)).lo > 0 && sub(inv(ap), dvd(mul(cc, add(pt(1), tau)), mul(LN2, Pu))).lo > 0, "");
    report("K3 near 1 a<=E9/100", sub(sub(pt(1), pt(0.01 * E9)), mul(cc, isq(dvd(hD(0.02 * E9), LN2)))).lo > 0, "");
    report("K3 1-b small x<=E9/100", sub(sub(pt(1), pt(0.01 * E9)), dvd(hD(0.01 * E9), ap)).lo > 0, "");
    report("K4 h(E9/100)<27/400", GT(rat(27, 400), hD(0.01 * E9)), "");
    report("K4 h(E9/5000)<21/10000", GT(rat(21, 10000), hD(E9 / 5000)), "");
    report("W1 small y<=E9/10", mul(rat(2, 5), isq(dvd(sc(2, ellI(0.1 * E9)), LN2))).lo > 1, "");
    report("W1 near 1 b<=E9/25", GT(dvd(rat(4, 25), LN2), hD(0.04 * E9)), "");
    report("W2 small y<=0.69*E9", GT(isq(dvd(sc(2, ellI(0.69 * E9)), LN2)), pt(0.032)) && 0.69 * E9 < rr.lo, "");
    l = ellI(0.001 * E9);
    I w3a = dvd(mul(pt(0.001 * E9), sqr(add(l, pt(1)))), isq(dvd(sc(2, l), LN2)));
    T x1 = tc(pt(0.001 * E9)); int o = ORD; ORD = 0; I ph001 = Phi(x1).v; ORD = o;
    report("W3 small x<=E9/1000", GT(rat(9, 100), w3a) && GT(rat(9, 100), mul(sqr(hr), ph001)), "");
    I h400 = hD(E9 / 400);
    report("W3 near 1 a<=E9/400", GT(rat(9, 100), mul(sL, mul(h400, isq(h400)))) && GT(rat(9, 100), mul(sqr(hr), isq(dvd(hD(E9 / 200), LN2)))), "");
    report("X0 x<=1e-10*E9", mul(k4, isq(dvd(sc(2, ellI(D0E * E9)), LN2))).lo >= 0.4, "");
    report("X1 a<=2e-7*E9", GT(dvd(sqr(mul(pt(0.032), k4)), LN2), hD(D1E * E9)), "");
    I beta = pt(4e-5 * E9);
    report("K11 beta=4e-5*E9", GT(dvd(mul(mul(k4, sub(pt(2), mul(add(pt(1), tau), beta))), ellI(8e-5 * E9)), LN2), sc(2, LN2)), "");
    report("Z1 z1=5e-5*E9", GT(dvd(sqr(dvd(qL, sc(11, rr))), LN2), hD(5e-5 * E9)), "");
    double d1 = D1E * E9;
    I Phi1 = isq(dvd(hD(2 * d1), LN2)); I la = ellI(1e-6 * E9);
    I hb = dvd(la, PhiU(la)), yb = isq(dvd(LN2, sc(2, la)));
    I c1 = sub(sub(pt(1), pt(d1)), dvd(Phi1, sc(5, ap)));
    I z3 = sub(sub(sub(mul(c1, hb), mul(LN2, Phi1)), mul(rat(1, 5), add(pt(1), dvd(hD(d1), ap)))), mul(rat(1, 5), mul(Phi1, yb)));
    report("Z3 delta1*E9, y_a*E9", z3.lo > 0 && c1.lo > 0, "");
    I l1 = ellI(1e-3 * E9); I rho1 = sub(pt(1), mul(rho0, isq(dvd(LN2, sc(2, l1))))); I beta1 = add(pt(1e-3 * E9), dvd(hD(1e-3 * E9), ap));
    report("Z4a y1*E9", sub(mul(rat(11, 50), rho1), dvd(beta1, pt(5))).lo > 0, "");
    I kap = rat(39, 100), dl = rat(8, 10000), dd = pt(E9 / 65536);
    report("C1 small x<=0.03*E9", GT(sub(kap, rat(1, 12)), sc(2, mul(pt(0.03 * E9), neg(ilog(pt(0.03 * E9)))))), "");
    report("C2 x<=2^-16*E9", sub(mul(dl, sub(rat(1, 12), dd)), mul(add(sub(kap, rat(1, 12)), dd), hD(E9 / 65536))).lo > 0, "");
    report("C2 x>=1-E9/64", sub(sub(mul(mul(kap, pt(1.5)), mul(hR(rat(2, 3)), sub(pt(1), pt(E9 / 64)))), mul(dl, rat(11, 12))), mul(add(pt(1), kap), hD(E9 / 64))).lo > 0, "");
}

/* S6: the scalar margin. All constants re-derived from the paper text. */
static I EPS;
static int TMIN_OVR = 0; static I TMIN_V;
static I M_of_v(I v) {
    I L = LN2, r = rr, m = mm, q = qq, p = pp, c = cc;
    I z = rat(77, 200), lam = rat(1232, 1000), Gam = pt(63), kF = rat(1, 20), kE = rat(11, 50);
    I Db = dvd(pt(11), sub(pt(1), inv(mul(pt(44), mul(q, L)))));
    I gam = dvd(hr, r), alpha = sub(pt(1), dvd(gam, p)), t0 = sub(r, m);
    I bm = neg(dvd(mul(hr, gam), mul(p, alpha))), D = sub(neg(dvd(hr, mul(p, alpha))), t0);
    I es = mul(gam, m), tmin = TMIN_OVR ? TMIN_V : dvd(sub(z, rat(27, 400)), rat(4, 5)), Cb = dvd(Db, sqr(tmin)), as = dvd(m, r);
    I ell = dvd(EPS, m);
    I k0 = dvd(add(inv(as), dvd(EPS, mul(mul(c, tmin), as))), es);
    I k1 = add(k0, inv(mul(mul(mul(neg(p), c), tmin), as)));
    I k2 = dvd(Cb, mul(mul(c, tmin), as));
    I BE = dvd(sqr(dvd(alpha, sub(pt(1), ell))), sc(2, kE));
    I Ds = (sub(L, es).lo > sub(es, z).hi) ? sub(L, es) : sub(es, z);
    I bY = rat(2, 5);
    I u = mul(sc(2, EPS), sqr(v)), su = mul(isq(sc(2, EPS)), v);
    I C0 = dvd(add(add(dvd(mul(gam, sub(EPS, sc(0.5, u))), z), dvd(mul(gam, EPS), es)), mul(mul(gam, Cb), u)), sub(pt(1), ell));
    I Ebar = add(add(C0, sc(0.5, mul(mul(BE, k1), u))), mul(mul(isq(BE), su), isq(add(add(add(k0, mul(k2, u)), mul(k1, C0)), sc(0.25, mul(mul(BE, sqr(k1)), u))))));
    I n = dvd(add(pt(1), Ebar), es);
    I Kf = add(add(k0, mul(k1, Ebar)), mul(k2, u));
    I A2 = add(mul(dvd(mul(gam, Ds), sub(pt(1), bY)), add(dvd(sub(EPS, sc(0.5, u)), z), mul(EPS, n))),
               mul(dvd(u, sub(pt(1), bY)), add(dvd(sqr(alpha), mul(mul(sc(8, bY), kE), tmin)), mul(mul(gam, sub(es, z)), Cb))));
    I B0 = add(add(mul(dvd(mul(q, L), mul(sqr(p), sqr(c))), A2), sc(2, EPS)), dvd(mul(mul(q, L), sqr(EPS)), mul(sqr(c), z)));
    I hr2 = hR(sqr(r));
    I astar = dvd(mul(m, hr2), mul(r, hr));
    I d0 = add(sub(sc(2, m), astar), dvd(mul(sc(2, mul(q, L)), as), mul(p, c)));
    I D0 = sub(lam, astar), rD = neg(dvd(d0, D0));
    I rho0_ = dvd(add(add(sub(pt(1), astar), B0), mul(Gam, sub(sc(2, EPS), u))), D0);
    I kap = rat(39, 100), C = rat(13, 25), a0 = rat(3, 8);
    I beta = sub(pt(1), dvd(D, kap));
    I b0 = iexp(mul(beta, ilog(a0)));
    I W = dvd(sub(kap, D), sub(pt(1), b0));
    I eta = mul(dvd(sub(kap, D), D), dvd(sub(b0, a0), sub(pt(1), b0)));
    I kc = sub(t0, dvd(mul(bm, EPS), es)), kh = mul(eta, kc), A0 = sc(0.5, W), Bc = mul(W, sub(a0, pt(0.5)));
    I Cm = add(mul(bm, Cb), dvd(add(pt(1), C), mul(sc(2, kF), tmin)));
    I Btot = add(add(mul(eta, add(kc, D)), Bc), mul(add(kh, m), rD));
    if (!(D0.lo > 0 && add(kh, m).lo > 0 && Btot.lo > 0 && Btot.hi < A0.lo && kc.lo > 0 && D.lo > 0 && D.hi < kap.lo && W.lo > 0 && eta.lo > 0
          && sub(sc(2, sub(m, EPS)), lam).lo > 0 && sub(Gam, sc(12.5, q)).lo > 0 && sub(es, z).lo > 0))
        return mk(-1, -1);
    I Qf = add(add(add(dvd(mul(mul(gam, add(A0, Btot)), sub(EPS, sc(0.5, u))), z), mul(mul(mul(gam, sub(A0, Btot)), EPS), n)),
                   mul(mul(mul(add(A0, Btot), alpha), su), isq(dvd(Kf, sc(2, kE))))), mul(mul(mul(sub(A0, Btot), gam), Cb), u));
    return sub(sub(sub(kh, mul(add(kh, m), rho0_)), Qf), mul(Cm, u));
}
static void group_S6(void) {
    printf("== S6: M(u)>0 on [0,2eps]%s (Program B)\n", RED ? " [RED eps=93e-6]" : "");
    EPS = RED ? rat(93, 1000000) : rat(92, 1000000);
    double st[200][2]; int sp = 1; st[0][0] = 0; st[0][1] = 1; long n = 0; double mn = INFINITY;
    while (sp > 0) {
        sp--; double a = st[sp][0], b = st[sp][1]; n++;
        I v = M_of_v(mk(a, b));
        if (v.lo > 0) { if (v.lo < mn) mn = v.lo; continue; }
        if (b - a < 1e-9 || sp > 190) { char buf[200]; snprintf(buf, 200, "v in [%g,%g] M in [%g,%g]", a, b, v.lo, v.hi); report("S6", 0, buf); return; }
        double c = 0.5 * (a + b); st[sp][0] = a; st[sp][1] = c; sp++; st[sp][0] = c; st[sp][1] = b; sp++;
    }
    char buf[200]; snprintf(buf, 200, "intervals=%ld min lower=%.4e (also all sign conditions)", n, mn); report("S6 M(u)>0", 1, buf);
}

static long S5EPS_PPM = 0;
static I PsiB(I e, I t, I s, I lam) {
    I q = qq, L = LN2, p = pp, c = cc;
    I kE = rat(11, 50), Gam = pt(63), Db = dvd(pt(11), sub(pt(1), inv(mul(pt(44), mul(q, L)))));
    I se = mul(s, e), qLt2 = mul(mul(q, L), sqr(t));
    I U0 = imax(mul(rat(809, 500), se), sub(sc(2, se), qLt2));
    I a1 = sub(add(U0, qLt2), sc(2, se));
    I a2 = mul(mul(sc(2, kE), t), sqr(sub(e, mul(hr, t))));
    I a3 = mul(dvd(sqr(t), Db), sub(sub(s, dvd(e, p)), mul(c, t)));
    I w = sub(mul(sL, t), mul(rat(7, 5), e));
    I a4 = w.hi <= 0 ? pt(0) : mul(rat(4, 25), sqr(imax(w, pt(0))));
    I mxv = imax(imax(a1, a2), imax(a3, a4));
    return sub(add(U0, mul(Gam, mxv)), mul(lam, e));
}
static void group_S5(void) {
    printf("== S5: Psi>=0 on the constrained box%s (Program B)\n", RED ? " [RED lambda=1.240]" : "");
    I eps = S5EPS_PPM ? rat(S5EPS_PPM, 1000000) : rat(92, 1000000), lam = rat(1232, 1000), lamU = (RED && !S5EPS_PPM) ? rat(1240, 1000) : lam;
    I e1 = dvd(sub(sc(2, sub(mm, eps)), lam), mul(sc(2, qq), rat(49, 25)));
    I smin = sub(mm, eps);
    typedef struct { double e0, e1, t0, t1, s0, s1; } Bx;
    long cap = 1 << 22; Bx *st = malloc(sizeof(Bx) * cap); long sp = 0, n = 0, disc = 0; double mn = INFINITY;
    Bx b0 = {e1.lo, 77.0 / 200, 0, 1, smin.lo, 1}; st[sp++] = b0;
    if (!(e1.lo <= e1.hi && smin.lo <= smin.hi)) { report("S5", 0, "bad constants"); return; }
    while (sp > 0) {
        Bx b = st[--sp]; n++;
        I E = mk(b.e0, b.e1), Tt = mk(b.t0, b.t1), S = mk(b.s0, b.s1);
        I c1 = sub(sub(pt(b.e0), mul(rat(4, 5), pt(b.t1))), rat(27, 400)), c2 = sub(sub(pt(b.e0), mul(rat(6, 5), pt(b.t1))), rat(21, 10000));
        if (c1.lo > 0 || c2.lo > 0) { disc++; continue; }
        I v = PsiB(E, Tt, S, lamU);
        if (!bad(v) && v.lo > 0) { if (v.lo < mn) mn = v.lo; continue; }
        if (b.e0 >= e1.hi) {
            I c1p = sub(sub(pt(b.e0), mul(rat(4, 5), pt(b.t1))), rat(27, 400)), c2p = sub(sub(pt(b.e0), mul(rat(6, 5), pt(b.t1))), rat(21, 10000));
            I pv = PsiB(pt(b.e0), pt(b.t1), pt(b.s0), lamU);
            if (c1p.hi <= 0 && c2p.hi <= 0 && pv.hi < 0) { char buf[300]; snprintf(buf, 300, "counterexample at e=%g t=%g s=%g Psi in [%g,%g]", b.e0, b.t1, b.s0, pv.lo, pv.hi); report("S5", 0, buf); free(st); return; }
        }
        if (n > 2000000000L || (b.e1 - b.e0 < 1e-14 && b.t1 - b.t0 < 1e-14)) {
            char buf[300]; snprintf(buf, 300, "stuck e[%g,%g] t[%g,%g] s[%g,%g] v[%g,%g]", b.e0, b.e1, b.t0, b.t1, b.s0, b.s1, v.lo, v.hi); report("S5", 0, buf); free(st); return; }
        if (sp + 2 >= cap) { cap *= 2; st = realloc(st, sizeof(Bx) * cap); }
        double we = (b.e1 - b.e0) / 0.4, wt = b.t1 - b.t0, ws = (b.s1 - b.s0) / 0.4;
        Bx l = b, r = b;
        if (we >= wt && we >= ws) { double md = 0.5 * (b.e0 + b.e1); l.e1 = md; r.e0 = md; }
        else if (wt >= ws) { double md = 0.5 * (b.t0 + b.t1); l.t1 = md; r.t0 = md; }
        else { double md = 0.5 * (b.s0 + b.s1); l.s1 = md; r.s0 = md; }
        st[sp++] = l; st[sp++] = r;
    }
    char buf[200]; snprintf(buf, 200, "boxes=%ld discarded=%ld min lower=%.3e", n, disc, mn); report("S5 Psi>0", 1, buf); free(st);
}

/* one-variable functions */
static T fR(T x) { return Rf(x); }
static T fS1b(T x) { return tadd(tsub(tmul(tsub(x, tc(rr)), th(x)), tsc(dvd(hr, pp), Vf(x))), Rf(x)); }
static T fB(T x) { return Bf(x); }
static T f1mB(T x) { return tsub(tc(pt(1)), Bf(x)); }
static T fH1(T x) { return tsub(tadd(tsc(rat(4, 5), Phi(x)), tc(rat(27, 400))), th(x)); }
static T fH2(T x) { return tsub(tadd(tsc(rat(6, 5), Phi(x)), tc(rat(21, 10000))), th(x)); }
static T fBop(T x) { return tsub(th(tmul(x, x)), tsc(rat(809, 500), tmul(x, th(x)))); }
static T fK6(T x) { return tsub(tsc(rat(7, 5), th(x)), tsc(sL, Phi(x))); }
static T fRred(T x) { return tsub(Rf(x), tc(rat(1, 1000))); }
static T fW1(T y) { T one = tc(pt(1)); T z = tsub(one, y); return tadd(tsc(rat(2, 5), Phi(y)), tmul(z, tlog(z))); }
static T fW1red(T y) { T one = tc(pt(1)); T z = tsub(one, y); return tadd(tsc(rat(1, 5), Phi(y)), tmul(z, tlog(z))); }
static T fW2(T y) { return tsub(Phi(y), tsc(pt(0.032), y)); }
static T fW3(T x) { T v = Vf(x); return tsub(tsc(rat(9, 100), Phi(x)), tmul(v, v)); }

static void group_K(void) {
    printf("== K: one-variable kernel inequalities (Program B)\n");
    if (RED) { bb1("RED K1 R-0.001", fRred, 0.15, 0.99, 2, 2); return; }
    conv1("H1 R''>0", fR, 0.66, 0.72); bb1("K1 R>0", fR, 0.15, 0.99, 0.66, 0.72);
    conv1("H2 S1b''>0", fS1b, 0.66, 0.72); bb1("K2 S1b>0", fS1b, 0.005, 0.995, 0.66, 0.72);
    conv1("H3 b''>0", fB, 0.66, 0.72); bb1("K3 b>0", fB, 0.002, 0.99, 0.66, 0.72);
    bb1("K3 1-b>0", f1mB, 0.01, 0.75, 2, 2);
    bb1("K4 (4/5)phi+27/400-h", fH1, 0.01, 0.99, 2, 2);
    bb1("K4 (6/5)phi+21/10000-h", fH2, 1.0 / 5000, 1 - 1.0 / 5000, 2, 2);
    bb1("K5 h(t^2)-(809/500)t h(t)", fBop, 0.001, 0.999, 2, 2);
    bb1("K6 7h/5 - sqrt(L)phi on [0.62,0.719]", fK6, 0.62, 0.719, 2, 2);
}
static void group_W(void) {
    printf("== W: helper inequalities (Program B)\n");
    if (RED) { bb1("RED W1 (1/5)phi-omega", fW1red, 0.1, 0.96, 2, 2); return; }
    bb1("W1 (2/5)phi-omega", fW1, 0.1, 0.96, 2, 2);
    bb1("W2 phi-0.032y", fW2, 0.69, 1 - 4e-5, 2, 2);
    bb1("W3 (9/100)phi-v^2", fW3, 0.001, 1 - 1.0 / 400, 2, 2);
}
static I Z2Phi1, Z2hd1;
static T fZ2(T y) {
    I ap = neg(pp);
    T a = tsub(th(y), tc(Z2hd1));
    T b = tsc(mul(LN2, Z2Phi1), Phi(y));
    T c = tsc(dvd(Z2Phi1, pt(5)), Bf(y));
    T d = tsc(dvd(add(pt(1), dvd(Z2hd1, ap)), pt(5)), Phi(y));
    T e = tsc(rat(11, 50), Rf(y));
    T r = tadd(tsub(tsub(tsub(a, b), c), d), e);
    if (RED) r = tsub(r, tc(rat(1, 1000)));
    return r;
}
static int Z4b(double y0, double y1) {
    I ap = neg(pp); I Y = mk(y0, y1); int o = ORD; ORD = 0;
    I ph = Phi(tc(Y)).v, bb = Bf(tc(Y)).v, RR = Rf(tc(Y)).v; ORD = o;
    I l0 = ellI(D0E), Pu0 = PhiU(l0);
    I Klo = add(pt(1), mul(mul(tau, sub(pt(1), pt(D0E))), sub(pt(1), Y)));
    I yK = mul(Y, Klo);
    I A = sub(yK, dvd(ph, sc(5, ap)));
    I P = sub(sub(rat(11, 50), mul(LN2, ph)), dvd(bb, pt(5)));
    I lyK = neg(ilog(yK));
    I base = sub(sub(mul(yK, add(l0, lyK)), mul(dvd(ph, pt(5)), add(pt(1), dvd(add(l0, pt(1)), ap)))), mul(rat(11, 50), rho0));
    I pm = mul(P, isq(dvd(sc(2, l0), LN2))), pu = mul(P, Pu0);
    double mnlo = fmin(pm.lo, pu.lo);
    I Lam = add(base, pt(mnlo));
    double Pneg = fmax(0, -P.lo);
    int beta_ok = Lam.lo > 0 && sub(A, dvd(mul(pt(Pneg), add(pt(1), tau)), mul(LN2, Pu0))).lo >= 0;
    I C0 = add(add(iabs(mul(yK, ilog(yK))), mul(dvd(ph, pt(5)), add(pt(1), inv(ap)))), mul(rat(11, 50), rho0));
    I need = mul(pt(D0E), add(add(mul(iabs(A), l0), C0), mul(iabs(P), Pu0)));
    int alpha_ok = sub(mul(rat(11, 50), RR), need).lo > 0;
    return beta_ok || alpha_ok;
}
static void group_Z(void) {
    printf("== Z: edge checks for G2 (Program B)\n");
    Z2Phi1 = isq(dvd(hD(2 * D1E), LN2)); Z2hd1 = hD(D1E);
    bb1("Z2 on [1e-6, 1-5e-5]", fZ2, 1e-6, 1 - 5e-5, 2, 2);
    if (RED) return;
    double st[4096][2]; int sp = 1; st[0][0] = 1e-3; st[0][1] = 1; long n = 0;
    while (sp > 0) {
        sp--; double a = st[sp][0], b = st[sp][1]; n++;
        if (Z4b(a, b)) continue;
        if (b - a < 1e-13 || sp > 4000) { report("Z4b", 0, "stuck"); return; }
        double c = 0.5 * (a + b); st[sp][0] = a; st[sp][1] = c; sp++; st[sp][0] = c; st[sp][1] = b; sp++;
    }
    char buf[100]; snprintf(buf, 100, "intervals=%ld", n); report("Z4b on [1e-3,1] (alpha or beta)", 1, buf);
}
static I KAPc;
static T fC1(T x) { T one = tc(pt(1)); T k = tsub(tmul(tsc(pt(2), tmul(x, tsub(one, x))), tlog(tmul(tsub(one, x), tinvT(x)))), x); return tsub(tc(sub(KAPc, rat(1, 12))), k); }
static T Gcov(T a, T x) {
    I dl = rat(8, 10000), C = rat(13, 25);
    T t1 = tmul(tmul(tsc(KAPc, th(a)), x), tinvT(a));
    T t2 = tsc(dl, tsub(a, x));
    T t3 = tmul(tadd(tsub(x, a), tc(KAPc)), th(x));
    return tadd(tsub(tadd(t1, t2), t3), tsc(C, Rf(x)));
}
static T Pkl(T a, T x) {
    T one = tc(pt(1));
    T KL = tsub(tsub(tsc(pt(-1), th(x)), tmul(x, tlog(a))), tmul(tsub(one, x), tlog(tsub(one, a))));
    return tsub(tsc(KAPc, KL), tmul(tsub(a, x), tsub(th(a), th(x))));
}
static void group_C(void) {
    printf("== C: covariance inequality (Program B)\n");
    KAPc = RED ? rat(37, 100) : rat(39, 100);
    bb1("C1 kappa-1/12-k(x)", fC1, 0.03, 0.5, 2, 2);
    bb2("C2 Gcov>0 on [1/12,2/3]x[2^-16,63/64]", Gcov, 1.0 / 12, 2.0 / 3, 1.0 / 65536, 63.0 / 64, 2, 2, 0, 2);
    /* C3 on {2/3<=a<=9/10, 0<=x<=1-a}: boxes with x0 >= 1-a0 are outside */
    typedef struct { double a0, a1, x0, x1; } Bx; Bx st[100000]; int sp = 0; long n = 0; double mn = INFINITY;
    Bx b0 = {2.0 / 3, 0.9, 0, 0.34}; st[sp++] = b0;   /* 18 Nov: 0.34 > 1/3 closes the sliver at a = 2/3, x = 1/3 */
    while (sp > 0) {
        Bx b = st[--sp]; n++;
        if (b.x0 >= 1 - b.a0) continue;
        I A = mk(b.a0, b.a1), X = mk(b.x0, b.x1);
        I v = val2(Pkl, A, X);
        if (!bad(v) && v.lo > 0) { if (v.lo < mn) mn = v.lo; continue; }
        I vc = val2(Pkl, pt(0.5 * (b.a0 + b.a1)), pt(0.5 * (b.x0 + b.x1)));
        I ga = der2(Pkl, A, X, 0), gx = der2(Pkl, A, X, 1);
        double am_ = 0.5 * (b.a0 + b.a1), xm_ = 0.5 * (b.x0 + b.x1);
        double ra_ = up(fmax(up(am_ - b.a0), up(b.a1 - am_))), rx_ = up(fmax(up(xm_ - b.x0), up(b.x1 - xm_)));
        I lb = sub(sub(vc, mul(iabs(ga), pt(ra_))), mul(iabs(gx), pt(rx_)));
        if (!bad(lb) && lb.lo > 0) { if (lb.lo < mn) mn = lb.lo; continue; }
        if (sp > 99990 || (b.a1 - b.a0 < 1e-13 && b.x1 - b.x0 < 1e-13)) { report("C3", 0, "stuck"); return; }
        Bx l = b, r = b;
        if (b.a1 - b.a0 > b.x1 - b.x0) { double md = 0.5 * (b.a0 + b.a1); l.a1 = md; r.a0 = md; } else { double md = 0.5 * (b.x0 + b.x1); l.x1 = md; r.x0 = md; }
        st[sp++] = l; st[sp++] = r;
    }
    char buf[100]; snprintf(buf, 100, "boxes=%ld min lower=%.3e", n, mn); report("C3 kappa KL - (a-x)(h(a)-h(x)) > 0", 1, buf);
}
static int has0(I a) { return !bad(a) && a.lo <= 0 && 0 <= a.hi; }
static void group_H(void) {
    printf("== RB: contact identities at the r-ball (Program B)\n");
    { F2 Gs[4] = {G1, G2, G3, G4}; const char *nm[4] = {"RB G1", "RB G2", "RB G3", "RB G4"};
      for (int k = 0; k < 4; k++) report(nm[k], has0(val2(Gs[k], rr, rr)) && has0(der2(Gs[k], rr, rr, 0)) && has0(der2(Gs[k], rr, rr, 1)), "value and gradient contain 0"); }
    report("RB R", has0(val1(fR, rr)) && has0(der1(fR, rr)), "value and derivative contain 0");
    report("RB S1b", has0(val1(fS1b, rr)) && has0(der1(fS1b, rr)), "value and derivative contain 0");
    report("RB b", has0(val1(fB, rr)) && has0(der1(fB, rr)), "value and derivative contain 0");
    printf("== H: contact Hessians (Program B)\n");
    conv2("H4 G1", G1, 0.675, 0.705); conv2("H5 G2", G2, 0.675, 0.705); conv2("H6 G3", G3, 0.675, 0.705); conv2("H7 G4", G4, 0.675, 0.705);
}
static void group_G(const char *w) {
    F2 G = !strcmp(w, "G1") ? G1 : !strcmp(w, "G2") ? G2 : !strcmp(w, "G3") ? G3 : G4;
    double corner = !strcmp(w, "G2") ? 5e-5 : 4e-5;
    printf("== %s on [1e-10,1-2e-7]^2 minus [1-%g,1]^2, contact [0.675,0.705]^2%s (Program B)\n", w, corner, RED ? " [RED]" : "");
    char nm[100]; snprintf(nm, 100, "%s > 0 off contact", w);
    bb2(nm, G, D0E, 1 - D1E, D0E, 1 - D1E, 0.675, 0.705, 1, 1 - corner);
}

static void group_DBG(void) {
    ORD = 0;
    double xs[][2] = {{0.15, 0.99}, {0.6, 0.8}, {0.69, 0.692}, {0.5, 0.5}};
    for (int i = 0; i < 4; i++) {
        I X = mk(xs[i][0], xs[i][1]);
        I hx = th(tc(X)).v, ph = Phi(tc(X)).v, R = Rf(tc(X)).v, hrx = th(tsc(rr, tc(X))).v;
        printf("x[%g,%g]: h[%g,%g] phi[%g,%g] h(rx)[%g,%g] R[%g,%g]\n", X.lo, X.hi, hx.lo, hx.hi, ph.lo, ph.hi, hrx.lo, hrx.hi, R.lo, R.hi);
    }
    I a = mk(0.15, 0.99), b = mk(0.1, 0.2); I s_ = sub(a, b); printf("sub [%g,%g]\n", s_.lo, s_.hi);
    printf("hr [%g,%g] qL [%g,%g] LN2 [%.17g,%.17g]\n", hr.lo, hr.hi, qL.lo, qL.hi, LN2.lo, LN2.hi);
}

/* Variant (16 Nov): t_min replaced by t_min'' = c0 + a(m - eps) + b z from the certificate phi >= c0 + a x + b h. */
static I VC0, VA, VB;
static T fTM(T x) { return tsub(tsub(tsub(Phi(x), tc(VC0)), tsc(VA, x)), tsc(VB, th(x))); }
static void group_V(long ppm) {
    VC0 = rat(-78949, 250000); VA = rat(19731, 62500); VB = rat(151707, 100000);
    long ppm_run = RED ? ppm + 2 : ppm;
    printf("== V: variant t_min'' at eps = %ld e-6%s (Program B)\n", ppm_run, RED ? " [RED: eps raised by 2e-6, t_min'' kept]" : "");
    report("a,b >= 0", VA.lo >= 0 && VB.lo >= 0, "");
    report("edge h(1/4) <= 1/(b^2 L)", GT(inv(mul(sqr(VB), LN2)), hR(rat(1, 4))) , "");
    report("edge -c0 - a > 0", sub(neg(VC0), VA).lo > 0, "");
    if (!RED) bb1("TM phi - c0 - a x - b h > 0 on [0,3/4]", fTM, 0.0, 0.75, 2, 2);
    TMIN_V = add(add(VC0, mul(VA, sub(mm, rat(ppm, 1000000)))), mul(VB, rat(77, 200)));
    TMIN_OVR = 1;
    { char b_[160]; snprintf(b_, 160, "t_min'' in [%.15f, %.15f]", TMIN_V.lo, TMIN_V.hi); report("t_min'' > paper t_min", TMIN_V.lo > dvd(sub(rat(77, 200), rat(27, 400)), rat(4, 5)).hi, b_); }
    EPS = rat(ppm_run, 1000000);
    report("lambda - 2 Gamma eps > 1", sub(rat(1232, 1000), mul(pt(126), EPS)).lo > 1, "");
    report("2(m-eps) > lambda", sub(sc(2, sub(mm, EPS)), rat(1232, 1000)).lo > 0, "");
    double st[400][2]; int sp = 1; st[0][0] = 0; st[0][1] = 1; long n = 0; double mn = INFINITY;
    while (sp > 0) {
        sp--; double a = st[sp][0], b = st[sp][1]; n++;
        I v = M_of_v(mk(a, b));
        if (v.lo > 0) { if (v.lo < mn) mn = v.lo; continue; }
        if (b - a < 1e-10 || sp > 390) { char buf[200]; snprintf(buf, 200, "v in [%g,%g] M in [%g,%g]", a, b, v.lo, v.hi); report("S6v", 0, buf); return; }
        double c = 0.5 * (a + b); st[sp][0] = a; st[sp][1] = c; sp++; st[sp][0] = c; st[sp][1] = b; sp++;
    }
    { char buf[200]; snprintf(buf, 200, "intervals=%ld min lower=%.4e (all sign conditions inside)", n, mn); report("S6v M(u)>0", 1, buf); }
    S5EPS_PPM = ppm_run; group_S5();
}
int main(int argc, char **argv) {
    mpfr_init2(MT, 53);
    const char *g = argc > 1 ? argv[1] : "SC";
    RED = argc > 2 && !strcmp(argv[2], "red");
    constants();
    time_t now = time(NULL);
    printf("Program B (binary64 intervals + MPFR %s) group %s red=%d start %s", mpfr_get_version(), g, RED, ctime(&now));
    if (!strcmp(g, "DBG")) group_DBG();
    else if (!strcmp(g, "V")) group_V(argc > 3 ? atol(argv[3]) : 98);
    else if (!strcmp(g, "SC")) group_SC();
    else if (!strcmp(g, "S6")) group_S6();
    else if (!strcmp(g, "S5")) group_S5();
    else if (!strcmp(g, "K")) group_K();
    else if (!strcmp(g, "W")) group_W();
    else if (!strcmp(g, "Z")) group_Z();
    else if (!strcmp(g, "C")) group_C();
    else if (!strcmp(g, "SCX")) group_SCX();
    else if (!strcmp(g, "ZX")) { D0E *= E9; D1E *= E9; group_Z(); }
    else if (!strcmp(g, "H")) group_H();
    else group_G(g);
    now = time(NULL);
    printf("END group %s red=%d %s %s", g, RED, FAILED ? "FAILED" : "ALL PASS", ctime(&now));
    return FAILED;
}
