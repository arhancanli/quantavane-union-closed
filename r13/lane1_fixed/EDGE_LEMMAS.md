# Quantavane replay of Sadhu-Ursin: our own reductions of every finite inequality to checked pieces

Written 16 Nov 2026 (Quantavane union-closed lane). These reductions are ours; the authors' edge arguments
exist only as comments in their C++ programs. Every numbered "check" (SCn, Kn, Wn, Zn, Cn, Hn, Gn, Sn) is run
by both Program A (`progA_*.py`, Arb balls) and Program B (`progB.c`, binary64 intervals + MPFR).

Notation: h(u) = -u log u - (1-u) log(1-u), L = log 2, ell(x) = log(1/x), r the root, m, q, p = h'(r) < 0,
h_r = h(r), c = r - h_r/p, tau = (1/2 - r^2)/(r^2 (1-r)^2), I(x,y) = xy[1 + tau(1-x)(1-y)],
phi(x) = sqrt(h(I(x,x))/L), R, v, b, F, J as in the paper. B(x,y) = h(xy) - y h(x) - x h(y), so F = B + qL phi phi.
omega(y) = -(1-y) log(1-y) >= 0. rho0 = -(1-r) log(1-r)/(qL). kappa4 = L(q - 4/25).

## Analytic facts (no computation)

A1. h is increasing on [0,1/2], decreasing on [1/2,1], h(1-u) = h(u), h <= L, h concave, h(0) = 0, so h is
    subadditive: h(u+v) <= h(u) + h(v) for u,v >= 0, u+v <= 1; and h(t u) >= t h(u) for t in [0,1].
    Also u log(1/u) <= h(u) <= u log(1/u) + u (because 0 <= -(1-u) log(1-u) <= u).
A2. x -> I(x,x) is increasing on [0,1]: derivative 2x[1 + tau(1-x)(1-2x)] and (1-x)(1-2x) >= -1/8, tau < 1/2.
    I(r,r) = r^2 + (1/2 - r^2) = 1/2. Hence I(x,x) <= 1/2 iff x <= r, phi(r) = 1, phi' (r) = 0
    (d/dx h(I(x,x)) = h'(1/2) * ... = 0), phi increasing on [0,r] and decreasing on [r,1].
A3. (small x) For 0 < x <= r: x^2 <= I(x,x) <= 1/2, so L phi(x)^2 = h(I(x,x)) >= h(x^2) >= 2 x^2 ell(x):
    phi(x) >= x sqrt(2 ell(x)/L).  Upper: if (1+tau)x^2 <= 1/2, L phi^2 <= h((1+tau)x^2)
    <= (1+tau)x^2 (2 ell(x) + 1 - log(1+tau)) =: x^2 L Phiu(ell(x))^2.
A4. (x near 1) Put nu(a) = 1 - I(1-a,1-a) = 2a - a^2 - tau a^2 (1-a)^2. For 0 <= a <= 1/4:
    a <= a(2 - (1+tau)a) <= nu(a) <= 2a <= 1/2, hence
    L phi(1-a)^2 = h(nu) >= h(a);  L phi(1-a)^2 >= a(2-(1+tau)a) log(1/(2a));  L phi(1-a)^2 <= h(2a).
A5. D(x) := h(rx) - r h(x) - x h_r has D'' = r(1-r)/((1-x)(1-rx)) > 0, D(0) = D(1) = 0, D'(0+) = (1-r)log(1-r).
    R = phi + D/(qL) (uses phi(r) = 1). Hence R <= phi and R(x) >= phi(x) - rho0 x.
    Near 1: h(r(1-a)) >= (1-a) h(r) (A1), so -D(1-a) <= r h(a) and R(1-a) >= phi(1-a) - (r/(qL)) h(a).
A6. B is convex in x (d^2/dx^2 B = (y/x)[1/(1-x) - 1/(1-xy)] >= 0), B(0,y) = 0, d/dx B(0+,y) = (1-y)log(1-y):
    B(x,y) >= -x omega(y).  By A1, h(xy) >= x h(y): B >= -y h(x) (and symmetric).
A7. (paper's supplement, re-proved) -B = I(X;Y|Z) <= H(X|Z) = P h(a/P) <= (a+b) h(a/(a+b)) <= 2L sqrt(ab),
    a = 1-x, b = 1-y, P = a+b-ab; last step from h(t) <= 2L sqrt(t(1-t)) (N(t) = 2t(1-t)h' - (1-2t)h is concave
    on (0,1/2) with N(0+) = N(1/2) = 0).
A8. b(x) <= x + h(x)/|p| (c phi >= 0). b(1-a) <= 1 + h(a)/|p|.
A9. Lemma 2.1 of the paper (J >= 0), re-derived line by line (NOTES.md).

## Contact identities (analytic) used with the Hessian checks H1-H7

At z* = (r,r): F = 0 (identity h(r^2) + qL = 2 r h_r), dF = 0 (r h'(r^2) - r p - h_r = r log(1+r) + (1-r) log(1-r) = 0
and phi'(r) = 0); J = 0 (h(1/2) = L phi(r)^2) and dJ = 0; R(r) = R'(r) = 0 (R = F(r,.)/(qL)); b(r) = 0 (c = r - h_r/p),
b'(r) = 1 - h'(r)/p - c phi'(r) = 0; v(r) = 0. Hence each of G1..G4 (and R, S1b, b in one variable) vanishes with zero
gradient at the contact. If its interval Hessian on a convex box K containing the r-enclosure is positive definite
(diagonal entries > 0 and lower(H11) lower(H22) > max|H12|^2), Taylor's formula with integral remainder along the
segment from z* gives G >= 0 on K. For G4: K lies in [0, 0.719]^2 where (d - 7h/5) < 0 (checked), so G4 = F on K.

## Endpoints and contact boxes (added after the referee report, 16 Nov)

All one-variable statements are ">= 0" on the closed interval; at x = 0 and x = 1 both sides are evaluated directly
(R(0) = R(1) = 0, S1b(0) = S1b(1) = 0, b(0) = 0, b(1) = 1, phi(0) = phi(1) = 0, h(0) = h(1) = 0), and every strict
"> 0" below refers to the open edge intervals. Contact boxes, the same in both programs: H1-H3 (R, S1b, b) on
K = [0.66, 0.72]; H4-H7 (G1-G4) on K = [0.675, 0.705]^2. The branch-and-bound runs skip exactly the boxes lying
inside K (parameter `contact` in progA_checks.py, `c0,c1` / `cx0,cx1` in progB.c) and certify everything else.
Both K contain the r-enclosure [0.6909077387254070, 0.6909077387254082] and lie in [0, 0.719].

## One-variable kernel inequalities (paper Props S1, S3)

K1 R >= 0. Edges: x <= 3/20: R >= x (sqrt(2 ell/L) - rho0) > 0 [SC: sqrt(2 log(20/3)/L) > rho0].
   x = 1-a, a <= 1/100: R >= sqrt(h(a)/L) - (r/qL) h(a) > 0 iff h(a) < (qL)^2/(r^2 L) [SC: h(1/100) < that].
   Middle [3/20, 99/100]: branch and bound; contact interval by R'' > 0 (H1).
K2 S1b(x) = (x - r + h_r/|p|) h(x) - (h_r^2/|p|) phi(x) + R(x) >= 0. With R >= phi - rho0 x:
   x <= 1/200: S1b/x >= (h_r/|p| - r) ell + (1 - h_r^2/|p|) sqrt(2 ell/L) - rho0, increasing in ell [SC].
   0 < a <= 1/200: S1b >= (1 - h_r^2/|p|) sqrt(h(a)/L) - kp(a) h(a) with kp(a) = r/qL - (1 - r + h_r/|p| - a) > 0.
   Divided by sqrt(h(a)) the bound is (1 - h_r^2/|p|)/sqrt(L) - kp(a) sqrt(h(a)), decreasing in a (kp increases,
   h increases on [0,1/2]); so it suffices that it is > 0 at a = 1/200, coefficient included [SC].
   Middle [1/200, 199/200]: B&B; contact by S1b'' > 0 (H2).
K3 b >= 0: x <= 1/500: b/x >= 1 + ell/|p| - c Phiu(ell), increasing for ell >= ell(1/500) [SC: value and
   derivative 1/|p| - c(1+tau)/(L Phiu) > 0 there]. a <= 1/100: b(1-a) >= 1 - a - c sqrt(h(2a)/L) > 0 [SC].
   Middle [1/500, 99/100]: B&B, contact by b'' > 0 (H3).
   1 - b >= 0: x <= 1/100: 1 - b >= 1 - x - h(x)/|p| >= 1 - 1/100 - h(1/100)/|p| > 0 [SC];
   a <= 1/4: 1 - b(1-a) = a - h(a)/|p| + c phi(1-a) >= a + sqrt(h(a)) (c/sqrt(L) - sqrt(h(a))/|p|) > 0 since
   sqrt(h(a)) <= sqrt(L) < c|p|/sqrt(L) [SC]. Middle [1/100, 3/4]: B&B.
K4 h <= (4/5)phi + 27/400: x or 1-x <= 1/100: h(x) <= h(1/100) < 27/400 [SC]. Middle [1/100, 99/100]: B&B.
   h <= (6/5)phi + 21/10000: x or 1-x <= 1/5000: h <= h(1/5000) < 21/10000 [SC]. Middle: B&B.
K5 h(t^2) >= (809/500) t h(t): Boppana's inequality h(x^2) >= phi_g x h(x), phi_g = (1+sqrt5)/2 > 809/500
   (Boppana, arXiv 2301.09664; also arXiv 2601.19327). Checked numerically by both programs as well (K5).

## Helper inequalities for the edges of the two-variable inequalities

W1 omega(y) <= (2/5) phi(y) on [0,1]. y <= 1/10: omega <= y and phi >= y sqrt(2 ell/L), (2/5) sqrt(2 log 10/L) > 1 [SC].
   1-y = b <= 1/25: omega = b log(1/b) <= h(b), phi >= sqrt(h(b)/L), and (2/5) sqrt(h/L) >= h iff h <= 4/(25L) [SC].
   Middle [1/10, 24/25]: B&B.
W2 phi(y) >= (32/1000) y on [0, 1 - 4e-5]. y <= 69/100 (< r): phi/y >= sqrt(2 ell/L) >= sqrt(2 log(100/69)/L) > 0.032
   [SC]. Middle [69/100, 1 - 4e-5]: B&B.
W3 v(x)^2 <= (9/100) phi(x) on [0,1]. v lies between -h_r phi and h, so v^2 <= max(h, h_r phi)^2.
   x <= 1/1000: h^2/phi <= x (ell+1)^2 / sqrt(2 ell/L) (increasing in x there) and (h_r phi)^2/phi = h_r^2 phi <=
   h_r^2 phi(1/1000) [SC both <= 9/100]. a <= 1/400: h(a)^2/phi(1-a) <= sqrt(L) h(a)^{3/2} and
   h_r^2 phi(1-a) <= h_r^2 sqrt(h(2a)/L) [SC]. Middle [1/1000, 399/400]: B&B.

## Two-variable inequalities G1 = S2a, G3 = F-quadratic, G4 = F-positive-part (paper Props S2, S3)

Reduction. R <= phi gives G1 >= B + (qL - 1/10) phi(x)phi(y); W3 gives G3 >= B + (qL - 99/2500) phi phi;
[d - 7h/5]_+ <= sqrt(L) phi gives G4 >= B + kappa4 phi phi. As kappa4 < qL - 1/10 < qL - 99/2500 [SC], all three
are >= Gamma(x,y) := B + kappa4 phi(x) phi(y). Gamma >= 0 on the following edge sets (and their mirror images):
X0  0 < x <= delta0 = 1e-10, all y: Gamma >= x phi(y)[kappa4 sqrt(2 ell(delta0)/L) - 2/5] >= 0 (A6, A3, W1) [SC].
X1  x = 1-a, a <= delta1 = 2e-7, y <= 1 - beta: Gamma >= kappa4 sqrt(h(a)/L) (0.032 y) - y h(a) >= 0 iff
    h(a) <= (0.032 kappa4)^2/L; h increasing [SC at a = delta1] (A6 second form, A4, W2).
K11 x, y >= 1 - beta, beta = 4e-5: Gamma >= sqrt(ab)[kappa4 (2-(1+tau)beta) log(1/(2beta))/L - 2L] >= 0 (A7, A4) [SC].
x = 0 or y = 0: every G vanishes or is F(0,y) = 0 >= RHS = 0. Remaining compact set
[1e-10, 1-2e-7]^2 minus [1-4e-5, 1]^2: branch and bound (G1, G3, G4 checks), contact boxes by H4, H6, H7.

## G2 = S2b (J inequality)

G2 = J - (1/5)(phi(x)b(y) + phi(y)b(x)) + (11/50)(R(x) + R(y)). Edge sets (and mirrors):
Z1 (pairing) If (11/50)R(x) >= phi(x)/5 and (11/50)R(y) >= phi(y)/5 (i.e. R >= (10/11)phi), then, as 0 <= b <= 1,
   G2 >= J + [(11/50)R(x) - phi(x)/5] + [(11/50)R(y) - phi(y)/5] >= J >= 0 (no division by phi; phi(1) = 0 is fine). By A5 and A4, R(1-a) >= (10/11)phi(1-a)
   when h(a) <= (qL/(11 r))^2 / L; at a <= z1 = 5e-5 [SC]. So G2 >= 0 on [1 - z1, 1]^2.
Z2 x = 1-a, a <= delta1 = 2e-7, y in [y_a, 1 - z1], y_a = 1e-6. Using h(I) >= h(y) - h(ay) (subadditivity,
   y - ay <= I <= y), phi(1-a) <= Phi1 := sqrt(h(2 delta1)/L) (A4), b(1-a) <= 1 + h(delta1)/|p|, R(1-a) >= 0, b(y) >= 0:
   G2 >= Z2(y) := h(y) - h(delta1) - L Phi1 phi(y) - (1/5)Phi1 b(y) - (1/5) phi(y)(1 + h(delta1)/|p|) + (11/50)R(y).
   Check Z2 > 0 on [y_a, 1 - z1] (one-variable B&B, check Z2).
Z3 x = 1-a, a <= delta1, 0 < y <= y_a: J >= h(I) - L phi phi with (1-delta1) y <= I <= y <= 1/2, so
   h(I) >= (1-delta1) h(y); b(y) <= y + h(y)/|p|; b(x) <= 1 + h(delta1)/|p|; R >= 0. Dividing by phi(y):
   G2/phi(y) >= (1 - delta1 - Phi1/(5|p|)) h(y)/phi(y) - L Phi1 - (1/5)(1 + h(delta1)/|p|) - (1/5)Phi1 y/phi(y),
   with h/phi >= ell/sqrt((1+tau)(2 ell + 1 - log(1+tau))/L) (increasing in ell; A1, A3) and y/phi <= sqrt(L/(2 ell))
   [SC at ell = log(1/y_a)]. (y = 0: G2(x,0) = (11/50) R(x) >= 0.)
Z4 0 < x <= delta0' = 1e-10, all y. (a) y <= y1 = 1e-3: J >= 0, R(x)/phi(x), R(y)/phi(y) >= 1 - rho0 sqrt(L/(2 ell(y1)))
   (A3, A5), b(x), b(y) <= y1 + h(y1)/|p| (A8): G2 >= phi(x)[(11/50) rho1 - beta1/5] + phi(y)[same] >= 0 [SC].
   (b) y in [y1, 1]: with K_lo = 1 + tau(1 - delta0')(1-y), I >= x y K_lo and I <= 3x/2 <= 1/e,
   h(I) >= I log(1/I) >= x y K_lo (ell_x + log(1/(y K_lo))); b(x) <= x(1 + (ell_x + 1)/|p|); R(x) >= phi(x) - rho0 x;
   phi(x) between x sqrt(2 ell_x/L) and x Phiu(ell_x). So G2 >= (11/50) R(y) + x Lambda(ell_x, y),
   Lambda(l,y) = y K_lo (l + log(1/(y K_lo))) - (1/5)phi(y)(1 + (l+1)/|p|) - (11/50) rho0 + min(P sqrt(2l/L), P Phiu(l)),
   P = P(y) = 11/50 - L phi(y) - b(y)/5. For each y-interval either
   (beta) Lambda(ell0, y) > 0 and A(y) - max(-P,0)(1+tau)/(L Phiu(ell0)) >= 0 with A = y K_lo - phi(y)/(5|p|)
          (then Lambda is nondecreasing in l >= ell0), or
   (alpha) (11/50) R(y) >= delta0' (|A| ell0 + C0 + |P| Phiu(ell0)), C0 = |y K_lo log(y K_lo)| + (1/5)phi(y)(1 + 1/|p|)
          + (11/50) rho0 (e^-l times each term decreases for l >= ell0 >= 2).
   Check Z4b covers [y1, 1]. (x = 0: G2 = (11/50)R(y) >= 0.)
Remaining compact set [1e-10, 1-2e-7]^2 minus [1 - 5e-5, 1]^2: branch and bound (check G2), contact by H5.

## Covariance inequality (paper Prop S4)

For X in [0,1] with mean a it suffices that some slope lam has psi_a(x) := (x-a+kappa)h(x) - C R(x)
<= kappa h(a) + lam (x - a) for all x (take expectations). Cov(X,h(X)) = E[(X-a)h(X)].
C1 a <= 1/12: x -> (x - a + kappa) h(x) is concave on (0,1) iff k(x) := 2x(1-x)h'(x) - x <= kappa - a; for x >= 1/2,
   k <= -x < 0; for x <= 3/100, k <= 2x ell(x) <= (6/100) log(100/3) [SC]; [3/100, 1/2]: B&B of kappa - 1/12 - k.
   Jensen then gives E[(X-a+kappa)h(X)] <= kappa h(a), and -C E R <= 0.
C2 1/12 <= a <= 2/3: lam = kappa h(a)/a - delta, delta = 8/10000. Gcov(a,x) = kappa h(a) x/a + delta(a-x)
   - (x + kappa - a) h(x) + C R(x) >= 0. x <= d = 2^-16: Gcov >= delta(1/12 - d) - (kappa - 1/12 + d) h(d) [SC].
   x >= 1 - 1/64: Gcov >= kappa (3/2) h(2/3)(63/64) - delta (11/12) - (1 + kappa) h(1/64) [SC] (h(a)/a decreasing).
   [1/12, 2/3] x [d, 63/64]: B&B (check C2).
C3 2/3 <= a <= 9/10: pointwise (x-a)(h(x)-h(a)) <= kappa KL(x||a) suffices (Cov = E[(X-a)(h(X)-h(a))],
   h(a) - E h(X) = E KL(X||a)); for x >= 1-a the left side is <= 0. B&B on {x <= 1-a} (check C3).
C4 a >= 9/10, x <= 1-a: left <= a h(a) <= h(a) and KL(x||a) >= KL(1-a||a) = (2a-1) log(a/(1-a)), increasing in a,
   while h decreases: need h(9/10) <= kappa (4/5) log 9 [SC].

## Psi (paper Prop S5) and the final scalar (paper Prop S6)

S5 Psi >= 0 on the stated constrained box (3-variable B&B, boxes violating a domain constraint on the whole box
   are discarded). S6 M(u) > 0 on u = 2 eps v^2, v in [0,1] (one-variable B&B), plus the sign conditions
   D0 > 0, k_h + m > 0, 0 < Btot < A0 and all constants' signs used in the text (check SC0).
