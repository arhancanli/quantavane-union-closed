# c_UC >= c* + 2.9e-4: one self-contained proof (v4: v3 plus moment-hull cuts for N3, Lemma 6.4 and N6)

Quantavane union-closed lane, round 13. Draft written 18 Nov 2026, 20:4x +04; updated 21:2x to the 2.3e-4 parameter set; v2 21:57 (all fixes of
the full fresh referee); v2.1 23:0x (all fixes of the short referee of the added parts); v3 19 Nov 00:1x (the new input (I9) and
the 2.7e-4 parameters; short fresh referee VALID WITH FIXES, applied; v3 kept as `SELF_CONTAINED_v3_e270.md`); v4 19 Nov 13:4x
(round 14: Lemma 6.4, the finite inequality N6, N3 on D_G cut by N6's cuts, and the 2.85e-4 and 2.9e-4 parameters; v2.1 kept as `SELF_CONTAINED_v2_1_e230.md`; Appendices A-B fold in every finite-check reduction, so this file stands alone; v1 is kept as
`SELF_CONTAINED_v1_refereed.md`). Public copy (deposit of 19 Nov 2026; this status sentence is the only change): not yet peer reviewed. Every lemma taken from Sadhu-Ursin (SU: A. Sadhu and T. Ursin, "An improved lower bound for union-closed
families", Zenodo 10.5281/zenodo.23322347, v1, 12 Oct 2026, unrefereed) is credited to them and re-proved here. SU equations are
cited by their LaTeX labels (paper.tex in `../../../union-closed/lit/zen_23322347/UC92_compact/`). The new parts (Lemma 6.1,
Theorem JB, Lemmas 6.2-6.3, Lemma 3.2(d), N1-N5, Appendix B; NOTEBOOK.md Sec. 13) are ours. Status of every ingredient: Section 9.

**Theorem.** Every finite union-closed family with a nonempty member has an element that belongs to more than a fraction
c* + 29/100000 = 0.38317525996679779... of its members. Hence c_UC >= c* + 2.9e-4. (With `params_e285.json` in place of
`params_e290.json` the same proof gives c* + 2.85e-4. v3 of this document, without Lemma 6.4 and
with `params_e270.json`, gives c* + 2.7e-4; v2.1, also without (I9) and with `params_e230.json`, gives c* + 2.3e-4 and was refereed
in full.) (If the family has one member, that member is
nonempty and each of its elements has abundance 1; the rest of the proof assumes |F| >= 2. The same proof with the parameter files
`params.json` and `params_e225.json` gives c* + 2e-4 and c* + 2.25e-4. Every parameter set is certified by two programs; the
2.25e-4, 2.3e-4 and 2.7e-4 sets were checked by short fresh referees; the 2.85e-4 and 2.9e-4 sets use the same proof plus Lemma 6.4, and that change
was checked by a short fresh referee (19 Nov 13:5x: VALID WITH FIXES, all applied).)
Here c* = 1 - m with m from (1.1) below; c* = 0.38288525996679779... is the ceiling of every single-coordinate (mean-weighted
product) entropy inequality (Costa-Sadhu, arXiv 2610.02295; Sadhu, Zenodo 10.5281/zenodo.22225967; see also Jiang, arXiv 2609.08291).

## 1. Notation and constants (SU Section 2)
h(x) = -x log x - (1-x) log(1-x) (nats; h(0) = h(1) = 0), L = log 2. Let r be the unique root in (0,1) of
-(1-r) log(1-r) - r log(1+r) = 0 (the left side is strictly concave, 0 at r = 0 with right derivative 1, and -L at 1), and
  m = 1/2 + (1-r) log(1+r)/(2L),  q = 2m - 1,  c* = 1 - m.                                                         (1.1)
Enclosure: 0.69090773872540775 < r < 0.69090773872540776 (certA/certB certify brackets of width < 1e-40 by the sign change;
lane 1's Program B one of width 1.1e-15, which suffices for its checks).
  tau = (1/2 - r^2)/(r^2 (1-r)^2) in (0, 1/2),  I(x,y) = x y [1 + tau (1-x)(1-y)],  phi(x) = sqrt(h(I(x,x))/L) in [0,1],
  F(x,y) = h(xy) + qL phi(x) phi(y) - y h(x) - x h(y),  J(x,y) = h(I(x,y)) - L phi(x) phi(y),
  R(x) = F(r,x)/(qL),  v(x) = h(x) - h_r phi(x),  b(x) = x - h(x)/p - c phi(x),
with h_r = h(r), p = h'(r) = log((1-r)/r) < 0, c = r - h_r/p > 0. Note I(r,r) = 1/2 exactly, so phi(r) = 1.
  gamma = h_r/r, alpha = 1 - gamma/p, t0 = r - m, b_m = -h_r gamma/(p alpha), D = -h_r/(p alpha) - t0,          (1.2)
  e* = gamma m, a* = m h(r^2)/(r h_r), a_s = m/r, d0 = 2m - a* + 2 q L a_s/(p c),
  k_E = 11/50, kappa_F = 1/20, D_b = 11/(1 - 1/(44 q L)), kappa = 39/100, C = 13/25.
Identities: c = r alpha (since gamma/p = 1 - alpha); b_m = gamma (t0 + D), hence b_m m = (t0 + D) e*. (SU eq:localconstants.)
Numerical values (both programs): m = 0.6171147400332022, q = 0.2342294801, tau = 0.4966, D = 0.290088, t0 = 0.073793,
e* = 0.552327, a* = 0.999719, d0 = -0.012509.

**Local quantities of a law (SU eq:gap).** For a law mu of a random variable X in [0,1] and a CIID coupling K of mu (an
auxiliary Z with law nu and conditional laws mu_z with int mu_z dnu = mu; K = int mu_z x mu_z dnu(z)):
  s = E X, e = E h(X), t = E phi(X), U = E_{mu x mu} h(XY), V = E_K h(I(X,Y)), g = U + q V - 2 s e,
and R-bar = E R(X), b-bar = E b(X), v-bar = E v(X) = e - h_r t, t_Z = E[phi(X) | Z], b_Z = E[b(X) | Z]. V, g and every
quantity built from them (beta below) depend on the coupling K as well as on the law; x_+ = max(x, 0), x_- = max(-x, 0).
EF := E_{mu x mu} F(X,Y) = U + q L t^2 - 2 s e.                                                                   (1.3)
From the definition of b: c t = s - e/p - b-bar (SU eq:ct).                                                       (1.4)

## 2. Certified inputs (SU Prop. 2.2, eq:hregular, eq:hproduct, eq:point2, eq:hsmall, supplement Props S1-S4)
Pointwise, for all x, y in [0,1] (I1-I7), and one statement about all laws (I8), each with its reduction to finite checks:
  (I1) R(x) >= 0;  0 <= b(x) <= 1;  (x - r) h(x) - (h_r/p) v(x) >= -R(x)   [SU eq:point1; supplement S1, and S3 for 0 <= b <= 1]
       finite checks: K1, K2, K3 (+ the 1 - b part), contact H1-H3, r-ball RB; reductions App. A.3 (A1-A5, A8), A.2.
  (I2) F(x,y) >= kappa_F [phi(x) R(y) + phi(y) R(x)]                                          [SU prop:input; S2a]
       checks: G1, contact H4, RB; edges via Gamma (App. A.5: X0, X1, K11 with A3, A4, A6, A7, W1, W2; the reduction to
       Gamma uses R <= phi, A5).
  (I3) F(x,y) >= k_E [v(x)^2 phi(y) + v(y)^2 phi(x)]                                          [SU eq:localF; S3]
       checks: G3, H6, RB; edges via Gamma (App. A.5, with W3).
  (I4) J(x,y) >= (1/5)[phi(x) b(y) + phi(y) b(x)] - (11/50)[R(x) + R(y)]                      [SU eq:localJ; S2b]
       checks: G2, H5, RB; edges App. A.6 (Z1-Z4 with A1, A3-A5, A8; they also use (I1) (R >= 0 from K1, 0 <= b <= 1 from K3)
       and J >= 0, Lemma 3.1).
  (I5) h(x) <= (4/5) phi(x) + 27/400  and  h(x) <= (6/5) phi(x) + 21/10000                    [SU eq:hregular, eq:hsmall]
       checks: K4 (App. A.3).
  (I6) h(xy) >= (809/1000)[y h(x) + x h(y)]                                                    [SU eq:hproduct]
       from h(t^2) >= (809/500) t h(t) on [0,1] by SU Lemma lem:prodred (Sec. 3.4); that inequality is certified on ALL of
       [0,1] in Appendix B (our programs bopA.py, bopB.py; it is also Boppana's inequality, arXiv 2301.09664, with phi > 809/500).
  (I7) F(x,y) >= (4/25) [sqrt(L) phi(x) - (7/5) h(x)]_+ [sqrt(L) phi(y) - (7/5) h(y)]_+       [SU eq:point2; S3]
       checks: G4, H7, RB, K6 (G4 = F on the contact box); edges via Gamma (App. A.5).
  (I8) for every random variable Y in [0,1] with mean a: Cov(Y, h(Y)) <= kappa (h(a) - E h(Y)) + C E R(Y)  [SU eq:covinput; S4]
       reduced to finite checks in App. A.7 (C1-C4; they also use R >= 0 from K1).
The r-ball check RB is a guard against a coding error in the contact identities, not a certificate of its own.
  (I9) F(x,y) >= (3/20)[v(x)^2 phi(y) + v(y)^2 phi(x)] + (51/2000)[phi(x) R(y) + phi(y) R(x)]   [ours, round 13]
       a joint form of (I2)-(I3) (it is NOT a convex combination of them: (51/2000) exceeds the chord value (1/20)(1 - (3/20)/(11/50)) =
       0.0159). Certificate: App. A.9 (edge reduction to Gamma via R <= phi and W3 with kappa_K = qL - 2(51/2000) - (9/50)(3/20) =
       0.08436 >= kappa4, so lane 1's edge sets X0, X1, K11 apply; contact Hessian; r-ball check; compact bisection; both programs).
The scalar conditions marked [SC] in Appendix A are check group SC (and the edge closures, App. A.8).

## 3. SU's local lemmas, re-proved (credit: Sadhu-Ursin, Section 2 and Lemma 5.1)
### 3.1 Lemma (SU Lemma lem:analyticJ): J >= 0 on [0,1]^2.
Proof (SU's, condensed; lane 1 re-derived it line by line, NOTES 16 Nov). Put a(u) = log h(u) on (0,1). For 0 < u < 1/2,
a'(u) > 0 and -u a''(u) >= (1 + u) a'(u): with omega = log((1-u)/u), l = -log(1-u), h = u omega + l, one computes
(h^2/u)(-u a'' - (1+u) a') = omega/(1-u) - u omega^2 + (l/u)(1/(1-u) - (1+u) omega); since 1 <= l/u <= 1/(1-u) this is at least
one of 1/(1-u) - u omega^2 and 1/(1-u)^2 - u omega^2 - u omega/(1-u), both positive because u log^2(1/u) <= 4/e^2,
u omega/(1-u) = omega e^{-omega} <= 1/e and 4/e^2 + 1/e < 1 (0 < omega < log(1/u)). With z = I(x,y):
I_x = y[1 + tau(1-y)(1-2x)], I_y = x[1 + tau(1-x)(1-2y)], I_xy = 1 + tau(1-2x)(1-2y), all positive (tau < 1/2), and
z I_xy - I_x I_y = tau x^2 y^2, so z I_xy/(I_x I_y) <= 1 + z (equivalent to A B C >= tau with A = 1 + tau(1-x)(1-y),
B = 1 + tau(1-y)(1-2x), C = 1 + tau(1-x)(1-2y): if both x, y <= 1/2 all factors are >= 1; if exactly one exceeds 1/2, at most one
of B, C is below 1 and it is >= 1 - tau > tau; if both exceed 1/2, A >= 1 and B, C >= 1 - tau/2, product >= 9/16 > tau).
Hence d^2/dx dy log h(I(x,y)) = a''(z) I_x I_y + a'(z) I_xy <= 0 (for z < 1/2 by the two displays; for z >= 1/2 since a' <= 0,
a'' < 0). Integrating over [x,y]^2 (x < y) and using symmetry: log h(I(x,x)) + log h(I(y,y)) <= 2 log h(I(x,y)), i.e.
h(I(x,y)) >= L phi(x) phi(y). Boundary cases by continuity. QED.

### 3.2 Lemma (SU Lemma lem:local, the parts used). For every law and every CIID coupling:
(a) g = EF + q E_K J + q L Var(t_Z) (SU eq:decomposition), and all three terms are >= 0; in particular g >= 0.
(b) EF >= 2 k_E t E v(X)^2 >= 2 k_E t v-bar^2, EF >= 2 kappa_F t R-bar, and (by (I9)) EF >= 2t((3/20) E v(X)^2 + (51/2000) R-bar)
    >= 2t((3/20) v-bar^2 + (51/2000) R-bar).
(c) E_K J >= (2/5) E(t_Z b_Z) - (11/25) R-bar   (a step inside SU's proof of Lemma lem:local).
(d) Jp := q E_K J + q L Var(t_Z) >= max(0, q b-bar (2t/5 - 1/(25L)) - (11 q/25) R-bar)   (ours).
(e) if t > 0: g >= t^2 b-bar / D_b.
Proof. (a) Both marginals of K equal mu, so E_{mu x mu} F = U + qLt^2 - 2se (independence) and E_K J = V - L E t_Z^2
(conditional independence given Z: E_K[phi(X) phi(Y)] = E t_Z^2); as E t_Z = t, g - EF - q E_K J = qL (E t_Z^2 - t^2) =
qL Var(t_Z). Nonnegativity: EF >= 0 by (I2) and R >= 0, phi >= 0; E_K J >= 0 by 3.1. (b) Integrate (I3), (I2) under
mu x mu; the right sides factor; then E v^2 >= v-bar^2 (Jensen). (c) Integrate (I4) under K; conditional independence gives
E_K[phi(X) b(Y)] = E(t_Z b_Z). (d) Jp >= 0 by (a). With B = b-bar and y = sqrt(Var t_Z): t B = E(t_Z b_Z) + E[(t - t_Z) b_Z]
<= E(t_Z b_Z) + y sqrt(E b_Z^2) <= E(t_Z b_Z) + y sqrt(B) (0 <= b_Z <= 1 by (I1)); so by (c)
Jp >= (2q/5)(t B - y sqrt(B)) - (11q/25) R-bar + qL y^2 >= (2q/5) t B - q B/(25 L) - (11q/25) R-bar, minimising over y >= 0.
(e) (SU's proof.) From (c), E(t_Z b_Z) <= (5/2) E_K J + (11/10) R-bar; with t <= 1, q > 5/22 and (b):
(t/11) E(t_Z b_Z) <= q E_K J + EF; with (a): g >= (t/11) E(t_Z b_Z) + qL y^2 >= (t^2 B)/11 - (t sqrt(B)/11) y + qL y^2
>= t^2 B (1 - 1/(44 qL))/11 = t^2 B/D_b, completing the square in y. QED.

### 3.3 Lemma (SU eq:signed and eq:identity, the part used). For every law: E[(X - m) h(X)] >= -D e + b_m s - b_m b-bar - R-bar.
Proof (SU's). Integrate the third inequality of (I1): E[(X - r) h(X)] >= (h_r/p) v-bar - R-bar, with v-bar = e - h_r t.
Then E[(X-m) h(X)] = t0 e + E[(X-r) h(X)] >= (t0 + h_r/p) e - (h_r^2/p) t - R-bar; substitute t = (s - e/p - b-bar)/c (1.4):
the coefficient of e is t0 + h_r/p + h_r^2/(p^2 c) = t0 + h_r/(p alpha) = -D (using c = r alpha, gamma = h_r/r), that of s is
-h_r^2/(p c) = b_m, that of b-bar is -b_m. QED.

### 3.4 Lemma (SU Lemma lem:prodred). If h(t^2) >= (809/500) t h(t) on [0,1], then (I6) holds.
Proof (SU's). G(u) = e^u h(e^{-u}) is concave on u > 0 (G'' = -[t/(1-t) + log(1-t)]/t with t = e^{-u}, and the bracket is
positive: 0 at t = 0, derivative t/(1-t)^2). For x = e^{-u}, y = e^{-v}: (y h(x) + x h(y))/(xy) = G(u) + G(v) <= 2 G((u+v)/2) =
2 h(sqrt(xy))/sqrt(xy), so y h(x) + x h(y) <= 2 sqrt(xy) h(sqrt(xy)) <= (1000/809) h(xy). Boundary cases by continuity. QED.

## 4. Families, budgets and the fibre inequality (SU Section 3, re-proved)
### 4.1 Setting (SU Section 3).
Let F be a finite union-closed family with |F| >= 2. Replacing every element by one representative of its membership pattern
(and intersecting members with the set of representatives) keeps distinct members distinct, keeps unions, and keeps the largest
abundance; so F is a family of subsets of [n]. A = (A_1..A_n) is the ABSENCE vector of a uniformly random member (A_i = 1 iff
i is absent); its support is closed under coordinatewise product; H = H(A) = log |F| > 0. For an order of [n],
Q_i = P(A_i = 1 | A_<i); s_i = E Q_i = P(A_i = 1) (order independent), e_i = E h(Q_i) = H(A_i | A_<i), and the local
quantities of Section 1 for the law of Q_i, with the coupling of 4.2.

### 4.2 Lemma (SU Lemma lem:protocols). For every coordinate order:
  sum_i U_i <= H, sum_i V_i <= H, and sum_i g_i + 2 sum_i (s_i - m) e_i <= 0.                                       (4.1)
If S is placed first and p_i = P(A_i = 1 | A_S) for i not in S, then sum_{i not in S} E[(p_i - m) h(Q_i)] <= 0.     (4.2)
Proof (SU's). Independent copies A^1, A^2: their coordinatewise product is in the support, so its entropy is <= H. By the
chain rule, H(A^1 A^2) = sum_i H(A^1_i A^2_i | A^1_<i A^2_<i) >= sum_i H(A^1_i A^2_i | A^1_<i, A^2_<i) = sum_i E h(Q^1_i Q^2_i)
= sum_i U_i (conditioning on the two separate prefixes is finer than on their product). Sign coupling: with independent fair
signs sigma_i (independent of everything before coordinate i), draw the two i-th bits independently with parameters
Q^j_i + sigma_i sqrt(tau) Q^j_i (1 - Q^j_i) (j = 1, 2), which lie in [0,1] and average to Q^j_i; by induction each input vector
has the law of A, so the product is again in the support. Given the earlier signs the two prefixes are i.i.d., so the pair
(Q^1_i, Q^2_i) is a CIID coupling of the law of Q_i (latent = earlier signs); averaging the current sign, the product bit is 1
with probability I(Q^1_i, Q^2_i). The chain rule conditioned on both prefixes and all earlier signs gives sum_i V_i <= H. Since
H = sum_i e_i and 1 + q = 2m: sum_i g_i <= (1 + q) H - 2 sum_i s_i e_i = -2 sum_i (s_i - m) e_i.
For (4.2): fix a with P(A_S = a) > 0. The remaining coordinates, given A_S = a, are uniform on a nonempty collection closed under
product (both factors agree on S). (4.1) for it, with its local gaps >= 0 (3.2(a)), gives sum_{i not in S} (p_i(a) - m)
H(A_i | A_S = a, A_{P_i}) <= 0, where P_i are the coordinates outside S before i; this conditional entropy is E[h(Q_i) | A_S = a]
since the full prefix of i is (A_S, A_{P_i}). Multiply by P(A_S = a) and sum. QED.
**Abundance hypothesis.** If every abundance is <= c* + eps, then s_i = 1 - (abundance of i) >= m - eps for every i, and
(4.1) gives sum_i beta_i <= 2 eps H in every order, where beta_i = g_i + 2 (s_i - m + eps) e_i (each beta_i >= 0).

### 4.3 Peeling (SU, after eq:budget). Fix z = 77/200. Starting from B = [n], remove i from B whenever H(A_i | A_{B \ i}) < z,
until no such i remains; G = the survivors, T = the removed ones listed in reverse removal order. Orders "G in any order, then
T in this order" are called compatible. In every compatible order e_i >= z for i in G (conditioning on a subset of G \ i can only
raise H(A_i | A_{G \ i}) >= z) and e_i < z for i in T (its prefix is exactly its conditioning set at removal); the laws of Q_i,
i in T, do not depend on the order of G (the prefix of i contains all of G). (SU eq:peel.)

## 5. The subset profile (SU Lemma lem:beta, re-proved for every a0 in (0,1))
Put d = |G|, H_G = H(A_G). For t in [0,1] let S_t be a random subset of G, each index in it independently with probability 1 - t
(independent of A); g_H(t) = E H(A_{S_t}), delta(t) = g_H(t) - (1 - t) H_G, C_H(t) = E TC(A_{G \ S_t} | A_{S_t}) with
TC(A_{G\S} | A_S) = sum_{i in G\S} H(A_i | A_S) - H(A_{G\S} | A_S). Fix a0 in (0,1), beta = 1 - D/kappa (> 0 since D < kappa),
b0 = a0^beta, density f(t) = (beta b0/(1 - b0)) t^{-beta-1} on [a0, 1], Theta ~ f, W = (kappa - D)/(1 - b0),
eta = E Theta = ((kappa - D)/D)(b0 - a0)/(1 - b0).
**Lemma.** (i) C_H(t) = delta(t) - t delta'(t); (ii) kappa E C_H(Theta) - D E delta(Theta) = W delta(a0);
(iii) delta(a0) <= (a0 - 1/2)(H_G - d e*) + (1/2) sum_{j=1}^d |a_j - e*|, where a_j is the average over all orders of G of the
j-th chain increment H(A_{pi(j)} | A_{pi(1..j-1)}).
Proof (SU's). (i) Differentiating the finite subset expansion, g_H'(t) = -sum_{i in G} E H(A_i | A_{S_{t,-i}}) (S_{t,-i}: the
same selection on G \ i), so -t g_H'(t) = E sum_{i not in S_t} H(A_i | A_{S_t}) (i is outside S_t with probability t, independently
of S_{t,-i}); and E H(A_{G\S_t} | A_{S_t}) = H_G - g_H(t). Subtract. (ii) f integrates to 1 ((b0/(1-b0))(a0^{-beta} - 1) = 1),
E Theta = (beta b0/((1-b0)(1-beta)))(1 - a0^{1-beta}) = eta (beta/(1-beta) = (kappa-D)/D). Since kappa (t f)' = -(kappa - D) f and
delta(1) = 0, integration by parts gives kappa E C_H - D E delta = kappa a0 f(a0) delta(a0) = (kappa beta/(1-b0)) delta(a0)
= W delta(a0). (iii) Let N ~ Bin(d, 1 - a0); given |S_{a0}| = k the set is uniform among k-subsets, so g_H(a0) = E H-bar_N =
sum_j a_j P(N >= j) and (1 - a0) H_G = (1 - a0) sum_j a_j; since sum_j P(N >= j) = E N = d(1 - a0), e* may be subtracted:
delta(a0) = sum_j (a_j - e*)[P(N >= j) - (1 - a0)]. The bracket lies in [a0 - 1, a0], within 1/2 of a0 - 1/2, and
sum_j (a_j - e*) = H_G - d e*. Nothing here uses a0 = 3/8. QED.

## 6. The new parts (ours; NOTEBOOK.md Sec. 13)
**Parameters** (`params_e290.json`). eps = 29/100000, a0 = 19/100, mu = 14/5, theta = 31, Psi = 1513/100000, kernel (I9), five cuts (Lemma 6.4)
(z = 77/200, R_max = 431/2000, e1 = 1/200, eight minorants). With beta, b0, W, eta from Sec. 5, put A0 = W/2,
B' = eta (t0 + D) + W (a0 - 1/2) - mu d0, Gamma0 = eta t0 - mu (1 - a*). Values: eta = 0.460654846537, W = 0.288324007081,
A0 = 0.144162003540, B' = 0.113268919745, A0 + B' = 0.257431 > 0, A0 - B' = 0.030893 > 0, Gamma0 = 0.0332055490313.
(`params_e285.json` differs only in eps = 57/200000. The 2.7e-4 set of v3, `params_e270.json`: eps = 27/100000, a0 = 21/100,
theta = 35, Psi = 1537/100000, no cuts; eta = 0.48121287198, W = 0.30317225585, B' = 0.12221008247, A0 + B' = 0.27380,
A0 - B' = 0.029376, Gamma0 = 0.0347225873757.)
Per law (Sec. 1): x = e - e*, sigma' = E[(X - m) h(X)] + D e - b_m m, beta(X) = g + 2 (s - m + eps) e,
  l_G(X) = (A0 + B') x_+ + (A0 - B') x_- + [C R-bar - sigma']_+ + mu (a* e + d0 x - U),
  l_T(X) = (Gamma0 + m + mu) e - (E sqrt(X h(X)))^2 - mu U.

**Lemma 6.1 (= NOTEBOOK Lemma 13.1).** For Q in [0,1] and any sigma-field S: E[E[Q|S] h(Q)] >= (E sqrt(Q h(Q)))^2.
Proof. E[E[Q|S] h(Q)] = E[E[Q|S] E[h(Q)|S]] >= E[(E[sqrt(Q h(Q)) | S])^2] (conditional Cauchy-Schwarz) >= (E sqrt(Q h(Q)))^2
(Jensen). QED.

**Theorem JB (NOTEBOOK 13.3).** Let F be a finite union-closed family, |F| >= 2, every abundance <= c* + eps. Suppose
  (P_G) l_G(X) <= Psi e + theta beta(X) for every law with s >= m - eps, e >= z and every CIID coupling;
  (P_T) l_T(X) <= Psi e + theta beta(X) for every law with s >= m - eps, e <= z and every CIID coupling;
  (M) Gamma0 - Psi - 2 eps theta > 0.
Then F does not exist.
Proof. (a) In every compatible order pi: sum_i U_i^pi <= H and sum_i beta_i^pi <= 2 eps H ((4.1) and the abundance hypothesis).
(b) Fibre inequality (4.2) for the order (S, G \ S, T), S a subset of G.
(c) Let pi be a uniformly random order of G and Theta ~ f, and given Theta let N ~ Bin(d, 1 - Theta); (Theta, N) is independent of
A and of pi, and pi of A. S = the first N elements of pi. Then S has the law of S_Theta, and (S, G\S, T) with both parts in the order of pi is pi.
(d) For i in G \ S: E[(p_i - m) h(Q_i)] = E[(Q_i - m) h(Q_i)] - E Cov(Q_i, h(Q_i) | A_S); by (I8) for the conditional laws given
A_S (means p_i), E Cov <= kappa (E h(p_i) - e_i^pi) + C R-bar_i^pi; and E[(Q_i - m) h(Q_i)] = -D e_i^pi + b_m m + sigma'_i^pi.
For i in T: Lemma 6.1 with the sigma-field of A_S (E[Q_i | A_S] = p_i since A_S is in the prefix of i):
E[(p_i - m) h(Q_i)] >= (E sqrt(Q_i h(Q_i)))^2 - m e_i. As sum_{G\S} (E h(p_i) - e_i^pi) = TC(A_{G\S} | A_S):
  0 >= sum_{G\S} (-D e_i^pi + b_m m) - kappa TC(A_{G\S}|A_S) - sum_{G\S} [C R-bar_i^pi - sigma'_i^pi]_+
       + sum_T [(E sqrt(Q_i h(Q_i)))^2 - m e_i].
(e) Average over (pi, Theta, N). Since sum_{G\S} e_i^pi = H(A_{G\S}|A_S) = H_G - H(A_S) and E|G \ S| = eta d,
E sum_{G\S}(-D e_i + b_m m) = -D(eta H_G - E delta(Theta)) + b_m m eta d, and with b_m m = (t0 + D) e* and d e* = H_G - Delta H
(Delta H := H_G - d e* = sum_{i in G} x_i^pi for every pi) this is eta t0 H_G - eta (t0 + D) Delta H + D E delta(Theta). With
Sec. 5 (ii): ... - kappa E C_H(Theta) = eta t0 H_G - eta(t0 + D) Delta H - W delta(a0). By Sec. 5 (iii) and the triangle inequality
(a_j - e* is the average over pi of x^pi at position j): delta(a0) <= (a0 - 1/2) Delta H + (1/2) E_pi sum_{i in G} |x_i^pi|.
The bracket terms are >= 0 and depend only on pi, so their average is <= E_pi sum_{i in G} [C R-bar_i^pi - sigma'_i^pi]_+.
(f) Add mu (sum_i U_i^pi - H) <= 0 (same pi), written as mu sum_i (U_i^pi - e_i^pi), with
U_i - e_i = -(1 - a*) e_i + d0 x_i - (a* e_i + d0 x_i - U_i) for i in G. Collect, using (W/2)|x| + B' x = (A0 + B') x_+ +
(A0 - B') x_- and H = H_G + sum_T e_i:  0 >= E_pi [Gamma0 H - sum_{i in G} l_G(X_i^pi) - sum_{i in T} l_T(X_i)].
(g) By (P_G), (P_T) and (a): for each pi, sum_G l_G + sum_T l_T <= Psi H + theta sum_i beta_i^pi <= (Psi + 2 eps theta) H. So
0 >= (Gamma0 - Psi - 2 eps theta) H > 0, a contradiction. QED.

**Lemma 6.2 (= 13.2, G side).** For a law with s >= m - eps, e >= z, and any CIID coupling, with t = (s - e/p - b-bar)/c (1.4),
v-bar = e - h_r t, x = e - e*: (s, e, b-bar, R-bar) is in
  D_G = {s in [m - eps, 1], e in [z, L], b-bar in [0,1], R-bar in [0, R_max], 0 < t <= 1, e <= (4/5) t + 27/400, e <= h(s)}
and l_G - theta beta <= Lambda_G, where EF_lo = max(2 k_E t v-bar^2, 2 kappa_F t R-bar, 2 t ((3/20) v-bar^2 + (51/2000) R-bar))
(the third entry is (I9) integrated against mu x mu, then Jensen E v^2 >= v-bar^2),
Jp_lo = max(0, q b-bar (2t/5 - 1/(25 L)) - (11q/25) R-bar), and
  Lambda_G = (A0 + B') x_+ + (A0 - B') x_- + [b_m b-bar + (1 + C) R-bar + b_m (m - s)]_+ + mu (a* e + d0 x + q L t^2 - 2 s e)
             - (mu + theta) EF_lo - theta Jp_lo - 2 theta (s - m + eps) e.
Proof. U = EF - qLt^2 + 2se (1.3), g = EF + Jp (3.2 a); EF >= EF_lo (3.2 b), Jp >= Jp_lo (3.2 d); C R-bar - sigma' <=
b_m b-bar + (1 + C) R-bar + b_m (m - s) (3.3). The coefficients of EF and Jp in l_G - theta beta are -(mu + theta) and -theta.
Domain: (I1) and N1 (R <= R_max); t = E phi in [0,1]; (I5) integrated; e <= h(s) by Jensen; t > 0 since e >= z > 27/400. QED.
**Lemma 6.3 (= 13.3, T side).** For a law with s >= m - eps, e <= z, and any CIID coupling, (e, t, s) is in
  D_T = {e in [0, z], t in [0,1], s in [m - eps, 1], e <= (4/5) t + 27/400, e <= (6/5) t + 21/10000, e <= h(s)}
and l_T - theta beta <= Lambda_T := (Gamma0 + m + mu) e - sig_lo^2 - mu U0 - theta (g_lo + 2 (s - m + eps) e), where
U0 = max((809/500) s e, 2 s e - qLt^2), g_lo = max(U0 + qLt^2 - 2se, 2 k_E t (e - h_r t)^2, (t^2/D_b)(s - e/p - c t)_+,
(4/25)((sqrt(L) t - 7e/5)_+)^2), sig_lo = max(0, max_k (a_k t + b_k e + c_k s + d_k)) (coefficients in params_e270.json, identical
in params_e230.json, params_e270.json, params_e285.json and params_e290.json; the kernel (I9) is identical in the last three).
Proof. U >= (809/500) s e by (I6) integrated under mu x mu; U >= 2se - qLt^2 since EF >= 0; g >= EF >= U0 + qLt^2 - 2se;
g >= 2 k_E t v-bar^2 (3.2 a, b); g >= t^2 b-bar/D_b (3.2 e; trivial if t = 0) with b-bar = s - e/p - ct (1.4); g >= EF >=
(4/25)(E[sqrt(L) phi - 7h/5]_+)^2 >= (4/25)((sqrt(L) t - 7e/5)_+)^2 ((I7) integrated, then convexity of (.)_+). E sqrt(X h(X))
>= sig_lo by N2 integrated (and sqrt(X h(X)) >= 0). Coefficients: -mu for U, -theta for g. Domain: (I5) and Jensen. QED.
**Lemma 6.4 (moment-hull cuts; ours, round 14).** Let g(x) = (1, x, h(x), phi(x), R(x)) and let alpha in Q^5 satisfy alpha . g(x) > 0
for every x in [0,1] (finite inequality N6). Then for every law of X in [0,1],
alpha_0 + alpha_1 E X + alpha_2 E h(X) + alpha_3 E phi(X) + alpha_4 E R(X) >= 0. In Lemma 6.2, E X = s, E h(X) = e,
E phi(X) = t = (s - e/p - b-bar)/c (1.4) and E R(X) = R-bar, so (s, e, b-bar, R-bar) also lies in
  D_G^cut = D_G intersected with {alpha_k . (1, s, e, t, R-bar) >= 0 for every cut k of the parameter file}.
Proof. g is continuous on [0,1] (h, phi and R are), so alpha . g(X) is a bounded random variable that is > 0; take expectations. QED.
(Why it matters: D_G treats (s, e, t, R-bar) as independent. The N3 maximisers that bind at eps = 2.8e-4 without cuts are the
corner e = z, R-bar ~ R_max, and the contact region s = m - eps, e ~ 0.61, b-bar ~ 0. A floating-point separation LP says neither is a
moment point of any law; for the corner this is now rigorous (cut 4, alpha0 = 17993/62500, certified by N6, excludes it), for the
contact region it remains a floating-point statement. The cuts were found by the LP in floating point; only N6 certifies them.)

## 7. The finite inequalities (ours) and the analytic small-e part
N1: R <= R_max = 431/2000 on [0,1]. N2: sqrt(x h(x)) >= a_k phi(x) + b_k h(x) + c_k x + d_k on [0,1], k = 1..8.
N3: Lambda_G - Psi e <= 0 on D_G^cut (Lemma 6.4). N4: Lambda_T - Psi e <= 0 on D_T. N5: (M).
N6: alpha_k . g(x) > 0 on [0,1] for each of the 5 cuts alpha_k of `params_e285.json` and `params_e290.json` (key 'cuts', identical; R as in N1, i.e.
R(x) = [h(r x) + qL phi(x) - x h_r - r h(x)]/(qL), which is F(r, x)/(qL) since phi(r) = 1).
Small e in N4 (analytic): put P_m = mu(2(m - eps) - 1) - m - Gamma0 + Psi. Using sig_lo^2 >= 0, U0 >= 2se - qLt^2, and
L t^2 <= (98/25) e^2 + (25/2) g_lo (from g_lo >= (4/25)((sqrt(L) t - 7e/5)_+)^2), with s = m - eps + delta:
Lambda_T - Psi e <= e [-P_m - 2(mu + theta) delta + (98/25) mu q e] + g_lo [(25/2) mu q - theta] <= 0 if theta >= (25/2) mu q
(8.198 <= 31) and e <= P_m/((98/25) mu q) (P_m = 0.0190282551214, bound 0.0074014 >= 1/200); N4 is then checked by bisection on
e in [1/200, z]. N5 margin: Gamma0 - Psi - 2 eps theta = 0.0000955490313377 (N5 would fail just above eps = 2.9154e-4).
(For the 2.85e-4 set: P_m = 0.0190562551214, bound 0.0074123, N5 margin 0.000405549031338.)
(For the 2.7e-4 set of v3: theta = 35, P_m = 0.0178632, bound 0.0069482, N5 margin 0.00045258737569.)
Proof of the Theorem: (P_G) from Lemmas 6.2 and 6.4 and N3 (with N1 and N6); (P_T) from Lemma 6.3, N4 and N2; (M) is N5; Theorem JB.

## 8. Certificates for N1-N6
**The 2.85e-4 set (v4).** `params_e285.json`: eps = 57/200000, a0 = 19/100, mu = 14/5, theta = 31, Psi = 1513/100000, z, R_max, e1,
minorants and the kernel (I9) as in `params_e270.json`, plus 5 cuts. Program A `certA_cut.py` and Program B `certB_cut.py` are
`certA_jk.py` / `certB_jk.py` plus (i) N6 (1-D bisection; Program A regroups alpha . g as alpha0 + (alpha1 - alpha4 h_r/(qL)) x +
(alpha2 - alpha4 r/(qL)) h(x) + (alpha3 + alpha4) phi(x) + (alpha4/(qL)) h(r x); Program B first tries the plain sum with R kept as one
enclosure, then the regrouped form) and (ii) in N3, a box is dropped when some cut has a negative upper bound on the whole box (nothing
else differs, apart from Program B's N3 box cap, raised from 60M to 250M). Program A (`certA_e285.log`): N5 margin 0.000405549031338,
N1 757, N2 38,814, N6 232,393, N4 3,534,487, N3 10,170,911 boxes, all PASS. Program B (run from a fresh clone at commit 86fcd4fe;
`fresh_clone/fresh_certB_e285.log`): N5 0.0004055490313376952, N1 757, N2 38,814, N6 232,393, N4 3,532,243, N3 54,270,139 boxes, all PASS.
Red controls at 2.85e-4 (each must fail): N3 with Psi = 0.0148 (below the float sup 0.0148324): POINT controls in both programs, +1.916e-5 at (s_min, 0.61414,
8.598e-4, 2.428e-3) with every cut holding there; N3 without the corner cut, and N3 without any cut: point controls in both,
+0.019788 at the corner (s_min, z, 0.4962, 0.215), every remaining cut holding; N6 with alpha0 of cut 1 set to -1e-6: FAIL in A
(located near x = 0.6538, where the modified cut crosses 0: B's value at that point is +2e-14), B ends at its box cap (8M boxes, no PASS; that run used certB_cut.py from before the N3-cap edit; that
earlier source is not recoverable, so this run is supporting evidence only) and B's own point value at x = 0, from the committed file, is -1e-6; N5 with eps = 59/200000 (2.95e-4): FAIL in both (-2.1445e-4). (A file with eps = 2.9e-4 was first listed as a red control;
the referee showed it passes N5, +9.555e-5, which led to the 2.9e-4 set.) Scripts `r14_redpoint.py`, `r14_redpointN6.py`; logs
`redpoint285.log`, `certAB_red285_quick.log`.
**The 2.9e-4 set.** `params_e290.json` (= `params_e285.json` with eps = 29/100000), from a second fresh clone at commit 1e2a13c2:
Program A (`certA_e290.log` in the working tree; `fresh_clone/fresh_certA_e290.log` from the clone, N3 11,092,699): N5
0.0000955490313377, N1 757, N2 38,814, N6 232,393, N4 3,540,481, N3 11,092,699, all PASS; Program B
(`fresh_clone/fresh_certB_e290.log`): N5 9.55490313376952e-05, N1 757, N2 38,814, N6 232,393, N4 3,538,273, N3 59,190,707, all PASS.
Red controls at 2.9e-4: Psi = 0.0148 (+3.338e-5), no corner cut and no cuts (+0.019802 at the corner): point controls in both
programs (`redpoint290.log`); N5 with eps = 2.95e-4 fails (above); N6 does not depend on eps (same cuts). Hashes of every
source, parameter file and log of both sets: `hashes_e285_e290.txt`.
**The 2.7e-4 set (v3).**
Program A `certA_jk.py` (python-flint/Arb) and Program B `certB_jk.py` (binary64 outward intervals + mpmath.iv), which are
`certA.py` / `certB.py` plus the kernels of `params['kernels']` in EF_lo (nothing else differs), reading `params_e270.json`: logs
`certA_e270.log`, `certB_e270.log`: N1 757 / 757, N2 38,814 / 38,814, N4 2,981,701 / 2,979,457, N3 2,928,205 / 15,068,425 boxes,
N5 margin 0.000452587375686; fresh-clone rerun at commit 4e1e968f with identical counts (Program B rebuilt from source;
`fresh_clone/fresh_*_e270.log`); hashes `hashes_e270_full.txt`. (I9): Appendix A.9.
Red controls at 2.7e-4 (each must fail): N5 with eps = 2.8e-4: FAIL in both (full runs); N4 without minorants: FAIL in both (full
runs, at SU's Lemma 5.2 corner); N3 without (I9): FAIL in A (located at s = m - eps, e = 0.652, R-bar = 0.0099), B point control
+2.0e-3; N3 with theta = 32: FAIL in A (located at the large-R-bar corner e = z, R-bar = R_max), B point control +2.2e-2 (at a
point of B's search domain, s = SMIN); N3 with Psi = 0.0148 (below the float sup 0.014925): A ended at its box cap (30M boxes,
"FAIL: box cap exceeded", i.e. no PASS) and A's own point control is +7.4e-5, B point control +7.6e-5. Scripts `redpointA270.py`,
`redpointB270.py`, `redpointB270_theta32.py`; logs `certA_red270.log`, `certB_red270.log`, `redpoint270_saved.log`.
Earlier sets (no (I9)): 2.3e-4 (`params_e230.json`, `certA.py`/`certB.py`, logs `certA_e230.log`, `certB_e230.log`, fresh clone at
2906d9c4, hashes `hashes_e225_e230.txt`), 2.25e-4 (`params_e225.json`) and 2e-4 (`params.json`, fresh clone at 7a5c7192f, hashes
`hashes.txt`).
Red controls (each must fail), stated exactly, for the earlier sets: for the 2.3e-4 set, N3 (Psi below the float sup; theta = 38), N4 (no minorants) and
N5 (eps = 2.4e-4) fail in Program A as full runs; in Program B, N4 and N5 fail as full runs and the two N3 controls fail as POINT
controls (B's own evaluator is strictly positive at a feasible point, so no B run can accept a box containing it). N1 and N2 do not
depend on the parameter set; their red controls were run once, for the 2e-4 files: R_max = 0.2148 fails N1 (A located; B at its box
cap) and a raised minorant fails N2 (A located; B by a point control). Logs `certA_red*.log`, `certB_red*.log`, `certB_redpoint.log`.

## 9. Status of every ingredient (what our lab has checked, and what not)
| ingredient | where proved / certified | checked by our lab |
|---|---|---|
| Lemmas 3.1, 3.2(a)(b)(e), 3.3, 3.4 (SU Sec. 2, Lemma 5.1) | SU, re-proved above | lane 1 line-by-line reading (16 Nov) and its referee of paper.tex; full fresh referee of this document (18 Nov 22:0x): each re-proof matches SU |
| Lemma 3.2(c) (inside SU's proof), 3.2(d) (ours) | above | full referee 18 Nov: (d) checked by hand |
| Lemma 4.2 (SU Lemma 3.1), peeling 4.3, Sec. 5 (SU Lemma 6.1, every a0) | SU, re-proved above | lane 1 reading + referee; JB referee (a0-independence); full referee 18 Nov |
| (I1) | S1, S3 | K1, K2, K3 (with 1 - b), H1-H3, RB; App. A.2-A.3 |
| (I2) | S2a | G1, H4, RB; App. A.5 (X0, X1, K11; W1, W2; A3, A4, A6, A7) |
| (I3) | S3 | G3, H6, RB; App. A.5 (with W3) |
| (I4) | S2b | G2, H5, RB; App. A.6 (Z1-Z4) |
| (I5) | S3 | K4 |
| (I6) | lem:prodred + h(t^2) >= (809/500) t h(t) | Appendix B: our two programs on ALL of [0,1] (analytic edges + bisection), red controls; the key identities of Boppana's proof checked symbolically |
| (I7) | S3 | G4, H7, K6, RB; App. A.5 |
| (I8) | S4 | C1-C4 (App. A.7) |
| all lane-1 checks above | lane 1's programs, patched copies in `lane1_fixed/` | code-audit referee (18 Nov 21:4x): all certified with fixes; fixes applied in `lane1_fixed/` in both programs (exact radii, edge closures SCX/ZX and closure.py, C3 on [0, 0.34], exact 4/5 in C4, exact C3 radius, r-ball check RB); full re-run of every group with the final sources; sources, binaries and logs hashed in `lane1_fixed/hashes_final.txt` |
| (I9) joint kernel | ours (App. A.9) | two programs (edge reduction, RB, contact Hessian, bisection), red controls (27/1000; a = 1 breaks the edge reduction; contact perturbation fails RB), fresh clone; short fresh referee of the 2.7e-4 change (19 Nov 00:0x): VALID |
| Lemma 6.1, Theorem JB, Lemmas 6.2, 6.3 | ours (NOTEBOOK 13) | JB referee 18 Nov 20:1x: VALID; full referee 18 Nov 22:0x: valid |
| Lemma 6.4, N6 (cuts) | ours (round 14) | two programs (N6 and N3 with the cuts), red controls (Sec. 8), fresh clone; short fresh referee of the change (19 Nov 13:5x): VALID WITH FIXES (wording, one red control; all applied); short referee of the 2.9e-4 delta (16:5x): VALID WITH FIXES (wording; applied) |
| N1-N5 | ours | two programs, red controls (stated exactly in Sec. 8), fresh clones (2e-4, 2.3e-4, 2.7e-4, 2.85e-4 and 2.9e-4 sets); short referee of the 2.7e-4 set (19 Nov 00:0x): VALID WITH FIXES (documentation only, applied); referee read both programs; short referee of the 2.25e-4 and 2.3e-4 sets |
NOT checked by our lab, or checked only partly: (1) The parts new in v2 (Section 2's citations, Appendix A's marked parts, Appendix B,
Sections 8-9) were refereed by a short fresh referee (18 Nov 22:3x-22:5x): mathematically correct, FIX-level issues only, all applied
in v2.1; v2.1 itself has not been re-read by a referee. (2) Trust base: python-flint/Arb, MPFR, mpmath, the
compilers; nothing is formally verified. (3) The region coverage of Appendix A (edge sets + bisection domains + contact boxes) is
argued in the text, not checked by a program. (4) Red controls for N1 and N2 exist only for the 2e-4 files (these checks do not depend
on the parameter set), and several Program B red controls are point controls (Sec. 8). For the Appendix A groups the red controls are
lane 1's 16 Nov runs (coarse, on the unpatched code), the code-audit referee's fine controls (18 Nov, on its patched copies; all
bite) and our contact-only and SCX controls on the final copies. SU's own programs are not used anywhere.
Unused parts of SU: Section 4 (t_min, the E-bound, A2, Prop. 4.1), Lemma 5.2 (Gamma, S5), Section 7 (S6) and the paper's final
margin. So the trust base is S1-S4 (Appendix A), (I6) (Appendix B), (I9) (A.9) and N1-N6.

## Appendix A. Reductions of the certified inputs I1-I5, I7, I8 to finite checks
Source: lane 1's EDGE_LEMMAS.md (16 Nov 2026), folded in verbatim in substance; the checks are run by lane 1's two programs
(patched copies in `lane1_fixed/`: Program A `progA_checks.py` with Arb balls and jets, Program B `progB.c` with binary64 intervals
and MPFR). Marks: [SC] a scalar condition checked by both programs; [added 16 Nov after the referee] and [added 18 Nov] parts that
postdate the 16 Nov referee of EDGE_LEMMAS.md.
Notation: ell(x) = log(1/x); B(x,y) = h(xy) - y h(x) - x h(y), so F = B + qL phi phi; omega(y) = -(1-y) log(1-y) >= 0;
rho0 = -(1-r) log(1-r)/(qL); kappa4 = L(q - 4/25).

### A.1 Analytic facts
A1. h is increasing on [0,1/2], decreasing on [1/2,1], h(1-u) = h(u), h <= L, h concave, h(0) = 0; so h(u+v) <= h(u) + h(v) for
    u, v >= 0, u + v <= 1, and h(t u) >= t h(u) for t in [0,1]. Also u log(1/u) <= h(u) <= u log(1/u) + u.
A2. x -> I(x,x) is increasing on [0,1] (derivative 2x[1 + tau(1-x)(1-2x)] and (1-x)(1-2x) >= -1/8, tau < 1/2); I(r,r) = 1/2,
    so I(x,x) <= 1/2 iff x <= r, phi(r) = 1, phi'(r) = 0, phi increasing on [0,r] and decreasing on [r,1].
A3. (small x) For 0 < x <= r: x^2 <= I(x,x) <= 1/2, so L phi(x)^2 = h(I(x,x)) >= h(x^2) >= 2 x^2 ell(x): phi(x) >= x sqrt(2 ell(x)/L).
    Upper: if (1+tau)x^2 <= 1/2, L phi^2 <= h((1+tau)x^2) <= (1+tau)x^2 (2 ell(x) + 1 - log(1+tau)) =: x^2 L Phiu(ell(x))^2.
A4. (x near 1) nu(a) = 1 - I(1-a,1-a) = 2a - a^2 - tau a^2 (1-a)^2. For 0 <= a <= 1/4: a <= a(2 - (1+tau)a) <= nu(a) <= 2a <= 1/2,
    hence L phi(1-a)^2 = h(nu) >= h(a); L phi(1-a)^2 >= a(2-(1+tau)a) log(1/(2a)); L phi(1-a)^2 <= h(2a).
A5. D(x) := h(rx) - r h(x) - x h_r has D'' = r(1-r)/((1-x)(1-rx)) > 0, D(0) = D(1) = 0, D'(0+) = (1-r) log(1-r). R = phi + D/(qL)
    (phi(r) = 1). Hence R <= phi and R(x) >= phi(x) - rho0 x. Near 1: h(r(1-a)) >= (1-a) h(r), so -D(1-a) <= r h(a) and
    R(1-a) >= phi(1-a) - (r/(qL)) h(a).
A6. B is convex in x (d^2/dx^2 B = (y/x)[1/(1-x) - 1/(1-xy)] >= 0), B(0,y) = 0, d/dx B(0+,y) = (1-y) log(1-y): B >= -x omega(y).
    By A1, h(xy) >= x h(y): B >= -y h(x) (and symmetrically).
A7. -B = I(X;Y|Z) <= H(X|Z) = P h(a/P) <= (a+b) h(a/(a+b)) <= 2L sqrt(ab), a = 1-x, b = 1-y, P = a+b-ab (SU's supplement,
    re-proved); last step from h(t) <= 2L sqrt(t(1-t)) (N(t) = 2t(1-t)h' - (1-2t)h is concave on (0,1/2), N(0+) = N(1/2) = 0).
A8. b(x) <= x + h(x)/|p| (c phi >= 0); b(1-a) <= 1 + h(a)/|p|.
A9. J >= 0 (Lemma 3.1).

### A.2 Contact identities and contact boxes
At (r,r): F = 0 (identity h(r^2) + qL = 2 r h_r), dF = 0 (r h'(r^2) - r p - h_r = r log(1+r) + (1-r) log(1-r) = 0 and phi'(r) = 0);
J = 0 (h(1/2) = L phi(r)^2) and dJ = 0; R(r) = R'(r) = 0 (R = F(r,.)/(qL)); b(r) = 0 (c = r - h_r/p), b'(r) = 1 - h'(r)/p - c phi'(r) = 0;
v(r) = 0. Hence each of G1-G4 (and R, S1b, b in one variable) vanishes with zero gradient at the contact. If its interval Hessian
on a convex box K containing the r-enclosure is positive definite (diagonal entries > 0, lower(H11) lower(H22) > max|H12|^2), Taylor's
formula with integral remainder from the contact gives G >= 0 on K. For G4: K = [0.675, 0.705]^2 lies in [0.62, 0.719]^2, where
sqrt(L) phi - 7h/5 < 0 (check K6), so G4 = F on K.
[added 18 Nov] Check RB (both programs): the value and gradient enclosures of G1-G4 at (r,r) and of R, S1b, b at r contain 0. This
guards the identities above against a wrong formula; a red control that perturbs only the contact (h_r + 1e-9 inside R) now FAILS.
[added 16 Nov after the referee] All one-variable statements are ">= 0" on the closed interval; at x = 0 and x = 1 both sides are
evaluated directly (R(0) = R(1) = 0, S1b(0) = S1b(1) = 0, b(0) = 0, b(1) = 1, phi(0) = phi(1) = 0, h(0) = h(1) = 0); every strict "> 0"
below refers to the open edge intervals. Contact boxes, the same in both programs: H1-H3 (R, S1b, b) on K = [0.66, 0.72]; H4-H7
(G1-G4) on K = [0.675, 0.705]^2. The bisection runs skip exactly the boxes inside K and certify everything else. Both K contain the
r-enclosure [0.6909077387254070, 0.6909077387254082] and lie in [0, 0.719].

### A.3 One-variable inputs (S1, S3)
K1 R >= 0. x <= 3/20: R >= x (sqrt(2 ell/L) - rho0) > 0 [SC: sqrt(2 log(20/3)/L) > rho0]. x = 1-a, a <= 1/100:
   R >= sqrt(h(a)/L) - (r/qL) h(a) > 0 iff h(a) < (qL)^2/(r^2 L) [SC at a = 1/100]. Middle [3/20, 99/100]: bisection; contact by H1.
K2 S1b(x) = (x - r + h_r/|p|) h(x) - (h_r^2/|p|) phi(x) + R(x) >= 0 (this is the third part of (I1)). With R >= phi - rho0 x:
   x <= 1/200: S1b/x >= (h_r/|p| - r) ell + (1 - h_r^2/|p|) sqrt(2 ell/L) - rho0, increasing in ell [SC]. 0 < a <= 1/200:
   S1b >= (1 - h_r^2/|p|) sqrt(h(a)/L) - kp(a) h(a), kp(a) = r/qL - (1 - r + h_r/|p| - a) > 0; divided by sqrt(h(a)) this is
   decreasing in a, so it suffices at a = 1/200 [SC]. Middle [1/200, 199/200]: bisection; contact by H2.
K3 b >= 0: x <= 1/500: b/x >= 1 + ell/|p| - c Phiu(ell), increasing for ell >= ell(1/500) [SC: value and derivative
   1/|p| - c(1+tau)/(L Phiu) > 0 there]. a <= 1/100: b(1-a) >= 1 - a - c sqrt(h(2a)/L) > 0 [SC]. Middle [1/500, 99/100]: bisection,
   contact by H3. 1 - b >= 0: x <= 1/100: 1 - b >= 1 - x - h(x)/|p| >= 1 - 1/100 - h(1/100)/|p| > 0 [SC]; a <= 1/4:
   1 - b(1-a) = a - h(a)/|p| + c phi(1-a) >= a + sqrt(h(a))(c/sqrt(L) - sqrt(h(a))/|p|) > 0 since sqrt(h(a)) <= sqrt(L) < c|p|/sqrt(L)
   [SC]. Middle [1/100, 3/4]: bisection.
K4 h <= (4/5) phi + 27/400: x or 1-x <= 1/100: h(x) <= h(1/100) < 27/400 [SC]; middle [1/100, 99/100]: bisection.
   h <= (6/5) phi + 21/10000: x or 1-x <= 1/5000: h <= h(1/5000) < 21/10000 [SC]; middle: bisection.
K6 7h/5 - sqrt(L) phi > 0 on [0.62, 0.719] (bisection; used only for H7, whose box [0.675, 0.705] lies inside).

### A.4 Helper inequalities
W1 omega(y) <= (2/5) phi(y) on [0,1]. y <= 1/10: omega <= y and phi >= y sqrt(2 ell/L), (2/5) sqrt(2 log 10/L) > 1 [SC].
   1 - y = b <= 1/25: omega = b log(1/b) <= h(b), phi >= sqrt(h(b)/L), and (2/5) sqrt(h/L) >= h iff h <= 4/(25L) [SC]. Middle
   [1/10, 24/25]: bisection.
W2 phi(y) >= (32/1000) y on [0, 1 - 4e-5]. y <= 69/100: phi/y >= sqrt(2 ell/L) >= sqrt(2 log(100/69)/L) > 0.032 [SC]. Middle
   [69/100, 1 - 4e-5]: bisection.
W3 v(x)^2 <= (9/100) phi(x) on [0,1]. v lies between -h_r phi and h, so v^2 <= max(h, h_r phi)^2. x <= 1/1000:
   h^2/phi <= x (ell+1)^2 / sqrt(2 ell/L) (increasing in x there) and h_r^2 phi <= h_r^2 phi(1/1000) [SC both <= 9/100]. a <= 1/400:
   h(a)^2/phi(1-a) <= sqrt(L) h(a)^{3/2} and h_r^2 phi(1-a) <= h_r^2 sqrt(h(2a)/L) [SC]. Middle [1/1000, 399/400]: bisection.

### A.5 Two-variable inputs (I2), (I3), (I7): checks G1, G3, G4
Reduction. R <= phi gives G1 >= B + (qL - 1/10) phi(x)phi(y); W3 gives G3 >= B + (qL - 99/2500) phi phi; [sqrt(L) phi - 7h/5]_+
<= sqrt(L) phi gives G4 >= B + kappa4 phi phi. As kappa4 < qL - 1/10 < qL - 99/2500 [SC], all three are >= Gamma(x,y) := B +
kappa4 phi(x) phi(y), and Gamma >= 0 on these edge sets (and mirrors):
X0  0 < x <= delta0 = 1e-10, all y: Gamma >= x phi(y)[kappa4 sqrt(2 ell(delta0)/L) - 2/5] >= 0 (A6, A3, W1) [SC].
X1  x = 1-a, a <= delta1 = 2e-7, y <= 1 - beta: Gamma >= kappa4 sqrt(h(a)/L)(0.032 y) - y h(a) >= 0 iff h(a) <= (0.032 kappa4)^2/L;
    h is increasing there, so [SC at a = delta1] suffices (A6 second form, A4, W2).
K11 x, y >= 1 - beta, beta = 4e-5: Gamma >= sqrt(ab)[kappa4 (2 - (1+tau) beta) log(1/(2 beta))/L - 2L] >= 0 (A7, A4) [SC].
x = 0 or y = 0: every G vanishes. Remaining compact set [1e-10, 1 - 2e-7]^2 minus [1 - 4e-5, 1]^2: bisection (G1, G3, G4), contact
boxes by H4, H6, H7.

### A.6 The J input (I4): check G2
G2 = J - (1/5)(phi(x) b(y) + phi(y) b(x)) + (11/50)(R(x) + R(y)). Edge sets (and mirrors):
Z1 If R >= (10/11) phi at both points, then, as 0 <= b <= 1, G2 >= J + [(11/50)R(x) - phi(x)/5] + [(11/50)R(y) - phi(y)/5] >= J >= 0.
   By A5 and A4, R(1-a) >= (10/11) phi(1-a) when h(a) <= (qL/(11 r))^2/L; at a <= z1 = 5e-5 [SC]. So G2 >= 0 on [1 - z1, 1]^2.
Z2 x = 1-a, a <= delta1 = 2e-7, y in [y_a, 1 - z1], y_a = 1e-6. Using h(I) >= h(y) - h(ay) (y - ay <= I <= y), phi(1-a) <= Phi1 :=
   sqrt(h(2 delta1)/L) (A4), b(1-a) <= 1 + h(delta1)/|p|, R(1-a) >= 0, b(y) >= 0: G2 >= Z2(y) := h(y) - h(delta1) - L Phi1 phi(y)
   - (1/5) Phi1 b(y) - (1/5) phi(y)(1 + h(delta1)/|p|) + (11/50) R(y); check Z2 > 0 on [y_a, 1 - z1] (bisection).
Z3 x = 1-a, a <= delta1, 0 < y <= y_a: J >= h(I) - L phi phi with (1 - delta1) y <= I <= y <= 1/2, so h(I) >= (1 - delta1) h(y);
   b(y) <= y + h(y)/|p|; b(x) <= 1 + h(delta1)/|p|; R >= 0. Dividing by phi(y): G2/phi(y) >= (1 - delta1 - Phi1/(5|p|)) h(y)/phi(y)
   - L Phi1 - (1/5)(1 + h(delta1)/|p|) - (1/5) Phi1 y/phi(y), with h/phi >= ell/sqrt((1+tau)(2 ell + 1 - log(1+tau))/L) and
   y/phi <= sqrt(L/(2 ell)); the bound increases in ell, so [SC at ell = log(1/y_a)] suffices. (y = 0: G2(x,0) = (11/50) R(x) >= 0.)
Z4 0 < x <= delta0' = 1e-10, all y. (a) y <= y1 = 1e-3: J >= 0, R(x)/phi(x), R(y)/phi(y) >= 1 - rho0 sqrt(L/(2 ell(y1))) (A3, A5),
   b(x), b(y) <= y1 + h(y1)/|p| (A8): G2 >= phi(x)[(11/50) rho1 - beta1/5] + phi(y)[same] >= 0 [SC]. (b) y in [y1, 1]: with
   K_lo = 1 + tau(1 - delta0')(1-y), I >= x y K_lo and I <= 3x/2 <= 1/e, h(I) >= I log(1/I) >= x y K_lo (ell_x + log(1/(y K_lo)));
   b(x) <= x(1 + (ell_x + 1)/|p|); R(x) >= phi(x) - rho0 x; phi(x) between x sqrt(2 ell_x/L) and x Phiu(ell_x). So
   G2 >= (11/50) R(y) + x Lambda(ell_x, y), Lambda(l,y) = y K_lo (l + log(1/(y K_lo))) - (1/5) phi(y)(1 + (l+1)/|p|) - (11/50) rho0
   + min(P sqrt(2l/L), P Phiu(l)), P = P(y) = 11/50 - L phi(y) - b(y)/5. For each y-interval either (beta) Lambda(ell0, y) > 0 and
   A(y) - max(-P,0)(1+tau)/(L Phiu(ell0)) >= 0 with A = y K_lo - phi(y)/(5|p|) (then Lambda is nondecreasing in l >= ell0), or
   (alpha) (11/50) R(y) >= delta0' (|A| ell0 + C0 + |P| Phiu(ell0)), C0 = |y K_lo log(y K_lo)| + (1/5) phi(y)(1 + 1/|p|) + (11/50) rho0
   (e^-l times each term decreases for l >= ell0 >= 2). Check Z4b covers [y1, 1]. (x = 0: G2 = (11/50) R(y) >= 0.)
Remaining compact set [1e-10, 1 - 2e-7]^2 minus [1 - 5e-5, 1]^2: bisection (G2), contact by H5.

### A.7 The covariance input (I8): checks C1-C4
For X in [0,1] with mean a it suffices that some slope lam has psi_a(x) := (x - a + kappa) h(x) - C R(x) <= kappa h(a) + lam (x - a)
for all x (take expectations; Cov(X, h(X)) = E[(X - a) h(X)]).
C1 a <= 1/12: x -> (x - a + kappa) h(x) is concave on (0,1) iff k(x) := 2x(1-x) h'(x) - x <= kappa - a; for x >= 1/2, k <= -x < 0;
   for x <= 3/100, k <= 2x ell(x) <= (6/100) log(100/3) [SC]; [3/100, 1/2]: bisection of kappa - 1/12 - k. Jensen then gives
   E[(X - a + kappa) h(X)] <= kappa h(a), and -C E R <= 0.
C2 1/12 <= a <= 2/3: lam = kappa h(a)/a - delta, delta = 8/10000; Gcov(a,x) = kappa h(a) x/a + delta(a - x) - (x + kappa - a) h(x)
   + C R(x) >= 0. x <= d = 2^-16: Gcov >= delta(1/12 - d) - (kappa - 1/12 + d) h(d) [SC]. x >= 1 - 1/64: Gcov >= kappa (3/2) h(2/3)(63/64)
   - delta (11/12) - (1 + kappa) h(1/64) [SC] (h(a)/a decreasing). [1/12, 2/3] x [d, 63/64]: bisection (C2).
C3 2/3 <= a <= 9/10: pointwise (x - a)(h(x) - h(a)) <= kappa KL(x||a) suffices (Cov = E[(X-a)(h(X) - h(a))], h(a) - E h(X) =
   E KL(X||a) since E X = a; C E R >= 0); for x >= 1 - a the left side is <= 0. Bisection on {x <= 1 - a} [added 18 Nov: run on
   x in [0, 0.34] in BOTH programs (skip test x0 >= 1 - a0 exact), closing a sliver at a = 2/3, x = 1/3; 91 boxes in each].
C4 a >= 9/10, x <= 1 - a: left side <= a h(a) <= h(a), and KL(x||a) >= KL(1-a||a) = (2a - 1) log(a/(1-a)), increasing in a, while h
   decreases: need h(9/10) <= kappa (4/5) log 9 [SC; added 18 Nov: Program B now uses the exact 4/5].

### A.8 Edge closures and exact radii [added 18 Nov, code audit]
The bisection ranges are binary doubles while the edge sets end at exact rationals, which left uncovered slivers of width
4e-27 to 6e-17. Closed in BOTH programs by re-checking every edge [SC] condition with its threshold enlarged by the factor 1 + 1e-9
(so each edge set overlaps its bisection range), each evaluated at the end of the enlarged set where the bound is weakest, with
delta0' and delta1 enlarged likewise in Z2, Z4b, X0, X1: Program A `lane1_fixed/closure.py` (25/25 PASS) and the enlarged
constants in `lane1_fixed/progA_checks.py`; Program B groups SCX (23/23 PASS) and ZX (Z2, Z4b with enlarged delta0', delta1;
PASS) in `lane1_fixed/progB.c` (red control: factor 2.5 instead of 1 + 1e-9 fails 15 of the 23 SCX conditions). Program A's
mean-value radius is now computed exactly from the box endpoints (it was up to one ulp too small); Program B's C3 radius likewise.
With these fixes every group passes; the box counts equal lane 1's 16 Nov counts except C3 (87 -> 91, from the range 0.34).


### A.9 The joint kernel input (I9) [added 18 Nov, round 13]: check GK2
GK2(x,y) = F(x,y) - (3/20)[v(x)^2 phi(y) + v(y)^2 phi(x)] - (51/2000)[phi(x) R(y) + phi(y) R(x)].
Reduction at the edges. By A5 (R <= phi) and W3 (v^2 <= (9/100) phi), the subtracted terms are at most
(2 (51/2000) + 2 (9/100)(3/20)) phi(x) phi(y), so GK2 >= B + kappa_K phi(x) phi(y) with
kappa_K = qL - 51/1000 - 27/1000 = 0.08436 > kappa4 = 0.05145 (interval check in both programs). Since phi phi >= 0,
B + kappa_K phi phi >= Gamma = B + kappa4 phi phi, which is >= 0 on lane 1's edge sets X0, X1, K11 (A.5, with the closures of A.8).
On x = 0 or y = 0, GK2 = 0 (F(0,y) = 0, v(0) = phi(0) = R(0) = 0). Contact: GK2 vanishes with zero gradient at (r,r) (A.2:
F, v, R and their first derivatives at r), r-ball check RB passes in both programs, and the interval Hessian is positive definite on
[0.675, 0.705]^2 (Program A: 24^2 cells, min det 1.75e-4; Program B: 2,096 cells, min det 5.1e-7). Remaining compact set
[1e-10, 1 - 2e-7]^2 minus [1 - 4e-5, 1]^2: bisection, Program A 71,299 boxes (exact mean-value radii), Program B 69,973 boxes.
Red control: coefficient 27/1000 in place of 51/2000 (above the grid maximum 0.026044 at a = 3/20) FAILS in both programs, with
counterexamples near (0.917, 0.917) (A) and (0.906, 0.906) (B). Files: `jk_kernels/progA_checks.py` (group GK2), `jk_kernels/progB.c`
(group GK2, binary `progB_jk`), logs `jk_kernels/logs/`.
Further red controls: a = 1 in place of 3/20 makes kappa_K = -0.0686 < kappa4 and the edge-reduction check FAILS in both programs
(the Hessian and the bisection fail too); the contact perturbation h_r + 1e-9 inside R makes RB FAIL in both programs while the
Hessian still passes (`jk_kernels/red_edge/`, `jk_kernels/red_contact/`). The progA_checks.py banner prints "python-flint 0.7.1";
the runs used python-flint 0.9.0 (cosmetic, kept so that the hashed sources stay identical).

## Appendix B. (I6) on all of [0,1]: h(t^2) >= c t h(t), c = 809/500
G(t) = h(t^2) - c t h(t). (i) t <= 1/100, ell = log(1/t): h(t^2) >= t^2 log(1/t^2) = 2 t^2 ell and t h(t) <= t^2 (ell + 1) (A1), so
G >= t^2 [(2 - c) ell - c], increasing in ell and positive at ell = log 100 ((2 - c) log 100 - c = 0.14117... > 0). (ii) u = 1 - t <=
1e-4, ell = log(1/u): t^2 = 1 - w with w = 2u - u^2 <= 2u, so h(t^2) = h(w) >= w log(1/w) >= (2u - u^2)(ell - log 2), and
t h(t) <= h(u) <= u (ell + 1); so G >= u [(2 - u - c) ell - (2 - u) log 2 - c] >= u [(2 - 1e-4 - c) log(1e4) - 2 log 2 - c] =
u (0.51313...) > 0. (iii) [1/100, 1 - 1e-4]: bisection with monotone enclosures of h (A1): Program A `boppana/bopA.py` (Arb) and
Program B `boppana/bopB.py` (binary64 outward + mpmath.iv) both PASS (7,607 intervals each); the two edge constants are checked
in both. Red controls: c = 1.618035 > phi fails near t = 0.6197 (A; B point control G(1/phi) <= -4.2e-7); c = 7/4 fails; c = 1.61803
< phi passes (the inequality is sharp at t = 1/phi: Boppana's inequality h(x^2) >= phi x h(x), arXiv 2301.09664, holds with equality
there; the key identities of its proof (f''' = p/(x(1-x^2)^2), the double roots at 0 and 1/phi, f''(1/phi) > 0) were checked
symbolically, `boppana/bop_sympy.py` and `bop_sympy.log`). The control c = 1.61803 (passes) was run in Program A only.
