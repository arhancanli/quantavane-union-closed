# Lemma F: forced coordinates stay nearly forced, on average, under random unions (21 November 2026)

**Status: proved below; checked by two internal referees (fresh context: VALID; a second referee of this public text: VALID
WITH FIXES, applied) and by two separately written numerical checks; not yet peer reviewed.** The statement has no finite part, so the checks are sanity checks, not certificates.

## Setting

Let F be a finite union-closed family of subsets of a finite ground set V with |F| >= 2, and let A, B be independent and uniform on
F. Fix an order of V. For a coordinate i and a member X, write X_<i for the pattern of X on the coordinates before i (the set
of those coordinates that X contains). For every pattern x realised by a member, let

    Q_i(x) = P(i not in A | A_<i = x),

the conditional probability that i is absent. Write h for the binary entropy (natural logarithms), d(p || q) for the binary
relative entropy (with 0 log 0 = 0), and put

    P_i = Q_i(A_<i) Q_i(B_<i),     c_i = A_<i u B_<i,     U_i = E h(P_i).

(U_i is the term of the entropy method for independent copies: H(A u B) >= sum_i U_i.)

## The identity

**Proposition.** log|F| - sum_i U_i = sum_i E d( P_i || Q_i(c_i) ), and every term on the right is finite and nonnegative.

*Proof.* Let mu be the uniform measure on F and S = A u B. Since F is union-closed, S lies in F, so -log mu(S) = log|F|.
Every prefix S_<i is a realised pattern, and by the chain rule

    -log mu(S) = sum_i -log mu(S_i | S_<i),

where mu(S_i | S_<i) is the conditional probability under mu (absence of i has probability Q_i(S_<i)). Given (A_<i, B_<i),
the copies' i-th bits are independent, absent with probabilities Q_i(A_<i) and Q_i(B_<i), so S_i is absent with probability
P_i, and

    E[ -log mu(S_i | S_<i) | A_<i, B_<i ] = -P_i log Q_i(c_i) - (1 - P_i) log(1 - Q_i(c_i)) = h(P_i) + d(P_i || Q_i(c_i)).

Taking expectations and summing over i gives log|F| = sum_i U_i + sum_i E d(P_i || Q_i(c_i)). Finiteness: if P_i > 0 then
members X, Y with patterns A_<i, B_<i and i absent exist, and X u Y has pattern c_i with i absent, so Q_i(c_i) > 0; if
P_i < 1 then one copy's pattern admits i present, hence so does c_i, and Q_i(c_i) < 1. QED.

The identity holds for the sum only: for a single coordinate, E d(P_i || Q_i(c_i)) and H(A_i | A_<i) - U_i can differ, and
H(A_i | A_<i) - U_i can be negative: for F = {empty set, {1}, {1, 2}} in the order 1, 2, H(A_i | A_<i) - U_i = 0.2877 and -0.0959,
while E d(P_i || Q_i(c_i)) = 0.1336 and 0.0581; both sums equal log 3 - sum_i U_i = 0.1918.

## Lemma F

**Lemma F.** With the notation above,

    sum_i E[ Q_i(c_i) 1{P_i = 0} ]  <=  sum_i E[ -log(1 - Q_i(c_i)) 1{P_i = 0} ]  <=  log|F| - sum_i U_i.

In particular this holds with the indicator 1{Q_i(A_<i) = 0} (one copy's prefix forces i present).

*Proof.* On {P_i = 0} the term of the Proposition is d(0 || q) = -log(1 - q) >= q with q = Q_i(c_i) < 1 (finite by the
Proposition, since P_i = 0 < 1). All other terms are nonnegative. QED.

**Meaning.** If one copy's prefix forces i present, then the union's prefix must, on average, also (nearly) force i present:
the total of the union-prefix absence probabilities on those events is at most the slack log|F| - sum_i U_i of the entropy
method with independent copies. Pinsker's inequality (d(0 || q) >= 2q^2) with Cauchy-Schwarz gives only
sqrt(n (log|F| - sum_i U_i)/2), n = |V|; the bound here is linear in the slack.

**Near a counterexample.** Notation of `../r13/SELF_CONTAINED.md` (c* = 0.38288525996679779..., m = 1 - c*, q = 2m - 1,
s_i = P(i not in A), e_i = H(A_i | A_<i), so sum_i e_i = H = log|F|). If every element of V has abundance at most c* + eps (so
s_i >= m - eps), then the per-law gaps g_i = U_i + q V_i - 2 s_i e_i >= 0 (Lemma 3.2(a) there, Sadhu and Ursin's Lemma lem:local,
which rests on the certified inputs (I1)-(I2) and Lemma 3.1) and the budget sum_i V_i <= H (Lemma 4.2 there) give
sum_i U_i >= 2 sum_i s_i e_i - q sum_i V_i >= (2(m - eps) - q) H = (1 - 2 eps) H, so the slack is at most 2 eps H. This step uses
the budget sum_i V_i <= H on its own; it does not use sum_i U_i <= H or the third inequality of (4.1). So in a near-counterexample a
forced prefix is "absorbing" under a random union, up to a total of 2 eps H.

Lemma F does not by itself improve any constant. It is the starting point for the open question of whether forced prefixes
of high-entropy coordinates are unforced by an independent random union often enough, which would bound the forced atoms
(the mass P(Q_i(A_<i) = 0) that the law of Q_i puts at 0) that limit the method.

## Checks

- `lemma_f_check_a.py`: the internal referee's own program (numpy, binary64), copied unchanged; exact enumeration of all
  ordered pairs (A, B) on random union-closed families (two seeds, 1,477 families each, n = 3 to 8, 2 to 76 members;
  `check.sh` runs seed 1): identity residual <= 1e-15; no violation of the one-copy form with 1{Q_i(A_<i) = 0} (both
  inequalities; "Lemma F" in the program's output) or of the second inequality of Lemma F with 1{P_i = 0} ("two-sided
  version"). It has no red control. Log `lemma_f_check_a.log`.
- `lemma_f_check_b.py`: written separately; Q_i as exact rationals, logarithms with mpmath at 40 digits; 300 random
  union-closed families (n = 3 to 6, 2 to 17 members, each containing the empty set; all ordered pairs): identity residual <= 2e-39, no violation: PASS. Red controls,
  each of which must fail: (1) the same computation on 151 families that are not union-closed (the union of two members can
  leave the family): FAIL on all (necessarily: in such a family some union leaves F, which makes a term infinite); (2) the union-prefix conditional Q_i(c_i) multiplied by 0.999 on the 300 union-closed
  families: the identity fails on all (residual up to 3.6e-3). Log `lemma_f_check_b.log`.

`sh ../check.sh` runs check B and its red controls (and check A if numpy is installed).
