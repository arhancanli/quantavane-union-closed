# Proposition P: a product-measure barrier for the method (19 November 2026)

Notation as in `../r13/SELF_CONTAINED.md` (Sections 1-4): h is the binary entropy, m = 1 - c*, q and I(s, t) as there,
and "SU" is Sadhu-Ursin (Zenodo 10.5281/zenodo.23322347). For s in (0, 1) put

    g(delta_s) = h(s^2) + q h(I(s,s)) - 2 s h(s),

and let eps_P be the root of g(delta_{m - eps}) = 2 eps h(m - eps):

    eps_P = 0.0012044432051...,   c* + eps_P = 0.3840897031719...

**Proposition P.** For every eps in [1.205e-3, m - 1/2], put s = m - eps. The product measure Bern(s)^n on absence
vectors (each coordinate's law is the point mass at s, in every order and in every coordinate fibre) has abundance
exactly 1 - s = c* + eps, and it satisfies:

1. SU's combined budget (4.1), g(delta_s) + 2 (s - m) h(s) <= 0, per coordinate;
2. the U-budget sum U <= H (h(s^2) - h(s) < 0 per coordinate; about -2.29e-3 at eps_P);
3. SU's fibre inequality (4.2), (s - m) h(s) < 0 (about -8.0e-4), in every fibre (fibres of a product are products);
4. every per-law inequality used in the proof of `../r13/SELF_CONTAINED.md` ((I1)-(I9), Lemma 6.1), since these hold
   for every law; it has zero chain-rule slack and no forced atom.

It violates the separate V-budget, h(I(s,s)) - h(s) <= 0 (about +9.77e-3 per coordinate), and the uniformity identity of
a uniform measure on a union-closed family.

**Consequence (scope).** An argument that uses only (i) sum U <= H and sum (U + qV) <= (1 + q) H on the measure and its
coordinate fibres in every order, (ii) SU's fibre inequality, (iii) per-law inequalities, and (iv) peeling and
subset-profile identities, all as hypotheses about an arbitrary measure, cannot prove c_UC > c* + 1.205e-3 = 0.38409:
the product measure satisfies every such hypothesis and has abundance c* + eps. This covers SU's argument and ours
(`../r13/`), with any parameters, kernels or cuts. It does not cover arguments that use uniformity (H = log|F|, the
per-coordinate identity) beyond forced atoms, nor forcing or slack lemmas whose conclusions fail on products. To go
past 0.38409, an argument must use sum V <= H with more weight than the combined budget gives it, or another input that
product measures violate.

## Certificate

- `d2_certP_A.py` (python-flint, Arb at 200 bits) and `d2_certP_B.py` (mpmath interval arithmetic at 50 digits), written
  separately: at eps = 1205/10^6 (and 121/10^5) both certify (C) < 0, (U) < 0, (F) < 0 and (V) > 0: PASS. Red controls:
  at eps = 12044/10^7 and 119/10^5, just below eps_P, both report FAIL, as they must. Log: `d2_certP.log`.
- Internal referee (fresh context, own code `referee/p_barrier.py`, Arb at 256 bits): reproduced every number; the root
  of C is unique on its bracket (dC/ds in 1.29293 +- 8e-6); C < 0 for every eps in [1.205e-3, m - 1/2] (135 boxes); the
  U-budget holds on the same range (37 boxes); also certified at eps = 12045/10^7. It wrote a direct proof that the
  argument of `../r13/` fails at every eps >= eps_P (`referee/p_meta.py`). Verdict: VALID WITH FIXES; the fixes (the
  scope statement above) are applied.
- `d1_product_barrier.py`: the floating-point computation (mpmath, 40 digits) that located eps_P; it proves nothing.

This is a statement about a method, not a bound on c_UC.
