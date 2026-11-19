# Union-closed families: c_UC ≥ c* + 2.9·10⁻⁴, a barrier for the method, and unions of translates

**Status: verified by two independent programs and an internal referee; not yet peer reviewed.**

Author: Arhan Canli (Quantavane, arhan@quantavane.com). On the use of AI, see
[Authorship and the role of AI](#authorship-and-the-role-of-ai).

## The main result

A finite family of sets is union-closed if the union of any two members is a member. Frankl's union-closed sets
conjecture says that every finite union-closed family with a nonempty member has an element in at least half of its
members. Let c_UC be the largest constant for which "at least a proportion c_UC" is true.

**Theorem.** c_UC ≥ c* + 29/100000 = 0.38317525996679779…, where c* = 0.38288525996679… is the limit of the
single-coordinate entropy method (Costa and Sadhu).

This is 2.9·10⁻⁴ above the bound 0.3828852549667978 of Costa and Sadhu (arXiv:2610.02295). The proof is
computer-assisted: `r13/SELF_CONTAINED.md` (v4) is the complete proof, with every finite inequality reduced to checks
that two programs certify.

## The earlier work this builds on

- **Method.** A. Sadhu and T. Ursin, *An improved lower bound for union-closed families*, Zenodo,
  [doi:10.5281/zenodo.23322347](https://doi.org/10.5281/zenodo.23322347) (2026): the fibre inequalities, peeling, the
  subset average and the local inputs (I1)–(I8). Their lemmas are re-proved in `r13/SELF_CONTAINED.md` with credit.
  They claim c* + 9.2·10⁻⁵ (not refereed).
- **The constant c* and the bound improved.** S. Costa and A. Sadhu, *Beyond Liu's 0.382709 threshold for the
  union-closed sets conjecture*, [arXiv:2610.02295](https://arxiv.org/abs/2610.02295) (2026); A. Sadhu, *The exact
  ceiling of the product relaxation*, Zenodo, doi:10.5281/zenodo.22225967 (2026).
- **(I6)** rests on R. B. Boppana's inequality (arXiv:2301.09664); `r13/boppana/` certifies the form we need directly.
- The entropy method: J. Gilmer (arXiv:2211.09055); R. Alweiss, B. Huang and M. Sellke; Z. Chase and S. Lovett;
  W. Sawin; S. Cambie; L. Yu; J. Liu; Y. Jiang (references in `r13/SELF_CONTAINED.md`).

New here: one shared gap budget with per-law losses in one random order (Theorem 6.2), the Cauchy–Schwarz fibre bound
for removed coordinates (Lemma 6.1), Lemma 3.2(d), the joint kernel input (I9) and its certificate, the moment-hull
cuts (Lemma 6.4) and the finite inequalities N1–N6.

## Checking it: one command

    PYTHON=/path/to/python sh check.sh          # Python with python-flint and mpmath; about 15-30 minutes
    PYTHON=/path/to/python sh check.sh full     # also Program B on N4 and N3, about 3 hours

- **Program A** (`r13/certA_cut.py`, python-flint/Arb) and **Program B** (`r13/certB_cut.py`, binary64 intervals
  widened outward plus mpmath.iv), written separately, certify N1–N6 for `r13/params_e290.json`. Recorded counts
  (fresh clone at the lab's commit 1e2a13c2): N5 margin 9.5549·10⁻⁵; N1 757; N2 38,814; N6 232,393; N4 3,540,481 (A),
  3,538,273 (B); N3 11,092,699 boxes (A), 59,190,707 (B); all PASS (`r13/certA_e290.log`,
  `r13/fresh_clone/fresh_certA_e290.log`, `r13/fresh_clone/fresh_certB_e290.log`; hashes `r13/hashes_e285_e290.txt`).
- **Red controls** (each must fail): eps = 2.95·10⁻⁴ fails N5 in both programs (run by `check.sh`); Ψ below its float
  supremum and the cut set removed are refused as point controls in both programs (`r13/redpoint290.log`); earlier
  parameter sets have their own controls (Section 8 of the proof).
- **The certified inputs** (I1)–(I8), by two programs each, after a code audit: `r13/lane1_fixed/` (Program A
  `progA_checks.py` groups SC, K, W, Z, C, H, G1–G4 and `closure.py`; Program B `progB.c`, C with MPFR), logs in
  `r13/lane1_fixed/logs/`. (I9): `r13/jk_kernels/` (group GK2, both programs, red controls). Boppana: `r13/boppana/`.
  These were run in the lab on 18–19 November 2026; their scripts contain the lab's paths.

## Also here

- **Proposition P** (`barrier/PROPOSITION_P.md`): the product measure Bern(m − ε)ⁿ satisfies every input of the
  Sadhu–Ursin method and of ours for ε ≥ 1.205·10⁻³, so no argument of that kind can go past c* + 1.205·10⁻³ = 0.38409.
  Certified by two programs, red controls, an internal referee (VALID WITH FIXES). It is a statement about the method.
- **Unions of translates** (`translates/note.pdf`): for a finite group G and R ⊆ G, the family of unions of left
  translates of R satisfies avg F_L(R) + avg F_L(R⁻¹) ≥ |G|, with equality exactly for right cosets; so the conjecture
  holds for F_L(R) whenever R⁻¹ = α(gRh) for an automorphism α. This extends Aaronson, Ellis and Leader's theorem for
  abelian groups and is a quantitative form of the duality of Johnson and Vaughan for these families. Internal referee
  VALID WITH FIXES; exhaustive checks in groups of order up to 27 (`translates/`). The general non-abelian case is open.

The few files edited for this public copy (status sentences, local paths) are listed with their hashes in
`PUBLIC_COPY_CHANGES.md`.

## Not yet done

- Outside review: none yet. A paper is in preparation.
- Lean 4: not started. The plan is a Lean check of the finite certificate layer; the analytic reduction stays a written
  proof.

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane, https://quantavane.com), who takes responsibility for the results.

The arguments, the certificates, the checkers and these documents were written and carried out by AI agents (Claude,
by Anthropic) working for Quantavane under Arhan Canli's direction. Programs described as written separately share no
code but were written with the same AI system. No mathematician outside Quantavane has reviewed this repository; the
proofs can be read and checked, and the computations rerun with the programs here or with your own.

## Licence and citation

Documents, parameter files, logs and text are under CC BY 4.0 (`LICENSE`); the code is under the MIT licence
(`LICENSE-CODE`). Please cite as: Arhan Canli, *Union-closed families: c_UC ≥ c* + 2.9·10⁻⁴*,
https://github.com/arhancanli/quantavane-union-closed (2026), together with Sadhu–Ursin and Costa–Sadhu.
