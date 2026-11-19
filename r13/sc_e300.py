# Round 15: SELF_CONTAINED.md v4 (2.9e-4) -> v5 at c* + 3.0e-4 (Lemma 6.5, N7). Run ONLY after Program B (fresh clone at 9ddf4556)
# passes and the short fresh referee of the round-15 delta passes; copy v4 to SELF_CONTAINED_v4_e290.md first.
# Usage: sc_e300.py <B N3> <B N4> <B N6> <B N7> <referee verdict text> [path]
import sys
B3, B4, B6, B7, REF = sys.argv[1:6]
p = sys.argv[6] if len(sys.argv) > 6 else 'research/tier2/union-closed/r13/SELF_CONTAINED.md'
s = open(p).read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old[:80]
    s = s.replace(old, new)
rep("# c_UC >= c* + 2.9e-4: one self-contained proof (v4: v3 plus moment-hull cuts for N3, Lemma 6.4 and N6)",
    "# c_UC >= c* + 3.0e-4: one self-contained proof (v5: v4 plus Jensen-free kernels, Lemma 6.5 and N7)")
rep("(round 14: Lemma 6.4, the finite inequality N6, N3 on D_G cut by N6's cuts, and the 2.85e-4 and 2.9e-4 parameters;",
    "(round 14: Lemma 6.4, the finite inequality N6, N3 on D_G cut by N6's cuts, and the 2.85e-4 and 2.9e-4 parameters; v4 kept as\n"
    "`SELF_CONTAINED_v4_e290.md`); v5 (round 15: Lemma 6.5, N7, w_lo in EF_lo, and the 3.0e-4 parameters;")
rep("c* + 29/100000 = 0.38317525996679779... of its members. Hence c_UC >= c* + 2.9e-4. (With `params_e285.json` in place of\n"
    "`params_e290.json` the same proof gives c* + 2.85e-4.",
    "c* + 3/10000 = 0.38318525996679779... of its members. Hence c_UC >= c* + 3.0e-4 (parameter file `params_r15_e300.json`).\n"
    "(Without Lemma 6.5, the same proof with `params_e290.json` or `params_e285.json` gives c* + 2.9e-4 or c* + 2.85e-4.")
rep("the 2.85e-4 and 2.9e-4 sets use the same proof plus Lemma 6.4, and that change\nwas checked by a short fresh referee (19 Nov 13:5x: VALID WITH FIXES, all applied).)",
    "the 2.85e-4 and 2.9e-4 sets use the same proof plus Lemma 6.4, and that change\nwas checked by a short fresh referee (19 Nov 13:5x: VALID WITH FIXES, all applied); the 3.0e-4 set adds Lemma 6.5, checked by a\nshort fresh referee (" + REF + ").)")
rep("""**Parameters** (`params_e290.json`). eps = 29/100000, a0 = 19/100, mu = 14/5, theta = 31, Psi = 1513/100000, kernel (I9), five cuts (Lemma 6.4)
(z = 77/200, R_max = 431/2000, e1 = 1/200, eight minorants). With beta, b0, W, eta from Sec. 5, put A0 = W/2,
B' = eta (t0 + D) + W (a0 - 1/2) - mu d0, Gamma0 = eta t0 - mu (1 - a*). Values: eta = 0.460654846537, W = 0.288324007081,
A0 = 0.144162003540, B' = 0.113268919745, A0 + B' = 0.257431 > 0, A0 - B' = 0.030893 > 0, Gamma0 = 0.0332055490313.
(`params_e285.json` differs only in eps = 57/200000.""",
"""**Parameters** (`params_r15_e300.json`). eps = 3/10000, a0 = 9/50, mu = 14/5, theta = 61/2, Psi = 703/50000, kernel (I9), five cuts
(Lemma 6.4) and 23 six-component cuts (Lemma 6.5) (z = 77/200, R_max = 431/2000, e1 = 1/200, eight minorants). With beta, b0, W, eta
from Sec. 5, put A0 = W/2, B' = eta (t0 + D) + W (a0 - 1/2) - mu d0, Gamma0 = eta t0 - mu (1 - a*). Values: eta = 0.449986952230,
W = 0.281033923979, A0 = 0.140516961990, B' = 0.108836662284, A0 + B' = 0.24935 > 0, A0 - B' = 0.031680 > 0, Gamma0 = 0.0324183331206.
(The 2.9e-4 set of v4, `params_e290.json`: eps = 29/100000, a0 = 19/100, theta = 31, Psi = 1513/100000, five cuts, no 6-cuts;
eta = 0.460654846537, W = 0.288324007081, B' = 0.113268919745, Gamma0 = 0.0332055490313; `params_e285.json` differs from it only in
eps = 57/200000.""")
rep("contact region it remains a floating-point statement. The cuts were found by the LP in floating point; only N6 certifies them.)\n",
    "contact region it remains a floating-point statement. The cuts were found by the LP in floating point; only N6 certifies them.)\n"
    "**Lemma 6.5 (no Jensen step; ours, round 15).** Let alpha in Q^6 with alpha . g6(x) > 0 for every x in [0,1], where\n"
    "g6(x) = (1, x, h(x), phi(x), R(x), v(x)^2) (finite inequality N7). Then every law of X in [0,1] satisfies\n"
    "alpha_0 + alpha_1 s + alpha_2 e + alpha_3 t + alpha_4 R-bar + alpha_5 E v(X)^2 >= 0. In Lemma 6.2 put w = E v(X)^2 and\n"
    "w_lo = max(v-bar^2, max over the 6-cuts with alpha_5 > 0 of -(alpha_0 + alpha_1 s + alpha_2 e + alpha_3 t + alpha_4 R-bar)/alpha_5).\n"
    "Then w >= w_lo, and Lemma 6.2 holds with v-bar^2 replaced by w_lo in EF_lo.\n"
    "Proof. g6 is continuous on [0,1]; take expectations; w >= v-bar^2 is Jensen. (I3) integrated against mu x mu gives EF >= 2 k_E t w\n"
    "and (I9) gives EF >= 2 t ((3/20) w + (51/2000) R-bar) (independent copies, no Jensen step); both increase in w since t >= 0, so\n"
    "they are >= the same expressions at w_lo. 6-cuts with alpha_5 = 0 are cuts of Lemma 6.4; alpha_5 < 0 is not used. QED.\n")
rep("""and l_G - theta beta <= Lambda_G, where EF_lo = max(2 k_E t v-bar^2, 2 kappa_F t R-bar, 2 t ((3/20) v-bar^2 + (51/2000) R-bar))
(the third entry is (I9) integrated against mu x mu, then Jensen E v^2 >= v-bar^2),""",
"""and l_G - theta beta <= Lambda_G, where EF_lo = max(2 k_E t w_lo, 2 kappa_F t R-bar, 2 t ((3/20) w_lo + (51/2000) R-bar)) with
w_lo = max(v-bar^2, max over the 6-cuts k with alpha_5^k > 0 of -(alpha_0^k + alpha_1^k s + alpha_2^k e + alpha_3^k t + alpha_4^k R-bar)/alpha_5^k)
(Lemma 6.5; with no 6-cuts, w_lo = v-bar^2 as in v4; the first and third entries are (I3) and (I9) integrated against mu x mu, which
give 2 k_E t E v^2 and 2 t ((3/20) E v^2 + (51/2000) R-bar), increasing in E v^2 >= w_lo),""")
rep("Proof. U = EF - qLt^2 + 2se (1.3), g = EF + Jp (3.2 a); EF >= EF_lo (3.2 b),",
    "Proof. U = EF - qLt^2 + 2se (1.3), g = EF + Jp (3.2 a); EF >= EF_lo (3.2 b without its last Jensen step, and Lemma 6.5),")
rep("R(x) = [h(r x) + qL phi(x) - x h_r - r h(x)]/(qL), which is F(r, x)/(qL) since phi(r) = 1).",
    "R(x) = [h(r x) + qL phi(x) - x h_r - r h(x)]/(qL), which is F(r, x)/(qL) since phi(r) = 1).\n"
    "N7: alpha_k . g6(x) > 0 on [0,1] for each of the 23 6-cuts of `params_r15_e300.json` (key 'cuts6'; v(x)^2 = (h(x) - h_r phi(x))^2).\n"
    "In the 3.0e-4 set, N3 uses EF_lo with w_lo (Lemma 6.5) in place of v-bar^2.")
rep("(8.198 <= 31) and e <= P_m/((98/25) mu q) (P_m = 0.0190282551214, bound 0.0074014 >= 1/200); N4 is then checked by bisection on\n"
    "e in [1/200, z]. N5 margin: Gamma0 - Psi - 2 eps theta = 0.0000955490313377 (N5 would fail just above eps = 2.9154e-4).",
    "(8.198 <= 30.5) and e <= P_m/((98/25) mu q) (P_m = 0.0186894710321, bound 0.0072696 >= 1/200); N4 is then checked by bisection on\n"
    "e in [1/200, z]. N5 margin: Gamma0 - Psi - 2 eps theta = 0.0000583331206310 (N5 would fail just above eps = 3.00956e-4).\n"
    "(For the 2.9e-4 set: theta = 31, P_m = 0.0190282551214, bound 0.0074014, N5 margin 0.0000955490313377.)")
rep("Proof of the Theorem: (P_G) from Lemmas 6.2 and 6.4 and N3 (with N1 and N6);",
    "Proof of the Theorem: (P_G) from Lemmas 6.2, 6.4 and 6.5 and N3 (with N1, N6 and N7);")
rep("## 8. Certificates for N1-N6\n",
    "## 8. Certificates for N1-N7\n"
    "**The 3.0e-4 set (v5).** `params_r15_e300.json`. Program A `certA_v2.py` and Program B `certB_v2.py` are `certA_cut.py` /\n"
    "`certB_cut.py` plus N7 and the w_lo rule in N3 (and larger N3 box caps: A 300M, B 1,000M). Program A: N5 0.0000583331206310,\n"
    "N1 757, N2 38,814, N6 232,393, N7 272,367, N4 3,734,909 (`certA_r15_e300_cap30M.log`, whose N3 stopped at the old 30M cap), N3\n"
    "43,249,271 boxes (`certA_r15_e300.log`, cap 300M), all PASS. Program B from a fresh clone at commit 9ddf4556\n"
    f"(`fresh_clone/fresh_certB_r15_e300.log`): N5 5.83331206310347e-05, N1 757, N2 38,814, N6 {B6}, N7 {B7}, N4 {B4}, N3 {B3}, all PASS.\n"
    "Red controls (`redpoint300_final.log`, committed code): a 6-cut with alpha_0 = -1e-6: N7 FAIL in A (located near x = 0.6518), B's own\n"
    "value at x = 0 is -1e-6; eps = 3.05e-4: N5 FAIL in both (-2.4667e-4); Psi = 0.0138: point control +6.344e-5 in both; no 6-cuts\n"
    "(Jensen step back): point control +2.087e-4 in both, at (s_min, 0.60623, 6.635e-4, 1.845e-3).\n")
rep("| Lemma 6.4, N6 (cuts) | ours (round 14) |",
    "| Lemma 6.5, N7 (Jensen-free kernels) | ours (round 15) | two programs (N7 and N3 with w_lo), red controls (Sec. 8), fresh clone; short fresh referee: " + REF + " |\n| Lemma 6.4, N6 (cuts) | ours (round 14) |")
rep("(I9) (A.9) and N1-N6.", "(I9) (A.9) and N1-N7.")
open(p, 'w').write(s); print("SELF_CONTAINED v5 (3.0e-4) written")
