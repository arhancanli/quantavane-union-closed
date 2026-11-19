# Round 14: finalise SELF_CONTAINED.md v4 after Program B. Usage: sc_final.py <B285 N3> <B290 N3 or NONE> <B290 N4 or NONE> <A290-fresh N3 or NONE>
import sys
B285, B290, B290_4, A290f = sys.argv[1:5]
top = B290 != 'NONE'
p = 'research/tier2/union-closed/r13/SELF_CONTAINED.md'; s = open(p).read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new)
if top:
    rep("# c_UC >= c* + 2.85e-4: one self-contained proof (v4: v3 plus moment-hull cuts for N3, Lemma 6.4 and N6)",
        "# c_UC >= c* + 2.9e-4: one self-contained proof (v4: v3 plus moment-hull cuts for N3, Lemma 6.4 and N6)")
    rep("c* + 57/200000 = 0.38317025996679779... of its members. Hence c_UC >= c* + 2.85e-4. (v3 of this document, without Lemma 6.4 and",
        "c* + 29/100000 = 0.38317525996679779... of its members. Hence c_UC >= c* + 2.9e-4. (With `params_e285.json` in place of\n"
        "`params_e290.json` the same proof gives c* + 2.85e-4. v3 of this document, without Lemma 6.4 and")
rep("`params.json` and `params_e225.json` gives c* + 2e-4 and c* + 2.25e-4. All four parameter sets are certified by two programs; the\n"
    "2.25e-4, 2.3e-4 and 2.7e-4 sets were checked by short fresh referees; the 2.85e-4 set uses the same proof plus Lemma 6.4.)",
    "`params.json` and `params_e225.json` gives c* + 2e-4 and c* + 2.25e-4. Every parameter set is certified by two programs; the\n"
    "2.25e-4, 2.3e-4 and 2.7e-4 sets were checked by short fresh referees; the 2.85e-4" + (" and 2.9e-4 sets use" if top else " set uses") +
    " the same proof plus Lemma 6.4, and that change\nwas checked by a short fresh referee (19 Nov 13:5x: VALID WITH FIXES, all applied).)")
pf = 'params_e290.json' if top else 'params_e285.json'
epsl = 'eps = 29/100000' if top else 'eps = 57/200000'
rep("""**Parameters** (`params_e270.json`). eps = 27/100000, a0 = 21/100, mu = 14/5, theta = 35, Psi = 1537/100000, kernel (I9) (z = 77/200,
R_max = 431/2000, e1 = 1/200, eight minorants). With beta, b0, W, eta from Sec. 5, put A0 = W/2,
B' = eta (t0 + D) + W (a0 - 1/2) - mu d0, Gamma0 = eta t0 - mu (1 - a*). Values: eta = 0.48121287198, W = 0.30317225585,
B' = 0.12221008247, A0 + B' = 0.27380 > 0, A0 - B' = 0.029376 > 0, Gamma0 = 0.0347225873757.""",
f"""**Parameters** (`{pf}`). {epsl}, a0 = 19/100, mu = 14/5, theta = 31, Psi = 1513/100000, kernel (I9), five cuts (Lemma 6.4)
(z = 77/200, R_max = 431/2000, e1 = 1/200, eight minorants). With beta, b0, W, eta from Sec. 5, put A0 = W/2,
B' = eta (t0 + D) + W (a0 - 1/2) - mu d0, Gamma0 = eta t0 - mu (1 - a*). Values: eta = 0.460654846537, W = 0.288324007081,
A0 = 0.144162003540, B' = 0.113268919745, A0 + B' = 0.257431 > 0, A0 - B' = 0.030893 > 0, Gamma0 = 0.0332055490313.
(`params_e285.json` differs only in eps = 57/200000. The 2.7e-4 set of v3, `params_e270.json`: eps = 27/100000, a0 = 21/100,
theta = 35, Psi = 1537/100000, no cuts; eta = 0.48121287198, W = 0.30317225585, B' = 0.12221008247, A0 + B' = 0.27380,
A0 - B' = 0.029376, Gamma0 = 0.0347225873757.)""")
rep("N6: alpha_k . g(x) > 0 on [0,1] for each of the 5 cuts alpha_k of `params_e285.json` (key 'cuts';",
    "N6: alpha_k . g(x) > 0 on [0,1] for each of the 5 cuts alpha_k of `params_e285.json` and `params_e290.json` (key 'cuts', identical;")
if top:
    rep("""(8.198 <= 31) and e <= P_m/((98/25) mu q) (P_m = 0.0190562551214, bound 0.0074123 >= 1/200); N4 is then checked by bisection on
e in [1/200, z]. N5 margin: Gamma0 - Psi - 2 eps theta = 0.000405549031338 (N5 would fail just above eps = 2.9154e-4).""",
"""(8.198 <= 31) and e <= P_m/((98/25) mu q) (P_m = 0.0190282551214, bound 0.0074014 >= 1/200); N4 is then checked by bisection on
e in [1/200, z]. N5 margin: Gamma0 - Psi - 2 eps theta = 0.0000955490313377 (N5 would fail just above eps = 2.9154e-4).
(For the 2.85e-4 set: P_m = 0.0190562551214, bound 0.0074123, N5 margin 0.000405549031338.)""")
rep("`fresh_clone/fresh_certB_e285.log`): RESULTS_B_PENDING.",
    f"`fresh_clone/fresh_certB_e285.log`): N5 0.0004055490313376952, N1 757, N2 38,814, N6 232,393, N4 3,532,243, N3 {B285} boxes, all PASS.")
red = ("N3 with Psi = 0.0148 (below the float sup 0.0148324): POINT controls in both programs, +1.916e-5 at (s_min, 0.61414,\n"
       "8.598e-4, 2.428e-3) with every cut holding there; N3 without the corner cut, and N3 without any cut: point controls in both,\n"
       "+0.019788 at the corner (s_min, z, 0.4962, 0.215), every remaining cut holding; N6 with alpha0 of cut 1 set to -1e-6: FAIL in A\n"
       "(located near x = 0.6538, where that cut crosses 0), B ends at its box cap (8M boxes, no PASS) and B's own point value at x = 0 is\n"
       "-1e-6; N5 with eps = 59/200000 (2.95e-4): FAIL in both (-2.1445e-4). (A file with eps = 2.9e-4 was first listed as a red control;\n"
       "the referee showed it passes N5, +9.555e-5, which led to the 2.9e-4 set.) Scripts `r14_redpoint.py`, `r14_redpointN6.py`; logs\n"
       "`redpoint285.log`, `certAB_red285_quick.log`.")
if top:
    red += (f"\n**The 2.9e-4 set.** `params_e290.json` (= `params_e285.json` with eps = 29/100000), from a second fresh clone at commit 1e2a13c2:\n"
            f"Program A (`certA_e290.log` in the working tree; `fresh_clone/fresh_certA_e290.log` from the clone, N3 {A290f}): N5\n"
            f"0.0000955490313377, N1 757, N2 38,814, N6 232,393, N4 3,540,481, N3 11,092,699, all PASS; Program B\n"
            f"(`fresh_clone/fresh_certB_e290.log`): N5 9.55490313376952e-05, N1 757, N2 38,814, N6 232,393, N4 {B290_4}, N3 {B290}, all PASS.\n"
            "Red controls at 2.9e-4: Psi = 0.0148 (+3.338e-5), no corner cut and no cuts (+0.019802 at the corner): point controls in both\n"
            "programs (`redpoint290.log`); N5 with eps = 2.95e-4 fails (above); N6 does not depend on eps (same cuts).")
rep("Red controls at 2.85e-4 (each must fail): RED_PENDING.", "Red controls at 2.85e-4 (each must fail): " + red)
rep("short fresh referee of the change: REFEREE_PENDING |", "short fresh referee of the change (19 Nov 13:5x): VALID WITH FIXES (wording, one red control; all applied) |")
rep("fresh clones (2e-4, 2.3e-4 and 2.7e-4 sets)", "fresh clones (2e-4, 2.3e-4, 2.7e-4, 2.85e-4" + (" and 2.9e-4" if top else "") + " sets)")
open(p, 'w').write(s); print("SELF_CONTAINED v4 finalised", "(2.9e-4)" if top else "(2.85e-4)")
