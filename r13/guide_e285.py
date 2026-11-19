# Round 14: update papers/union-closed-0383155/AUTHOR_GUIDE.md to draft v3 (c* + 2.85e-4). Run only after the referee passes.
P = 'papers/union-closed-0383155/AUTHOR_GUIDE.md'
s = open(P).read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new)
rep('(draft v2, c* + 2.7e-4)', '(draft v3, c* + 2.85e-4)')
rep("updated 19 Nov 2026 (v2: constant c* + 2.7e-4; v1 at c* + 2.3e-4 kept as\n`main_v1_e230.tex`)",
    "updated 19 Nov 2026 (v2: constant c* + 2.7e-4, kept as `main_v2_e270.tex`; v3, 19 Nov afternoon: c* + 2.85e-4 with the\nmoment cuts of Lemma 7.3; v1 at c* + 2.3e-4 kept as\n`main_v1_e230.tex`)")
rep("`research/tier2/union-closed/r13/SELF_CONTAINED.md` (v3, with Appendices A-B).",
    "`research/tier2/union-closed/r13/SELF_CONTAINED.md` (v4, with Appendices A-B; Lemma 6.4 there is Lemma 7.3 here).")
rep("c* + 2.7e-4 = 0.38315525996679779... of its members, so c_UC >= c* + 2.7e-4. This is 2.7e-4 above the best bound proved in a",
    "c* + 2.85e-4 = 0.38317025996679779... of its members, so c_UC >= c* + 2.85e-4. This is 2.85e-4 above the best bound proved in a")
rep("and 2.9 times the increment claimed", "and 3.1 times the increment claimed")
rep("the joint kernel input (I9) and its certificate, Lemmas 7.1-7.2\nand the finite inequalities N1-N5.",
    "the joint kernel input (I9) and its certificate, Lemmas 7.1-7.3 (7.3: the moment cuts)\nand the finite inequalities N1-N6.")
rep("* N1-N5: two programs (Arb; binary64 outward + mpmath), red controls, fresh-clone reproduction (at three commits, including the\n  2.7e-4 set at 4e1e968f). (I9): two programs, red controls, fresh clone.",
    "* N1-N6: two programs (Arb; binary64 outward + mpmath), red controls, fresh-clone reproduction (at four commits, including the\n  2.7e-4 set at 4e1e968f and the 2.85e-4 set at 86fcd4fe). (I9): two programs, red controls, fresh clone.")
rep("the new parameters): VALID WITH FIXES (documentation only, applied in v3).",
    "the new parameters): VALID WITH FIXES (documentation only, applied in v3); short referee of the 2.85e-4 change (the cuts,\n  N6, the new parameters): REFEREE_VERDICT.")
rep("7. Lemmas 7.1, 7.2 (45 min):", "6c. Lemma 7.3, the moment cuts (15 min): the one-line proof, and why t = E phi(X) and R-bar = E R(X) in Lemma 7.1.\n7. Lemmas 7.1, 7.2 (45 min):")
rep("`research/tier2/union-closed/r13/certA_jk.py params_e270.json` and `certB_jk.py params_e270.json`\n   (Program B's N3 takes about 25 minutes);",
    "`research/tier2/union-closed/r13/certA_cut.py params_e285.json` and `certB_cut.py params_e285.json`\n   (Program A takes about 15 minutes; Program B's N4 and N3 together take B_TIME);")
rep('"c_UC >= c* + 2.7e-4 = 0.38315525996679779...",', '"c_UC >= c* + 2.85e-4 = 0.38317025996679779...",')
rep("(including the sentence on the new input (I9))", "(including the sentences on the new input (I9) and on the moment cuts)")
open(P, 'w').write(s); print("AUTHOR_GUIDE updated to v3")
