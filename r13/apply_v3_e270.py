# Apply the 2.7e-4 update to SELF_CONTAINED.md (v3). Run ONLY after both programs pass N1-N5 at params_e270.json and the
# kernel K2 certificate (jk_kernels/) passes in both programs with red controls.
p = 'SELF_CONTAINED.md'; s = open(p).read(); R = []
def rep(a, b):
    assert s.count(a) == 1, ("MISSING", a[:80]); R.append((a, b))
rep("# c_UC >= c* + 2.3e-4: one self-contained proof (v2.1, after the full referee and the short referee of the added parts)",
    "# c_UC >= c* + 2.7e-4: one self-contained proof (v3: v2.1 plus one new certified kernel input (I9))")
rep("**Theorem.** Every finite union-closed family with a nonempty member has an element that belongs to more than a fraction\nc* + 23/10^5 = 0.38311525996679779... of its members. Hence c_UC >= c* + 2.3e-4.",
    "**Theorem.** Every finite union-closed family with a nonempty member has an element that belongs to more than a fraction\nc* + 27/10^5 = 0.38315525996679779... of its members. Hence c_UC >= c* + 2.7e-4. (v2.1 of this document, with the same proof\nwithout (I9) and the parameter file `params_e230.json`, gives c* + 2.3e-4; that version was refereed in full.)")
rep("The scalar conditions marked [SC] in Appendix A are check group SC (and the edge closures, App. A.8).",
    "  (I9) F(x,y) >= (3/20)[v(x)^2 phi(y) + v(y)^2 phi(x)] + (51/2000)[phi(x) R(y) + phi(y) R(x)]   [ours, round 13]\n"
    "       a joint form of (I2)-(I3) (it is NOT a convex combination of them: (51/2000) exceeds the chord value (1/20)(1 - (3/20)/(11/50)) =\n"
    "       0.0159). Certificate: App. A.9 (edge reduction to Gamma via R <= phi and W3 with kappa_K = qL - 2(51/2000) - (9/50)(3/20) =\n"
    "       0.08436 >= kappa4, so lane 1's edge sets X0, X1, K11 apply; contact Hessian; r-ball check; compact bisection; both programs).\n"
    "The scalar conditions marked [SC] in Appendix A are check group SC (and the edge closures, App. A.8).")
rep("and l_G - theta beta <= Lambda_G, where EF_lo = max(2 k_E t v-bar^2, 2 kappa_F t R-bar),",
    "and l_G - theta beta <= Lambda_G, where EF_lo = max(2 k_E t v-bar^2, 2 kappa_F t R-bar, 2 t ((3/20) v-bar^2 + (51/2000) R-bar))\n(the third entry is (I9) integrated against mu x mu, then Jensen E v^2 >= v-bar^2),")
rep("**Parameters** (`params_e230.json`). eps = 23/100000, a0 = 21/100, mu = 14/5, theta = 42, Psi = 743/50000 (z = 77/200,",
    "**Parameters** (`params_e270.json`). eps = 27/100000, a0 = 21/100, mu = 14/5, theta = 35, Psi = 1537/100000, kernel (I9) (z = 77/200,")
open(p, 'w').write(''.join([s]) if False else s)
for a, b in R: s = s.replace(a, b)
open(p, 'w').write(s); print('applied', len(R), '- now update Secs. 7, 8, 9 numbers and add App. A.9 by hand')
