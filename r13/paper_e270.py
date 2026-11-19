# Update papers/union-closed-0383155/main.tex from c* + 2.3e-4 to c* + 2.7e-4 (run only after the referee of the 2.7e-4 changes).
import sys
p = sys.argv[1]; s = open(p).read(); R = []
def rep(a, b):
    assert s.count(a) == 1, ("MISSING", a[:90]); R.append((a, b))
rep(r" \cuc\ \ge\ \cs+2.3\cdot10^{-4}=0.38311525996679\ldots,", r" \cuc\ \ge\ \cs+2.7\cdot10^{-4}=0.38315525996679\ldots,")
rep(r"single coordinate, all losses are placed in one uniformly random coordinate order and paid from one shared gap budget, and the" + "\n" +
    r"fibre terms of low-entropy coordinates are bounded by a Cauchy--Schwarz inequality instead of by zero.",
    r"single coordinate, all losses are placed in one uniformly random coordinate order and paid from one shared gap budget, the" + "\n" +
    r"fibre terms of low-entropy coordinates are bounded by a Cauchy--Schwarz inequality instead of by zero, and one new two-variable" + "\n" +
    r"kernel inequality (certified here) lets two of their local estimates share the same quantity.")
rep(r"$\cs+23\cdot10^{-5}=0.38311525996679779\ldots$ of its members. Hence $\cuc\ge\cs+2.3\cdot10^{-4}$.",
    r"$\cs+27\cdot10^{-5}=0.38315525996679779\ldots$ of its members. Hence $\cuc\ge\cs+2.7\cdot10^{-4}$.")
rep(r"We also re-tune the free parameter $a_0$ of their subset average. Lemmas~\ref{lem:J}--\ref{lem:beta} are theirs; we give their",
    r"Finally, their two kernel inequalities (I2) and (I3) bound the same quantity $\mathrm{EF}$ from below by the two losses separately;" + "\n" +
    r"we certify a joint form (I9), which is not a convex combination of them. We also re-tune the free parameter $a_0$ of their" + "\n" +
    r"subset average. Lemmas~\ref{lem:J}--\ref{lem:beta} are theirs; we give their")
rep(r"\item $\Cov(Y,h(Y))\le\kappa\,(h(a)-\E h(Y))+C\,\E R(Y)$." + "\n" + r"\end{enumerate}",
    r"\item $\Cov(Y,h(Y))\le\kappa\,(h(a)-\E h(Y))+C\,\E R(Y)$;" + "\n" +
    r"\item $F(x,y)\ge\frac3{20}[v(x)^2\phi(y)+v(y)^2\phi(x)]+\frac{51}{2000}[\phi(x)R(y)+\phi(y)R(x)]$." + "\n" + r"\end{enumerate}")
rep(r"lemmas of Appendix~\ref{app:edges}; (I6) follows from Lemma~\ref{lem:prodred} and the inequality of Appendix~\ref{app:boppana}.",
    r"lemmas of Appendix~\ref{app:edges}; (I6) follows from Lemma~\ref{lem:prodred} and the inequality of Appendix~\ref{app:boppana}." + "\n" +
    r"(I9) is new; it is certified in Appendix~\ref{app:I9}. It is not a convex combination of (I2) and (I3): at $a=\frac3{20}$ the chord" + "\n" +
    r"between them allows only $\frac1{20}(1-\frac{3/20}{11/50})=0.0159$ in place of $\frac{51}{2000}$.")
rep(r"$\ell_G-\theta\beta\le\Lambda_G$, where $\mathrm{EF}_{\rm lo}=\max(2k_Et\vb^2,2\kappa_Ft\Rb)$,",
    r"$\ell_G-\theta\beta\le\Lambda_G$, where $\mathrm{EF}_{\rm lo}=\max(2k_Et\vb^2,\ 2\kappa_Ft\Rb,\ 2t(\frac3{20}\vb^2+\frac{51}{2000}\Rb))$,")
rep(r"$U=\mathrm{EF}-qLt^2+2se$ and $g=\mathrm{EF}+J_p$; use Lemma~\ref{lem:local}(b),(d) and Lemma~\ref{lem:signed}; the",
    r"$U=\mathrm{EF}-qLt^2+2se$ and $g=\mathrm{EF}+J_p$; use Lemma~\ref{lem:local}(b),(d), (I9) integrated against $\mu\otimes\mu$ with" + "\n" +
    r"Jensen, and Lemma~\ref{lem:signed}; the")
rep(r"We take $\varepsilon=23/10^5$, $a_0=21/100$, $\mu_U=14/5$, $\theta=42$, $\Psi=743/50000$, $R_{\max}=431/2000$, and the eight" + "\n" +
    r"coefficient vectors listed in the ancillary file \texttt{params\_e230.json}.",
    r"We take $\varepsilon=27/10^5$, $a_0=21/100$, $\mu_U=14/5$, $\theta=35$, $\Psi=1537/10^5$, $R_{\max}=431/2000$, and the eight" + "\n" +
    r"coefficient vectors listed in the ancillary file \texttt{params\_e270.json}.")
rep(r"\item $\Gamma_0-\Psi-2\varepsilon\theta>0$ (its value is $5.4259\cdot10^{-4}$).",
    r"\item $\Gamma_0-\Psi-2\varepsilon\theta>0$ (its value is $4.5259\cdot10^{-4}$).")
rep(r"g_{\rm lo}[\frac{25}2\mu_Uq-\theta]$ with $P_m=\mu_U(2(m-\varepsilon)-1)-m-\Gamma_0+\Psi=0.0175772$; both brackets are negative for" + "\n" +
    r"$e\le e_1$ ($\frac{25}2\mu_Uq=8.198\le\theta$ and $P_m/(\frac{98}{25}\mu_Uq)=0.006837\ge e_1$).",
    r"g_{\rm lo}[\frac{25}2\mu_Uq-\theta]$ with $P_m=\mu_U(2(m-\varepsilon)-1)-m-\Gamma_0+\Psi=0.0178632$; both brackets are negative for" + "\n" +
    r"$e\le e_1$ ($\frac{25}2\mu_Uq=8.198\le\theta$ and $P_m/(\frac{98}{25}\mu_Uq)=0.006948\ge e_1$).")
rep(r"throughout. Box counts (A/B): (N1) $757/757$, (N2) $38{,}814/38{,}814$, (N3) $730{,}237/3{,}562{,}679$," + "\n" +
    r"(N4) $2{,}094{,}473/2{,}092{,}447$; (N5) margin $5.4259\cdot10^{-4}$ in both. A fresh clone of the repository reproduces every count.",
    r"throughout. Box counts (A/B): (N1) $757/757$, (N2) $38{,}814/38{,}814$, (N3) $2{,}928{,}205/15{,}068{,}425$," + "\n" +
    r"(N4) $2{,}981{,}701/2{,}979{,}457$; (N5) margin $4.5259\cdot10^{-4}$ in both; (I9) $71{,}299/69{,}973$. A fresh clone of the" + "\n" +
    r"repository reproduces every count.")
rep(r"$\varepsilon=2.4\cdot10^{-4}$ fails (N5); a smaller $R_{\max}$ fails (N1); a raised minorant fails (N2).",
    r"$\varepsilon=2.8\cdot10^{-4}$ fails (N5); dropping (I9) fails (N3); a smaller $R_{\max}$ fails (N1); a raised minorant fails (N2);" + "\n" +
    r"a larger coefficient $\frac{27}{1000}$ in (I9) is rejected.")
rep(r"The same proof with $(\varepsilon,a_0,\theta,\Psi)=(2\cdot10^{-4},6/25,47,163/10^4)$ and $(2.25\cdot10^{-4},21/100,44,71/5000)$" + "\n" +
    r"is also certified. Numerically, the relaxation of Lemma~\ref{lem:G} stops near $\varepsilon\approx2.4\cdot10^{-4}$; an ideal" + "\n" +
    r"version that evaluated the per-law quantities exactly would reach about $3.7\cdot10^{-4}$. These two numbers are floating-point" + "\n" +
    r"estimates, not theorems.",
    r"Without (I9) the same proof is certified with $(\varepsilon,a_0,\theta,\Psi)=(2\cdot10^{-4},6/25,47,163/10^4)$," + "\n" +
    r"$(2.25\cdot10^{-4},21/100,44,71/5000)$ and $(2.3\cdot10^{-4},21/100,42,743/50000)$. Numerically, more joint kernels of the" + "\n" +
    r"form (I9) would move the relaxation of Lemma~\ref{lem:G} to about $\varepsilon\approx2.8\cdot10^{-4}$, and an ideal version that" + "\n" +
    r"evaluated the per-law quantities exactly would reach about $3.7\cdot10^{-4}$. These numbers are floating-point estimates, not" + "\n" +
    r"theorems.")
rep(r"\section{The inequality $h(t^2)\ge\frac{809}{500}t\,h(t)$}\label{app:boppana}",
    r"\subsection{The joint input (I9)}\label{app:I9}" + "\n" +
    r"Let $G_9$ be the difference of the two sides of (I9). By $R\le\phi$ (A5) and $v^2\le\frac9{100}\phi$ (W3),"
    r" $G_9\ge B+\kappa_K\phi(x)\phi(y)$ with" + "\n" +
    r"$\kappa_K=qL-\frac{51}{1000}-\frac9{50}\cdot\frac3{20}=0.08436>\kappa_4$, so $G_9\ge\Gamma\ge0$ on the edge sets of the inputs (I2), (I3)," + "\n" +
    r"(I7) above; $G_9=0$ on $x=0$ and $y=0$; at the contact $G_9$ vanishes with zero gradient and its interval Hessian is positive" + "\n" +
    r"definite on $[0.675,0.705]^2$; the rest of $[0,1]^2$ is covered by bisection in both programs." + "\n\n" +
    r"\section{The inequality $h(t^2)\ge\frac{809}{500}t\,h(t)$}\label{app:boppana}")
for a, b in R: s = s.replace(a, b)
open(p, 'w').write(s); print('applied', len(R))
