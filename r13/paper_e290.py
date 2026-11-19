# Round 14: update papers/union-closed-0383155/main.tex (v2, 2.7e-4) to v3 at c* + 2.9e-4 (moment-hull cuts, Lemma lem:cuts, N6).
# Run ONLY after the short fresh referee of the change passes. Usage: paper_e285.py <B N3 count> <B N4 count>
import sys
B3, B4 = sys.argv[1], sys.argv[2]
P = 'papers/union-closed-0383155/main.tex'
s = open(P).read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new)
rep(r" \cuc\ \ge\ \cs+2.7\cdot10^{-4}=0.38315525996679\ldots,", r" \cuc\ \ge\ \cs+2.9\cdot10^{-4}=0.38317525996679\ldots,")
rep(r"kernel inequality (certified here) lets two of their local estimates share the same quantity.",
    r"kernel inequality (certified here) lets two of their local estimates share the same quantity, and five certified linear"
    "\n" r"inequalities on the curve $x\mapsto(x,h(x),\phi(x),R(x))$ remove points of the relaxation that are not moments of any law.")
rep(r"$\cs+27\cdot10^{-5}=0.38315525996679779\ldots$ of its members. Hence $\cuc\ge\cs+2.7\cdot10^{-4}$.",
    r"$\cs+29\cdot10^{-5}=0.38317525996679779\ldots$ of its members. Hence $\cuc\ge\cs+2.9\cdot10^{-4}$.")
rep(r"subset average. Lemmas~\ref{lem:J}--\ref{lem:beta} are theirs;",
    r"subset average. The reduction of Lemma~\ref{lem:G} treats four moments of a law as independent; Lemma~\ref{lem:cuts} adds"
    "\n" r"linear inequalities between them that hold for every law, which remove the points where the relaxation was weakest."
    "\n" r"Lemmas~\ref{lem:J}--\ref{lem:beta} are theirs;")
rep(r"the five" "\n" r"inequalities (N1)--(N5) of Section~\ref{sec:finite}.", r"the six" "\n" r"inequalities (N1)--(N6) of Section~\ref{sec:finite}.")
rep(r"""$\E\sqrt{Xh(X)}\ge\sigma_{\rm lo}$ by (N2). The coefficients of $U$ and $g$ are $-\mu_U$ and $-\theta$.
\end{proof}
""", r"""$\E\sqrt{Xh(X)}\ge\sigma_{\rm lo}$ by (N2). The coefficients of $U$ and $g$ are $-\mu_U$ and $-\theta$.
\end{proof}

\begin{lemma}[moment cuts]\label{lem:cuts}
Let $\alpha\in\mathbb Q^5$ with $\alpha_0+\alpha_1x+\alpha_2h(x)+\alpha_3\phi(x)+\alpha_4R(x)>0$ for every $x\in[0,1]$. Then every
random variable $X$ in $[0,1]$ satisfies $\alpha_0+\alpha_1\E X+\alpha_2\E h(X)+\alpha_3\E\phi(X)+\alpha_4\E R(X)\ge0$. In
particular, in Lemma~\ref{lem:G}, where $\E X=s$, $\E h(X)=e$, $\E\phi(X)=t=(s-e/p-\bb)/c$ and $\E R(X)=\Rb$, the point
$(s,e,\bb,\Rb)$ satisfies $\alpha_0+\alpha_1s+\alpha_2e+\alpha_3t+\alpha_4\Rb\ge0$.
\end{lemma}
\begin{proof}
The function $\alpha_0+\alpha_1x+\alpha_2h(x)+\alpha_3\phi(x)+\alpha_4R(x)$ is continuous and positive on $[0,1]$; take expectations.
\end{proof}
Let $D_G^{\rm cut}$ be the set of points of $D_G$ that satisfy the conclusion of Lemma~\ref{lem:cuts} for each of the five vectors
$\alpha^{(1)},\dots,\alpha^{(5)}$ of the ancillary file. They were found by a separation linear program in floating point and are
certified by (N6) below; one of them removes the corner $e=z$, $\Rb\approx R_{\max}$, and four are tangent to the curve at $x=0$ and
near the contact law.
""")
rep(r"""We take $\varepsilon=27/10^5$, $a_0=21/100$, $\mu_U=14/5$, $\theta=35$, $\Psi=1537/10^5$, $R_{\max}=431/2000$, and the eight
coefficient vectors listed in the ancillary file \texttt{params\_e270.json}. Then $\eta=0.48121287$, $W=0.30317226$,
$B'=0.12221008$ and $\Gamma_0=0.03472259$.""",
r"""We take $\varepsilon=29/10^5$, $a_0=19/100$, $\mu_U=14/5$, $\theta=31$, $\Psi=1513/10^5$, $R_{\max}=431/2000$, and the eight
coefficient vectors and five cut vectors listed in the ancillary file \texttt{params\_e290.json}. Then $\eta=0.46065485$,
$W=0.28832401$, $B'=0.11326892$ and $\Gamma_0=0.03320555$.""")
rep(r"\item $\Lambda_G-\Psi e\le0$ on $D_G$;", r"\item $\Lambda_G-\Psi e\le0$ on $D_G^{\rm cut}$;")
rep(r"\item $\Gamma_0-\Psi-2\varepsilon\theta>0$ (its value is $4.5259\cdot10^{-4}$).",
    r"\item $\Gamma_0-\Psi-2\varepsilon\theta>0$ (its value is $9.555\cdot10^{-5}$);" "\n"
    r"\item $\alpha^{(k)}_0+\alpha^{(k)}_1x+\alpha^{(k)}_2h(x)+\alpha^{(k)}_3\phi(x)+\alpha^{(k)}_4R(x)>0$ on $[0,1]$ for $k=1,\dots,5$.")
rep(r"$P_m=\mu_U(2(m-\varepsilon)-1)-m-\Gamma_0+\Psi=0.0178632$", r"$P_m=\mu_U(2(m-\varepsilon)-1)-m-\Gamma_0+\Psi=0.0190283$")
rep(r"$P_m/(\frac{98}{25}\mu_Uq)=0.006948\ge e_1$", r"$P_m/(\frac{98}{25}\mu_Uq)=0.007401\ge e_1$")
rep(r"By Lemmas~\ref{lem:G}, \ref{lem:T} and (N1)--(N4), (P\textsubscript{1})",
    r"By Lemmas~\ref{lem:G}, \ref{lem:T}, \ref{lem:cuts} and (N1)--(N4), (N6), (P\textsubscript{1})")
rep(r"""$D_G\subset\mathbb R^4$ and $D_T\subset\mathbb R^3$, discarding a box only when it violates a defining constraint of the domain
throughout. Box counts (A/B): (N1) $757/757$, (N2) $38{,}814/38{,}814$, (N3) $2{,}928{,}205/15{,}068{,}425$,
(N4) $2{,}981{,}701/2{,}979{,}457$; (N5) margin $4.5259\cdot10^{-4}$ in both; (I9) $71{,}299/69{,}973$. A fresh clone of the
repository reproduces every count.""",
r"""$D_G\subset\mathbb R^4$ and $D_T\subset\mathbb R^3$, discarding a box only when it violates a defining constraint of the domain
(for $D_G^{\rm cut}$, including the five cuts) throughout. Box counts (A/B): (N1) $757/757$, (N2) $38{,}814/38{,}814$,
(N3) $11{,}092{,}699/""" + B3 + r"""$, (N4) $3{,}540{,}481/""" + B4 + r"""$, (N6) $232{,}393/232{,}393$; (N5) margin $9.555\cdot10^{-5}$
in both; (I9) $71{,}299/69{,}973$. A fresh clone of the repository reproduces every count.""")
rep(r"""Red controls, each required to fail: $\Psi$ just below the numerically found maximum of $\Lambda_G/e$ and a smaller $\theta$ fail
(N3);""", r"""Red controls, each required to fail: $\Psi$ just below the numerically found maximum of $\Lambda_G/e$, removing the corner cut,
and removing all cuts fail (N3); making one cut negative at $x=0$ ($\alpha_0=-10^{-6}$) fails (N6); $\varepsilon=2.95\cdot10^{-4}$
fails (N5). For the earlier parameter set without cuts, a smaller $\theta$ fails (N3);""")
rep(r"$\varepsilon=2.8\cdot10^{-4}$ fails (N5); dropping (I9) fails (N3);", r"$\varepsilon=2.8\cdot10^{-4}$ fails (N5) there; dropping (I9) fails (N3);")
rep(r"""Without (I9), the same proof is also certified for the parameter sets""",
    r"""Without the cuts, the same proof is certified for $(\varepsilon,a_0,\theta,\Psi)=(2.7\cdot10^{-4},21/100,35,1537/10^5)$. Without the
cuts and without (I9), it is also certified for the parameter sets""")
rep(r"""Numerically, more joint kernels of the
form (I9) would move the relaxation of Lemma~\ref{lem:G} to about $\varepsilon\approx2.8\cdot10^{-4}$, and an ideal version that
evaluated the per-law quantities exactly would reach about $3.7\cdot10^{-4}$.""",
r"""Numerically, with the cuts the
maximiser of the relaxation is already the moment point of a law, so further gains need stronger kernel inequalities than (I2),
(I3) and (I9); an ideal version that evaluated the per-law quantities exactly would reach about $3.7\cdot10^{-4}$.""")
open(P, 'w').write(s); print("main.tex updated to v3 (2.85e-4)")
