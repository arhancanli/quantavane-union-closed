# Round 15: update papers/union-closed-0383155/main.tex (v3, 2.9e-4) to v4 at c* + 3.0e-4 (Jensen-free kernels, Lemma lem:jensen,
# N7). Run ONLY after Program B (fresh clone) passes and the short fresh referee of the round-15 delta passes; keep v3 as
# main_v3_e290.tex first. Usage: paper_e300.py <B N3> <B N4> <B N6> <B N7> [path to main.tex]   (counts in TeX form, e.g. 216{,}000{,}000)
import sys
B3, B4, B6, B7 = sys.argv[1:5]
P = sys.argv[5] if len(sys.argv) > 5 else 'papers/union-closed-0383155/main.tex'
s = open(P).read()
def rep(old, new):
    global s
    assert s.count(old) == 1, old[:80]
    s = s.replace(old, new)
rep(r" \cuc\ \ge\ \cs+2.9\cdot10^{-4}=0.38317525996679\ldots,", r" \cuc\ \ge\ \cs+3\cdot10^{-4}=0.38318525996679\ldots,")
rep(r"$\cs+29\cdot10^{-5}=0.38317525996679779\ldots$ of its members. Hence $\cuc\ge\cs+2.9\cdot10^{-4}$.",
    r"$\cs+3\cdot10^{-4}=0.38318525996679779\ldots$ of its members. Hence $\cuc\ge\cs+3\cdot10^{-4}$.")
rep(r"kernel inequality (certified here) lets two of their local estimates share the same quantity, and five certified linear",
    r"kernel inequality (certified here) lets two of their local estimates share the same quantity, and certified linear")
rep(r"subset average. The reduction of Lemma~\ref{lem:G} treats four moments of a law as independent; Lemma~\ref{lem:cuts} adds",
    r"subset average. The reduction of Lemma~\ref{lem:G} treats four moments of a law as independent and used Jensen's inequality for"
    "\n" r"$\E v(X)^2$; Lemmas~\ref{lem:cuts} and~\ref{lem:jensen} add")
rep(r"the six" "\n" r"inequalities (N1)--(N6) of Section~\ref{sec:finite}.", r"the seven" "\n" r"inequalities (N1)--(N7) of Section~\ref{sec:finite}.")
rep(r"""certified by (N6) below; one of them removes the corner $e=z$, $\Rb\approx R_{\max}$, and four are tangent to the curve at $x=0$ and
near the contact law.
""", r"""certified by (N6) below; one of them removes the corner $e=z$, $\Rb\approx R_{\max}$, and four are tangent to the curve at $x=0$ and
near the contact law.

\begin{lemma}[no Jensen step]\label{lem:jensen}
Let $\alpha\in\mathbb Q^6$ with $\alpha_0+\alpha_1x+\alpha_2h(x)+\alpha_3\phi(x)+\alpha_4R(x)+\alpha_5v(x)^2>0$ for every $x\in[0,1]$ and
$\alpha_5>0$. In Lemma~\ref{lem:G}, put $w=\E v(X)^2$ and $w_{\rm lo}=\max(\vb^2,\max_k-(\alpha^{(k)}_0+\alpha^{(k)}_1s+\alpha^{(k)}_2e+
\alpha^{(k)}_3t+\alpha^{(k)}_4\Rb)/\alpha^{(k)}_5)$ over such vectors $\alpha^{(k)}$. Then $w\ge w_{\rm lo}$, and the conclusion of
Lemma~\ref{lem:G} holds with $\vb^2$ replaced by $w_{\rm lo}$ in $\mathrm{EF}_{\rm lo}$.
\end{lemma}
\begin{proof}
$w\ge\vb^2$ by Jensen, and $w\ge-(\dots)/\alpha^{(k)}_5$ by taking expectations as in Lemma~\ref{lem:cuts}. Integrated against
$\mu\otimes\mu$, (I3) gives $\mathrm{EF}\ge2k_Etw$ and (I9) gives $\mathrm{EF}\ge2t(\frac3{20}w+\frac{51}{2000}\Rb)$; both increase in
$w$ since $t\ge0$.
\end{proof}
The ancillary file lists $23$ such vectors (all with $\alpha_5=1$), found by the same separation method and certified by (N7).
""")
rep(r"""We take $\varepsilon=29/10^5$, $a_0=19/100$, $\mu_U=14/5$, $\theta=31$, $\Psi=1513/10^5$, $R_{\max}=431/2000$, and the eight
coefficient vectors and five cut vectors listed in the ancillary file \texttt{params\_e290.json}. Then $\eta=0.46065485$,
$W=0.28832401$, $B'=0.11326892$ and $\Gamma_0=0.03320555$.""",
r"""We take $\varepsilon=3/10^4$, $a_0=9/50$, $\mu_U=14/5$, $\theta=61/2$, $\Psi=703/50000$, $R_{\max}=431/2000$, and the eight
coefficient vectors, five cut vectors and $23$ vectors of Lemma~\ref{lem:jensen} listed in the ancillary file
\texttt{params\_r15\_e300.json}; $\mathrm{EF}_{\rm lo}$ uses $w_{\rm lo}$. Then $\eta=0.44998695$, $W=0.28103392$, $B'=0.10883666$ and
$\Gamma_0=0.03241833$.""")
rep(r"\item $\Gamma_0-\Psi-2\varepsilon\theta>0$ (its value is $9.555\cdot10^{-5}$);",
    r"\item $\Gamma_0-\Psi-2\varepsilon\theta>0$ (its value is $5.833\cdot10^{-5}$);")
rep(r"on $[0,1]$ for $k=1,\dots,5$.",
    r"on $[0,1]$ for $k=1,\dots,5$;" "\n" r"\item the inequality of Lemma~\ref{lem:jensen} on $[0,1]$ for each of the $23$ vectors.")
rep(r"-m-\Gamma_0+\Psi=0.0190283$", r"-m-\Gamma_0+\Psi=0.0186895$")
rep(r"\mu_Uq)=0.007401\ge e_1$", r"\mu_Uq)=0.007270\ge e_1$")
rep(r"By Lemmas~\ref{lem:G}, \ref{lem:T}, \ref{lem:cuts} and (N1)--(N4), (N6), (P\textsubscript{1})",
    r"By Lemmas~\ref{lem:G}, \ref{lem:T}, \ref{lem:cuts}, \ref{lem:jensen} and (N1)--(N4), (N6), (N7), (P\textsubscript{1})")
rep(r"""(N3) $11{,}092{,}699/59{,}190{,}707$, (N4) $3{,}540{,}481/3{,}538{,}273$, (N6) $232{,}393/232{,}393$; (N5) margin $9.555\cdot10^{-5}$
in both;""", r"""(N3) $43{,}249{,}271/""" + B3 + r"""$, (N4) $3{,}734{,}909/""" + B4 + r"""$, (N6) $232{,}393/""" + B6 + r"""$, (N7) $272{,}367/""" + B7 + r"""$;
(N5) margin $5.833\cdot10^{-5}$ in both;""")
rep(r"""and removing all cuts fail (N3); making one cut negative at $x=0$ ($\alpha_0=-10^{-6}$) fails (N6); $\varepsilon=2.95\cdot10^{-4}$
fails (N5).""", r"""and removing the vectors of Lemma~\ref{lem:jensen} fail (N3); making one vector of Lemma~\ref{lem:jensen} negative at $x=0$
($\alpha_0=-10^{-6}$) fails (N7); $\varepsilon=3.05\cdot10^{-4}$ fails (N5). For the $2.9\cdot10^{-4}$ set: removing the corner cut or all
cuts fails (N3), and making one cut negative at $x=0$ fails (N6).""")
rep(r"""With $\varepsilon=57/200000$ (that is, $\cs+2.85\cdot10^{-4}$) and the other parameters unchanged, the margin of (N5) is
$4.0555\cdot10^{-4}$ (box counts of (N3) $10{,}170{,}911/54{,}270{,}139$).""",
r"""Without Lemma~\ref{lem:jensen}, the proof is certified for $(\varepsilon,a_0,\theta,\Psi)=(29/10^5,19/100,31,1513/10^5)$ ((N5) margin
$9.555\cdot10^{-5}$) and for $\varepsilon=57/200000$ with the same other parameters ((N5) margin $4.0555\cdot10^{-4}$).""")
rep(r"""Numerically, with the cuts the
maximiser of the relaxation is already the moment point of a law, so further gains need stronger kernel inequalities than (I2),
(I3) and (I9);""", r"""Numerically, further gains need
stronger kernel inequalities than (I2), (I3) and (I9);""")
open(P, 'w').write(s); print("main.tex updated to v4 (3.0e-4)")
