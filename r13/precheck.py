# Round 13 pre-check (FLOAT, not a certificate): lane 1's scalar model of SU's margin M(u) (upside.py), with three knobs:
# fQ3 (Q3 multiplier), Gam (Lemma 5.2 constant), split=True: Q3 and the C_m-type loss draw on ONE gap budget u
# (adversary picks the split u1 + u2 = u) instead of both charging the full u. t_min'' = 0.463 (Theorem V).
import numpy as np
from math import log, sqrt
L=log(2)
def h(x): return 0.0 if x<=0 or x>=1 else -x*log(x)-(1-x)*log(1-x)
f=lambda r: -(1-r)*np.log1p(-r)-r*np.log1p(r)
lo,hi=0.5,0.9
for _ in range(200):
    mid=(lo+hi)/2; lo,hi=(mid,hi) if f(mid)>0 else (lo,mid)
r=(lo+hi)/2; m=0.5+(1-r)*np.log1p(r)/(2*L); q=2*m-1; hr=h(r); p=log((1-r)/r); c=r-hr/p
kF=1/20; kE=11/50; Db=11/(1-1/(44*q*L)); gam=hr/r; alpha=1-gam/p; t0=r-m; bm=-hr*gam/(p*alpha); D=-hr/(p*alpha)-t0
es=gam*m; a_s=m/r; astar=m*h(r*r)/(r*hr)
def minM(eps,fQ3=1.0,Gam=63,z=77/200,lam=1.232,a0=3/8,kap=0.39,C=13/25,bY=0.4,n=240,split=False,te=0.463066,fCm=1.0):
    Cb=Db/te**2; ell=eps/m
    k0=(1/a_s+eps/(c*te*a_s))/es; k1=k0+1/(-p*c*te*a_s); k2=Cb/(c*te*a_s); BE=(alpha/(1-ell))**2/(2*kE)
    Ds=max(L-es,es-z); d0=2*m-astar+2*q*L*a_s/(p*c); D0=lam-astar; rD=-d0/D0
    beta=1-D/kap; b0=a0**beta; W=(kap-D)/(1-b0); eta=(kap-D)/D*(b0-a0)/(1-b0)
    kc=t0-bm*eps/es; kh=eta*kc; A0=W/2; B=W*(a0-0.5); Cm=fCm*(bm*Cb+(1+C)/(2*kF*te)); Btot=eta*(kc+D)+B+(kh+m)*rD
    def M_at(u,u1,u2):
        C0=(gam*(eps-u/2)/z+gam*eps/es+gam*Cb*u)/(1-ell)
        Eb=C0+BE*k1*u/2+sqrt(BE*u)*sqrt(k0+k2*u+k1*C0+BE*k1**2*u/4)
        nn=(1+Eb)/es; K=k0+k1*Eb+k2*u
        A2=gam*Ds/(1-bY)*((eps-u/2)/z+eps*nn)+u/(1-bY)*(alpha**2/(8*bY*kE*te)+gam*(es-z)*Cb)
        B0=q*L/(p*p*c*c)*A2+2*eps+q*L*eps**2/(c*c*z)
        rho0=(1-astar+B0+Gam*(2*eps-u))/D0
        Q=gam*(A0+Btot)*(eps-u/2)/z+gam*(A0-Btot)*eps*nn+fQ3*(A0+Btot)*alpha*sqrt(u1*K/(2*kE))+(A0-Btot)*gam*Cb*u2
        return kh-(kh+m)*rho0-Q-Cm*u2
    best=1e9
    for v in np.linspace(0,1,n):
        u=2*eps*v*v
        if split:
            val=min(M_at(u,u*w,u*(1-w)) for w in np.linspace(0,1,41))
        else:
            val=M_at(u,u,u)
        best=min(best,val)
    return best
def maxeps(**kw):
    lo,hi=0.0,3e-3
    if minM(1e-9,**kw)<=0: return 0.0
    for _ in range(36):
        mid=(lo+hi)/2
        if minM(mid,**kw)>0: lo=mid
        else: hi=mid
    return lo*1e6
print("Theorem V model (t_min''=0.463):", round(maxeps(),1),"ppm")
for split in [False, True]:
  for fQ3 in [1.0, 0.9, 0.83]:
    row=[maxeps(fQ3=fQ3,Gam=G,split=split) for G in [63,35,20,10,0]]
    print(f"split={split!s:5} fQ3={fQ3}: Gamma 63/35/20/10/0 ->", " ".join(f"{x:6.1f}" for x in row))
print("split, fQ3=0.9, Cm x0.5, Gamma 20/0:", round(maxeps(fQ3=0.9,Gam=20,split=True,fCm=0.5),1), round(maxeps(fQ3=0.9,Gam=0,split=True,fCm=0.5),1))
