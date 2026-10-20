"""Two rigor additions: (11) launch-error contraction, (12) MC convergence + CIs."""
import os, numpy as np, matplotlib.pyplot as plt
from tt_engine import TableTennisSim
FIG="figures"; os.makedirs(FIG,exist_ok=True)
sim=TableTennisSim(); SEED=20260202

# ---------- FIGURE 11: launch-error sensitivity contracts with spin ----------
v0,th0,h=12.0,8.0,0.15; dv=dth=0.15
spins=np.arange(0,6001,200)
def xl(v,th,r): return sim.simulate(v,th,r,h=h,x0=0.15)[3]
dXdv=np.array([abs((xl(v0+dv,th0,r)-xl(v0-dv,th0,r))/(2*dv)) for r in spins])
dXdth=np.array([abs((xl(v0,th0+dth,r)-xl(v0,th0-dth,r))/(2*dth)) for r in spins])
with plt.style.context('bmh'):
    fig,ax=plt.subplots(figsize=(10,5))
    ax.plot(spins,dXdv,color='crimson',lw=2.5,marker='o',ms=3,
            label=r'$|\partial x_{land}/\partial v_0|$  (m per m/s of launch-speed error)')
    ax.plot(spins,dXdth,color='navy',lw=2.5,marker='s',ms=3,
            label=r'$|\partial x_{land}/\partial \theta_0|$  (m per degree of launch-angle error)')
    ax.annotate(f"-{(1-dXdv[-1]/dXdv[0])*100:.0f}%",xy=(6000,dXdv[-1]),
                xytext=(-70,28),textcoords='offset points',
                arrowprops=dict(arrowstyle='->'),color='crimson',fontsize=10)
    ax.set_title("Topspin Contracts the Sensitivity of Landing Position to Launch Error\n"
                 "(nominal shot 12 m/s, 8 deg; deterministic Jacobian)")
    ax.set_xlabel("Topspin (RPM)"); ax.set_ylabel("Landing-position sensitivity")
    ax.set_ylim(0,None); ax.legend()
    fig.savefig(f"{FIG}/fig11_sensitivity.png",dpi=300,bbox_inches='tight'); plt.close(fig)
print(f"[F11] |dX/dv| {dXdv[0]:.3f} -> {dXdv[-1]:.3f} m/(m/s) (-{(1-dXdv[-1]/dXdv[0])*100:.0f}%); "
      f"|dX/dth| {dXdth[0]:.3f} -> {dXdth[-1]:.3f} m/deg (-{(1-dXdth[-1]/dXdth[0])*100:.0f}%)")

# ---------- FIGURE 12: Monte Carlo convergence + 95% CI ----------
def success(rpm,bv,ba,br,lam,N,rng):
    s=0
    for _ in range(N):
        v=bv+rng.normal(0,lam*1.5); a=ba+rng.normal(0,lam*1.5)
        r=max(0.0,br+rng.normal(0,lam*400))
        if sim.simulate(v,a,r,h=0.15,x0=0.15)[2]=='safe': s+=1
    return s/N
def wilson(p,n,z=1.96):
    d=1+z*z/n; c=(p+z*z/(2*n))/d; hw=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-hw,c+hw
Ns=[25,50,100,200,400,800,1600]
top=[];flat=[];topci=[];flatci=[]
for N in Ns:
    rng=np.random.default_rng(SEED)
    pt=success(3500,12,8,3500,0.5,N,rng); pf=success(0,9,5,0,0.5,N,rng)
    top.append(pt*100); flat.append(pf*100)
    lo,hi=wilson(pt,N); topci.append((pt*100-lo*100,hi*100-pt*100))
    lo,hi=wilson(pf,N); flatci.append((pf*100-lo*100,hi*100-pf*100))
topci=np.array(topci).T; flatci=np.array(flatci).T
with plt.style.context('bmh'):
    fig,ax=plt.subplots(figsize=(10,5))
    ax.errorbar(Ns,top,yerr=topci,color='red',marker='o',capsize=4,lw=2,label='Topspin loop (3500 RPM)')
    ax.errorbar(Ns,flat,yerr=flatci,color='blue',marker='s',capsize=4,lw=2,label='Flat drive (0 RPM)')
    ax.set_xscale('log'); ax.set_xticks(Ns); ax.set_xticklabels([str(n) for n in Ns])
    ax.set_title("Monte Carlo Convergence with 95% Confidence Intervals\n"
                 "(safe-landing rate at execution-error scale lambda=0.5)")
    ax.set_xlabel("Number of Monte Carlo shots N (log scale)")
    ax.set_ylabel("Probability of safe landing (%)"); ax.legend()
    fig.savefig(f"{FIG}/fig12_convergence.png",dpi=300,bbox_inches='tight'); plt.close(fig)
# final CIs at N=500 for the text
rng=np.random.default_rng(SEED)
pt=success(3500,12,8,3500,0.5,500,rng); pf=success(0,9,5,0,0.5,500,rng)
lt=wilson(pt,500); lf=wilson(pf,500)
print(f"[F12] N=500 lambda=0.5: topspin {pt*100:.0f}% (95% CI {lt[0]*100:.0f}-{lt[1]*100:.0f}), "
      f"flat {pf*100:.0f}% (95% CI {lf[0]*100:.0f}-{lf[1]*100:.0f})")
print("DONE")
