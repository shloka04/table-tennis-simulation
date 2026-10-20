import os,numpy as np,matplotlib.pyplot as plt
from tt_engine import TableTennisSim
FIG="figures"; os.makedirs(FIG,exist_ok=True)
s=TableTennisSim(); ANG,H=8.0,0.15; LAM=0.4; N=300; SEED=20260202
def safe_frac(v,rpm,rng):
    c=0
    for _ in range(N):
        vv=v+rng.normal(0,LAM*1.5); aa=ANG+rng.normal(0,LAM*1.5); rr=max(0,rpm+rng.normal(0,LAM*400))
        if s.simulate(vv,aa,rr,h=H,x0=0.15)[2]=='safe': c+=1
    return c/N

spins=np.arange(0,8001,400)
# panel A: safe-fraction vs spin for a few speeds (shows interior optimum)
fig,(a1,a2)=plt.subplots(1,2,figsize=(14,5))
for V,c in [(10,'#2ecc71'),(12,'#f39c12'),(14,'#e74c3c')]:
    rng=np.random.default_rng(SEED)
    fr=np.array([safe_frac(V,r,rng) for r in spins])*100
    a1.plot(spins,fr,'o-',color=c,ms=3,lw=2,label=f"{V} m/s")
    a1.plot(spins[np.argmax(fr)],fr.max(),'*',color=c,ms=16,mec='k')
a1.set_xlabel("Topspin (RPM)"); a1.set_ylabel("Safe-landing probability (%)")
a1.set_title("Error tolerance vs spin (stars = optimum)\nover-spinning a slow ball is counter-productive")
a1.legend(title="Launch speed"); a1.grid(alpha=0.3)

# panel B: optimal spin vs launch speed
Vs=np.arange(9.5,15.1,0.5); optn=[]
for V in Vs:
    rng=np.random.default_rng(SEED)
    fr=np.array([safe_frac(V,r,rng) for r in spins])
    optn.append(spins[int(np.argmax(fr))])
optn=np.array(optn)
a2.plot(Vs,optn,'o-',color='navy',lw=2.5,ms=5)
a2.set_xlabel("Launch speed (m/s)"); a2.set_ylabel("Error-tolerance-optimal topspin (RPM)")
a2.set_title("Optimal topspin rises with pace\n(spin that maximises safe-landing rate under execution error)")
a2.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(f"{FIG}/fig14_optimal_spin.png",dpi=300,bbox_inches='tight'); plt.close(fig)
print("[F14] optimal topspin (max error tolerance) vs speed:")
for V,o in zip(Vs,optn): print(f"   {V:.1f} m/s -> {o} RPM")
