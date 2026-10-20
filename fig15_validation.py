"""
Parameter-recovery validation of the fitting pipeline on SYNTHETIC data with
known ground truth. This validates the METHOD (can we recover launch speed,
angle, and spin from noisy position samples?), not the physics. It is exactly
the pipeline to run on real high-speed-camera footage once available.
"""
import os,numpy as np,matplotlib.pyplot as plt
from scipy.optimize import minimize
from tt_engine import TableTennisSim
FIG="figures"; os.makedirs(FIG,exist_ok=True)
s=TableTennisSim(); DT=0.001

def sample(v0,ang,rpm,fps=120,h=0.30,x0=0.0):
    xs,ys,o,xl,d=s.simulate(v0,ang,rpm,h=h,x0=x0)
    step=int((1.0/fps)/DT)
    idx=np.arange(0,len(xs),step)
    return np.array(xs)[idx], np.array(ys)[idx]

def fit(xo,yo,fps=120,h=0.30,x0=0.0):
    step=int((1.0/fps)/DT)
    def cost(p):
        v0,ang,rpm=p
        if v0<1 or v0>30 or abs(ang)>60: return 1e6
        xs,ys,*_=s.simulate(v0,ang,rpm,h=h,x0=x0)
        xs=np.array(xs); ys=np.array(ys)
        idx=np.arange(0,len(xs),step)[:len(xo)]
        if len(idx)<len(xo): return 1e6
        return np.sum((xs[idx]-xo)**2+(ys[idx]-yo)**2)
    # initial guess from first two samples
    dtc=step*DT
    v0g=np.hypot(xo[1]-xo[0],yo[1]-yo[0])/dtc
    angg=np.degrees(np.arctan2(yo[1]-yo[0],xo[1]-xo[0]))
    r=minimize(cost,[v0g,angg,1000],method='Nelder-Mead',
               options={'xatol':1e-2,'fatol':1e-8,'maxiter':2000})
    return r.x

# ---- recovery experiment ----
TRUE=(11.0,7.0,3000.0); rng=np.random.default_rng(20260202)
print(f"[F15] ground truth: v0={TRUE[0]} m/s, angle={TRUE[1]} deg, spin={TRUE[2]} RPM")
for noise_mm in [2,5,10]:
    errs=[]
    for _ in range(40):
        xo,yo=sample(*TRUE)
        xo=xo+rng.normal(0,noise_mm/1000,len(xo)); yo=yo+rng.normal(0,noise_mm/1000,len(yo))
        v,a,r=fit(xo,yo)
        errs.append((v-TRUE[0],a-TRUE[1],r-TRUE[2]))
    errs=np.array(errs)
    print(f"  noise +-{noise_mm} mm (n=40): "
          f"v0 err {errs[:,0].mean():+.2f}+-{errs[:,0].std():.2f} m/s | "
          f"angle {errs[:,1].mean():+.2f}+-{errs[:,1].std():.2f} deg | "
          f"spin {errs[:,2].mean():+.0f}+-{errs[:,2].std():.0f} RPM")

# ---- figure: one noisy realisation + true + recovered ----
xo,yo=sample(*TRUE); xn=xo+rng.normal(0,0.005,len(xo)); yn=yo+rng.normal(0,0.005,len(yo))
vf,af,rf=fit(xn,yn)
xs_t,ys_t,*_=s.simulate(*TRUE,h=0.30,x0=0.0)
xs_f,ys_f,*_=s.simulate(vf,af,rf,h=0.30,x0=0.0)
# trajectory RMSE (fitted vs true, noiseless)
L=min(len(xs_t),len(xs_f)); rmse=np.sqrt(np.mean((np.array(xs_t[:L])-np.array(xs_f[:L]))**2+(np.array(ys_t[:L])-np.array(ys_f[:L]))**2))
fig,ax=plt.subplots(figsize=(10,5))
ax.plot(xs_t,ys_t,'-',color='green',lw=2,label='True trajectory')
ax.plot(xs_f,ys_f,'--',color='red',lw=2,label=f'Recovered fit (v={vf:.1f}, a={af:.1f}, spin={rf:.0f})')
ax.plot(xn,yn,'k.',ms=8,label='Noisy samples (120 fps, +-5 mm)')
ax.set_xlabel("Horizontal distance (m)"); ax.set_ylabel("Height above table (m)")
ax.set_title(f"Parameter-Recovery Validation on Synthetic Data\n"
             f"true (11 m/s, 7 deg, 3000 RPM); fit RMSE vs truth = {rmse*1000:.1f} mm")
ax.legend(); ax.grid(alpha=0.3)
fig.savefig(f"{FIG}/fig15_validation.png",dpi=300,bbox_inches='tight'); plt.close(fig)
print(f"[F15] single-fit trajectory RMSE vs ground truth = {rmse*1000:.1f} mm")
