import numpy as np, glob, os, argparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import validate_real_data as V
def A(): return argparse.Namespace(tcol='t',framecol=None,xcol='x',ycol='y',spincol='spin_rpm',
    fps=142.9,units='m',y_from_floor=False,fit_only=False)
MED=37.8
shots=[]
for f in sorted(glob.glob('shots/*.csv')):
    t,x,y,spin=V.load(f,A())
    if len(t)<10: continue
    h,x0=y[0],x[0]; v3,a3=V.launch_from_early(t,x,y)
    mx,my=V.model_at(v3,a3,spin,h,x0,t); rmse=V.rmse_r2(x,y,mx,my)[0]*1000
    shots.append((f,spin,rmse,t,x,y,v3,a3,h,x0))
# representative topspin & backspin: RMSE closest to the median
top=[sh for sh in shots if sh[1]>200]; bak=[sh for sh in shots if sh[1]<-200]
tp=min(top,key=lambda z:abs(z[2]-MED)); bk=min(bak,key=lambda z:abs(z[2]-MED))
fig,axes=plt.subplots(1,2,figsize=(9.2,3.6))
for ax,(f,spin,rmse,t,x,y,v3,a3,h,x0),ttl in [(axes[0],tp,'Topspin'),(axes[1],bk,'Backspin')]:
    tm,xm,ym=V.model_traj(v3,a3,spin,h,x0,t.max())
    ax.axhline(0,color='0.6',lw=1.2,zorder=0)
    ax.plot(xm,ym,'-',color='0.15',lw=1.8,label='Model (blind prediction)',zorder=2)
    ax.plot(x,y,'o',mfc='none',mec='0.15',ms=4.5,mew=1.0,label='Recorded (SpinDOE)',zorder=3)
    ax.set_xlim(left=0); ax.set_ylim(bottom=min(-0.02,ym.min()-0.02))
    ax.set_xlabel('Horizontal distance (m)'); 
    ax.set_title(f'{ttl}: {abs(spin):.0f} RPM   (RMSE {rmse:.0f} mm)',fontsize=10)
    ax.grid(alpha=0.25)
axes[0].set_ylabel('Height above table (m)')
axes[0].legend(fontsize=8,loc='upper right',framealpha=0.9)
plt.tight_layout()
os.makedirs('figures',exist_ok=True)
plt.savefig('figures/fig16_real_validation.png',dpi=200,bbox_inches='tight')
print(f"topspin shot {os.path.basename(tp[0])} rmse {tp[2]:.1f} mm; backspin {os.path.basename(bk[0])} rmse {bk[2]:.1f} mm")
print("saved figures/fig16_real_validation.png")
