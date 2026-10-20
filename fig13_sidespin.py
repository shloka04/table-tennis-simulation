import os,numpy as np,matplotlib.pyplot as plt
from tt_engine3d import TableTennisSim3D
FIG="figures"; os.makedirs(FIG,exist_ok=True)
s=TableTennisSim3D()
V,ANG,H=7.0,12.0,0.15
cases=[(0,'0 RPM (straight)','#888888'),(2000,'2000 RPM','#2ecc71'),
       (4000,'4000 RPM','#f39c12'),(6000,'6000 RPM','#e74c3c')]
fig,ax=plt.subplots(figsize=(11,5))
# table outline (top-down): x in [0,2.74], z in [-0.7625,0.7625]
ax.add_patch(plt.Rectangle((0,-s.half_width),s.L,2*s.half_width,fill=False,ec='saddlebrown',lw=2))
ax.plot([s.net_x,s.net_x],[-s.half_width,s.half_width],'k--',lw=1.5,label='Net')
defl={}
for rs,lab,c in cases:
    xs,ys,zs,o,xl,zl=s.simulate(V,ANG,0,rs,h=H,x0=0.15)
    ax.plot(xs,zs,color=c,lw=2.5,label=f"{lab}  ({o}, {zl*100:+.0f} cm)")
    defl[rs]=zl
ax.set_xlabel("Forward distance x (m)"); ax.set_ylabel("Lateral deflection z (m)")
ax.set_title("Sidespin Serve: Top-Down View of Lateral Curve\n"
             f"(launch {V:.0f} m/s, {ANG:.0f} deg; pure sidespin; 3D engine)")
ax.set_xlim(0,2.9); ax.set_ylim(-0.85,0.85); ax.legend(loc='upper left',fontsize=9); ax.grid(alpha=0.3)
fig.savefig(f"{FIG}/fig13_sidespin.png",dpi=300,bbox_inches='tight'); plt.close(fig)
print("[F13] lateral deflection at landing:")
for rs in [0,2000,4000,6000]:
    print(f"   {rs} RPM sidespin -> {defl[rs]*100:+.1f} cm")

# deflection vs sidespin (saturation curve)
rss=np.arange(0,6001,500); dd=[abs(s.simulate(V,ANG,0,rr,h=H,x0=0.15)[5])*100 for rr in rss]
fig,ax=plt.subplots(figsize=(9,4))
ax.plot(rss,dd,'o-',color='purple',lw=2); ax.set_xlabel("Sidespin (RPM)")
ax.set_ylabel("Lateral deflection at landing (cm)")
ax.set_title("Lateral Deflection vs Sidespin (saturating, as C_L saturates)"); ax.grid(alpha=0.3)
fig.savefig(f"{FIG}/fig13b_sidespin_curve.png",dpi=300,bbox_inches='tight'); plt.close(fig)
print(f"[F13b] deflection saturates: 6000 RPM gives {dd[-1]:.0f} cm vs 2000 RPM {dd[4]:.0f} cm "
      f"(ratio {dd[-1]/dd[4]:.2f}x for 3x the spin)")
