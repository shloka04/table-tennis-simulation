"""Convert SpinDOE trajectories to validate_real_data columns: t,x,y,spin_rpm.
Frame: Z=0 table surface (height->y, no floor offset); Y=flight axis, X=lateral.
Only the X-axis spin component is topspin/backspin; keep spin-in-plane shots and
trim each to its pre-first-bounce arc. spin sign left as SpinDOE's (verified separately)."""
import os, csv, numpy as np, pandas as pd
RPM = 60/(2*np.pi)
def first_arc_end(Z):
    n=len(Z); end=n-1
    c=np.where(Z<=0.005)[0]; c=c[c>=3]
    if len(c): end=min(end,int(c[0]))
    for i in range(2,n-1):
        if Z[i]<Z[i-1] and Z[i]<0.10 and (Z[i:].max()-Z[i])>0.03:
            end=min(end,i); break
    return end
def convert(src, out, flip=False, topfrac_min=0.6, side_max_rpm=1500, spin_cap_rpm=9000):
    idx=pd.read_csv(os.path.join(src,"index.csv")); os.makedirs(out,exist_ok=True)
    kept=0; rpms=[]
    for _,r in idx.iterrows():
        f=os.path.join(src,"%03d.csv"%int(r['traj_file']))
        if not os.path.exists(f): continue
        a=np.genfromtxt(f,delimiter=';')
        if a.ndim<2 or len(a)<10: continue
        t,X,Y,Z=a[:,0],a[:,1],a[:,2],a[:,3]
        top=r['x_spin']*RPM*(-1 if flip else 1)
        side=np.hypot(r['y_spin'],r['z_spin'])*RPM
        mag=np.sqrt(r['x_spin']**2+r['y_spin']**2+r['z_spin']**2)+1e-9
        if abs(r['x_spin'])/mag<topfrac_min: continue
        if side>side_max_rpm: continue
        if abs(top)>spin_cap_rpm: continue
        e=first_arc_end(Z)
        t,X,Y,Z=t[:e+1],X[:e+1],Y[:e+1],Z[:e+1]
        if len(t)<10: continue
        step=np.hypot(np.diff(X),np.diff(Y)); horiz=np.concatenate([[0.0],np.cumsum(step)])
        t=t-t[0]
        with open(os.path.join(out,"%03d.csv"%int(r['traj_file'])),'w',newline='') as fo:
            w=csv.writer(fo); w.writerow(['t','x','y','spin_rpm'])
            for i in range(len(t)): w.writerow([f"{t[i]:.5f}",f"{horiz[i]:.5f}",f"{Z[i]:.5f}",f"{top:.1f}"])
        kept+=1; rpms.append(top)
    rpms=np.array(rpms)
    print(f"wrote {kept} shots to {out}  |  topspin(+) {int((rpms>200).sum())}, backspin(-) {int((rpms<-200).sum())}, rpm range {rpms.min():.0f}..{rpms.max():.0f}")
    return kept
if __name__=='__main__':
    import sys
    src=sys.argv[1] if len(sys.argv)>1 else 'spindoe/data/trajectories/_t'
    convert(src,'shots',flip=False)
