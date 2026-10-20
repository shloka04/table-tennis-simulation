"""
validate_real_data.py -- validate the flight model against real tracked shots,
from ANY source: a public dataset, tracked broadcast/phone video, or robot logs.

The tool only needs, per shot, the ball's position over time. Spin is optional
(measured spin unlocks the strongest test; without it you still get Tiers 1-2).

--------------------------------------------------------------------------------
QUICK START
  # self-check on synthetic data (no files needed):
  python validate_real_data.py

  # a dataset/CSV with your own column names, positions in millimetres,
  # height measured from the FLOOR (so subtract table height 0.76 m):
  python validate_real_data.py shot.csv --tcol time --xcol X --ycol Z \
      --units mm --y-from-floor

  # only frame numbers, no time column, 240 fps footage:
  python validate_real_data.py shot.csv --framecol frame --fps 240

  # broadcast video with no spin info -> fit-only (Tiers 1-2, honest labelling):
  python validate_real_data.py *.csv --fit-only
--------------------------------------------------------------------------------
"""
import sys, os, csv, argparse
import numpy as np
from scipy.optimize import minimize
from tt_engine import TableTennisSim
s = TableTennisSim(); DT = 0.001; TABLE_H = 0.76

# ---------------- model helpers ----------------
def model_traj(v0, ang, rpm, h, x0, t_end):
    xs, ys, *_ = s.simulate(v0, ang, rpm, h=h, x0=x0, t_max=t_end + 0.05)
    return np.arange(len(xs)) * DT, np.array(xs), np.array(ys)

def model_at(v0, ang, rpm, h, x0, times):
    t, xs, ys = model_traj(v0, ang, rpm, h, x0, times.max())
    return np.interp(times, t, xs), np.interp(times, t, ys)

def landing_x(v0, ang, rpm, h, x0):
    t, xs, ys = model_traj(v0, ang, rpm, h, x0, 3.0)
    for i in range(1, len(ys)):
        if ys[i] <= 0 < ys[i-1]:
            f = ys[i-1]/(ys[i-1]-ys[i]); return xs[i-1]+f*(xs[i]-xs[i-1])
    return xs[-1]

def launch_from_early(t, x, y, k=6):
    k = min(k, len(t)); tt = t[:k] - t[0]
    if k < 3: return np.hypot((x[1]-x[0])/(t[1]-t[0]), (y[1]-y[0])/(t[1]-t[0])), \
                     np.degrees(np.arctan2(y[1]-y[0], x[1]-x[0]))
    vx0 = np.polyfit(tt, x[:k], 2)[1]; vy0 = np.polyfit(tt, y[:k], 2)[1]
    return np.hypot(vx0, vy0), np.degrees(np.arctan2(vy0, vx0))

def fit(times, xo, yo, h, x0):
    v0g, angg = launch_from_early(times, xo, yo)
    def cost(p):
        v0, ang, rpm = p
        if v0 < 1 or v0 > 40 or abs(ang) > 70: return 1e9
        mx, my = model_at(v0, ang, rpm, h, x0, times)
        return np.sum((mx-xo)**2 + (my-yo)**2)
    r = minimize(cost, [v0g, angg, 2000], method='Nelder-Mead',
                 options={'xatol':1e-3,'fatol':1e-10,'maxiter':3000})
    return r.x

def rmse_r2(xo, yo, mx, my):
    rmse = np.sqrt(np.mean((mx-xo)**2 + (my-yo)**2))
    obs = np.concatenate([xo, yo]); pred = np.concatenate([mx, my])
    ss = np.sum((obs-pred)**2); tot = np.sum((obs-obs.mean())**2)
    return rmse, (1 - ss/tot if tot > 0 else float('nan'))

# ---------------- flexible loader ----------------
def load(path, a):
    rows = list(csv.DictReader(open(path)))
    cols = rows[0].keys()
    def col(name): return [float(r[name]) for r in rows]
    if a.tcol and a.tcol in cols:       t = np.array(col(a.tcol))
    elif a.framecol and a.framecol in cols: t = np.array(col(a.framecol)) / a.fps
    else:                                t = np.arange(len(rows)) / a.fps
    x = np.array(col(a.xcol)); y = np.array(col(a.ycol))
    scale = {'m':1.0, 'cm':0.01, 'mm':0.001}[a.units]
    x *= scale; y *= scale
    if a.y_from_floor: y = y - TABLE_H          # convert floor-height to table-height
    x = x - x[0]                                # put contact point at x=0
    spin = None
    if a.spincol and a.spincol in cols and rows[0].get(a.spincol):
        spin = float(rows[0][a.spincol])
    # keep only the airborne, descending-into-table portion (y>=0)
    m = y >= -0.02
    return t[m], x[m], y[m], spin

# ---------------- analysis ----------------
def analyse(path, a):
    t, x, y, spin = load(path, a)
    if len(t) < 6:
        print(f"\n=== {os.path.basename(path)}: only {len(t)} usable points, skipped ==="); return None
    h, x0 = y[0], x[0]
    print(f"\n=== {os.path.basename(path)}  ({len(t)} frames, {(t[-1]-t[0])*1000:.0f} ms) ===")
    v1, a1, r1 = fit(t, x, y, h, x0)
    mx, my = model_at(v1, a1, r1, h, x0, t); rmse, r2 = rmse_r2(x, y, mx, my)
    print(f"[Tier 1 descriptive] v0={v1:.2f} m/s, angle={a1:.2f} deg, spin={r1:.0f} RPM "
          f"| RMSE {rmse*1000:.1f} mm, R^2 {r2:.4f}")
    cut = t <= (t[0] + 0.4*(t[-1]-t[0]))
    lx_pred = np.nan
    if cut.sum() >= 4:
        v2, a2, r2s = fit(t[cut], x[cut], y[cut], h, x0)
        mxh, myh = model_at(v2, a2, r2s, h, x0, t[~cut])
        prmse, _ = rmse_r2(x[~cut], y[~cut], mxh, myh)
        lx_pred = landing_x(v2, a2, r2s, h, x0)
        print(f"[Tier 2 hold-out]    fit first 40%, predict rest | pred RMSE {prmse*1000:.1f} mm "
              f"| landing pred {lx_pred:.3f} vs last-tracked {x[-1]:.3f} m ({abs(lx_pred-x[-1])*100:.1f} cm)")
    if spin is not None and not a.fit_only:
        v3, a3 = launch_from_early(t, x, y)
        mx3, my3 = model_at(v3, a3, spin, h, x0, t); rmse3, r23 = rmse_r2(x, y, mx3, my3)
        print(f"[Tier 3 no-fit]      MEASURED spin={spin:.0f} RPM, v0/angle from early frames, "
              f"predict all | RMSE {rmse3*1000:.1f} mm, R^2 {r23:.4f}")
    elif a.fit_only or spin is None:
        print("[Tier 3 skipped]     no measured spin -> descriptive/hold-out only "
              "(cannot claim prediction)")
    return dict(rmse=rmse, r2=r2, landing_err=abs(lx_pred-x[-1])*100 if not np.isnan(lx_pred) else np.nan)

def make_synthetic(name, v0, ang, rpm):
    xs, ys, *_ = s.simulate(v0, ang, rpm, h=0.30, x0=0.0)
    step = int((1/120)/DT); idx = np.arange(0, len(xs), step); rng = np.random.default_rng(1)
    with open(name, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['t','x','y','spin_rpm'])
        for i in idx: w.writerow([i*DT, xs[i]+rng.normal(0,0.004), ys[i]+rng.normal(0,0.004), rpm])

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="Validate the flight model against tracked shots.")
    ap.add_argument('files', nargs='*')
    ap.add_argument('--tcol', default='t'); ap.add_argument('--framecol', default=None)
    ap.add_argument('--xcol', default='x'); ap.add_argument('--ycol', default='y')
    ap.add_argument('--spincol', default='spin_rpm')
    ap.add_argument('--fps', type=float, default=120)
    ap.add_argument('--units', choices=['m','cm','mm'], default='m')
    ap.add_argument('--y-from-floor', dest='y_from_floor', action='store_true',
                    help="subtract table height (0.76 m) if y is measured from the floor")
    ap.add_argument('--fit-only', dest='fit_only', action='store_true',
                    help="no spin available: report Tiers 1-2 only, never claim prediction")
    a = ap.parse_args()
    if not a.files:
        print("No files -> self-check on synthetic data.")
        make_synthetic('_syn_loop.csv', 11, 7, 3000); make_synthetic('_syn_flat.csv', 9, 5, 0)
        a.files = ['_syn_loop.csv', '_syn_flat.csv']
    res = [r for r in (analyse(f, a) for f in a.files) if r]
    if res:
        import numpy as np
        print(f"\n=== SUMMARY ({len(res)} shots) ===")
        print(f"  median Tier-1 RMSE : {np.nanmedian([r['rmse'] for r in res])*1000:.1f} mm")
        print(f"  median Tier-1 R^2  : {np.nanmedian([r['r2'] for r in res]):.4f}")
        le=[r['landing_err'] for r in res if not np.isnan(r['landing_err'])]
        if le: print(f"  median landing err : {np.median(le):.1f} cm")
