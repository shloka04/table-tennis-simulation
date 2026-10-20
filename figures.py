"""Generate all ten figures from the single unified engine (tt_engine)."""
import os, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from tt_engine import TableTennisSim

FIG = "figures"; os.makedirs(FIG, exist_ok=True)
SEED = 20260202
sim  = TableTennisSim()
DT   = 0.001
SIG  = {"calm (0.016)": 0.5*DT**0.5, "club (0.063)": 2.0*DT**0.5,
        "adverse (0.158)": 5.0*DT**0.5}

def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=300, bbox_inches="tight")
    plt.show(); plt.close(fig)

# ============================ FIGURES 1 & 2 ============================
spins = np.linspace(0, 6000, 40); vels = np.linspace(3, 20, 40)
Z = np.zeros((len(vels), len(spins))); crit = []
for i, rpm in enumerate(spins):
    msafe = 0.0
    for j, v in enumerate(vels):
        o = sim.simulate(v, 8.0, rpm, h=0.15, x0=0.15)[2]
        Z[j, i] = {"net": 0, "safe": 1, "long": 2}[o]
        if o == "safe": msafe = max(msafe, v)
    crit.append(msafe)
crit = np.array(crit)

cmap = ListedColormap(['#e74c3c', '#2ecc71', '#3498db'])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
fig, ax = plt.subplots(figsize=(10, 6))
mesh = ax.pcolormesh(*np.meshgrid(spins, vels), Z, cmap=cmap, norm=norm, shading='auto')
ax.legend(handles=[Patch(color='#e74c3c', label='Net collision'),
                   Patch(color='#2ecc71', label='Safe landing'),
                   Patch(color='#3498db', label='Long (overshoots)')],
          loc='upper left', fontsize=11, framealpha=0.85)
ax.set_title("Landing Outcome vs Topspin and Launch Speed\n"
             "(launch angle 8\u00b0, contact height 0.15 m; deterministic RK4)", fontsize=13)
ax.set_xlabel("Topspin (RPM)", fontweight='bold')
ax.set_ylabel("Launch speed (m/s)", fontweight='bold')
cb = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2]); cb.ax.set_yticklabels(['Net', 'Safe', 'Long'])
cb.set_label('Landing outcome')
save(fig, "fig1_phase_diagram.png")

crit_spins = np.linspace(0, 6000, 40)
crit_v = np.array([sim.max_safe_speed(8.0, rpm, h=0.15, x0=0.15) or np.nan
                   for rpm in crit_spins])
with plt.style.context('bmh'):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(crit_spins, crit_v, color='purple', lw=2.5)
    ax.fill_between(crit_spins, crit_v, color='purple', alpha=0.2)
    v0, v6 = crit_v[0], crit_v[-1]
    ax.annotate(f"{v0:.1f} m/s", (0, v0), textcoords="offset points", xytext=(8, 6))
    ax.annotate(f"{v6:.1f} m/s  (+{(v6/v0-1)*100:.0f}% vs no spin)", (6000, v6),
                textcoords="offset points", xytext=(-165, -16))
    ax.set_title("Maximum Safe Launch Speed vs Topspin  (8\u00b0 launch)", fontsize=13)
    ax.set_xlabel("Topspin (RPM)"); ax.set_ylabel("Maximum safe launch speed (m/s)")
    ax.set_ylim(0, None)
    save(fig, "fig2_critical_velocity.png")
v3 = crit_v[np.argmin(abs(crit_spins-3000))]
print(f"[F1/2] N=1600 shots. critical v 0/3000/6000 = {v0:.2f}/{v3:.2f}/{v6:.2f} m/s "
      f"(+{(v3/v0-1)*100:.0f}% / +{(v6/v0-1)*100:.0f}%)")

# ============================ FIGURE 3 ============================
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot([0, sim.L], [0, 0], color='orange', lw=5, label='Table')
ax.plot([sim.net_x, sim.net_x], [0, sim.net_h], color='black', lw=2, label='Net')
runs = [("Topspin", +3000, 'red'), ("No spin", 0, 'blue'), ("Backspin", -3000, 'green')]
res = {}
for name, rpm, c in runs:
    xs, ys, o, xl, d = sim.simulate(6.0, 6.5, rpm, h=0.39, x0=0.0)
    res[name] = (xl, d); ax.plot(xs, ys, color=c, lw=2.5, zorder=4)
d_no = res["No spin"][1]; handles = []
for name, rpm, c in runs:
    xl, d = res[name]
    if rpm == 0:
        t = f"No spin, 0 RPM  -  descent {d:.1f} deg, lands {xl:.2f} m"
    else:
        dd = d - d_no
        t = (f"{name}, {rpm:+d} RPM  -  descent {d:.1f} deg "
             f"({dd:+.1f} deg {'steeper' if dd>0 else 'flatter'}), lands {xl:.2f} m")
    handles.append(Line2D([0], [0], color=c, lw=2.5, label=t))
th, _ = ax.get_legend_handles_labels()
ax.legend(handles=th[:2] + handles, fontsize=9, loc='upper right')
ax.set_xlim(0, 3.0); ax.set_ylim(0, 1.0)
ax.set_xlabel("Horizontal Distance (m)", fontweight='bold')
ax.set_ylabel("Height above table (m)", fontweight='bold')
ax.set_title("How Spin Direction Shapes the Trajectory\n"
             "(all launched at 6.0 m/s, 6.5\u00b0, 0.39 m contact height)",
             fontsize=13, fontweight='bold')
save(fig, "fig3_spin_direction.png")
print(f"[F3] descent: topspin {res['Topspin'][1]-d_no:+.2f} deg, backspin {res['Backspin'][1]-d_no:+.2f} deg")

# ============================ FIGURES 4 & 5 ============================
rng = np.random.default_rng(SEED)
spins2 = np.linspace(0, 6000, 25); ITER = 600
def batch(rpm, sg, keep):
    out = []
    for _ in range(ITER):
        _, _, o, xl, _ = sim.simulate(12, 8, rpm, h=0.15, x0=0.15, sigma=sg, rng=rng)
        if o in keep: out.append(xl)
    return np.array(out)

with plt.style.context('bmh'):
    fig, ax = plt.subplots(figsize=(10, 5))
    stds = {}
    for label, sg in SIG.items():
        col = []
        for rpm in spins2:
            L = batch(rpm, sg, ('safe',))
            col.append(np.std(L) if len(L) > 5 else np.nan)
        stds[label] = col
        ax.plot(spins2, col, marker='o', ms=4, lw=2,
                label=f"{label.split()[0]}:  \u03c3 = {sg:.3f} m\u00b7s\u207b\u00b9\u141f\u00b2")
    ax.set_title("Landing-Position Spread vs Topspin, Three Turbulence Levels\n"
                 "(launch 12 m/s, 8\u00b0; 600 Monte-Carlo shots per point)")
    ax.set_xlabel("Topspin (RPM)"); ax.set_ylabel("Landing-position std. dev. (m)")
    ax.legend(title="Turbulence intensity (velocity-diffusion)")
    ax.text(0.015, 0.03,
            "Curves begin where \u22655 of 600 shots land in-bounds;\n"
            "at 12 m/s the low-spin shots overshoot and are excluded.",
            transform=ax.transAxes, fontsize=8, va='bottom',
            bbox=dict(boxstyle='round', fc='white', alpha=0.75))
    save(fig, "fig4_landing_variance_vs_spin.png")

    msd = []
    for rpm in spins2:
        L = batch(rpm, SIG["club (0.063)"], ('safe', 'long'))
        msd.append(np.mean((L - L.mean())**2) if len(L) > 5 else np.nan)
    m0 = next(m for m in msd if not np.isnan(m)); mN = msd[-1]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(spins2, msd, color='darkred', marker='s', ms=4, lw=2)
    ax.annotate(f"\u2212{(1-mN/m0)*100:.0f}% from 0 \u2192 6000 RPM",
                xy=(6000, mN), xytext=(-185, 34), textcoords='offset points',
                arrowprops=dict(arrowstyle='->'), fontsize=10)
    ax.set_title("Mean-Squared Landing Displacement vs Topspin\n"
                 "(club turbulence \u03c3 = 0.063 m\u00b7s\u207b\u00b9\u141f\u00b2; 600 shots per point)")
    ax.set_xlabel("Topspin (RPM)"); ax.set_ylabel("MSD of landing position (m\u00b2)")
    save(fig, "fig5_msd_vs_spin.png")
print(f"[F5] MSD {m0:.6f} -> {mN:.6f} m^2  ({(1-mN/m0)*100:.0f}% reduction)")

# ============================ FIGURE 6 ============================
rng = np.random.default_rng(SEED)
spins3 = np.linspace(0, 6000, 35); IT6 = 600; ent = []; sd = []; lands_at = {}
for rpm in spins3:
    L = []
    for _ in range(IT6):
        _, _, o, xl, _ = sim.simulate(12, 8, rpm, h=0.15, x0=0.15, sigma=SIG["club (0.063)"], rng=rng)
        if o != 'net': L.append(xl)
    L = np.array(L)
    if len(L) < 20: ent.append(np.nan); sd.append(np.nan); continue
    lands_at[rpm] = L
    c = L - L.mean()
    counts, _ = np.histogram(c, bins=np.linspace(-0.15, 0.15, 30))
    p = counts / counts.sum(); p = p[p > 0]
    ent.append(-np.sum(p * np.log2(p))); sd.append(np.std(L))
ent = np.array(ent)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5))
a1.plot(spins3, ent, color='darkblue', lw=2.5, marker='o', ms=4)
a1.fill_between(spins3, ent, alpha=0.12, color='blue')
a1.set_title("Landing-Distribution Entropy vs Topspin"); a1.set_xlabel("Topspin (RPM)")
a1.set_ylabel("Shannon entropy H (bits)"); a1.grid(alpha=0.4)
# fresh distributions at EXACT round spin rates for clean labels
for rpm, col in zip([0, 2000, 5000], ['#e74c3c', '#f39c12', '#2ecc71']):
    Lr = []
    for _ in range(IT6):
        _, _, o, xl, _ = sim.simulate(12, 8, rpm, h=0.15, x0=0.15,
                                      sigma=SIG["club (0.063)"], rng=rng)
        if o != 'net':
            Lr.append(xl)
    Lr = np.array(Lr); c = Lr - Lr.mean()
    cnt, _ = np.histogram(c, bins=np.linspace(-0.15, 0.15, 30))
    p = cnt / cnt.sum(); p = p[p > 0]; Hr = -np.sum(p * np.log2(p))
    a2.hist(Lr, bins=30, alpha=0.5, color=col, density=True,
            label=f"{rpm} RPM  (H = {Hr:.2f} bits)")
a2.axvline(sim.net_x, color='green', ls='--', label='Net (1.37 m)')
a2.axvline(sim.L, color='black', ls='--', label='Far edge (2.74 m)')
a2.set_title("Landing Distributions at Three Topspin Rates"); a2.set_xlabel("Landing position (m)")
a2.set_ylabel("Probability density"); a2.legend(fontsize=9); a2.grid(alpha=0.4)
fig.suptitle("Information-Theoretic View: Topspin Suppresses Landing Entropy\n"
             "(entropy from a fixed 30-bin grid, so bin width is constant across spins)",
             fontsize=12, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.93])
save(fig, "fig6_shannon_entropy.png")
H0 = next(h for h in ent if not np.isnan(h))
print(f"[F6] entropy {H0:.2f} -> {ent[-1]:.2f} bits ({(1-ent[-1]/H0)*100:.0f}% reduction)")

# ============================ FIGURES 7 & 8 ============================
rng = np.random.default_rng(SEED)
def mc(rpm, runs=300):
    trs, lands = [], []
    for _ in range(runs):
        v = 12.5 + rng.normal(0, 0.4); a = 8 + rng.normal(0, 0.4)
        xs, ys, o, xl, _ = sim.simulate(v, a, rpm, h=0.15, x0=0.15,
                                        sigma=SIG["club (0.063)"], rng=rng)
        trs.append((xs, ys)); lands.append(xl)
    return trs, np.array(lands)
tp, lp = mc(3500); ta, la = mc(0)
on_tp = np.mean((lp > sim.net_x) & (lp <= sim.L)) * 100
on_ta = np.mean((la > sim.net_x) & (la <= sim.L)) * 100
vr = np.var(la) / np.var(lp)
with plt.style.context('bmh'):
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot([0, sim.L], [0, 0], color='orange', lw=6, zorder=3, label='Table')
    ax.plot([sim.net_x, sim.net_x], [0, sim.net_h], color='black', lw=3, zorder=3, label='Net')
    for xs, ys in tp: ax.plot(xs, ys, color='red', alpha=0.04, lw=0.8)
    for xs, ys in ta: ax.plot(xs, ys, color='blue', alpha=0.04, lw=0.8)
    ax.plot([], [], 'r-', lw=2.5, label=f'Topspin 3500 RPM - {on_tp:.0f}% on table')
    ax.plot([], [], 'b-', lw=2.5, label=f'No spin 0 RPM - {on_ta:.0f}% on table')
    ax.set_xlim(0, 3.3); ax.set_ylim(0, 0.9)
    ax.set_xlabel("Horizontal Distance (m)"); ax.set_ylabel("Height above table (m)")
    ax.set_title("Monte-Carlo Trajectory Bundles: Topspin vs No Spin\n"
                 "(identical launch 12.5 m/s, 8\u00b0; club turbulence; 300 shots each)")
    ax.legend()
    save(fig, "fig7_stochastic_spread.png")

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.hist(lp, bins=40, alpha=0.6, color='red', label=f"Topspin (std={np.std(lp):.3f} m)")
    ax.hist(la, bins=40, alpha=0.6, color='blue', label=f"No spin (std={np.std(la):.3f} m)")
    ax.axvline(sim.net_x, color='green', ls='--', lw=2, label='Net')
    ax.axvline(sim.L, color='black', ls='--', lw=2, label='Far edge')
    ax.set_xlabel("Landing Position (m)"); ax.set_ylabel("Frequency")
    ax.set_title("Landing-Position Distribution: Topspin vs No Spin\n"
                 "(launch 12.5 m/s, 8\u00b0; club turbulence; 300 shots each)")
    ax.legend()
    save(fig, "fig8_landing_distribution.png")
print(f"[F7/8] topspin std {np.std(lp):.3f}, flat std {np.std(la):.3f}, "
      f"variance ratio {vr:.2f}x, on-table {on_tp:.0f}% vs {on_ta:.0f}%")

# ============================ FIGURE 9 ============================
spins4 = np.linspace(0, 6000, 40); spins4b = np.linspace(0, 3000, 40)
def curve(angle, h, grid):
    return [sim.max_safe_speed(angle, rpm, h=h, x0=0.15) or np.nan for rpm in grid]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5))
a1.plot(spins4, curve(7, 0.55, spins4), 'r-', lw=2, label='High contact point (0.55 m)')
a1.plot(spins4, curve(7, 0.20, spins4), 'b-', lw=2, label='Low contact point (0.20 m)')
a1.set_title("Steeper launch angle (7\u00b0)"); a1.set_xlabel("Topspin (RPM)")
a1.set_ylabel("Maximum safe launch speed (m/s)"); a1.legend(); a1.grid(alpha=0.3)
a2.plot(spins4b, curve(1, 0.55, spins4b), 'r-', lw=2, label='High contact point (0.55 m)')
a2.plot(spins4b, curve(1, 0.20, spins4b), 'b-', lw=2, label='Low contact point (0.20 m)')
a2.set_title("Near-flat launch angle (1\u00b0)"); a2.set_xlabel("Topspin (RPM)")
a2.set_ylabel("Maximum safe launch speed (m/s)"); a2.legend(); a2.grid(alpha=0.3)
fig.suptitle("Contact Height \u00d7 Topspin: a Lower Contact Point Benefits More From Spin\n"
             "(both panels vary topspin at a fixed launch angle; \u201chigh/low\u201d = contact height above table)",
             fontsize=12, fontweight='bold')
fig.tight_layout(); save(fig, "fig9_biomechanics.png")
print("[F9] biomechanics generated")

# ============================ FIGURE 10 ============================
rng = np.random.default_rng(42)
lam = np.linspace(0, 1.5, 30); ITER10 = 500
profiles = [("Flat Shot (0 RPM)", 'blue', 9, 5, 0),
            ("Topspin Shot (3500 RPM)", 'red', 12, 8, 3500)]
res10 = {}
for name, c, bv, ba, br in profiles:
    rates = []
    for e in lam:
        s = 0
        for _ in range(ITER10):
            v = bv + rng.normal(0, e*1.5); a = ba + rng.normal(0, e*1.5)
            r = max(0.0, br + rng.normal(0, e*400))
            if sim.simulate(v, a, r, h=0.15, x0=0.15)[2] == 'safe': s += 1
        rates.append(s/ITER10*100)
    res10[name] = rates
with plt.style.context('bmh'):
    fig, ax = plt.subplots(figsize=(10, 6))
    for name, c, *_ in profiles: ax.plot(lam, res10[name], label=name, color=c, marker='o', ms=4, lw=2.5)
    ax.axhline(50, color='gray', ls='--', alpha=0.6, label='50% threshold')
    ax.set_title("Execution-Error Tolerance: Probability of a Safe Landing\n"
                 "(500 Monte-Carlo shots per point; each shot integrated with RK4)")
    ax.set_xlabel("Execution-error scale \u03bb   "
                  "(per-shot noise: speed \u00b11.5\u03bb m/s, angle \u00b11.5\u03bb\u00b0, spin \u00b1400\u03bb RPM)")
    ax.set_ylabel("Probability of safe landing (%)")
    ax.set_ylim(-5, 105); ax.legend()
    save(fig, "fig10_human_error.png")
i = np.argmin(abs(lam - 0.5))
print(f"[F10] at lambda=0.5: topspin {res10['Topspin Shot (3500 RPM)'][i]:.0f}% "
      f"vs flat {res10['Flat Shot (0 RPM)'][i]:.0f}%")
print("\nALL FIGURES DONE.")
