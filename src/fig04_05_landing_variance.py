"""
Figure 4: Standard deviation of landing position vs. spin, at three turbulence
          levels (sigma = 0.5, 2.0, 5.0 m s^-1/2).
Figure 5: Mean squared displacement (MSD) of landing position vs. spin.

Note on the two different filters:
  Figure 4 keeps only shots that land ON the table (outcome == 'safe'), which is
  why its curves start near 1,750-2,000 RPM - below that the shot overshoots at
  12 m/s and there are too few on-table landings to compute a std dev.
  Figure 5 keeps every landing that reaches the far side ('safe' or 'long'), so
  the MSD curve is defined across the whole 0-6,000 RPM range.
"""

import numpy as np
import matplotlib.pyplot as plt

from engine import TableTennisPhysics

import os

FIGURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       os.pardir, "figures")
os.makedirs(FIGURES, exist_ok=True)


eng = TableTennisPhysics()

spins             = np.linspace(0, 6000, 25)
turbulence_levels = [0.5, 2.0, 5.0]
iterations        = 150
fixed_v           = 12.0
fixed_angle       = 8.0

plt.style.use('bmh')

# ---------------------------------------------------------------- Figure 4
print("Running landing variance simulation...")
fig, ax = plt.subplots(figsize=(10, 5))
all_variances = {}

for turb in turbulence_levels:
    variances = []
    for rpm in spins:
        landings = []
        for _ in range(iterations):
            xs, ys, outcome = eng.simulate_trajectory(
                fixed_v, fixed_angle, rpm, turbulence_scale=turb)
            if outcome == 'safe':
                landings.append(xs[-1])
        variances.append(np.std(landings) if len(landings) > 5 else np.nan)
    all_variances[turb] = variances
    ax.plot(spins, variances, label=f'Turbulence \u03c3 = {turb}',
            marker='o', markersize=4, linewidth=2)

ax.set_title("Landing Position Std Dev vs Spin under Turbulence\n"
             f"({fixed_v} m/s, {fixed_angle}\u00b0 launch, "
             f"{iterations} runs each)", fontsize=13)
ax.set_xlabel("Spin (RPM)", fontsize=12)
ax.set_ylabel("Landing Position Standard Deviation (m)", fontsize=12)
ax.legend(fontsize=10)
fig.savefig(os.path.join(FIGURES, "fig4_landing_variance_vs_spin.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig)

print("\n--- Variance suppression by spin ---")
for turb in turbulence_levels:
    v = all_variances[turb]
    valid = [(spins[i], v[i]) for i in range(len(v)) if not np.isnan(v[i])]
    if len(valid) >= 2 and valid[0][1] > 0:
        (rpm0, std0), (rpm_hi, std_hi) = valid[0], valid[-1]
        print(f"  Turbulence {turb}: std dev {std0:.4f} m @ {rpm0:.0f} RPM "
              f"\u2192 {std_hi:.4f} m @ {rpm_hi:.0f} RPM "
              f"({(1 - std_hi/std0)*100:.1f}% reduction)")

# ---------------------------------------------------------------- Figure 5
print("\nGenerating MSD plot...")
msd_values = []
turb_fixed = 2.0
for rpm in spins:
    landings = []
    for _ in range(iterations):
        xs, ys, outcome = eng.simulate_trajectory(
            fixed_v, fixed_angle, rpm, turbulence_scale=turb_fixed)
        if outcome in ('safe', 'long'):
            landings.append(xs[-1])
    if len(landings) > 5:
        landings = np.array(landings)
        msd_values.append(np.mean((landings - landings.mean()) ** 2))
    else:
        msd_values.append(np.nan)

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(spins, msd_values, color='darkred', marker='s', markersize=4,
        linewidth=2, label=f'MSD (turbulence \u03c3={turb_fixed})')
ax.set_title("Mean Squared Displacement of Landing Position vs Spin\n"
             "(Suppression of Brownian Trajectory Diffusion)", fontsize=13)
ax.set_xlabel("Spin (RPM)", fontsize=12)
ax.set_ylabel("Mean Squared Displacement (m\u00b2)", fontsize=12)
ax.legend(fontsize=10)
fig.savefig(os.path.join(FIGURES, "fig5_msd_vs_spin.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

valid_msd = [(spins[i], msd_values[i]) for i in range(len(msd_values))
             if not np.isnan(msd_values[i])]
if len(valid_msd) >= 2 and valid_msd[0][1] > 0:
    (rpm0, msd0), (rpmN, msdN) = valid_msd[0], valid_msd[-1]
    print(f"\nMSD at {rpm0:.0f} RPM : {msd0:.6f} m\u00b2")
    print(f"MSD at {rpmN:.0f} RPM : {msdN:.6f} m\u00b2")
    print(f"MSD reduction       : {(1 - msdN/msd0)*100:.1f}%")
