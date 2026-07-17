"""
Figure 1: Phase diagram of spin rate vs. launch velocity (net / safe / long).
Figure 2: Critical velocity threshold (maximum safe launch velocity) vs. spin.

1,600 simulated shots = 40 spin values x 40 launch velocities.
This grid reproduces the published thresholds 9.28 / 12.85 / 14.77 m/s exactly.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm

from engine import TableTennisPhysics

import os

FIGURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       os.pardir, "figures")
os.makedirs(FIGURES, exist_ok=True)


eng = TableTennisPhysics()

LAUNCH_ANGLE = 8.0
spins      = np.linspace(0, 6000, 40)
velocities = np.linspace(3, 20, 40)
X, Y = np.meshgrid(spins, velocities)
Z = np.zeros_like(X)
critical_v = []

print(f"Running phase diagram sweep "
      f"({len(spins)} x {len(velocities)} = {len(spins)*len(velocities)} shots)...")
for i, rpm in enumerate(spins):
    max_safe_v = 0.0
    for j, v in enumerate(velocities):
        _, _, outcome = eng.simulate_trajectory(v, LAUNCH_ANGLE, rpm)
        if outcome == 'net':
            Z[j, i] = 0
        elif outcome == 'safe':
            Z[j, i] = 1
            max_safe_v = max(max_safe_v, v)
        else:  # 'long'
            Z[j, i] = 2
    critical_v.append(max_safe_v)
critical_v = np.array(critical_v)
print("Done.")

# ---------------------------------------------------------------- Figure 1
cmap = ListedColormap(['#e74c3c', '#2ecc71', '#3498db'])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

fig, ax = plt.subplots(figsize=(10, 6))
mesh = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=norm, shading='auto')
ax.legend(handles=[Patch(color='#e74c3c', label='Net collision'),
                   Patch(color='#2ecc71', label='Safe landing'),
                   Patch(color='#3498db', label='Long (overshoots)')],
          loc='upper left', fontsize=11, framealpha=0.85)
ax.set_title("Phase Diagram: Spin vs Launch Velocity\n"
             f"(launch angle {LAUNCH_ANGLE}\u00b0, player 0.5 m behind table, "
             "corrected Watts-Ferrer Magnus)", fontsize=13)
ax.set_xlabel("Spin (RPM)", fontsize=12, fontweight='bold')
ax.set_ylabel("Launch Velocity (m/s)", fontsize=12, fontweight='bold')
cbar = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2])
cbar.ax.set_yticklabels(['Net', 'Safe', 'Long'])
cbar.set_label('Outcome', fontsize=11)
fig.savefig(os.path.join(FIGURES, "fig1_phase_diagram.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
window = 5
padded = np.pad(critical_v, window // 2, mode='edge')
smoothed_v = np.convolve(padded, np.ones(window) / window,
                         mode='valid')[:len(critical_v)]

with plt.style.context('bmh'):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(spins, smoothed_v, color='purple', linewidth=2.5,
            label='Max safe velocity')
    ax.fill_between(spins, smoothed_v, color='purple', alpha=0.2)
    ax.set_title("Critical Velocity Threshold vs Spin", fontsize=13)
    ax.set_xlabel("Spin (RPM)", fontsize=12)
    ax.set_ylabel("Maximum Safe Launch Velocity (m/s)", fontsize=12)
    ax.set_ylim(0, None)
    ax.legend(fontsize=10, loc='upper left')
    fig.savefig(os.path.join(FIGURES, "fig2_critical_velocity.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)

v0    = smoothed_v[0]
v3000 = smoothed_v[np.argmin(abs(spins - 3000))]
v6000 = smoothed_v[-1]
print(f"Critical velocity at 0 RPM    : {v0:.2f} m/s")
print(f"Critical velocity at 3000 RPM : {v3000:.2f} m/s  "
      f"(+{(v3000/v0 - 1)*100:.0f}%)")
print(f"Critical velocity at 6000 RPM : {v6000:.2f} m/s  "
      f"(+{(v6000/v0 - 1)*100:.0f}%)")
