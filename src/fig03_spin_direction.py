"""
Figure 3: Deterministic trajectories for a topspin loop, a flat drive and a
          backspin chop launched under identical conditions.

Launch: x0 = 0.0, launch_height = 1.15 m, 6.5 degrees, 6.0 m/s.
Uses engine_odeint, with positive rpm = topspin (Magnus force directed down).

Each curve is labelled by its physics:
    topspin  +3000  lands ~1.64 m, steepest descent
    no spin      0  lands ~1.91 m
    backspin -3000  lands ~2.33 m, shallowest descent
"""

import os
import matplotlib.pyplot as plt

from engine_odeint import simulate, descent_angle, draw_table

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          os.pardir, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

LAUNCH_V      = 6.0
LAUNCH_ANG    = 6.5
LAUNCH_HEIGHT = 1.15
SPIN          = 3000

runs = [
    ("Topspin",  +SPIN, 'red'),
    ("No Spin",       0, 'blue'),
    ("Backspin", -SPIN, 'green'),
]

fig, ax = plt.subplots(figsize=(12, 5))
draw_table(ax)

results = {}
for name, rpm, colour in runs:
    xs, zs, outcome = simulate(LAUNCH_V, LAUNCH_ANG, rpm,
                               launch_height=LAUNCH_HEIGHT)
    results[name] = (xs[-1], descent_angle(xs, zs))
    label = f"{name} | {rpm:+d} RPM [{outcome}]" if rpm else \
            f"{name} | 0 RPM [{outcome}]"
    ax.plot(xs, zs, color=colour, lw=2.5, label=label, zorder=4)

ax.set_xlim(0.5, 3.5)
ax.set_ylim(0.0, 2.0)
ax.set_xlabel("Horizontal Distance (m)", fontsize=12, fontweight='bold')
ax.set_ylabel("Height (m)", fontsize=12, fontweight='bold')
ax.set_title("Topspin vs No Spin vs Backspin\n"
             "How spin direction shapes the trajectory",
             fontsize=13, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'fig3_spin_direction.png'), dpi=300,
            bbox_inches='tight')
plt.close(fig)

x_no, d_no = results["No Spin"]
print(f"{'shot':10s} {'landing (m)':>12s} {'descent (deg)':>14s} "
      f"{'vs no-spin':>12s}")
for name, _, _ in runs:
    x, d = results[name]
    print(f"{name:10s} {x:12.2f} {d:14.2f} {d - d_no:+12.2f}")
