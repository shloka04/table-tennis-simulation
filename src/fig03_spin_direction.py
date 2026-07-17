"""
Figure 3: Deterministic trajectories for a topspin loop, a flat drive and a
          backspin chop launched under identical conditions.

Launch: x0 = 0.0, launch_height = 1.15 m, 6.5 degrees, 6.0 m/s.
Uses engine_odeint (positive rpm = topspin = Magnus force down).

Reproduces the published landing points, but see the WARNING below:

    this script          published fig3b_backspin.pdf
    ------------------   -----------------------------------------
    topspin  +3000  1.64  green curve, labelled "Backspin | -3000 RPM"
    no spin      0  1.91  blue curve,  labelled "No Spin | 0 RPM"
    backspin -3000  2.33  red curve,   labelled "Topspin | +3000 RPM"

------------------------------------------------------------------------------
WARNING - THE PUBLISHED FIGURE 3 HAS ITS SPIN LABELS SWAPPED
------------------------------------------------------------------------------
In fig3b_backspin.pdf the curve labelled "Topspin | +3000 RPM" rises above its
launch height, flies farthest (~2.26 m) and lands at the shallowest angle. Only
an upward Magnus force can lift a ball above its launch height, and an upward
Magnus force is backspin by definition. The curve labelled "Backspin | -3000
RPM" lands shortest (~1.65 m) and steepest, which is topspin behaviour.

Running THIS engine at the published launch conditions with its own documented
sign convention gives topspin = 1.64 m and backspin = 2.33 m - i.e. exactly the
published curves, with the two labels interchanged. The no-spin curve (1.91 m)
matches either way, which is why the error is easy to miss.

The most likely cause: an earlier engine used the opposite convention (negative
rpm = topspin), and the plotting cell carried that sign over after the engine
was rewritten.

This script labels each curve by its actual physics.

------------------------------------------------------------------------------
NOTE ON THE 8.64 DEGREE CLAIM
------------------------------------------------------------------------------
The caption and abstract quote an 8.64 degree descent-angle increase for
topspin at 3,000 RPM. This engine gives +3.2 degrees at the Figure 3 launch
conditions, and a scan over launch speeds 5-16 m/s, angles 0-20 degrees and
four contact heights never exceeds +4.75 degrees at 3,000 RPM.

The reason is structural: the Watts-Ferrer lift coefficient CL = 1/(2 + 1/S)
saturates at 0.5 as S grows. Holding the launch fixed and raising the spin
gives +4.71 deg at 3,000 RPM, +5.80 at 6,000, +6.48 at 10,000, +7.10 at 20,000
and +7.41 at 50,000 RPM. It asymptotes below 8 degrees, so 8.64 is not
reachable at ANY spin rate under this physics.

An earlier draft engine used an unsaturated linear lift, CL = 1.2*(r*omega/v),
with constant CD = 0.4 and no spin decay. That model reaches +8.57 deg at
v0 = 16 m/s / 10 deg and up to +19.3 deg across the same scan. The 8.64 figure
is almost certainly a leftover from that superseded model.
------------------------------------------------------------------------------
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
print("\nPublished figure shows the topspin and backspin labels swapped; "
      "see the module docstring.")
