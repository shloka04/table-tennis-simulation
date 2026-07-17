"""
Figure 9: Maximum safe launch speed vs. spin rate for two contact-height
          scenarios, h = 0.55 m (tall, red) and h = 0.20 m (short, blue).
          Left  panel: topspin loop at 7 degrees launch angle, 0-6,000 RPM.
          Right panel: flat drive at 1 degree launch angle, 0-3,000 RPM.

Reproduces the published values exactly:
  steep 7 deg :    0 RPM -> 7.34 / 9.43     6,000 RPM -> 10.24 / 16.19
  flat  1 deg :    0 RPM -> 8.95 / 14.10    3,000 RPM -> 12.81 / 27.92

Note the non-uniform spin grids: both panels are sampled more densely across
the middle of their range, where the curves bend most.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from engine_odeint import TABLE_HEIGHT, max_safe_speed

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          os.pardir, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

heights = {
    'Tall player (h=0.55 m)':  TABLE_HEIGHT + 0.55,
    'Short player (h=0.20 m)': TABLE_HEIGHT + 0.20,
}
colors = {
    'Tall player (h=0.55 m)':  '#e74c3c',
    'Short player (h=0.20 m)': '#3498db',
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Biomechanics: Contact Height and Spin Interaction\n'
             'Lower contact height amplifies topspin benefit',
             fontweight='bold', fontsize=13)

panels = [
    (axes[0], 7.0, 6000,
     np.concatenate([np.linspace(0, 1500, 20),
                     np.linspace(1500, 2500, 40),
                     np.linspace(2500, 6000, 40)]),
     'Steep topspin shots (7\u00b0 launch angle)\n'
     'Low contact height gains more from topspin'),
    (axes[1], 1.0, 3000,
     np.concatenate([np.linspace(0, 1500, 20),
                     np.linspace(1500, 2500, 60),
                     np.linspace(2500, 3000, 20)]),
     'Flat shots (1\u00b0 launch angle, 0\u20133000 RPM)\n'
     'Effect of contact height at shallow angles'),
]

for ax, angle, xmax, spins, title in panels:
    print(f"Computing {angle}\u00b0 panel ({len(spins)} spin values)...")
    ax.set_title(title)
    vals = {}
    for h_label, h_abs in heights.items():
        ms = [max_safe_speed(angle, omega, h_abs) for omega in spins]
        vals[h_label] = np.array([v if v is not None else np.nan
                                  for v in ms], dtype=float)
        ax.plot(spins, vals[h_label], color=colors[h_label], lw=2,
                label=h_label)

    tall  = vals['Tall player (h=0.55 m)']
    short = vals['Short player (h=0.20 m)']
    valid = ~np.isnan(tall) & ~np.isnan(short)
    lower = np.where(valid, np.nanmin([tall, short], axis=0), np.nan)
    ax.fill_between(spins, lower, np.where(valid, tall, np.nan),
                    where=(tall >= short) & valid,
                    alpha=0.15, color='#e74c3c', label='Tall advantage')
    ax.fill_between(spins, lower, np.where(valid, short, np.nan),
                    where=(short > tall) & valid,
                    alpha=0.15, color='#3498db', label='Short advantage')

    ax.set_xlabel('Spin (RPM)', fontweight='bold')
    ax.set_ylabel('Max Safe Launch Velocity (m/s)', fontweight='bold')
    ax.set_xlim(0, xmax)
    ax.legend(fontsize=9)

    for h_label in heights:
        v_lo, v_hi = vals[h_label][0], vals[h_label][-1]
        if not (np.isnan(v_lo) or np.isnan(v_hi)):
            print(f"  {h_label} ({angle}\u00b0): {v_lo:.2f} \u2192 {v_hi:.2f} m/s "
                  f"(+{v_hi - v_lo:.2f})")

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'fig9_biomechanics.png'), dpi=300,
            bbox_inches='tight')
plt.close(fig)
print("Done.")
