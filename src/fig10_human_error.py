"""
Figure 10: Shot success rate vs. execution error magnitude lambda, comparing a
           topspin loop (12 m/s, 8 deg, 3,500 RPM) with a flat drive
           (9 m/s, 5 deg, 0 RPM).

The error magnitude lambda scales the standard deviation of three independent
Gaussian perturbations applied to the launch condition: speed, launch angle and
spin rate. lambda = 0 is a perfect stroke; lambda = 1.5 is a very sloppy one.
"""

import numpy as np
import matplotlib.pyplot as plt

from engine import TableTennisPhysics

import os

FIGURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       os.pardir, "figures")
os.makedirs(FIGURES, exist_ok=True)


eng = TableTennisPhysics()

# Both baselines must land safely with zero error, or the sweep is meaningless.
_, _, o1 = eng.simulate_trajectory(9.0,  5.0, 0)
_, _, o2 = eng.simulate_trajectory(12.0, 8.0, 3500)
print(f"Flat baseline    : {o1}   (must be 'safe')")
print(f"Topspin baseline : {o2}   (must be 'safe')")
assert o1 == 'safe' and o2 == 'safe', "baseline shots must land safely"

ITERATIONS     = 500
N_ERROR_STEPS  = 30
MAX_ERROR      = 1.5
VELOCITY_NOISE = 1.5
ANGLE_NOISE    = 1.5
SPIN_NOISE     = 400.0

PROFILES = [
    {"label": "Flat Shot (0 RPM)", "color": "darkorange",
     "base_v": 9.0, "base_ang": 5.0, "base_rpm": 0},
    {"label": "Topspin Shot (3500 RPM)", "color": "royalblue",
     "base_v": 12.0, "base_ang": 8.0, "base_rpm": 3500},
]

error_magnitudes = np.linspace(0, MAX_ERROR, N_ERROR_STEPS)
rng = np.random.default_rng(seed=42)
results = {}

print("\nRunning simulation...")
for profile in PROFILES:
    print(f"  {profile['label']} ...")
    rates = []
    for err in error_magnitudes:
        safe = 0
        for _ in range(ITERATIONS):
            v   = profile["base_v"] + rng.normal(0, err * VELOCITY_NOISE)
            ang = profile["base_ang"] + rng.normal(0, err * ANGLE_NOISE)
            rpm = max(0.0, profile["base_rpm"]
                      + rng.normal(0, err * SPIN_NOISE))
            _, _, outcome = eng.simulate_trajectory(v, ang, rpm)
            if outcome == 'safe':
                safe += 1
        rates.append(safe / ITERATIONS * 100)
    results[profile["label"]] = rates

plt.style.use('bmh')
fig, ax = plt.subplots(figsize=(10, 6))
for profile in PROFILES:
    ax.plot(error_magnitudes, results[profile["label"]],
            label=profile["label"], color=profile["color"],
            marker='o', markersize=4, linewidth=2.5)
ax.axhline(50, color='gray', linewidth=1.0, linestyle='--', alpha=0.6,
           label='50% threshold')
ax.set_title("Human Error Tolerance: Probability of Safe Landing\n"
             "(Monte Carlo \u00b7 n=500 per point, "
             "corrected Watts-Ferrer physics)", fontsize=13)
ax.set_xlabel("Error Magnitude  (0 = perfect, 1.5 = high clumsiness)",
              fontsize=12)
ax.set_ylabel("Probability of Safe Landing (%)", fontsize=12)
ax.set_ylim(-5, 105)
ax.legend(fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "fig10_human_error.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

print("\n--- Key results ---")
for target in [0.0, 0.5, 1.0, 1.5]:
    i = int(np.argmin(abs(error_magnitudes - target)))
    lam = error_magnitudes[i]
    top = results["Topspin Shot (3500 RPM)"][i]
    flat = results["Flat Shot (0 RPM)"][i]
    print(f"  \u03bb \u2248 {target:.1f} (grid {lam:.3f}) | "
          f"topspin {top:5.1f}%  vs  flat {flat:5.1f}%")
