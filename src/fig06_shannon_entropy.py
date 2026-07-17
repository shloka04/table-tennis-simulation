"""
Figure 6: Shannon entropy of the landing distribution vs. spin rate (left) and
          landing histograms at three representative spin rates (right).

The three spin rates shown on the right are the grid points nearest to
0 / 2,000 / 5,000 RPM, which on a 35-point grid over 0-6,000 RPM land on
0, 1,941 and 4,941 RPM - the values that appear in the published legend.
"""

import numpy as np
import matplotlib.pyplot as plt

from engine import TableTennisPhysics

import os

FIGURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       os.pardir, "figures")
os.makedirs(FIGURES, exist_ok=True)


eng = TableTennisPhysics()

spins      = np.linspace(0, 6000, 35)
iterations = 600
turbulence = 2.0
fixed_v    = 12.0
fixed_ang  = 8.0

print("Computing Shannon entropy vs spin rate...")
print(f"({iterations} runs per spin rate, turbulence \u03c3={turbulence})\n")

entropies, std_devs = [], []
all_landings = {}

for rpm in spins:
    landings = []
    for _ in range(iterations):
        xs, ys, outcome = eng.simulate_trajectory(
            fixed_v, fixed_ang, rpm, turbulence_scale=turbulence)
        if outcome != 'net':
            landings.append(xs[-1])

    if len(landings) < 20:
        entropies.append(np.nan)
        std_devs.append(np.nan)
        continue

    landings = np.array(landings)
    all_landings[rpm] = landings

    # Entropy of the *centred* distribution: measures spread, not mean shift.
    centred   = landings - np.mean(landings)
    bin_edges = np.linspace(-0.15, 0.15, 30)
    counts, _ = np.histogram(centred, bins=bin_edges)
    probs     = counts / counts.sum()
    probs     = probs[probs > 0]
    entropies.append(-np.sum(probs * np.log2(probs)) if len(probs) else 0.0)
    std_devs.append(np.std(landings))

entropies = np.array(entropies)
std_devs  = np.array(std_devs)

print("--- Shannon entropy at key spin rates ---")
for target_rpm in [0, 1000, 2000, 3000, 4000, 6000]:
    idx = np.argmin(abs(spins - target_rpm))
    if not np.isnan(entropies[idx]):
        print(f"  {target_rpm:4d} RPM | H = {entropies[idx]:.4f} bits | "
              f"\u03c3 = {std_devs[idx]:.4f} m")

valid = [(spins[i], entropies[i]) for i in range(len(entropies))
         if not np.isnan(entropies[i])]
H0, HN = valid[0][1], valid[-1][1]
print(f"\nEntropy at 0 RPM   : {H0:.4f} bits")
print(f"Entropy at 6000 RPM: {HN:.4f} bits")
print(f"Entropy reduction  : {(1 - HN/H0)*100:.1f}%")

# ---------------------------------------------------------------- Figure 6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(spins, entropies, color='darkblue', linewidth=2.5,
         marker='o', markersize=4)
ax1.fill_between(spins, entropies, alpha=0.12, color='blue')
ax1.annotate(f"H = {H0:.2f} bits\n(0 RPM)", xy=(spins[0], H0),
             xytext=(500, H0 - 0.3), fontsize=9, color='darkblue',
             arrowprops=dict(arrowstyle='->', color='darkblue'))
ax1.annotate(f"H = {HN:.2f} bits\n(6000 RPM)", xy=(spins[-1], HN),
             xytext=(4000, HN + 0.3), fontsize=9, color='darkblue',
             arrowprops=dict(arrowstyle='->', color='darkblue'))
ax1.set_title("Shannon Entropy of Landing Distribution vs Spin\n"
              "Topspin as an Entropy-Reducing Mechanism", fontsize=12)
ax1.set_xlabel("Spin (RPM)", fontsize=11)
ax1.set_ylabel("Shannon Entropy H (bits)", fontsize=11)
ax1.grid(True, alpha=0.4)

for rpm, color in zip([0, 2000, 5000], ['#e74c3c', '#f39c12', '#2ecc71']):
    idx = np.argmin(abs(spins - rpm))
    actual_rpm = spins[idx]
    if actual_rpm in all_landings:
        ax2.hist(all_landings[actual_rpm], bins=30, alpha=0.5, color=color,
                 density=True,
                 label=f"{int(actual_rpm)} RPM (H={entropies[idx]:.2f} bits)")
ax2.axvline(1.0,  color='green', linestyle='--', linewidth=1.5,
            label='Table start')
ax2.axvline(3.74, color='black', linestyle='--', linewidth=1.5,
            label='Table end')
ax2.set_title("Landing Distributions at Three Spin Rates\n"
              "Narrowing distribution = decreasing entropy", fontsize=12)
ax2.set_xlabel("Landing Position (m)", fontsize=11)
ax2.set_ylabel("Probability Density", fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.4)

fig.suptitle("Information-Theoretic Analysis: Spin as Entropy Suppression",
             fontsize=13, fontweight='bold')
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "fig6_shannon_entropy.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig)
