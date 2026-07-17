"""
Figure 7: Monte Carlo trajectory bundles for a topspin loop and a flat drive
          launched under identical conditions (14 m/s, 8 degrees).
Figure 8: Histogram of the resulting landing positions.

This script uses its own local integrator rather than engine.TableTennisPhysics
because the published figures are drawn in the "table surface at y = 0.76 m"
frame with launch-condition scatter (dvx, dvy) applied on top of the
environmental noise term.
"""

import numpy as np
import matplotlib.pyplot as plt

import os

FIGURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       os.pardir, "figures")
os.makedirs(FIGURES, exist_ok=True)


TABLE_HEIGHT = 0.76
TABLE_LENGTH = 2.74
TABLE_START  = 1.0
NET_POSITION = TABLE_START + TABLE_LENGTH / 2
NET_HEIGHT   = 0.1525
K_SPIN_DECAY = 0.08

RHO = 1.225
R   = 0.02
M   = 0.0027
A   = np.pi * R ** 2
G   = 9.81
MU  = 1.81e-5
D   = 2 * R


def spin_parameter(omega_mag, v_mag):
    if v_mag < 1e-6 or omega_mag < 1e-9:
        return 0.0
    return (omega_mag * R) / v_mag


def CL_watts_ferrer(S):
    if S < 1e-9:
        return 0.0
    return 1.0 / (2.0 + 1.0 / S)


def CD_reynolds(v):
    if v < 1e-6:
        return 0.47
    Re = (RHO * v * D) / MU
    return 0.55 - 0.08 / (1.0 + np.exp(1.5e-4 * (Re - 28000)))


def decay_spin(omega, dt):
    return omega * np.exp(-K_SPIN_DECAY * dt)


def acceleration_2d(vx, vy, omega_scalar):
    """omega_scalar < 0 => topspin (Magnus force directed downwards)."""
    v_mag = np.sqrt(vx ** 2 + vy ** 2)
    if v_mag < 1e-8:
        return 0.0, -G
    S    = spin_parameter(abs(omega_scalar), v_mag)
    CL   = CL_watts_ferrer(S)
    sign = np.sign(omega_scalar)
    mag  = sign * (np.pi * RHO * R ** 2 * CL * v_mag) / (2 * M)
    drag = (np.pi * RHO * R ** 2 * CD_reynolds(v_mag) * v_mag) / (2 * M)
    ax = -vy * mag - vx * drag
    ay = vx * mag - G - vy * drag
    return ax, ay


def simulate_stochastic(spin_rpm, vx0, vy0, sigma_env=0.2, dt=0.001,
                        x0=TABLE_START, y0=TABLE_HEIGHT + 0.15):
    pos = np.array([x0, y0])
    vel = np.array([vx0, vy0])
    omega = -(spin_rpm * 2 * np.pi / 60)   # positive rpm == topspin
    traj = [pos.copy()]
    for _ in range(3000):
        omega = decay_spin(omega, dt)
        ax, ay = acceleration_2d(vel[0], vel[1], omega)
        noise = sigma_env * np.sqrt(dt) * np.random.normal(0, 1, 2)
        vel += np.array([ax, ay]) * dt + noise
        pos += vel * dt
        traj.append(pos.copy())
        if pos[1] < TABLE_HEIGHT:
            break
    return np.array(traj)


def monte_carlo(spin_rpm, vx0, vy0, runs=300, sigma=0.2):
    all_traj, landings = [], []
    for _ in range(runs):
        traj = simulate_stochastic(spin_rpm, vx0 + np.random.normal(0, 0.3),
                                   vy0 + np.random.normal(0, 0.1), sigma)
        all_traj.append(traj)
        landings.append(traj[-1, 0])
    return all_traj, np.array(landings)


speed_ms  = 14.0
angle_deg = 8.0
vx0 = speed_ms * np.cos(np.radians(angle_deg))
vy0 = speed_ms * np.sin(np.radians(angle_deg))

print("Running Monte Carlo simulations...")
trajs_pro, land_pro = monte_carlo(3500, vx0, vy0, runs=300)
trajs_am,  land_am  = monte_carlo(0,    vx0, vy0, runs=300)
print("Done.\n")

print(f"Topspin landing std dev : {np.std(land_pro):.4f} m")
print(f"Flat    landing std dev : {np.std(land_am):.4f} m")
print(f"Variance ratio (flat/topspin): "
      f"{np.var(land_am)/np.var(land_pro):.2f}x")

table_end = TABLE_START + TABLE_LENGTH
pro_on_table = sum(1 for t in trajs_pro if TABLE_START < t[-1, 0] < table_end)
am_on_table  = sum(1 for t in trajs_am  if TABLE_START < t[-1, 0] < table_end)
print(f"\nTopspin shots landing on table : {pro_on_table}/300 "
      f"({pro_on_table/3:.1f}%)")
print(f"Flat shots landing on table    : {am_on_table}/300 "
      f"({am_on_table/3:.1f}%)")

mean_pro = simulate_stochastic(3500, vx0, vy0, sigma_env=0)
mean_am  = simulate_stochastic(0,    vx0, vy0, sigma_env=0)

plt.style.use('bmh')

# ---------------------------------------------------------------- Figure 7
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot([TABLE_START, TABLE_START + TABLE_LENGTH],
        [TABLE_HEIGHT, TABLE_HEIGHT], color='darkorange', linewidth=6,
        zorder=3, label="Table")
ax.plot([NET_POSITION, NET_POSITION],
        [TABLE_HEIGHT, TABLE_HEIGHT + NET_HEIGHT], color='black',
        linewidth=3, zorder=3, label="Net")
for t in trajs_pro:
    ax.plot(t[:, 0], t[:, 1], color='red', alpha=0.04, linewidth=0.8)
for t in trajs_am:
    ax.plot(t[:, 0], t[:, 1], color='blue', alpha=0.04, linewidth=0.8)
ax.plot(mean_pro[:, 0], mean_pro[:, 1], 'r-', linewidth=2.5,
        label=f'Topspin 3500 RPM \u2014 {pro_on_table/3:.0f}% on table')
ax.plot(mean_am[:, 0], mean_am[:, 1], 'b--', linewidth=2.5,
        label=f'No Spin 0 RPM \u2014 {am_on_table/3:.0f}% on table')
ax.set_xlim(TABLE_START - 0.3, TABLE_START + TABLE_LENGTH + 0.3)
ax.set_ylim(TABLE_HEIGHT - 0.05, TABLE_HEIGHT + 0.55)
ax.set_xlabel("Horizontal Distance (m)", fontsize=12)
ax.set_ylabel("Height (m)", fontsize=12)
ax.set_title(f"Topspin vs No Spin \u2014 Same Launch Conditions "
             f"({speed_ms} m/s, {angle_deg}\u00b0)", fontsize=13)
ax.legend(fontsize=10)
fig.savefig(os.path.join(FIGURES, "fig7_stochastic_spread.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 8
fig, ax = plt.subplots(figsize=(9, 4))
ax.hist(land_pro, bins=40, alpha=0.6, color='red',
        label=f"Topspin 3500 RPM (\u03c3={np.std(land_pro):.3f} m)")
ax.hist(land_am, bins=40, alpha=0.6, color='blue',
        label=f"No Spin 0 RPM (\u03c3={np.std(land_am):.3f} m)")
ax.axvline(TABLE_START, color='green', linestyle='--', linewidth=2,
           label="Table start")
ax.axvline(TABLE_START + TABLE_LENGTH, color='black', linestyle='--',
           linewidth=2, label="Table end")
ax.set_xlabel("Landing Position (m)", fontsize=12)
ax.set_ylabel("Frequency", fontsize=12)
ax.set_title("Landing Distribution \u2014 Topspin vs No Spin", fontsize=13)
ax.legend(fontsize=10)
fig.savefig(os.path.join(FIGURES, "fig8_landing_distribution.png"), dpi=300,
            bbox_inches="tight")
plt.close(fig)
