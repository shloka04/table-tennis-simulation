"""
engine_odeint.py
================
Deterministic + stochastic trajectory engine used for Figures 3 and 9.

This is a separate engine from `engine.py`. Both implement the same physics
(Watts-Ferrer saturated Magnus lift, Reynolds-dependent drag, exponential spin
decay) but differ in frame and integrator:

  engine.py         Euler integration, table surface at y = 0, launch from
                    x = 0.15. Used for Figures 1, 2, 4, 5, 6, 10.
  engine_odeint.py  scipy.odeint (LSODA, rtol=1e-8), table surface at
                    z = 0.76 m, launch from x = 0.0. Used for Figures 3, 9.

Sign convention
---------------
  omega0 positive = topspin  (Magnus force pushes the ball DOWN)
  omega0 negative = backspin (Magnus force pushes the ball UP)

This is the OPPOSITE of the internal storage in engine.py, which negates the
rpm on the way in and so reaches the same external convention (positive rpm
= topspin) by a different route. Both modules take positive rpm to mean
topspin at their public interface.

Outcome strings
---------------
  'net'   hits or clips the net, or lands on the server's own half
  'table' lands on the opponent's half
  'long'  overshoots the far edge, or times out

Note `classify_shot` returns 'safe' rather than 'table' for an on-table
landing, matching the outcome vocabulary of engine.py.
"""

import numpy as np
from scipy.integrate import odeint

# ---------------------------------------------------------------- constants
TABLE_LENGTH = 2.74
TABLE_HEIGHT = 0.76
NET_HEIGHT   = 0.1525
NET_X        = TABLE_LENGTH / 2      # 1.37 m

M   = 0.0027          # ball mass, kg
R   = 0.020           # ball radius, m
D   = 0.040           # ball diameter, m
A   = np.pi * R ** 2  # cross-sectional area, m^2
RHO = 1.225           # air density, kg/m^3
MU  = 1.81e-5         # dynamic viscosity, Pa s
G   = 9.81
KS  = 0.08            # spin decay constant, s^-1


# ------------------------------------------------------------- aerodynamics
def drag_coeff(v):
    """Reynolds-dependent drag with a sigmoid drag-crisis transition."""
    Re = RHO * abs(v) * D / MU
    return 0.55 - 0.08 / (1.0 + np.exp(1.5e-4 * (Re - 28000)))


def lift_coeff(S):
    """Watts-Ferrer saturated lift coefficient. Note this SATURATES at 0.5."""
    if S <= 0:
        return 0.0
    return 1.0 / (2.0 + 1.0 / S)


# --------------------------------------------------------- equations of motion
def equations_of_motion(state, t, omega0):
    x, z, vx, vz = state
    v_mag = np.sqrt(vx ** 2 + vz ** 2)

    omega = omega0 * np.exp(-KS * t)
    spin_sign = np.sign(omega) if abs(omega) > 1e-6 else 0.0

    CD = drag_coeff(v_mag)
    drag_f = 0.5 * RHO * A * CD * v_mag
    Fdx = -drag_f * vx
    Fdz = -drag_f * vz

    S = (abs(omega) * R / v_mag) if v_mag > 1e-6 else 0.0
    CL = lift_coeff(S)
    Fm = 0.5 * RHO * A * CL * v_mag ** 2
    if v_mag > 1e-6:
        Fmx = spin_sign * Fm * (vz / v_mag)
        Fmz = spin_sign * Fm * (-vx / v_mag)
    else:
        Fmx = Fmz = 0.0

    ax = (Fdx + Fmx) / M
    az = (Fdz + Fmz) / M - G
    return [vx, vz, ax, az]


# -------------------------------------------------------------- simulation
def simulate(v0, angle_deg, omega0_rpm, launch_height=1.15,
             dt=0.001, t_max=2.5):
    """Deterministic shot. Returns (xs, zs, outcome).

    `omega0_rpm` positive = topspin. `launch_height` is ABSOLUTE height above
    the floor, so a 0.20 m contact height is launch_height = 0.76 + 0.20.
    """
    omega0 = omega0_rpm * 2 * np.pi / 60.0
    angle = np.radians(angle_deg)
    t = np.arange(0, t_max, dt)
    state0 = [0.0, launch_height, v0 * np.cos(angle), v0 * np.sin(angle)]
    sol = odeint(equations_of_motion, state0, t, args=(omega0,),
                 rtol=1e-8, atol=1e-10)
    xs, zs = sol[:, 0], sol[:, 1]

    net_top = TABLE_HEIGHT + NET_HEIGHT
    for i in range(1, len(xs)):
        if xs[i] >= NET_X and xs[i - 1] < NET_X:
            frac = (NET_X - xs[i - 1]) / (xs[i] - xs[i - 1])
            z_at_net = zs[i - 1] + frac * (zs[i] - zs[i - 1])
            if z_at_net <= net_top:
                return xs[:i + 1], zs[:i + 1], 'net'
        if xs[i] > NET_X and zs[i] <= TABLE_HEIGHT < zs[i - 1]:
            return xs[:i + 1], zs[:i + 1], 'table'
        if xs[i] < NET_X and zs[i] <= TABLE_HEIGHT < zs[i - 1]:
            return xs[:i + 1], zs[:i + 1], 'net'
        if xs[i] < NET_X and zs[i] < TABLE_HEIGHT - 0.5:
            return xs[:i + 1], zs[:i + 1], 'net'
    return xs, zs, 'long'


def classify_shot(v0, angle_deg, omega_rpm, launch_height):
    """Returns 'net' | 'safe' | 'long', interpolating the landing point."""
    xs, zs, status = simulate(v0, angle_deg, omega_rpm, launch_height)
    if status == 'net':
        return 'net'
    for i in range(1, len(xs)):
        if xs[i] > NET_X and zs[i] <= TABLE_HEIGHT < zs[i - 1]:
            frac = (TABLE_HEIGHT - zs[i - 1]) / (zs[i] - zs[i - 1])
            x_land = xs[i - 1] + frac * (xs[i] - xs[i - 1])
            return 'safe' if NET_X < x_land <= TABLE_LENGTH else 'long'
    return 'long'


def max_safe_speed(angle_deg, omega_rpm, launch_height,
                   speed_min=3.0, speed_max=35.0, n=200):
    """Largest speed on the search grid that still lands on the table."""
    best = None
    for v in np.linspace(speed_min, speed_max, n):
        if classify_shot(v, angle_deg, omega_rpm, launch_height) == 'safe':
            best = v
    return best


def descent_angle(xs, zs):
    """Angle below horizontal at the final integration step, in degrees."""
    return np.degrees(np.arctan2(-(zs[-1] - zs[-2]), xs[-1] - xs[-2]))


def draw_table(ax):
    ax.plot([0, TABLE_LENGTH], [TABLE_HEIGHT, TABLE_HEIGHT],
            color='orange', lw=4, zorder=2, label='Table')
    ax.plot([NET_X, NET_X], [TABLE_HEIGHT, TABLE_HEIGHT + NET_HEIGHT],
            color='black', lw=2, zorder=3, label='Net')
