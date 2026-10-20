"""
tt_engine.py  --  Unified table-tennis flight engine (single coordinate frame).

One engine, one frame, derivable forces. Deterministic trajectories use RK4;
stochastic (turbulent) trajectories use Euler-Maruyama, because a Wiener noise
term is only first-order integrable. Every equation below is stated so it can be
re-derived on paper.

COORDINATE FRAME
    x : horizontal, along the table.  Table surface at y = 0.
    Near edge x = 0, far edge x = L = 2.74 m, net at x = L/2 = 1.37 m.
    Ball is launched from (x0, h) with contact height h above the table.

SIGN CONVENTION
    rpm > 0  ->  topspin  ->  Magnus force points DOWN-and-back.
    rpm < 0  ->  backspin ->  Magnus force points UP-and-forward.

FORCES  (2D, per unit mass after dividing by m)
    Drag    F_D = -0.5 rho A C_D(Re) |v| v          (opposes velocity)
    Magnus  F_M = +0.5 rho A C_L(S) |v| * s * (vy, -vx)
            where s = sign(rpm). Derived from F = 0.5 rho A C_L v^2 n_hat with
            n_hat = (w x v)/|w x v|; for topspin (w = -|w| z_hat) this evaluates
            to (vy, -vx)/|v|, i.e. down-and-back. (See derivation in README.)
    Gravity (0, -g)

AERODYNAMIC COEFFICIENTS
    Re     = rho |v| d / mu
    C_D(Re)= 0.55 - 0.08 / (1 + exp(1.5e-4 (Re - 28000)))     ~0.47..0.55
    S      = |w| r / |v|                                       (spin parameter)
    C_L(S) = 1 / (2 + 1/S)                                     saturates -> 0.5
    w(t)   = w0 exp(-k_s t)                                    exponential decay
"""

import math
import numpy as np


class TableTennisSim:
    # -- physical constants (ITTF ball; ISA dry air 15 C sea level) -----------
    g   = 9.81            # m s^-2
    m   = 0.0027          # kg
    r   = 0.020           # m
    d   = 0.040           # m
    rho = 1.225           # kg m^-3
    mu  = 1.81e-5         # Pa s
    k_s = 0.08            # s^-1   spin-decay constant
    L   = 2.74            # m      table length
    net_x = 1.37          # m
    net_h = 0.1525        # m

    def __init__(self, rho=None):
        self.A = math.pi * self.r ** 2
        if rho is not None:            # allow rho=0 for the vacuum self-test
            self.rho = rho

    # -- coefficients ---------------------------------------------------------
    def C_D(self, v):
        if v < 1e-9:
            return 0.47
        Re = self.rho * v * self.d / self.mu
        return 0.55 - 0.08 / (1.0 + math.exp(1.5e-4 * (Re - 28000.0)))

    def C_L(self, S):
        if S < 1e-12:
            return 0.0
        return 1.0 / (2.0 + 1.0 / S)

    # -- equations of motion (scalar, fast) -----------------------------------
    def _accel(self, x, y, vx, vy, t, omega0, s):
        """Return (dx, dy, dvx, dvy) = time-derivative of the state."""
        v = math.hypot(vx, vy)
        if v < 1e-9:
            return vx, vy, 0.0, -self.g
        omega = omega0 * math.exp(-self.k_s * t)      # decayed spin magnitude
        S  = omega * self.r / v
        cd = self.C_D(v)
        cl = self.C_L(S)
        kD = 0.5 * self.rho * self.A * cd * v / self.m
        kM = 0.5 * self.rho * self.A * cl * v / self.m
        ax = -kD * vx + s * kM * vy
        ay = -kD * vy - s * kM * vx - self.g
        return vx, vy, ax, ay

    # keep the array form for external callers / clarity
    def _deriv(self, state, t, omega0, s):
        x, y, vx, vy = state
        dx, dy, ax, ay = self._accel(x, y, vx, vy, t, abs(omega0), s)
        return np.array([dx, dy, ax, ay])

    # -- integrators ----------------------------------------------------------
    def simulate(self, v0, angle_deg, rpm, h=0.15, x0=0.15,
                 sigma=0.0, rng=None, dt=0.001, t_max=3.0):
        """Return (xs, ys, outcome, x_land, descent_deg).

        sigma = 0 -> deterministic RK4.  sigma > 0 -> Euler-Maruyama with a
        velocity Wiener increment of intensity sigma (m s^-1/2).
        outcome in {'net', 'safe', 'long'}.
        """
        if rng is None:
            rng = np.random
        s      = 0.0 if rpm == 0 else math.copysign(1.0, rpm)
        omega0 = abs(rpm) * 2 * math.pi / 60.0
        th     = math.radians(angle_deg)
        x, y   = x0, h
        vx, vy = v0 * math.cos(th), v0 * math.sin(th)
        xs, ys = [x], [y]
        t = 0.0
        sqrt_dt = math.sqrt(dt)
        acc = self._accel
        n = int(t_max / dt)

        for _ in range(n):
            x_prev, y_prev = x, y
            if sigma <= 0:
                # classical RK4 (scalar)
                k1x, k1y, k1vx, k1vy = acc(x, y, vx, vy, t, omega0, s)
                k2x, k2y, k2vx, k2vy = acc(x+0.5*dt*k1x, y+0.5*dt*k1y,
                                           vx+0.5*dt*k1vx, vy+0.5*dt*k1vy, t+0.5*dt, omega0, s)
                k3x, k3y, k3vx, k3vy = acc(x+0.5*dt*k2x, y+0.5*dt*k2y,
                                           vx+0.5*dt*k2vx, vy+0.5*dt*k2vy, t+0.5*dt, omega0, s)
                k4x, k4y, k4vx, k4vy = acc(x+dt*k3x, y+dt*k3y,
                                           vx+dt*k3vx, vy+dt*k3vy, t+dt, omega0, s)
                x  += dt/6.0*(k1x + 2*k2x + 2*k3x + k4x)
                y  += dt/6.0*(k1y + 2*k2y + 2*k3y + k4y)
                vx += dt/6.0*(k1vx + 2*k2vx + 2*k3vx + k4vx)
                vy += dt/6.0*(k1vy + 2*k2vy + 2*k3vy + k4vy)
            else:
                # semi-implicit Euler-Maruyama
                _, _, ax, ay = acc(x, y, vx, vy, t, omega0, s)
                vx += ax*dt + sigma*sqrt_dt*rng.normal()
                vy += ay*dt + sigma*sqrt_dt*rng.normal()
                x  += vx*dt
                y  += vy*dt
            t += dt
            xs.append(x); ys.append(y)

            if x_prev < self.net_x <= x and y < self.net_h:
                return xs, ys, 'net', x, math.degrees(math.atan2(-vy, vx))
            if y_prev > 0.0 >= y:
                frac = y_prev / (y_prev - y)
                x_land = x_prev + frac * (x - x_prev)
                dsc = math.degrees(math.atan2(-vy, vx))
                if x_land <= self.net_x:
                    return xs, ys, 'net', x_land, dsc
                elif x_land <= self.L:
                    return xs, ys, 'safe', x_land, dsc
                else:
                    return xs, ys, 'long', x_land, dsc
            if x > 6.0:
                break
        return xs, ys, 'long', x, math.degrees(math.atan2(-vy, vx))

    def _descent(self, st):
        vx, vy = st[2], st[3]
        return math.degrees(math.atan2(-vy, vx))   # angle below horizontal

    def max_safe_speed(self, angle_deg, rpm, h=0.15, x0=0.15,
                       vmin=3.0, vmax=35.0, n=120):
        best = None
        for v in np.linspace(vmin, vmax, n):
            if self.simulate(v, angle_deg, rpm, h=h, x0=x0)[2] == 'safe':
                best = v
        return best
