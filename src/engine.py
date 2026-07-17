"""
Core deterministic/stochastic flight engine for the table tennis simulations.

Sign convention
---------------
`rpm` is passed as a POSITIVE number for topspin.  Internally the scalar
angular velocity is stored as ``omega = -(rpm * 2*pi / 60)`` so that a positive
`rpm` produces a Magnus force directed DOWNWARDS (the defining property of
topspin).  A negative `rpm` therefore represents backspin (Magnus force up).
"""

import numpy as np
import math


class TableTennisPhysics:
    def __init__(self):
        self.g            = 9.81       # m s^-2
        self.m            = 0.0027     # kg   (ITTF ball)
        self.r            = 0.02       # m    (ITTF 40 mm ball)
        self.A            = np.pi * self.r ** 2
        self.rho          = 1.225      # kg m^-3
        self.mu           = 1.81e-5    # Pa s
        self.D            = 2 * self.r
        self.table_length = 2.74       # m
        self.net_x        = 1.37       # m
        self.net_h        = 0.1525     # m
        self.k_decay      = 0.08       # s^-1  spin decay constant

    # ------------------------------------------------------------------
    def spin_parameter(self, omega, v):
        """Non-dimensional spin parameter S = omega*r / v."""
        if v < 1e-6 or omega < 1e-9:
            return 0.0
        return (omega * self.r) / v

    def CL(self, S):
        """Watts-Ferrer saturated lift coefficient."""
        if S < 1e-9:
            return 0.0
        return 1.0 / (2.0 + 1.0 / S)

    def CD(self, v):
        """Reynolds-dependent drag coefficient (sigmoid transition)."""
        if v < 1e-6:
            return 0.47
        Re = (self.rho * v * self.D) / self.mu
        return 0.55 - 0.08 / (1.0 + math.exp(1.5e-4 * (Re - 28000)))

    def decay_spin(self, omega, dt):
        return omega * math.exp(-self.k_decay * dt)

    # ------------------------------------------------------------------
    def simulate_trajectory(self, v0, theta0_deg, rpm,
                            turbulence_scale=0.0, y0=0.15):
        """Integrate one trajectory. Returns (xs, ys, outcome).

        outcome in {'net', 'safe', 'long'}.
        Coordinates: table surface at y = 0, table spans x in [0, 2.74],
        net at x = 1.37.
        """
        theta = np.radians(theta0_deg)
        vx    = v0 * math.cos(theta)
        vy    = v0 * math.sin(theta)
        x, y  = 0.15, y0
        dt    = 0.001
        omega = -(rpm * 2 * math.pi / 60)
        xs, ys = [x], [y]

        while y >= 0 and x < 5.0:
            v = math.sqrt(vx ** 2 + vy ** 2)
            if v < 1e-6:
                break

            omega     = self.decay_spin(omega, dt)
            omega_mag = abs(omega)
            S         = self.spin_parameter(omega_mag, v)
            cl        = self.CL(S)
            sign      = math.copysign(1.0, omega)
            Fm        = sign * 0.5 * self.rho * self.A * cl * v ** 2
            cd        = self.CD(v)
            Fd        = 0.5 * self.rho * self.A * cd * v ** 2

            ax = (-Fd * (vx / v) - Fm * (vy / v)) / self.m
            ay = (-self.g - Fd * (vy / v) / self.m
                  + Fm * (vx / v) / self.m)

            if turbulence_scale > 0:
                ax += np.random.normal(0, turbulence_scale)
                ay += np.random.normal(0, turbulence_scale)

            vx += ax * dt
            vy += ay * dt
            x  += vx * dt
            y  += vy * dt
            xs.append(x)
            ys.append(y)

            if len(xs) > 1 and xs[-2] < self.net_x <= x:
                if y < self.net_h:
                    return xs, ys, 'net'

        if x <= self.net_x:
            # Dropped short of the net (never reached the opponent's half).
            return xs, ys, 'net'
        elif x <= self.table_length:
            return xs, ys, 'safe'
        else:
            return xs, ys, 'long'
