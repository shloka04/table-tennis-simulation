"""
SUPERSEDED - NOT USED BY ANY FIGURE. Kept for provenance only.

The pre-Watts-Ferrer engine: unsaturated linear lift CL = 1.2*(r*omega/v),
constant CD = 0.4, no spin decay.

This is retained solely because it documents the origin of the 8.64 degree
descent-angle figure quoted in the paper. Because its lift is linear and
unbounded, it reaches +8.57 deg at v0 = 16 m/s / 10 deg and up to +19.3 deg
across a launch scan. The Watts-Ferrer engine actually described in the paper
saturates below +7.5 deg at any spin rate and cannot produce 8.64.

Do not use this for new work.
"""

import numpy as np


class TableTennisPhysics:
    def __init__(self):
        self.g = 9.81
        self.m = 0.0027
        self.r = 0.02
        self.A = np.pi * self.r ** 2
        self.rho = 1.225
        self.table_length = 2.74
        self.net_x = 1.37
        self.net_h = 0.1525

    def simulate_trajectory(self, v0, theta0_deg, rpm,
                            turbulence_scale=0.0, y0=0.15):
        theta = np.radians(theta0_deg)
        vx = v0 * np.cos(theta)
        vy = v0 * np.sin(theta)
        x, y = 0.15, y0
        dt = 0.001
        omega = rpm * (2 * np.pi / 60)
        xs, ys = [x], [y]

        while y >= 0 and x < 5.0:
            v = np.sqrt(vx ** 2 + vy ** 2)
            if v == 0:
                break
            if turbulence_scale > 0:
                vx += np.random.normal(0, turbulence_scale * dt)
                vy += np.random.normal(0, turbulence_scale * dt)
                v = np.sqrt(vx ** 2 + vy ** 2)

            Cd = 0.4
            Cl = 1.2 * (self.r * omega / v) if v > 0 else 0
            Fd = 0.5 * self.rho * self.A * Cd * v ** 2
            Fm = 0.5 * self.rho * self.A * Cl * v ** 2

            ax = (-Fd * (vx / v) + Fm * (vy / v)) / self.m
            ay = (-self.g * self.m - Fd * (vy / v) - Fm * (vx / v)) / self.m

            vx += ax * dt
            vy += ay * dt
            x += vx * dt
            y += vy * dt
            xs.append(x)
            ys.append(y)

            if len(xs) > 1 and xs[-2] < self.net_x <= x:
                if y < self.net_h:
                    return xs, ys, 'net'

        if self.net_x < x <= self.table_length:
            return xs, ys, 'safe'
        return xs, ys, 'long'
