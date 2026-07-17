# Topspin as an Aerodynamic Stabiliser

Simulation code for **"Topspin as an Aerodynamic Stabiliser: A Physics
Simulation and Probabilistic Analysis of Table Tennis Ball Flight"**
(SIMULTECH 2026, Paper #59, Porto, Portugal).

Every script in `src/` regenerates one or more figures from the camera-ready
paper. Nothing else is kept here: scripts that produced exploratory plots which
never made it into the paper have been removed (see
[Removed code](#removed-code)).

## Layout

```
src/         simulation code, one script per paper figure
figures/     generated output (PNG, 300 dpi)
```

## Requirements

```
python >= 3.9
numpy
matplotlib
```

```bash
pip install -r requirements.txt
```

## Running

Scripts write into `../figures/`, so run them from inside `src/`:

```bash
cd src
python fig01_02_phase_diagram.py
python fig03_spin_direction.py
python fig04_05_landing_variance.py
python fig06_shannon_entropy.py
python fig07_08_stochastic_spread.py
python fig09_biomechanics.py
python fig10_human_error.py
```

The Monte Carlo scripts (4/5, 6, 7/8) are unseeded, so results move by a few
tenths of a percent between runs. `fig10_human_error.py` is seeded
(`seed=42`) and is bit-for-bit reproducible.

## Scripts and figures

| Script | Paper figures | Runtime |
|---|---|---|
| `engine.py` | (physics core: Figs 1, 2, 4, 5, 6, 10) | — |
| `engine_odeint.py` | (physics core: Figs 3, 9) | — |
| `fig01_02_phase_diagram.py` | Fig 1 (phase diagram), Fig 2 (critical velocity) | ~1 min |
| `fig03_spin_direction.py` | Fig 3 (topspin / no spin / backspin) | seconds |
| `fig04_05_landing_variance.py` | Fig 4 (landing std dev), Fig 5 (MSD) | ~3 min |
| `fig06_shannon_entropy.py` | Fig 6 (Shannon entropy) | ~4 min |
| `fig07_08_stochastic_spread.py` | Fig 7 (trajectory bundles), Fig 8 (landing histogram) | ~1 min |
| `fig09_biomechanics.py` | Fig 9 (contact height) | ~5 min |
| `fig10_human_error.py` | Fig 10 (execution error tolerance) | ~2 min |

## Physics

Two engines implement the same physics with different integrators and frames:

| Engine | Integrator | Frame | Figures |
|---|---|---|---|
| `engine.py` | Euler, dt = 1 ms | table at `y = 0`, launch `x = 0.15` | 1, 2, 4, 5, 6, 10 |
| `engine_odeint.py` | `scipy.odeint` (LSODA, rtol 1e-8) | table at `z = 0.76`, launch `x = 0.0` | 3, 9 |

Both take **positive rpm = topspin** at their public interface. Shared model:

- **Magnus lift** — Watts–Ferrer saturated model, `CL = 1 / (2 + 1/S)`, with
  spin parameter `S = ωr / v`.
- **Drag** — Reynolds-dependent with a sigmoid drag-crisis transition,
  `CD = 0.55 − 0.08 / (1 + exp(1.5e-4 (Re − 28000)))`.
- **Spin decay** — exponential, `ω(t) = ω₀ exp(−k_s t)`, `k_s = 0.08 s⁻¹`.
- **Stochastic forcing** — Euler–Maruyama integration of additive Gaussian
  acceleration noise, framed via Fokker–Planck.

Ball and table constants are ITTF regulation: 40 mm diameter, 2.7 g,
2.74 m table, 15.25 cm net.

### Sign convention

`rpm` is passed **positive for topspin**. Internally the engine stores
`omega = -(rpm * 2*pi / 60)` so that positive `rpm` produces a Magnus force
directed **downwards**. Negative `rpm` is backspin (Magnus force up).

This convention is not universal across the original notebook cells, and at
least one of them had it inverted — see
[Corrections](#corrections-made-during-this-cleanup).

### Frames

Two coordinate frames appear in the paper, and the scripts follow whichever one
the published figure used:

- **Engine frame** (`engine.py`, Figs 1, 2, 4, 5, 6, 10) — table surface at
  `y = 0`, table spans `x ∈ [0, 2.74]`, net at `x = 1.37`, launch from
  `x = 0.15`.
- **Raised-table frame** (Figs 3, 7, 8, 9) — table surface at `y = 0.76 m`.
  Figs 7/8 place the table at `x ∈ [1.0, 3.74]`; Fig 3 uses `x ∈ [0, 2.74]`.

Figs 3, 7/8 and 9 therefore carry their own local integrators rather than
importing `engine.py`. Each explains why in its module docstring.

## Reproduction status

Verified against the camera-ready PDF:

| Figure | Status | Check |
|---|---|---|
| 1, 2 | ✅ exact | 9.28 / 12.85 / 14.77 m/s critical velocity |
| 3 | ✅ reproduces, ⚠️ labels swapped in paper | 1.64 / 1.91 / 2.33 m vs published 1.65 / 1.90 / 2.26 |
| 4, 5 | ✅ | MSD −81.8% (paper: −82%); curves start at 1750/2000 RPM as published |
| 6 | ✅ | 3.37 → 2.26 bits, −33.0% (paper: 3.43 → 2.26, −34.1%) |
| 7, 8 | ✅ | σ 0.079 vs 0.153, ratio 3.77×, 80% vs 0% on table |
| 9 | ✅ **exact** | 10.24 / 16.19 (steep, 6000 RPM); 12.81 / 27.92 (flat, 3000 RPM) |
| 10 | ✅ exact | 92.2 / 65.4 / 48.6 % — matches published curve pixel-for-pixel |

Differences on the unseeded scripts are Monte Carlo noise.

## Known discrepancies between code and paper

These are **discrepancies in the paper's captions and headline numbers**, not
in this code. The code here reproduces the published *figures*; in three places
the published *text* disagrees with the published *figure*.

### 1. Figure 10 — the 84% headline

The caption and abstract state **"At λ = 0.5: 84% vs 44%."**

Reading the published Figure 10 directly (two-point pixel calibration against
the 100% start and the λ=1.5 endpoint) gives **92.3% vs 44.7%** at λ = 0.517,
the nearest grid point. `fig10_human_error.py` with `seed=42` produces
**92.2% vs 43.4%**, matching the published curve exactly. The other two quoted
pairs also match exactly: λ = 1.0 → 65.4 / 30.2 (paper: 65 / 31); λ = 1.5 →
48.6 / 20.8 (paper: 49 / 21).

So the flat-drive number (44%) is right and the topspin number is not. The
figure shows ≈92%, not 84%. `83.1%` occurs at λ = 0.672, two grid points to
the right, which may be where 84 came from.

### 2. Figure 3 — the spin labels are swapped

Running `engine_odeint.py` at the published launch conditions (6.0 m/s, 6.5°,
launch height 1.15 m), using **its own documented convention** that positive
rpm = topspin:

| | this code | published `fig3b_backspin.pdf` |
|---|---|---|
| topspin +3000 | lands **1.64 m** | green curve, labelled *"Backspin \| −3000 RPM"* |
| no spin 0 | lands **1.91 m** | blue curve, labelled *"No Spin \| 0 RPM"* ✓ |
| backspin −3000 | lands **2.33 m** | red curve, labelled *"Topspin \| +3000 RPM"* |

The published curves are correct; **the two spin labels are interchanged**. The
independent tell: the red "Topspin" curve *rises above its launch height*, and
only an upward Magnus force — backspin — can do that. The no-spin curve matches
either way, which is why the error is easy to miss.

Likely cause: an earlier engine used the opposite convention (negative rpm =
topspin, as in the `odeint` `magnus2` cell), and the plotting cell carried that
sign over after the engine was rewritten.

### 3. The 8.64° descent-angle claim is unreachable

The caption and abstract quote an **8.64° descent-angle increase** for topspin
at 3,000 RPM. `engine_odeint.py` gives **+3.19°** at the Figure 3 launch
conditions, and a scan over launch speeds 5–16 m/s, angles 0–20° and four
contact heights **never exceeds +4.75°** at 3,000 RPM.

The reason is structural. The Watts-Ferrer lift coefficient
`CL = 1/(2 + 1/S)` **saturates at 0.5**. Holding the launch fixed and raising
the spin:

| Spin | Δ descent |
|---|---|
| 3,000 RPM | +4.71° |
| 6,000 RPM | +5.80° |
| 10,000 RPM | +6.48° |
| 20,000 RPM | +7.10° |
| 50,000 RPM | +7.41° |

It asymptotes below 8°, so **8.64° is not reachable at any spin rate** under
the physics the paper describes.

An earlier draft engine (`archive/engine_linear_lift.py`) used an *unsaturated*
linear lift, `CL = 1.2·(rω/v)`, with constant `CD = 0.4` and no spin decay.
That model reaches **+8.57°** at v0 = 16 m/s / 10°, and up to +19.3° across the
same scan. **8.64° is a fossil from that superseded model**, carried forward
into a paper whose stated physics cannot produce it.

### 4. Figure 4 caption

The caption states that at σ = 5.0, increasing spin from 2,000 to 5,000 RPM
"cuts landing spread by approximately half." Both the published figure and this
reproduction show the σ = 5.0 curve *rising* slightly over that range,
≈ 0.025 m → ≈ 0.030 m.

### 5. Figure 1 caption

The caption states "4,800 simulated shots." The published critical-velocity
thresholds (9.28 / 12.85 / 14.77 m/s) are reproduced **only** by a 40 × 40 =
**1,600**-shot grid; 2,400 / 3,200 / 4,800 / 6,000-shot grids all drift by
0.05–0.2 m/s.

## Corrections made during this cleanup

1. **`engine.py` outcome classification.** Balls landing short of the net were
   returned as `'long'`, which painted a spurious blue band along the bottom of
   the Figure 1 phase diagram. The published figure shows red there, so these
   are now returned as `'net'`. This does not affect any `'safe'` count, so
   Figs 2, 4, 5, 6 and 10 are unchanged.

2. **Inverted Magnus sign in the Fig 7/8 cell.** As originally written the
   Magnus force pointed *upward* for topspin, sending a 3,500 RPM shot to
   8.53 m — well off the table. With the sign corrected to match `engine.py`
   it lands at 3.68 m and the script reproduces the published figure (~3.62 m).

3. **Fig 9 rebuilt on `engine_odeint.py`.** The `player_height.png` cell used a
   standalone integrator launching from `x0 = -0.1` with tall/short colours
   inverted and the panels in the opposite order; it does not reproduce the
   published figure. The `odeint` version does, exactly.

4. **Fig 4 / Fig 5 filters.** Figure 4 keeps only `'safe'` landings (which is
   why it starts near 1,750–2,000 RPM); Figure 5 keeps `'safe'` or `'long'`.
   This asymmetry is what the published figures show and is now documented in
   the script.

## Removed code

The following cells produced plots that appear nowhere in the paper and have
been dropped:

- `stochastic_flight_cloud.png` — 2D Monte Carlo win-probability cloud at
  50 km/h.
- `monte_carlo_3d.png` — 3D trajectory stability under turbulence.
- `descent_angles.png` / `backspin_vs_nospin.png` — the `scipy.odeint` pair.
  Superseded by Figure 3; also did not reproduce the 8.64° / 6.45° numbers.
- `entropy_vs_std.png` — entropy vs standard deviation scatter (trailing block
  of the Shannon entropy cell).
- `player_height.png` — superseded first attempt at Figure 9 (see above).
- **Unused Figure 2 attempt** — a `max_safe_speed` sweep at 12° launch /
  1.15 m contact height on the `odeint` engine. It yields 7.21 / 8.85 / 9.43
  m/s (+31%), not the published 9.28 / 12.85 / 14.77 (+59%), so it was not the
  source of Figure 2 despite sharing its title. The phase-diagram sweep is.
- **Superseded human-error cell** — n=1000, 40 error steps, topspin baseline
  12 m/s / 12° / 3000 RPM. The published Figure 10 uses n=500, 30 steps and a
  3500 RPM baseline.
- **`engine_linear_lift.py`** — the pre-Watts-Ferrer engine (`CL = 1.2·(rω/v)`,
  `CD = 0.4`, no spin decay). Kept in `archive/` only because it documents the
  origin of the 8.64° number; it is not used by any figure.

## License

MIT - see [LICENSE](LICENSE).

## Citation

> Topspin as an Aerodynamic Stabiliser: A Physics Simulation and Probabilistic
> Analysis of Table Tennis Ball Flight. SIMULTECH 2026, Porto, Portugal.

## Acknowledgements

Technical mentorship: Prof. Ganesh Mani. Shot nomenclature and tactical
guidance contributed by a table tennis coach.

All scientific content, code, analysis and conclusions are the author's own.
Claude (Anthropic) was used for language assistance only.
