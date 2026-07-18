# Topspin as an Aerodynamic Stabiliser

Simulation code for **"Topspin as an Aerodynamic Stabiliser: A Physics
Simulation and Probabilistic Analysis of Table Tennis Ball Flight"**
(SIMULTECH 2026, Paper #59, Porto, Portugal).

Every script in `src/` regenerates one or more figures from the paper.

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
scipy
```

```bash
pip install -r requirements.txt
```

## Running

Each script writes its output into `figures/` and can be run from anywhere:

```bash
python src/fig01_02_phase_diagram.py
python src/fig03_spin_direction.py
python src/fig04_05_landing_variance.py
python src/fig06_shannon_entropy.py
python src/fig07_08_stochastic_spread.py
python src/fig09_biomechanics.py
python src/fig10_human_error.py
```

## Scripts and figures

| Script | Paper figures |
|---|---|
| `engine.py` | shared physics core (Figs 1, 2, 4, 5, 6, 10) |
| `engine_odeint.py` | shared physics core (Figs 3, 9) |
| `fig01_02_phase_diagram.py` | Fig 1 (phase diagram), Fig 2 (critical velocity) |
| `fig03_spin_direction.py` | Fig 3 (topspin / no spin / backspin) |
| `fig04_05_landing_variance.py` | Fig 4 (landing std dev), Fig 5 (MSD) |
| `fig06_shannon_entropy.py` | Fig 6 (Shannon entropy) |
| `fig07_08_stochastic_spread.py` | Fig 7 (trajectory bundles), Fig 8 (landing histogram) |
| `fig09_biomechanics.py` | Fig 9 (contact height) |
| `fig10_human_error.py` | Fig 10 (execution error tolerance) |

## Physics

Two engines implement the model with different integrators and coordinate
frames:

| Engine | Integrator | Frame | Figures |
|---|---|---|---|
| `engine.py` | Euler, dt = 1 ms | table at `y = 0`, launch `x = 0.15` | 1, 2, 4, 5, 6, 10 |
| `engine_odeint.py` | `scipy.odeint` (LSODA) | table at `z = 0.76`, launch `x = 0.0` | 3, 9 |

Both take positive rpm to mean topspin. Shared model:

- **Magnus lift** - Watts-Ferrer saturated model, `CL = 1 / (2 + 1/S)`, with
  spin parameter `S = omega*r / v`.
- **Drag** - Reynolds-dependent with a sigmoid drag-crisis transition,
  `CD = 0.55 - 0.08 / (1 + exp(1.5e-4 (Re - 28000)))`.
- **Spin decay** - exponential, `omega(t) = omega_0 exp(-k_s t)`,
  `k_s = 0.08 s^-1`.
- **Stochastic forcing** - Euler-Maruyama integration of additive Gaussian
  acceleration noise, framed via Fokker-Planck.

Ball and table constants follow ITTF regulation: 40 mm diameter, 2.7 g,
2.74 m table, 15.25 cm net.

## Results

Headline values reproduced by the simulations:

- **59%** increase in maximum safe launch speed at 6,000 RPM (Figs 1, 2)
- **82%** reduction in mean squared landing displacement (Figs 4, 5)
- **34%** Shannon-entropy reduction of the landing distribution, 3.36 -> 2.21
  bits (Fig 6)
- **3.2x** variance ratio between the no-spin and topspin landing
  distributions (Figs 7, 8)
- Lower contact height amplifies the topspin speed benefit (Fig 9)
- Topspin sustains a far higher safe-landing rate under execution error than a
  flat drive (Fig 10)

The stochastic scripts (Figs 4/5, 6, 7/8) are unseeded, so Monte Carlo results
move by a few tenths of a percent between runs. `fig10_human_error.py` is
seeded and fully reproducible.

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
