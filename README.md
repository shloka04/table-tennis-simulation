# Topspin as an Aerodynamic Stabiliser

Coaches tell you to loop the ball rather than hit it flat. The usual explanation is that topspin brings the ball down so you can hit harder, which is true and incomplete.

What interested me was what topspin does once you accept that nobody executes a shot perfectly. Bat speed varies. Contact angle varies. The air is not still. If you model a shot as a single trajectory you cannot see any of that — so I modelled it as a distribution instead, and topspin stopped looking like a way to hit harder and started looking like a stabiliser.

It shrinks how far your mistakes move the ball.

Simulation code for *"Topspin as an Aerodynamic Stabiliser: A Physics Simulation and Probabilistic Analysis of Table Tennis Ball Flight"*.

An earlier version of the engine was published at **SIMULTECH 2026** in Porto — *Proceedings of the 16th International Conference on Simulation and Modeling Methodologies, Technologies and Applications*, pages 114–125, [doi:10.5220/0015182700004094](https://doi.org/10.5220/0015182700004094). Both that paper and the extended version are in [`paper/`](paper/), and every figure and number below regenerates from the scripts here.

The figures quoted throughout are from the extended version. The published SIMULTECH paper reports the earlier two-integrator engine, where the same effects came out at 8.64° of descent steepening, +59% safe speed and 82% displacement reduction.

## The model

Two layers.

A **deterministic** layer integrates the flight with drag, Magnus lift and spin decay, using RK4 at a 1 ms timestep. One coordinate frame, every force derivable on paper.

- **Drag** — Reynolds-dependent with a sigmoid drag-crisis transition, `C_D = 0.55 − 0.08 / (1 + exp(1.5e−4 (Re − 28000)))`
- **Magnus lift** — the Watts–Ferrer saturating model, `C_L = 1 / (2 + 1/S)` with spin parameter `S = ωr / v`
- **Spin decay** — exponential, `ω(t) = ω₀ exp(−k_s t)`, `k_s = 0.08 s⁻¹`

A **stochastic** layer adds Gaussian air turbulence and player execution error, turning every shot into a distribution of landing points. These runs use Euler–Maruyama rather than RK4, because a Wiener noise term is only first-order integrable — using RK4 there would be wrong, not more accurate.

Ball and table constants follow ITTF regulation: 40 mm diameter, 2.7 g, 2.74 m table, 15.25 cm net.

## What the simulations show

- Descent angle steepens by **3.2°**. Maximum safe launch speed rises from 9.18 to 14.83 m/s — **61%** at 6,000 RPM, and already 41% by 3,000.
- Mean squared landing displacement drops **83%**. Shannon entropy of the landing distribution falls **33%**, from 3.31 to 2.21 bits.
- Under moderate execution error a topspin loop lands safely **91%** of the time. A flat drive manages **45%** from the identical error distribution.
- At 12.5 m/s into club-level turbulence, 95% of topspin shots land on the table and essentially none of the flat ones do — a flat drive at that pace simply sails long. Landing spread is 3.5 times tighter.
- These are not separate effects. Topspin contracts the sensitivity of landing position to launch error by **62%**, and everything else follows from that.
- More spin is not always better. The optimum is interior and rises with pace: about 1,600 RPM at 10 m/s, climbing to the search ceiling by 13.5 m/s. Over-spinning a slow ball makes it land short.
- The 3D extension reproduces curved sidespin serves, deflecting the ball 20 to 27 cm across the table.

## Does it match real balls?

Everything above is simulated, so the real test is whether the same engine predicts shots it has never seen.

I used the public [SpinDOE dataset](https://github.com/cogsys-tuebingen/spindoe) (Gossard et al., 2023), which records ball trajectories alongside spin measured accurately by a printed-dot method. From it I kept the 133 trajectories with a measured spin and a long enough tracked arc.

The test is deliberately strict. For each shot the launch speed and angle come from the first few frames only, the spin is taken straight from the dataset, and the model predicts the **entire remaining flight** with nothing fitted to the later points. No parameter was tuned to these trajectories.

Across all 133 shots, spanning −2,900 to +4,700 RPM, the median error is **37.8 mm** over the whole arc, with a median R² of **0.994**. On a typical shot the predicted path stays within about four centimetres of the real ball across a flight of two to three metres, with every physical parameter fixed in advance.

## Running it

```bash
pip install -r requirements.txt
python run_all.py
```

That writes every figure except the real-data one into `figures/`. That last one needs the SpinDOE dataset, which isn't redistributed here:

```bash
git clone --depth 1 https://github.com/cogsys-tuebingen/spindoe.git
cd spindoe/data/trajectories && unzip -o trajectories.zip -d _t && cd ../../..
python convert_spindoe.py spindoe/data/trajectories/_t
python make_fig16.py
```

`run_all.py` picks the dataset up automatically if it is already present.

You can also point the validator at your own tracked shots. It needs only ball position over time, and is flexible about column names, units and frame rates:

```bash
python validate_real_data.py            # self-check on synthetic data, no files needed
python validate_real_data.py myshot.csv --framecol frame --fps 240 --units mm
```

Without a measured spin it reports descriptive and hold-out fits only, and says so, rather than claiming a prediction it cannot support.

## Files

About 1,000 lines of Python producing 17 figures.

| Script | Produces |
|---|---|
| `tt_engine.py` | the 2D engine; nearly every result in the paper comes from this one |
| `tt_engine3d.py` | 3D sidespin extension; reduces exactly to the 2D engine when sidespin is zero |
| `figures.py` | paper figures 1–10 |
| `extras.py` | paper figures 11 and 12 |
| `fig13_sidespin.py` | paper figure 13, plus a supplementary saturation curve |
| `fig14_optimal.py` | paper figure 14 |
| `fig15_validation.py` | parameter-recovery check on synthetic data; a method test, not a figure in the paper |
| `convert_spindoe.py` | turns SpinDOE trajectories into the validator's format |
| `validate_real_data.py` | validates the model against tracked shots from any source |
| `make_fig16.py` | **paper figure 15**, the real-data validation |

Script names run one ahead of the paper from `fig15` onward, because the synthetic parameter-recovery check did not make the final version. `make_fig16.py` is the one that produces the paper's Figure 15.

The 3D engine is verified to reduce exactly to the 2D engine at zero sidespin — identical landing positions to five decimal places across topspin, backspin and zero spin, with no lateral drift.

## What it doesn't do

The core analysis is two-dimensional. The 3D engine handles sidespin but the main results do not use it, and the larger errors on the real dataset came mostly from shots with strong sidespin, which a model in a single vertical plane cannot fully capture.

Air density is constant through the flight and the bounce is excluded, so this covers the airborne phase only. The lift law saturates smoothly, so it will not reproduce the non-monotonic "lift crisis" reported for table tennis balls at intermediate spin (Miyazaki et al., 2017), nor the inverse Magnus force seen at very low spin. The spin-decay rate is a single representative value within the published range, not one fitted to this particular ball and surface.

The real-data test uses a marked ball from a robot dataset. It has not been tried on broadcast or phone footage of live play, where spin must be inferred from the trajectory rather than read from a printed dot.

## Reproducibility

The Monte Carlo scripts are seeded, so repeated runs reproduce the same figures. Figure 12 reseeds per sample size deliberately, so its confidence intervals reflect sampling error at each N rather than one shared draw.

The numbers above are what these scripts print today. Two have shifted slightly since the engine was unified from the earlier two-integrator version: landing-displacement reduction is now 83% where it was 82%, and the entropy drop is 33% where it was 34%. Everything else is unchanged.

## Citation

> Maruthi, S. (2026). Topspin as an Aerodynamic Stabiliser: A Physics Simulation and Probabilistic Analysis of Table Tennis Ball Flight. In *Proceedings of the 16th International Conference on Simulation and Modeling Methodologies, Technologies and Applications (SIMULTECH 2026)*, pages 114–125. SCITEPRESS. doi:[10.5220/0015182700004094](https://doi.org/10.5220/0015182700004094)

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

The published paper is redistributed here under CC BY-NC-ND 4.0, the licence it was released under. The code in this repository is MIT.

## Acknowledgements

Research mentorship from Prof. Ganesh Mani (Carnegie Mellon University). Shot nomenclature and tactical interpretation came from conversations with an experienced table tennis coach. The real-trajectory test uses the SpinDOE dataset (Gossard et al., 2023).

I used Claude (Anthropic) for code debugging, language polishing and grammar checking. All scientific content, simulation code, analysis, figures and conclusions are my own.

## Licence

MIT — see [LICENSE](LICENSE).
