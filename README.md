# Topspin as an Aerodynamic Stabiliser

Simulation code for **"Topspin as an Aerodynamic Stabiliser: A Physics Simulation and Probabilistic Analysis of Table Tennis Ball Flight"**. An earlier baseline version of the engine was presented at SIMULTECH 2026 in Porto; this repository carries the extended work, with the unified engine, the three-dimensional sidespin extension and validation against recorded trajectories.

Both papers are in [`paper/`](paper/): the extended version, and the shorter SIMULTECH conference paper that preceded it. Every figure and number below can be regenerated from the scripts here.

## The question

Coaches tell you to loop the ball rather than hit it flat. The usual explanation is that topspin brings the ball down so you can hit harder, which is true but incomplete. What I wanted to know was what topspin does once you accept that nobody executes a shot perfectly. Your bat speed varies, your angle varies, your contact point varies, and the air is not still either.

So the model has two layers. The first integrates the flight deterministically with drag, Magnus lift and spin decay. The second adds Gaussian air turbulence and player execution error on top, so each shot becomes a distribution of possible landing points instead of one trajectory.

Looked at that way, topspin is less a trick for hitting harder than a stabiliser. It shrinks how far your mistakes move the ball.

## What the simulations show

- The descent angle steepens by **3.2°**. The maximum safe launch speed rises from 9.18 to 14.83 m/s, which is **61%** at 6,000 RPM and already 41% by 3,000.
- Mean squared landing displacement drops **83%**. The Shannon entropy of the landing distribution falls **33%**, from 3.31 to 2.21 bits.
- Under moderate execution error a topspin loop lands safely **91%** of the time. A flat drive manages **45%** from the same error distribution.
- Launched at 12.5 m/s into club-level turbulence, 95% of topspin shots land on the table and essentially none of the flat ones do, because a flat drive at that pace simply sails long. The landing spread is 3.5 times tighter.
- These are not separate effects. Topspin contracts the sensitivity of landing position to launch error by **62%**, and the rest follows from that.
- More spin is not always better. The optimum is interior and rises with pace: about 1,600 RPM at 10 m/s, climbing to the search ceiling by 13.5 m/s. Over-spinning a slow ball makes it land short.
- The 3D extension reproduces curved sidespin serves, deflecting the ball 20 to 27 cm across the table.

## Does it match real balls?

Everything above is simulated, so the real test is whether the same engine predicts shots it has never seen.

I used the public [SpinDOE dataset](https://github.com/cogsys-tuebingen/spindoe) (Gossard et al., 2023), which records ball trajectories alongside an accurately measured spin from a printed-dot method. From it I kept the 133 trajectories with a measured spin and a long enough tracked arc.

The test is deliberately strict. For each shot I read launch speed and angle from the first few frames only, take the spin straight from the dataset, and let the model predict the entire rest of the flight with nothing fitted to the later points. Across all 133 shots, spanning −2,900 to +4,700 RPM, the median error is **37.8 mm** over the whole arc, with a median R² of **0.994**.

On a typical shot, then, the predicted path stays within about four centimetres of the real ball across a flight of two to three metres, with every physical parameter fixed in advance.

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

`run_all.py` picks the dataset up automatically if it's already there.

You can also point the validator at your own tracked shots. It only needs ball position over time, and it's flexible about column names, units and frame rates:

```bash
python validate_real_data.py            # self-check on synthetic data, no files needed
python validate_real_data.py myshot.csv --framecol frame --fps 240 --units mm
```

Without a measured spin it reports descriptive and hold-out fits only, and says so, rather than claiming a prediction it can't support.

## Files

| Script | Produces |
|---|---|
| `tt_engine.py` | the 2D engine; nearly every result in the paper comes from this one |
| `tt_engine3d.py` | 3D extension with sidespin; reduces exactly to the 2D engine when sidespin is zero |
| `figures.py` | paper figures 1–10 |
| `extras.py` | paper figures 11 and 12 |
| `fig13_sidespin.py` | paper figure 13, plus a supplementary saturation curve |
| `fig14_optimal.py` | paper figure 14 |
| `fig15_validation.py` | parameter-recovery check on synthetic data; a method test, not a figure in the paper |
| `convert_spindoe.py` | turns SpinDOE trajectories into the validator's format |
| `validate_real_data.py` | validates the model against tracked shots from any source |
| `make_fig16.py` | **paper figure 15**, the real-data validation |

The script names run one ahead of the paper from `fig15` onward, because the synthetic parameter-recovery check didn't make the final version. `make_fig16.py` is the one that produces the paper's Figure 15.

## The physics

One engine, one coordinate frame, every force derivable on paper. `x` runs along the table with the surface at `y = 0`; the near edge is `x = 0`, the far edge `x = 2.74 m`, the net at `x = 1.37 m`. Positive RPM means topspin.

- **Drag** — Reynolds-dependent with a sigmoid drag-crisis transition, `C_D = 0.55 − 0.08 / (1 + exp(1.5e−4 (Re − 28000)))`.
- **Magnus lift** — the Watts–Ferrer saturating model, `C_L = 1 / (2 + 1/S)` with spin parameter `S = ωr / v`.
- **Spin decay** — exponential, `ω(t) = ω₀ exp(−k_s t)` with `k_s = 0.08 s⁻¹`.
- **Turbulence** — additive Gaussian velocity noise, integrated with Euler–Maruyama.

Deterministic runs use RK4. Stochastic runs use Euler–Maruyama instead, because a Wiener noise term is only first-order integrable, so using RK4 there would be wrong rather than more accurate.

Ball and table constants follow ITTF regulation: 40 mm diameter, 2.7 g, 2.74 m table, 15.25 cm net.

## What it doesn't do

The core analysis is two-dimensional. The 3D engine handles sidespin but the main results don't use it, and the larger errors on the real dataset came mostly from shots with strong sidespin, which a model living in a single vertical plane can't fully capture.

Air density is constant through the flight, and the bounce is excluded, so this covers the airborne phase only. The lift law saturates smoothly, so it won't reproduce the non-monotonic "lift crisis" reported for table tennis balls at intermediate spin (Miyazaki et al., 2017), nor the inverse Magnus force seen at very low spin. The spin-decay rate is a single representative value within the published range, not one fitted to the particular ball and surface.

The real-data test uses a marked ball from a robot dataset. It hasn't been tried on ordinary broadcast or phone footage of live play, where spin has to be inferred from the trajectory rather than read off a printed dot.

## Reproducibility

The Monte Carlo scripts are seeded, so repeated runs reproduce the same figures. Figure 12 reseeds per sample size on purpose, so its confidence intervals reflect sampling error at each N rather than one shared draw.

The numbers quoted above are what these scripts print today, and they're what you'll get if you run them. Two have shifted slightly as the engine was unified from the earlier two-integrator version: the landing-displacement reduction is now 83% where it used to come out at 82%, and the entropy drop is 33% where it was 34%. Everything else is unchanged.

## Citation

> Maruthi, S. (2026). *Topspin as an Aerodynamic Stabiliser: A Physics Simulation and Probabilistic Analysis of Table Tennis Ball Flight.* SIMULTECH 2026, Porto, Portugal.

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

## Licence

MIT — see [LICENSE](LICENSE).

## Acknowledgements

Research mentorship from Prof. Ganesh Mani (Carnegie Mellon University). Shot nomenclature and tactical interpretation came from conversations with an experienced table tennis coach. The real-trajectory test uses the SpinDOE dataset (Gossard et al., 2023).

I used Claude (Anthropic) for code debugging, language polishing and grammar checking. All scientific content, simulation code, analysis, figures and conclusions are my own.
