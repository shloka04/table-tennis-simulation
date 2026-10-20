"""Run every figure script in order. Produces figures/fig1..fig16.
Figs 1-15 need only tt_engine.py / tt_engine3d.py. Fig 16 needs the SpinDOE
dataset in spindoe/data/trajectories/_t/ (clone + unzip first)."""
import os, sys, subprocess
os.environ["MPLBACKEND"] = "Agg"          # headless-safe plotting
os.makedirs("figures", exist_ok=True)
def run(*cmd):
    print(f"\n===== {' '.join(cmd)} =====", flush=True)
    subprocess.run([sys.executable, *cmd], check=True, env={**os.environ})
for s in ["figures.py", "extras.py", "fig13_sidespin.py", "fig14_optimal.py", "fig15_validation.py"]:
    run(s)
# figure 16 (real-data validation) -- only if the SpinDOE data is present
if os.path.exists("spindoe/data/trajectories/_t/index.csv"):
    run("convert_spindoe.py", "spindoe/data/trajectories/_t")
    run("make_fig16.py")
else:
    print("\n[fig16] SpinDOE data not found. To make figure 16, first run:\n"
          "  git clone --depth 1 https://github.com/cogsys-tuebingen/spindoe.git\n"
          "  cd spindoe/data/trajectories && unzip -o trajectories.zip -d _t && cd ../../..\n"
          "  then re-run this script (or run convert_spindoe.py + make_fig16.py).")
pngs = sorted(f for f in os.listdir("figures") if f.endswith(".png"))
print(f"\nALL DONE -> {len(pngs)} figures in figures/:")
for p in pngs: print("   ", p)
