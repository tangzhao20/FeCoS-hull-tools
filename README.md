# FeCoS Hull & Plots
Minimal, script-first pipeline to compute **formation energy** and **energy above the convex hull** from `raw_data.csv` and make quick plots for the Fe–Co–S study.

## Quickstart
```bash
# 1) Create an env (conda or venv) and install deps
conda create -n fecos-tools -y python=3.10
conda activate fecos-tools
pip install pandas numpy matplotlib pymatgen monty python-ternary

# 2) Compute energies & hull (reads raw_data.csv, writes processed_data.csv)
python compute_energies.py raw_data.csv -o processed_data.csv

# 3) Correlation plots (Ehull/Eform vs Js)
python plot_correlations.py processed_data.csv --outdir plots

# 4) Ternary phase diagrams (color by Ehull and by Js)
python plot_phase_diagram.py processed_data.csv --outdir plots
```

## Files
- `compute_energies.py` – reads `raw_data.csv` (positional arg), writes `processed_data.csv` (default) or `-o/--output_csv`.
  - Uses `pymatgen` PhaseDiagram. Assumes elemental reference rows named `Fe_ref`, `Co_ref`, `S_ref` in the CSV.
- `plot_correlations.py` – **positional** `processed_csv`, saves to `--outdir` (default `.`). Creates:
  - `Eform_vs_Js.png`, `Ehull_vs_Js.png`, `Ehull_vs_Js_ylim2.png`.
- `plot_phase_diagram.py` – **positional** `processed_csv`, saves to `--outdir` (default `.`). Creates:
  - `hull_Ehull.png`, `hull_Js.png`.
- `raw_data.csv` – example input (columns: `index,F,E0,mag,vol,nFe,nCo,nS` with `*_ref` rows).
- `processed_data.csv` – example output produced by `compute_energies.py`.

## CLI reference
- `compute_energies.py`:
  - `input_csv` (positional) – e.g., `raw_data.csv`
  - `-o, --output_csv` (default: `processed_data.csv`)

- `plot_correlations.py`:
  - `processed_csv` (positional)
  - `--outdir` (default: `.`)

- `plot_phase_diagram.py`:
  - `processed_csv` (positional)
  - `--outdir` (default: `.`)

## Notes
- **Units**: energies in eV; counts `nFe,nCo,nS` are integers per cell.
- **Refs**: include `Fe_ref`, `Co_ref`, `S_ref` rows in `raw_data.csv` for per-atom elemental energies.
- **Outputs** land in `--outdir` (create e.g. `plots/` first: `mkdir -p plots`).

## Repo layout
```
FeCoS-hull-tools/
├─ compute_energies.py
├─ plot_correlations.py
├─ plot_phase_diagram.py
├─ raw_data.csv           # example input
├─ processed_data.csv     # example output
├─ requirements.txt
├─ .gitignore
└─ README.md
```

## License
MIT (see `LICENSE`).
