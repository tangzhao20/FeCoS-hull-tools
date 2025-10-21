# FeCoS Hull & Plots
Minimal, script-first pipeline to compute **formation energy** and **energy above the convex hull** from `raw_data.csv` and make quick plots for the Fe–Co–S study.

## Quickstart
```bash
# 1) Create an env (conda or venv) and install deps
conda create -n fecos-tools -y python=3.10
conda activate fecos-tools
pip install pandas numpy matplotlib pymatgen monty

# 2) Compute energies & hull (reads raw_data.csv, writes processed_data.csv)
python compute_energies.py raw_data.csv -o processed_data.csv

# 3) Correlation plots (Ehull/Eform vs Js)
python plot_correlations.py processed_data.csv --outdir plots

# 4) Ternary phase diagrams (Ehull-colored and Js-colored) + Gibbs triangles listing
#   (These scripts DO NOT depend on pymatgen.)
python ternary_ehull_plot.py --in processed_data.csv --out plots/FeCoS_PD.png --facet-tol 1e-4
python ternary_js_plot.py   --in processed_data.csv --out plots/ternary_Js.png --js-min 0 --js-max 2.5 --facet-tol 1e-4
python list_gibbs_triangles.py --in processed_data.csv --out plots/triangles.txt

Notes
	•	Headers expected in processed_data.csv: at least nFe,nCo,nS,Eform/E_form,Ehull/E_hull (and Js for the Js plot).
	•	Stable = Ehull ≤ zero_tol (default 1e-6) → plotted as black dots with labels.
	•	Tie-lines are derived from the lower convex hull in (fFe, fCo, Eform) using all stable vertices (elements + binaries + ternaries).
	•	For near-stable visualization, increase tolerance: --zero-tol 2e-3.

Files
	•	compute_energies.py – reads raw_data.csv, writes processed_data.csv. Uses pymatgen PhaseDiagram.
	•	plot_correlations.py – scatter plots of Eform/Ehull vs Js into --outdir.
	•	ternary_ehull_plot.py – ternary (Fe top, Co BL, S BR), color by Ehull, stable in black, Gibbs tie-lines, dashed grid. No pymatgen.
	•	ternary_js_plot.py – same ternary, color by Js. No pymatgen.
	•	list_gibbs_triangles.py – prints/writes all Gibbs triangles (facet triplets) by formulas. No pymatgen.
	•	raw_data.csv – example input.
	•	processed_data.csv – example output produced by compute_energies.py.

CLI reference
	•	compute_energies.py
	•	input_csv (positional), -o/--output_csv (default processed_data.csv)
	•	plot_correlations.py
	•	processed_csv (positional), --outdir (default .)
	•	ternary_ehull_plot.py
	•	--in, --out, --max-ehull (default 0.5), --zero-tol (default 1e-6), --facet-tol (default 1e-4), --dpi
	•	ternary_js_plot.py
	•	--in, --out, --js-min (default 0), --js-max (default 2.5), --zero-tol, --facet-tol, --dpi
	•	list_gibbs_triangles.py
	•	--in, --out (optional), --zero-tol (default 1e-6), --facet-tol (default 1e-4)

FeCoS-hull-tools/
├─ compute_energies.py
├─ plot_correlations.py
├─ ternary_ehull_plot.py
├─ ternary_js_plot.py
├─ list_gibbs_triangles.py
├─ raw_data.csv             # example input
├─ processed_data.csv       # example output
├─ plots/
│  ├─ Eform_vs_Js.png       # example Eform vs Js (optional)
│  ├─ Ehull_vs_Js.png       # example Ehull vs Js (optional)
│  ├─ Ehull_vs_Js_ylim2.png # example Ehull vs Js zoomed in (optional)
│  ├─ FeCoS_PD.png          # example Ehull ternary (optional)
│  └─ ternary_Js.png        # example Js ternary (optional)
├─ requirements.txt
├─ .gitignore
└─ README.md


