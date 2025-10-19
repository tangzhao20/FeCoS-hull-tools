#!/usr/bin/env python3
import inspect

# ───────────────────────────────────────────────────────────────────────────────
# Monkey-patch monty.dev.deprecated to drop unsupported 'deadline'
try:
    import monty.dev
    _orig = monty.dev.deprecated
    sig   = inspect.signature(_orig)
    if "deadline" not in sig.parameters:
        def _patched(*a, **kw):
            kw.pop("deadline", None)
            return _orig(*a, **kw)
        monty.dev.deprecated = _patched
except ImportError:
    pass
# ───────────────────────────────────────────────────────────────────────────────

import argparse
import pandas as pd
import numpy as np
import ternary
import matplotlib.pyplot as plt
import matplotlib as mpl
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram

def main():
    p = argparse.ArgumentParser(
        description="Plot Fe–Co–S convex hull colored by Ehull and Js"
    )
    p.add_argument("processed_csv", help="CSV with nFe,nCo,nS,E0,Ehull,Js,index")
    p.add_argument("--outdir",    default=".", help="where to save images")
    args = p.parse_args()

    df = pd.read_csv(args.processed_csv)
    needed = {"nFe","nCo","nS","E0","Ehull","Js","index"}
    if missing := needed - set(df.columns):
        raise KeyError(f"Missing columns: {missing}")

    # build entries for stability check
    entries = []
    for _, r in df.iterrows():
        comp = Composition({"Fe": int(r["nFe"]),
                            "Co": int(r["nCo"]),
                            "S" : int(r["nS"])})
        entries.append(ComputedEntry(comp,
                                    float(r["E0"]),
                                    entry_id=str(r["index"]),
                                    parameters={"Ehull": float(r["Ehull"]),
                                                "Js":   float(r["Js"])}))

    pdia = PhaseDiagram(entries)
    stable_pts = [(e.composition["Fe"],
                   e.composition["Co"],
                   e.composition["S"])
                  for e in pdia.stable_entries]

    # coords + data arrays
    pts   = [(int(r["nFe"]), int(r["nCo"]), int(r["nS"]))
             for _, r in df.iterrows()]
    ehull = df["Ehull"].to_numpy()
    js    = df["Js"].to_numpy()
    scale = max(sum(pt) for pt in pts)

    # prepare RGBA arrays
    norm_e = mpl.colors.Normalize(vmin=ehull.min(), vmax=ehull.max())
    cmap_e = mpl.colormaps["viridis"]
    colors_e = cmap_e(norm_e(ehull))

    norm_j = mpl.colors.Normalize(vmin=js.min(), vmax=js.max())
    cmap_j = mpl.colormaps["plasma"]
    colors_j = cmap_j(norm_j(js))

    # ─── Ehull plot ───
    fig, tax = ternary.figure(scale=scale)
    tax.boundary(); tax.gridlines(color="grey", multiple=1)
    tax.scatter(pts, marker="o", color=colors_e, label="Ehull (eV/at)")
    tax.scatter(stable_pts, marker="D", color="red", label="stable")
    tax.set_title("Fe–Co–S Hull (colored by Ehull)")
    tax.clear_matplotlib_ticks()
    fig.savefig(f"{args.outdir}/hull_Ehull.png", dpi=300)
    print("Saved hull_Ehull.png")

    # ─── Js plot ───
    fig2, tax2 = ternary.figure(scale=scale)
    tax2.boundary(); tax2.gridlines(color="grey", multiple=1)
    tax2.scatter(pts, marker="o", color=colors_j, label="Js (T)")
    tax2.scatter(stable_pts, marker="D", color="red", label="stable")
    tax2.set_title("Fe–Co–S Hull (colored by Js)")
    tax2.clear_matplotlib_ticks()
    fig2.savefig(f"{args.outdir}/hull_Js.png", dpi=300)
    print("Saved hull_Js.png")

if __name__ == "__main__":
    main()

