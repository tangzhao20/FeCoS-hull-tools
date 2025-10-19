#!/usr/bin/env python3
import argparse, math, inspect

# ─── Monkey-patch monty.dev.deprecated if needed ───
try:
    import monty.dev
    _orig = monty.dev.deprecated
    sig = inspect.signature(_orig)
    if "deadline" not in sig.parameters:
        def _patched(*a, **k):
            k.pop("deadline", None)
            return _orig(*a, **k)
        monty.dev.deprecated = _patched
except ImportError:
    pass
# ────────────────────────────────────────────────────

import pandas as pd
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram

# Constants for Js calculation
angm = 1e-10           # 1 Å in m
eV2J = 1.60217653e-19  # J per eV
mu_B = 5.7883818e-11   # MeV/T
mu_0 = math.pi*4e-7    # T·m/A

def compute_js(mag, vol):
    ms = mag * mu_B * eV2J / (vol * angm**3)  # MA/m
    return ms * 1e6 * mu_0                     # Tesla

def main():
    parser = argparse.ArgumentParser(
        description="Compute Eform, Ehull, and Js from a single raw_data.csv"
    )
    parser.add_argument("input_csv", help="raw_data.csv")
    parser.add_argument("-o","--output_csv", default="processed_data.csv")
    args = parser.parse_args()

    # 1) Read everything
    df = pd.read_csv(args.input_csv, sep=r"[\t,]+", engine="python")

    # 2) Find the elemental reference rows
    ref_df = df[df["index"].astype(str).str.endswith("_ref")]
    per_atom_ref = {}
    for _, row in ref_df.iterrows():
        # Build a Composition from whatever atom counts you gave
        comp = Composition({
            "Fe": int(row["nFe"]),
            "Co": int(row["nCo"]),
            "S" : int(row["nS"])
        })
        E0 = float(row["E0"])  # this is the TOTAL energy of that cell
        per_atom_ref[ list(comp.as_dict().keys())[0] ] = E0 / comp.num_atoms

    # sanity check
    for el in ("Fe","Co","S"):
        if el not in per_atom_ref:
            raise RuntimeError(f"Missing reference entry for {el}_ref in your CSV")

    # 3) Loop all rows and build entries + compute Eform & Js
    entries, eform_list, js_list = [], [], []
    for _, row in df.iterrows():
        idx  = str(row["index"])
        E0   = float(row["E0"])
        mag  = float(row["mag"])
        vol  = float(row["vol"])
        nFe  = int(row["nFe"])
        nCo  = int(row["nCo"])
        nS   = int(row["nS"])
        comp = Composition({"Fe": nFe, "Co": nCo, "S": nS})

        # formation energy per atom
        ref_total = (
            nFe*per_atom_ref["Fe"] +
            nCo*per_atom_ref["Co"] +
            nS *per_atom_ref["S"]
        )
        Ef = (E0 - ref_total) / comp.num_atoms
        eform_list.append(Ef)

        # compute Js from mag & volume
        Js = compute_js(mag, vol)
        js_list.append(Js)

        entries.append(
            ComputedEntry(comp, E0, entry_id=idx,
                          parameters={"Ef": Ef, "Js": Js})
        )

    # 4) Build the convex hull (now includes your Fe_ref, Co_ref, S_ref entries)
    pdia = PhaseDiagram(entries)
    ehull_list = [pdia.get_e_above_hull(e) for e in entries]

    # 5) Write back to CSV
    df["Eform"] = eform_list
    df["Ehull"] = ehull_list
    df["Js"]    = js_list
    df.to_csv(args.output_csv, index=False, float_format="%.6f")
    print(f"Wrote {len(df)} rows → {args.output_csv}")

if __name__ == "__main__":
    main()

