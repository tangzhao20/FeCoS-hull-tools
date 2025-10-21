#!/usr/bin/env python3
"""
list_gibbs_triangles.py — enumerate Gibbs triangles (lower-hull facets) by formula

CSV expected columns (your header): index,F,E0,mag,vol,nFe,nCo,nS,Eform,Ehull,Js

What it does
------------
- Keeps the lowest-Ehull entry per *reduced* composition (divide by gcd).
- Takes as "stable vertices" those with Ehull <= zero_tol,
  and also includes pure Fe, Co, S with Eform = 0.
- Builds the LOWER convex hull in (x=fFe, y=fCo, z=Eform)
  and lists each triangular facet by its 3 formulas (order-insensitive).

Usage
-----
python list_gibbs_triangles.py --in processed_data.csv
python list_gibbs_triangles.py --in processed_data.csv --zero-tol 1e-6 --facet-tol 1e-4 --out triangles.txt
"""
import argparse
import math
from itertools import combinations
from typing import Tuple, Dict, List, Set

import numpy as np
import pandas as pd


# ---------- helpers ----------
def gcd3(a: int, b: int, c: int) -> int:
    from math import gcd
    return gcd(gcd(abs(a), abs(b)), abs(c))

def reduce_composition(nFe: int, nCo: int, nS: int) -> Tuple[int, int, int]:
    g = gcd3(max(nFe,0), max(nCo,0), max(nS,0)) or 1
    return (nFe // g, nCo // g, nS // g)

def formula_from_counts(nFe: int, nCo: int, nS: int) -> str:
    parts = []
    for el, n in (("Fe", nFe), ("Co", nCo), ("S", nS)):
        if n > 0:
            parts.append(f"{el}{'' if n == 1 else n}")
    return "".join(parts) if parts else "—"

def fractions(nFe: int, nCo: int, nS: int) -> Tuple[float, float, float]:
    tot = nFe + nCo + nS
    if tot <= 0:
        return (0.0, 0.0, 0.0)
    return (nFe / tot, nCo / tot, nS / tot)

def fit_plane(p1, p2, p3, tol_area: float = 1e-12):
    """Given three points (x,y,z), return (a,b,c) for z = a x + b y + c. None if nearly collinear."""
    (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) = p1, p2, p3
    area2 = abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
    if area2 < tol_area:
        return None
    A = np.array([[x1, y1, 1.0],
                  [x2, y2, 1.0],
                  [x3, y3, 1.0]], dtype=float)
    b = np.array([z1, z2, z3], dtype=float)
    a, beta, c = np.linalg.solve(A, b)
    return a, beta, c


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="List Gibbs triangles (lower-hull facets) by their 3 endpoint formulas.")
    ap.add_argument("--in", dest="input_csv", required=True)
    ap.add_argument("--zero-tol", type=float, default=1e-6, help="stable if Ehull <= this")
    ap.add_argument("--facet-tol", type=float, default=1e-4, help="lower-hull inequality tolerance")
    ap.add_argument("--out", dest="out_txt", default="", help="optional output text file")
    args = ap.parse_args()

    df = pd.read_csv(args.input_csv)
    for col in ("nFe", "nCo", "nS", "Eform", "Ehull"):
        if col not in df.columns:
            raise ValueError(f"CSV missing column: {col}")

    # Keep lowest-Ehull per reduced composition
    best: Dict[Tuple[int,int,int], Dict] = {}
    for _, r in df.iterrows():
        key = reduce_composition(int(r["nFe"]), int(r["nCo"]), int(r["nS"]))
        eh = float(r["Ehull"])
        if key not in best or eh < best[key]["Ehull"]:
            best[key] = {
                "nFe": key[0], "nCo": key[1], "nS": key[2],
                "Eform": float(r["Eform"]),
                "Ehull": eh,
                "formula": r.get("formula", formula_from_counts(*key))
            }

    # Build stable vertices (x=fFe, y=fCo, z=Eform)
    verts: List[Dict] = []
    # Pure elements (Eform = 0 per atom)
    verts.append({"x": 1.0, "y": 0.0, "z": 0.0, "formula": "Fe"})
    verts.append({"x": 0.0, "y": 1.0, "z": 0.0, "formula": "Co"})
    verts.append({"x": 0.0, "y": 0.0, "z": 0.0, "formula": "S"})

    for rec in best.values():
        if rec["Ehull"] <= args.zero_tol:
            nFe, nCo, nS = rec["nFe"], rec["nCo"], rec["nS"]
            fFe, fCo, _ = fractions(nFe, nCo, nS)
            verts.append({"x": fFe, "y": fCo, "z": rec["Eform"], "formula": rec["formula"]})

    if len(verts) < 3:
        print("No facets: need at least 3 stable vertices (including elements).")
        return

    # Lower convex hull enumeration
    triangles: Set[Tuple[str, str, str]] = set()
    N = len(verts)
    for i, j, k in combinations(range(N), 3):
        p1 = (verts[i]["x"], verts[i]["y"], verts[i]["z"])
        p2 = (verts[j]["x"], verts[j]["y"], verts[j]["z"])
        p3 = (verts[k]["x"], verts[k]["y"], verts[k]["z"])
        plane = fit_plane(p1, p2, p3)
        if plane is None:
            continue
        a, b, c = plane
        # Keep if all vertices satisfy z >= a*x + b*y + c - facet_tol (lower hull)
        for v in verts:
            if v["z"] < (a*v["x"] + b*v["y"] + c) - args.facet_tol:
                break
        else:
            tri = tuple(sorted([verts[i]["formula"], verts[j]["formula"], verts[k]["formula"]]))
            triangles.add(tri)

    tri_list = sorted(triangles)

    # Print (and optionally write) results
    for idx, tri in enumerate(tri_list, 1):
        print(f"Triangle {idx}: {tri[0]}, {tri[1]}, {tri[2]}")

    if args.out_txt:
        with open(args.out_txt, "w") as f:
            for idx, tri in enumerate(tri_list, 1):
                f.write(f"Triangle {idx}: {tri[0]}, {tri[1]}, {tri[2]}\n")
        print(f"Wrote {len(tri_list)} triangles to {args.out_txt}")

if __name__ == "__main__":
    main()
