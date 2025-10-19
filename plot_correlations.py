#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    p = argparse.ArgumentParser()
    p.add_argument("processed_csv", help="CSV with Eform, Ehull, Js columns")
    p.add_argument("--outdir", default=".", help="output directory")
    args = p.parse_args()

    df = pd.read_csv(args.processed_csv)

    # 1) Formation energy vs Js (Js on x, Eform on y)
    plt.figure()
    plt.scatter(
        df["Js"], df["Eform"],
        facecolors='none', edgecolors='blue', linewidths=0.8
    )
    plt.axvline(1.0, color="grey", linestyle="--", label="Js = 1 T")
    plt.axhline(0.0, color="grey", linestyle="--", label="Eform = 0")
    plt.xlabel("Js (T)")
    plt.ylabel("Formation Energy per atom (eV)")
    plt.title("Eform vs Js")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/Eform_vs_Js.png", dpi=300)
    plt.close()
    print("Saved Eform_vs_Js.png")

    # 2a) Energy above hull vs Js (full scale)
    plt.figure()
    plt.scatter(
        df["Js"], df["Ehull"],
        facecolors='none', edgecolors='green', linewidths=0.8
    )
    plt.axvline(1.0, color="grey", linestyle="--", label="Js = 1 T")
    plt.axhline(0.1, color="grey", linestyle="--", label="Ehull = 0.1 eV/at")
    plt.xlabel("Js (T)")
    plt.ylabel("Energy above hull (eV/atom)")
    plt.title("Ehull vs Js")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/Ehull_vs_Js.png", dpi=300)
    plt.close()
    print("Saved Ehull_vs_Js.png")

    # 2b) Energy above hull vs Js (y-limit = 2 eV/atom)
    plt.figure()
    plt.scatter(
        df["Js"], df["Ehull"],
        facecolors='none', edgecolors='green', linewidths=0.8
    )
    plt.ylim(0, 2.0)
    plt.axvline(1.0, color="grey", linestyle="--", label="Js = 1 T")
    plt.axhline(0.1, color="grey", linestyle="--", label="Ehull = 0.1 eV/at")
    plt.xlabel("Js (T)")
    plt.ylabel("Energy above hull (eV/atom)")
    plt.title("Ehull vs Js (y ≤ 2 eV/at)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/Ehull_vs_Js_ylim2.png", dpi=300)
    plt.close()
    print("Saved Ehull_vs_Js_ylim2.png")

if __name__ == "__main__":
    main()

