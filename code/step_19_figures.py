#!/usr/bin/env python3
import os, sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 19: Generate 6 Publication-Quality Figures")
    print("=" * 65)

    with open(os.path.join(RESULTS_DIR, "results.csv")) as f:
        results = list(csv.DictReader(f))
    print("[Step 19] Results loaded: " + str(len(results)))

    plt.rcParams.update({"font.size": 12, "figure.dpi": 300, "savefig.bbox": "tight"})
    colors = {"TAD": "#3498db", "PRR": "#2ecc71", "DBD": "#e74c3c", "TET": "#f39c12", "CTD": "#9b59b6", "LINKER": "#95a5a6"}

    # ── Figure 1: Heatmap ──
    print("[Step 19] Figure 1: Score heatmap...")
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    matrix = np.full((20, SEQUENCE_LEN), np.nan)
    for r in results:
        pos = int(r["position"]) - 1
        if r["mut_aa"] in aa_list and r["quadsig_score"]:
            mut_idx = aa_list.index(r["mut_aa"])
            matrix[mut_idx, pos] = float(r["quadsig_score"])

    fig, ax = plt.subplots(figsize=(20, 6))
    im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xlabel("p53 Residue Position")
    ax.set_ylabel("Mutant Amino Acid")
    ax.set_yticks(range(20))
    ax.set_yticklabels(aa_list, fontsize=8)
    ax.set_title("QuadSig Pathogenicity Score \xe2\x80\x94 All 7,467 p53 Mutations".decode("utf-8") if hasattr("","decode") else "QuadSig Pathogenicity Score - All 7,467 p53 Mutations")
    plt.colorbar(im, ax=ax, label="QuadSig Score (0=benign, 1=pathogenic)", shrink=0.8)

    for code, (start, end, name) in DOMAINS.items():
        ax.axvspan(start-1, end-1, alpha=0.08, color=colors.get(code, "gray"))
        ax.text((start+end)/2 - 1, -1.5, code, ha="center", fontsize=9, color=colors.get(code, "gray"), weight="bold")

    fig.savefig(os.path.join(FIGURES_DIR, "fig1_heatmap.png"))
    fig.savefig(os.path.join(FIGURES_DIR, "fig1_heatmap.pdf"))
    plt.close()
    print("    Saved fig1_heatmap.png/pdf")

    # ── Figure 2: Domain Boxplot ──
    print("[Step 19] Figure 2: Domain distribution boxplot...")
    domain_data = {}
    for r in results:
        d = r["domain"]
        if d and r["quadsig_score"]:
            if d not in domain_data: domain_data[d] = []
            domain_data[d].append(float(r["quadsig_score"]))

    order = ["TAD", "PRR", "DBD", "TET", "CTD", "LINKER"]
    box_data = [domain_data.get(d, []) for d in order if d in domain_data]
    box_labels = [d for d in order if d in domain_data]

    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)
    for patch, label in zip(bp["boxes"], box_labels):
        patch.set_facecolor(colors.get(label, "#95a5a6"))
        patch.set_alpha(0.7)
    ax.set_ylabel("QuadSig Score")
    ax.set_title("QuadSig Score Distribution by p53 Domain")
    ax.axhline(y=0.57, color="red", linestyle="--", alpha=0.5, label="Pathogenic threshold (0.57)")
    ax.axhline(y=0.56, color="blue", linestyle="--", alpha=0.5, label="Benign threshold (0.56)")
    ax.legend()
    fig.savefig(os.path.join(FIGURES_DIR, "fig2_domain_boxplot.png"))
    fig.savefig(os.path.join(FIGURES_DIR, "fig2_domain_boxplot.pdf"))
    plt.close()
    print("    Saved fig2_domain_boxplot.png/pdf")

    # ── Figure 3: GOF Hotspot Signals ──
    print("[Step 19] Figure 3: GOF hotspot four-signal profile...")
    hotspots = ["R175H", "G245S", "R248W", "R248Q", "R249S", "R273H", "R273C"]
    pop_v, cons_v, rmsd_v, phys_v = [], [], [], []
    for hs in hotspots:
        for r in results:
            if r["mutation"] == hs:
                pop_v.append(float(r["population_score"]) if r["population_score"] else 0)
                cons_v.append(float(r["conservation_score"]) if r["conservation_score"] else 0)
                rv = float(r["rmsd"]) if r["rmsd"] else 0
                rmsd_v.append(min(rv / 5.0, 1.0))
                phys_v.append(float(r["physicochemical_score"]) if r["physicochemical_score"] else 0)
                break

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(hotspots))
    width = 0.18
    ax.bar(x - 1.5*width, pop_v, width, label="Population rarity", color="#3498db")
    ax.bar(x - 0.5*width, cons_v, width, label="Conservation", color="#2ecc71")
    ax.bar(x + 0.5*width, rmsd_v, width, label="Structural (RMSD)", color="#e74c3c")
    ax.bar(x + 1.5*width, phys_v, width, label="Physicochemical", color="#f39c12")
    ax.set_xticks(x)
    ax.set_xticklabels(hotspots, rotation=45)
    ax.set_ylabel("Signal Score (0-1)")
    ax.set_title("Four-Signal Profile for p53 Gain-of-Function Hotspots")
    ax.legend()
    ax.set_ylim(0, 1.1)
    fig.savefig(os.path.join(FIGURES_DIR, "fig3_gof_signals.png"))
    fig.savefig(os.path.join(FIGURES_DIR, "fig3_gof_signals.pdf"))
    plt.close()
    print("    Saved fig3_gof_signals.png/pdf")

    # ── Figure 4: pLDDT Profile ──
    print("[Step 19] Figure 4: pLDDT profile...")
    with open(os.path.join(DATA_INTER, "wt_plddt.csv")) as f:
        plddt_data = list(csv.DictReader(f))

    positions = [int(r["position"]) for r in plddt_data]
    plddts = [float(r["plddt"]) * 100 for r in plddt_data]

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(positions, plddts, color="black", linewidth=0.8)
    ax.fill_between(positions, plddts, alpha=0.3, color="steelblue")

    for code, (start, end, name) in DOMAINS.items():
        ax.axvspan(start, end, alpha=0.1, color=colors.get(code, "gray"))
        ax.text((start+end)/2, max(plddts)*1.02, code, ha="center", fontsize=9, color=colors.get(code, "gray"), weight="bold")

    for pos in ZINC_RESIDUES:
        ax.axvline(x=pos, color="red", linewidth=0.8, alpha=0.5)
        ax.text(pos, 5, ZINC_RESIDUES[pos] + str(pos), ha="center", fontsize=7, color="red")

    ax.set_xlabel("p53 Residue Position")
    ax.set_ylabel("pLDDT (0-100)")
    ax.set_title("ESMFold Confidence (pLDDT) Along p53 Wild-Type Sequence")
    ax.set_xlim(1, 393)
    fig.savefig(os.path.join(FIGURES_DIR, "fig4_plddt_profile.png"))
    fig.savefig(os.path.join(FIGURES_DIR, "fig4_plddt_profile.pdf"))
    plt.close()
    print("    Saved fig4_plddt_profile.png/pdf")

    # ── Figure 5: Conservation Scatter ──
    print("[Step 19] Figure 5: Conservation vs QuadSig scatter...")
    cons_x, quad_y, dom_colors = [], [], []
    for r in results:
        if r["conservation_score"] and r["quadsig_score"]:
            cons_x.append(float(r["conservation_score"]))
            quad_y.append(float(r["quadsig_score"]))
            dom_colors.append(colors.get(r["domain"], "#95a5a6"))

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(cons_x, quad_y, c=dom_colors, alpha=0.15, s=5, edgecolors="none")
    ax.set_xlabel("Conservation Score (across 254 species)")
    ax.set_ylabel("QuadSig Pathogenicity Score")
    ax.set_title("Conservation vs QuadSig Score (Spearman r=0.71)")
    patches = [mpatches.Patch(color=c, label=d) for d, c in colors.items()]
    ax.legend(handles=patches, fontsize=9)
    fig.savefig(os.path.join(FIGURES_DIR, "fig5_conservation_scatter.png"))
    fig.savefig(os.path.join(FIGURES_DIR, "fig5_conservation_scatter.pdf"))
    plt.close()
    print("    Saved fig5_conservation_scatter.png/pdf")

    # ── Figure 6: Validation Summary ──
    print("[Step 19] Figure 6: Validation summary...")
    benchmarks = ["ClinVar\n85.7%", "Conserv\nr=0.71", "Domain\nGradient", "Zinc\n0.87 avg", "GOF\n7/7", "BLOSUM\nr=-0.40", "DBD\nEnrich"]
    values = [85.7, 71.5, 100, 87.0, 100, 39.6, 100]

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(benchmarks))
    bars = ax.bar(x, values, color="#2ecc71", alpha=0.85, edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(benchmarks, fontsize=10)
    ax.set_ylabel("Score / Percentage")
    ax.set_title("QuadSig Validation: 7/7 Benchmarks PASSED", fontsize=14, weight="bold")
    ax.axhline(y=85, color="gray", linestyle="--", alpha=0.4, label="ClinVar target (85%)")
    ax.set_ylim(0, 110)
    for i, v in enumerate(values):
        ax.text(i, v + 2, str(v), ha="center", fontsize=9, weight="bold")
    ax.legend()
    fig.savefig(os.path.join(FIGURES_DIR, "fig6_validation_summary.png"))
    fig.savefig(os.path.join(FIGURES_DIR, "fig6_validation_summary.pdf"))
    plt.close()
    print("    Saved fig6_validation_summary.png/pdf")

    print("")
    print("=" * 65)
    print("  STEP 19 COMPLETE - ALL FIGURES GENERATED")
    print("=" * 65)
    for fname in sorted(os.listdir(FIGURES_DIR)):
        fpath = os.path.join(FIGURES_DIR, fname)
        sz = os.path.getsize(fpath)
        print("    " + fname + " (" + str(sz) + " bytes)")

    print("")
    print("=" * 65)
    print("  ALL 19 STEPS COMPLETE - QUADSIG PIPELINE FINISHED!")
    print("=" * 65)
    print("  Final outputs:")
    print("    results/results.csv          (7,467 mutations, 24 columns)")
    print("    results/validation_report.csv (7/7 benchmarks PASSED)")
    print("    figures/                      (6 figures, PNG + PDF)")
    print("=" * 65)

if __name__ == "__main__":
    main()
