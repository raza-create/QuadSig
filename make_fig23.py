#!/usr/bin/env python3
"""Figures 2 and 3 for the QuadSig v2 manuscript, redrawn at the template text width
(408 pt = 5.65 in) so that no text prints below 6 pt.

Reads ONLY results_v2/results.csv and results_v2/parameters.json (written by rerun_v2.py).
Writes results_v2/final_figs/{fig2,fig3}.{pdf,png} and figure_numbers.txt, which lists
every number drawn in the figures so they can be checked against results_v2/summary.txt.
Usage (from ~/QuadSig):  python3 make_fig23.py            [or: python3 make_fig23.py ROOT]
"""
import json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
from scipy.stats import spearmanr

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results_v2")
OUT = os.path.join(RES, "final_figs")
os.makedirs(OUT, exist_ok=True)

have = {f.name for f in font_manager.fontManager.ttflist}
FONT = ["Liberation Sans", "DejaVu Sans"] if "Liberation Sans" in have else ["DejaVu Sans"]
plt.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.3,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 600})
W = 5.65                                       # inches = \textwidth of interact.cls
LP, AMB, LB = "#c8413b", "#9a9a9a", "#3b6fb6"
DOMC = {"DBD": "#3b6fb6", "TET": "#e07b39", "TAD": "#d4a017", "LINKER": "#8c8c8c",
        "CTD": "#8e5ea2", "PRR": "#2a9d8f"}

d = pd.read_csv(os.path.join(RES, "results.csv"))
P = json.load(open(os.path.join(RES, "parameters.json")))
lo, hi = float(P["threshold_low"]), float(P["threshold_high"])
cv_acc = float(P["cv"]["accuracy"])
n = len(d)
log = []


def note(s):
    log.append(s); print(s)


def letter(ax, s, title=""):
    ax.set_title(f"$\\bf{{{s}}}$   {title}" if title else f"$\\bf{{{s}}}$", loc="left", fontsize=8, pad=4)


def layout_check(fig, name):
    """No two visible texts may overlap; every text must lie inside the figure; fonts >= 6 pt."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    skip = set()                       # tick labels outside the axis limits are never drawn
    for a_ in fig.axes:
        for axis in (a_.xaxis, a_.yaxis):
            v0, v1 = sorted(axis.get_view_interval())
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                if not (v0 - 1e-9 <= tk.get_loc() <= v1 + 1e-9):
                    skip.update({id(tk.label1), id(tk.label2)})
    items = []
    for t in fig.findobj(Text):
        if id(t) in skip or not t.get_visible() or not t.get_text().strip():
            continue
        bb = t.get_window_extent(r)
        if bb.width < 1:
            continue
        items.append((t, bb))
    bad = 0
    for t, bb in items:
        if bb.x0 < fb.x0 - 1 or bb.x1 > fb.x1 + 1 or bb.y0 < fb.y0 - 1 or bb.y1 > fb.y1 + 1:
            bad += 1; print(f"  [{name}] OUTSIDE FIGURE: {t.get_text()!r}")
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            if a.x0 < b.x1 - 0.5 and b.x0 < a.x1 - 0.5 and a.y0 < b.y1 - 0.5 and b.y0 < a.y1 - 0.5:
                bad += 1; print(f"  [{name}] OVERLAP: {items[i][0].get_text()!r} / {items[j][0].get_text()!r}")
    fmin = min(t.get_fontsize() for t, _ in items)
    note(f"  [{name}] layout check: {len(items)} texts, {bad} problems, smallest font {fmin:g} pt")
    return bad


# ================================ Figure 2 ===============================================
note(f"Figure 2  (n = {n:,}; thresholds {lo:.4f} / {hi:.4f})")
fig = plt.figure(figsize=(W, 4.05))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.08], height_ratios=[1.0, 0.30],
                      left=0.085, right=0.985, top=0.94, bottom=0.085, wspace=0.30, hspace=0.62)
CATS = [("likely_benign", "Likely benign", LB), ("ambiguous", "Ambiguous", AMB),
        ("likely_pathogenic", "Likely pathogenic", LP)]
cnt = {c: int((d.category == c).sum()) for c, _, _ in CATS}

# a: score histogram stacked by class (class ranges in the legend)
ax = fig.add_subplot(gs[0, 0])
bins = np.arange(0, 1.0001, 0.02)
ax.hist([d.loc[d.category == c, "quadsig_score"] for c, _, _ in CATS], bins=bins, stacked=True,
        color=[col for _, _, col in CATS], edgecolor="white", linewidth=0.25, zorder=2,
        label=[f"{lab} {rng} (n = {cnt[c]:,})" for (c, lab, _), rng in
               zip(CATS, [f"< {lo:.3f}", f"{lo:.3f}–{hi:.3f}", f"≥ {hi:.3f}"])])
ax.set_xlim(0, 1); ax.set_xlabel("QuadSig score"); ax.set_ylabel("Number of variants")
ax.set_ylim(0, ax.get_ylim()[1] * 1.42)
ax.legend(frameon=False, loc="upper left", handlelength=1.0, handleheight=0.8,
          borderaxespad=0.1, labelspacing=0.3)
letter(ax, "a", "Score distribution")
for c, lab, _ in CATS:
    note(f"  a: {lab} n = {cnt[c]:,}")

# b: score by domain (violin + box), sorted by mean; dashed line = likely-pathogenic threshold
ax = fig.add_subplot(gs[0, 1])
order = d.groupby("domain_code").quadsig_score.mean().sort_values(ascending=False).index.tolist()
data = [d.loc[d.domain_code == g, "quadsig_score"].to_numpy() for g in order]
vp = ax.violinplot(data, showextrema=False, widths=0.82)
for body, g in zip(vp["bodies"], order):
    body.set_facecolor(DOMC.get(g, "#777777")); body.set_edgecolor(DOMC.get(g, "#777777"))
    body.set_alpha(0.35); body.set_linewidth(0.5)
ax.boxplot(data, widths=0.16, showfliers=False, medianprops=dict(color="black", lw=0.8),
           boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
ax.axhline(hi, ls="--", lw=0.6, color="black", zorder=0)
for i, v in enumerate(data, 1):
    ax.text(i, 1.035, f"{v.mean():.2f}", ha="center", va="bottom", fontsize=6.3, fontweight="bold")
    note(f"  b: {order[i-1]} n = {len(v):,}, mean {v.mean():.2f}")
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels([f"{g}\nn = {len(v):,}" for g, v in zip(order, data)], fontsize=6.2)
ax.set_xlim(0.45, len(order) + 0.55); ax.set_ylim(0, 1.12); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_ylabel("QuadSig score")
letter(ax, "b", f"Score by domain (dashed: {hi:.3f})")

# c: share of the 7,467 variants in each class
ax = fig.add_subplot(gs[1, :])
left = 0.0
ax.set_xlim(0, 100); ax.set_ylim(-0.38, 0.62)   # fix the axes first so label widths can be measured
for c, lab, col in [CATS[2], CATS[1], CATS[0]]:
    w = 100 * cnt[c] / n
    ax.barh(0, w, left=left, color=col, edgecolor="white", height=0.62)
    t = ax.text(left + w / 2, 0, f"{lab}\n{cnt[c]:,} ({w:.1f}%)", ha="center", va="center",
                color="#222222" if c == "ambiguous" else "white", fontsize=6.4, fontweight="bold",
                linespacing=1.15)
    fig.canvas.draw()
    tw = t.get_window_extent(fig.canvas.get_renderer()).width
    sw = ax.transData.transform((left + w, 0))[0] - ax.transData.transform((left, 0))[0]
    if tw > sw - 10:                    # does not fit: put it above the segment instead
        t.set_position((left + w / 2, 0.34)); t.set_va("bottom"); t.set_color("#222222")
        t.set_text(f"{lab} {cnt[c]:,} ({w:.1f}%)")
    note(f"  c: {lab} {cnt[c]:,} ({w:.1f}%)")
    left += w
ax.set_yticks([])
ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xlabel(f"Percentage of all {n:,} variants")
ax.spines["left"].set_visible(False)
letter(ax, "c", "Final classification")
bad2 = layout_check(fig, "fig2")
for e in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"fig2.{e}"))
plt.close(fig)

# ================================ Figure 3 ===============================================
rho_c = spearmanr(d.conservation_score, d.quadsig_score).correlation
rho_b = spearmanr(d.blosum62, d.quadsig_score).correlation
note(f"Figure 3  (Spearman conservation {rho_c:.2f}; BLOSUM62 {rho_b:.2f}; CV accuracy {100*cv_acc:.1f}%)")
fig, axes = plt.subplots(1, 3, figsize=(W, 2.35), gridspec_kw=dict(wspace=0.42))
fig.subplots_adjust(left=0.075, right=0.975, top=0.80, bottom=0.17)
M = lambda v: f"{v:.2f}".replace("-", "\u2212")

ax = axes[0]
for g in order:
    x = d[d.domain_code == g]
    ax.scatter(x.conservation_score, x.quadsig_score, s=1.6, alpha=0.45, color=DOMC.get(g, "#777777"),
               lw=0, rasterized=True, label=g)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("Conservation score"); ax.set_ylabel("QuadSig score")
leg = ax.legend(frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(-0.02, 1.0), markerscale=3.5,
                handletextpad=0.1, columnspacing=0.6, borderaxespad=0.1, fontsize=6.0)
for h_ in (getattr(leg, "legend_handles", None) or leg.legendHandles):
    h_.set_alpha(1)                    # solid legend dots
ax.set_xticks([0, 0.5, 1.0]); ax.set_xticklabels(["0", "0.5", "1.0"])
ax.set_title(f"$\\bf{{a}}$   Spearman ρ = {M(rho_c)}", loc="left", fontsize=7, pad=24)

ax = axes[1]
vals = sorted(d.blosum62.unique())
bx = ax.boxplot([d.loc[d.blosum62 == v, "quadsig_score"] for v in vals], showfliers=False,
                patch_artist=True, widths=0.6, medianprops=dict(color="black", lw=0.8),
                boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
cm = plt.get_cmap("RdYlGn")
for i, p in enumerate(bx["boxes"]):
    p.set_facecolor(cm(i / max(len(vals) - 1, 1)))
ax.set_xticks(range(1, len(vals) + 1))
ax.set_xlabel("BLOSUM62 score"); ax.set_ylabel("QuadSig score"); ax.set_ylim(0, 1)
ax.set_xticklabels([str(v).replace("-", "\u2212") for v in vals])
ax.set_title(f"$\\bf{{b}}$   Spearman ρ = {M(rho_b)}", loc="left", fontsize=7, pad=24)

ax = axes[2]
bins = np.arange(0, 1.0001, 0.025)
nb = int((d.clinvar_label == "benign").sum()); npat = int((d.clinvar_label == "pathogenic").sum())
ax.hist(d.loc[d.clinvar_label == "benign", "quadsig_score"], bins=bins, color=LB, alpha=0.6,
        label=f"Benign (n = {nb})")
ax.hist(d.loc[d.clinvar_label == "pathogenic", "quadsig_score"], bins=bins, color=LP, alpha=0.6,
        label=f"Pathogenic (n = {npat})")
ax.axvline(hi, ymax=1 / 1.45, ls="--", lw=0.6, color="black", label=f"Threshold {hi:.3f}")
ax.set_xlim(0, 1); ax.set_xlabel("QuadSig score"); ax.set_ylabel("ClinVar variants")
ax.set_ylim(0, ax.get_ylim()[1] * 1.45)
ax.legend(frameon=False, loc="upper left", ncol=1, handlelength=1.2, borderaxespad=0.2,
          labelspacing=0.2, fontsize=6.0)
ax.set_xticks([0, 0.5, 1.0]); ax.set_xticklabels(["0", "0.5", "1.0"])
ax.set_title(f"$\\bf{{c}}$   CV accuracy {100*cv_acc:.1f}%", loc="left", fontsize=7, pad=24)
note(f"  c: ClinVar benign n = {nb}, pathogenic n = {npat}, threshold {hi:.3f}")
bad3 = layout_check(fig, "fig3")
for e in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"fig3.{e}"))
plt.close(fig)

open(os.path.join(OUT, "figure_numbers.txt"), "w").write("\n".join(log) + "\n")
print(f"\nwrote {OUT}/fig2.pdf, fig2.png, fig3.pdf, fig3.png, figure_numbers.txt"
      f"  ({'OK' if bad2 + bad3 == 0 else 'LAYOUT PROBLEMS - see above'}; font {FONT[0]})")
