#!/usr/bin/env python
"""Figures for the Isopedia isoform-usage validation.
  data/figures/isopedia_scatter.png   predicted vs measured, model and promoter-only baseline
  data/figures/isopedia_bars.png      per-isoform measured / model / promoter-only, grouped by gene
CSVs behind them: isopedia_validation.csv (both figures), isopedia_validation_summary.csv, isopedia_per_gene.csv
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
F = os.path.join(ROOT, "data", "figures")
SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#b4b2ab"
C = {"measured_n": "#2a78d6", "pred_model_n": "#eb6834", "pred_promoter_only_n": "#1baf7a"}
plt.rcParams.update({"font.size": 8, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6,
                     "font.family": "DejaVu Sans"})
d = pd.read_csv(os.path.join(F, "isopedia_validation.csv"))
S = pd.read_csv(os.path.join(F, "isopedia_validation_summary.csv")).set_index("predictor")

# ---------------- scatter
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.9), constrained_layout=True)
for ax, key, name in zip(axes, ["pred_model_n", "pred_promoter_only_n"],
                         ["kinetic splicing model", "CAGE promoter weights only (no splicing)"]):
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls="--", zorder=1, label="_")
    ax.scatter(d["measured_n"], d[key], s=34, color=C[key], edgecolor=SURF, linewidth=2, zorder=3)
    for r in d.itertuples():   # label only the points that carry the story; the CSV has them all
        if max(r.measured_n, getattr(r, key)) < 0.2:
            continue
        ax.annotate(r.gene, (r.measured_n, getattr(r, key)), textcoords="offset points",
                    xytext=(5, 3), fontsize=6, color=INK2)
    k = "model" if key == "pred_model_n" else "promoter_only"
    ax.set_title(f"{name}\nPearson r = {S.loc[k,'pearson_r']:.2f}   Spearman = {S.loc[k,'spearman_r']:.2f}   "
                 f"MAE = {S.loc[k,'mae']:.2f}\nwithin-gene pair accuracy = {S.loc[k,'pair_accuracy']:.2f} "
                 f"({S.loc[k,'pairs_correct']:.1f}/{int(S.loc[k,'pairs_n'])} pairs, chance = 0.50)",
                 fontsize=7.5, color=INK, loc="left", linespacing=1.6)
    ax.set_xlabel("measured isoform usage (Isopedia, renormalised within gene)")
    ax.set_ylabel("predicted fraction (renormalised)")
    ax.set_xlim(-0.03, 1.08); ax.set_ylim(-0.05, 1.08)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis="both", color=MUTED, lw=0.4, alpha=0.5); ax.set_axisbelow(True)
fig.suptitle(f"Kinetic Gene predicted vs Isopedia measured isoform usage — "
             f"{len(d)} isoforms, {d.gene.nunique()} genes, dashed line = perfect",
             fontsize=9.5, color=INK, x=0.012, ha="left")
fig.savefig(os.path.join(F, "isopedia_scatter.png"), dpi=200)
plt.close(fig)

# ---------------- grouped bars, one group per compared isoform, blocked by gene
d2 = d.sort_values(["gene", "measured_n"], ascending=[True, False]).reset_index(drop=True)
lab = [f"{g}  {n.split(';')[0]}" for g, n in zip(d2.gene, d2.app_names)]
y = np.arange(len(d2))[::-1].astype(float)
h = 0.26
fig, ax = plt.subplots(figsize=(7.6, 0.265 * len(d2) + 1.5), constrained_layout=True)
series = [("measured_n", "measured (Isopedia)"), ("pred_model_n", "kinetic splicing model"),
          ("pred_promoter_only_n", "CAGE promoter weights only")]
for i, (key, name) in enumerate(series):
    ax.barh(y + (1 - i) * h, d2[key], height=h - 0.04,
            color=C[key], label=name, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=6.5)
ax.set_xlabel("isoform fraction, renormalised within gene over the compared isoforms\n(the measured share of ALL gene reads is in isopedia_validation.csv / isopedia_per_gene.csv)")
ax.set_xlim(0, 1.02)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.grid(axis="x", color=MUTED, lw=0.4, alpha=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.012), ncol=3)
ax.set_title("Per-isoform comparison, grouped by gene", fontsize=9.5, color=INK, loc="left")
fig.savefig(os.path.join(F, "isopedia_bars.png"), dpi=200)
plt.close(fig)
print("wrote", os.path.join(F, "isopedia_scatter.png"), "and isopedia_bars.png")
