"""Figures for data/discovery.json -> data/figures/discovery_*.png (+ discovery_*.json with plotted numbers).
Run with any python that has matplotlib (utrrace has none; e.g. ~/miniconda3/envs/fungi/bin/python -I).
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = os.path.join(ROOT, "pipeline/discovery/work")
FIG = os.path.join(ROOT, "data/figures")
os.makedirs(FIG, exist_ok=True)
doc = json.load(open(os.path.join(ROOT, "data/discovery.json")))
iso = json.load(open(os.path.join(W, "isoforms.json")))
mo = json.load(open(os.path.join(W, "model_out.json")))["iso"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
GREY, BLUE, ORANGE, RED = "#b8b8b8", "#2a6fb0", "#e08a1e", "#c0392b"

# 1. all pairs: calm vs stress log2(A/B) at default params; robust calls coloured
xs, ys = [], []
for g in iso["groups"]:
    I = g["isoforms"]
    for a in range(len(I)):
        for b in range(a + 1, len(I)):
            ca, cb = mo[I[a]["id"]]["calm"]["global"][0], mo[I[b]["id"]]["calm"]["global"][0]
            sa, sb = mo[I[a]["id"]]["stress"]["global"][0], mo[I[b]["id"]]["stress"]["global"][0]
            if min(ca, cb, sa, sb) <= 1e-9:
                continue
            xs.append(math.log2(ca / cb)); ys.append(math.log2(sa / sb))
fig, ax = plt.subplots(figsize=(4.6, 4.2))
ax.scatter(xs, ys, s=3, c=GREY, lw=0, alpha=.5, label=f"all same-CDS pairs (n={len(xs)})", rasterized=True)
for k, col, lab in (("stress_flip", RED, "robust ISR flip"), ("isr_escape", ORANGE, "robust ISR escape")):
    pts = [(c["default"]["log2_calm_A_over_B"], c["default"]["log2_stress_A_over_B"]) for c in doc["calls"][k]]
    if pts:
        p = np.array(pts)
        ax.scatter(p[:, 0], p[:, 1], s=10, c=col, lw=0, label=f"{lab} (n={doc['meta']['counts'][k]})")
lim = max(1, np.percentile(np.abs(xs + ys), 99.5)) if xs else 1
ax.plot([-lim, lim], [-lim, lim], c="k", lw=.5, ls="--")
ax.axhline(0, c="k", lw=.4); ax.axvline(0, c="k", lw=.4)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("log2 output A/B, no stress"); ax.set_ylabel("log2 output A/B, stress (TC = 0.15)")
ax.set_title("Same protein, different 5' leader: stress reshuffles leaders", fontsize=9)
ax.legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "discovery_pairs_calm_vs_stress.png"), dpi=200); plt.close(fig)
json.dump({"x_log2_calm": [round(v, 3) for v in xs], "y_log2_stress": [round(v, 3) for v in ys]},
          open(os.path.join(FIG, "discovery_pairs_calm_vs_stress.json"), "w"))

# 2. validation vs Cao 2021 TE
val = doc["validation_cao2021"]
fig, axs = plt.subplots(1, 3, figsize=(9, 3.1))
for ax, (cell, v) in zip(axs, val.items()):
    p = np.array(v["points"]) if v["points"] else np.zeros((0, 2))
    if len(p):
        ax.scatter(p[:, 0], p[:, 1], s=8, c=BLUE, lw=0, alpha=.7)
    ax.axhline(0, c="k", lw=.4); ax.axvline(0, c="k", lw=.4)
    rho = v["spearman_log2"]
    ax.set_title(f"{cell}: n={v['pairs']}, sign agree {v['frac']}, rho={rho}", fontsize=8)
    ax.set_xlabel("model log2 A/B (no stress)")
axs[0].set_ylabel("measured log2 TE A/B (Cao 2021)")
fig.suptitle("Check: model vs measured within-gene TE (unambiguous, same CDS)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "discovery_validation_cao2021.png"), dpi=200); plt.close(fig)

# 3. top robust examples: output calm vs stress for both leaders, bars = default, whiskers = range over points x folds
ex = []
byname = {i["name"]: i for g in iso["groups"] for i in g["isoforms"]}
if "ATF4-206" in byname and "ATF4-202" in byname:  # MANE vs alternative leader (TSL1), same protein
    ex.append(("MANE vs alt leader", {"gene": "ATF4", "A": byname["ATF4-206"], "B": byname["ATF4-202"]}))
for k in ("stress_flip", "isr_escape"):
    hc = [c for c in doc["calls"][k] if c["both_high_conf"]] or doc["calls"][k]
    ex += [(k, c) for c in hc[:3]]
if ex:
    fig, axs = plt.subplots(1, len(ex), figsize=(2.3 * len(ex), 3.0), squeeze=False)
    dump = []
    for ax, (k, c) in zip(axs[0], ex):
        vals = []
        for side in ("A", "B"):
            r = mo[c[side]["id"]]
            for cond in ("calm", "stress"):
                a = np.array([r[cond]["global"], r[cond]["local"]])
                vals.append((f"{c[side]['name']}\n{cond}", a[0, 0], a.min(), a.max()))
        xs_ = np.arange(4)
        cols = [BLUE, RED, BLUE, RED]
        ax.bar(xs_, [v[1] for v in vals], color=cols)
        for x, v in zip(xs_, vals):
            ax.plot([x, x], [v[2], v[3]], c="k", lw=.8)
        ax.set_xticks(xs_); ax.set_xticklabels([v[0] for v in vals], fontsize=6, rotation=45, ha="right")
        ax.set_title(f"{c['gene']} ({k.replace('_', ' ')})", fontsize=8)
        dump.append({"gene": c["gene"], "call": k, "bars": [{"label": v[0], "default": v[1], "min": v[2], "max": v[3]} for v in vals]})
    axs[0][0].set_ylabel("protein / min (model)\nblue = no stress, red = stress")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "discovery_examples.png"), dpi=200); plt.close(fig)
    json.dump(dump, open(os.path.join(FIG, "discovery_examples.json"), "w"), indent=1)
print("figures written")
