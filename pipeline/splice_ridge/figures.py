"""Figures for the isoform-usage ridge: coefficient forest, predicted-vs-measured, ranking bars.
Palette: the dataviz default categorical slots (blue/orange/aqua/yellow/magenta/green/violet/red)
on the light chart surface #fcfcfb. Each figure writes its own CSV alongside the PNG/PDF."""
import os, csv, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
W = os.path.join(HERE, "work"); FIG = os.path.join(ROOT, "data", "figures")
os.makedirs(FIG, exist_ok=True)
SURF = "#fcfcfb"; INK = "#1f1f1e"; MUTED = "#6b6b68"; GRID = "#e4e4e0"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})

def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig); print("wrote", os.path.join(FIG, stem + ".png"))

# ---- 1. coefficient forest -------------------------------------------------
rows = list(csv.DictReader(open(os.path.join(W, "coefs.csv"))))
rows = [r for r in rows if abs(float(r["coef"])) > 1e-6][:18][::-1]
fig, ax = plt.subplots(figsize=(7.2, 0.33 * len(rows) + 1.0))
for i, r in enumerate(rows):
    c, lo, hi = float(r["coef"]), float(r["ci_lo"]), float(r["ci_hi"])
    sig = (lo > 0) or (hi < 0)
    col = C[0] if c > 0 else C[1]
    ax.plot([lo, hi], [i, i], color=col, lw=2, alpha=1.0 if sig else 0.35,
            solid_capstyle="round")
    ax.plot([c], [i], "o", ms=8, color=col, mec=SURF, mew=2, alpha=1.0 if sig else 0.4)
ax.axvline(0, color=MUTED, lw=1, ls=(0, (3, 3)))
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r["label"] for r in rows])
ax.set_xlabel("ridge coefficient (usage fraction per 1 SD of feature)")
ax.set_title("What predicts a transcript's share of its gene\n"
             "ridge on standardised features; bars = 95% CI from 200 bootstraps over genes; "
             "faded = CI spans 0", loc="left", fontsize=9.5, color=INK)
ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
h = [plt.Line2D([], [], color=C[0], lw=2, marker="o", ms=7, mec=SURF, mew=1.5, label="raises usage"),
     plt.Line2D([], [], color=C[1], lw=2, marker="o", ms=7, mec=SURF, mew=1.5, label="lowers usage")]
ax.legend(handles=h, frameon=False, loc="lower right", fontsize=8)
save(fig, "splice_ridge_coefficients")
import shutil; shutil.copy(os.path.join(W, "coefs.csv"), os.path.join(FIG, "splice_ridge_coefficients.csv"))

# ---- 2. predicted vs measured ---------------------------------------------
pr = list(csv.DictReader(open(os.path.join(W, "preds.tsv")), delimiter="\t"))
m = np.array([float(r["measured"]) for r in pr]); p = np.array([float(r["ridge"]) for r in pr])
core = json.load(open(os.path.join(W, "splice_ridge_core.json")))["models"]["usage"]
fig, ax = plt.subplots(figsize=(5.0, 5.0))
ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
ax.scatter(m, p, s=7, color=C[0], alpha=0.16, linewidths=0, zorder=2, rasterized=True)
b = np.linspace(0, 1, 11); c = np.clip(np.digitize(m, b) - 1, 0, 9)
bx = [m[c == k].mean() for k in range(10) if (c == k).sum() > 20]
by = [p[c == k].mean() for k in range(10) if (c == k).sum() > 20]
ax.plot(bx, by, "-o", color=C[1], lw=2, ms=7, mec=SURF, mew=1.5, zorder=3, label="binned mean")
ax.set_xlabel("measured usage fraction (Isopedia long-read panel)")
ax.set_ylabel("ridge prediction (held-out, gene-grouped CV)")
ax.set_title(f"Held-out isoform usage\nr$^2$={core['r2_cv']:.3f}  Pearson r={core['pearson']:.3f}  "
             f"Spearman={core['spearman']:.3f}  n={core['n']} transcripts / {core['n_genes']} genes",
             loc="left", fontsize=9.5)
ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.1, 1.03)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
save(fig, "splice_ridge_scatter")
with open(os.path.join(FIG, "splice_ridge_scatter.csv"), "w") as fh:
    fh.write("transcript_id,gene_name,measured,ridge_heldout\n")
    for r in pr: fh.write(f"{r['transcript_id']},{r['gene_name']},{r['measured']},{r['ridge']}\n")

# ---- 3. within-gene ranking ------------------------------------------------
wg = {r["method"]: r for r in csv.DictReader(open(os.path.join(W, "within_gene.csv")))}
LAB = {"ridge": "ridge (held-out)", "ridge_logit_softmax": "ridge on logit, softmax per gene",
       "ridge_splice_only": "ridge, splice scores only",
       "uniform": "uniform fractions", "longest_wins": "longest transcript wins",
       "dominant_tss_only": "dominant TSS only", "confound_tx_length": "CONFOUND: transcript length",
       "confound_n_exons": "CONFOUND: exon count", "confound_read_support": "CONFOUND: read support"}
order = ["ridge_logit_softmax", "ridge", "ridge_splice_only", "longest_wins", "dominant_tss_only", "uniform",
         "confound_tx_length", "confound_n_exons", "confound_read_support"]
order = [k for k in order if k in wg]
vals = [float(wg[k]["accuracy"]) for k in order]
ses = [float(wg[k]["se"]) for k in order]
cols = [C[0] if not k.startswith("confound") else C[1] for k in order]
fig, ax = plt.subplots(figsize=(7.4, 0.42 * len(order) + 1.4))
yy = np.arange(len(order))[::-1]
ax.barh(yy, vals, height=0.62, color=cols, zorder=2, edgecolor=SURF, linewidth=2)
ax.errorbar(vals, yy, xerr=np.array(ses) * 1.96, fmt="none", ecolor=INK, elinewidth=1.2, capsize=3, zorder=3)
ax.axvline(0.5, color=MUTED, lw=1.4, ls=(0, (3, 3)), zorder=4)
ax.text(0.5, len(order) - 0.35, " chance (0.5)", color=MUTED, fontsize=8, va="center")
for v, yv in zip(vals, yy):
    ax.text(v + 0.006, yv, f"{v:.3f}", va="center", fontsize=8.5, color=INK)
ax.set_yticks(yy); ax.set_yticklabels([LAB.get(k, k) for k in order])
ax.set_xlim(0.4, max(vals) + 0.09); ax.set_xlabel("fraction of within-gene isoform pairs ordered correctly")
ax.set_title("Within-gene ranking: which isoform of a gene is used more?\n"
             f"{wg[order[0]]['n_pairs']} isoform pairs; error bars = 95% CI; "
             "orange = trivial confound, not a model", loc="left", fontsize=9.5)
ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
save(fig, "splice_ridge_within_gene")
shutil.copy(os.path.join(W, "within_gene.csv"), os.path.join(FIG, "splice_ridge_within_gene.csv"))

# ---- 4. the 12 app genes: mechanistic model vs ridge vs baselines ----------
a12 = os.path.join(W, "app12.json")
if os.path.exists(a12):
    A = json.load(open(a12))
    L = {"mech": "mechanistic model (app defaults)",
         "ridge_compact": "ridge, like-for-like features (genes held out)",
         "ridge": "ridge, all features (genes held out)",
         "mes_donor_min": "MaxEntScan weakest donor", "mes_acc_min": "MaxEntScan weakest acceptor",
         "mes_donor_mean": "MaxEntScan mean donor", "mes_acc_mean": "MaxEntScan mean acceptor",
         "cage": "CAGE promoter weight", "is_dom_cage": "uses the dominant CAGE promoter",
         "uniform": "uniform fractions", "tx_len": "CONFOUND: transcript length",
         "n_exons": "CONFOUND: exon count", "read_support": "CONFOUND: read support"}
    ks = [k for k in L if k in A["metrics"] and math.isfinite(A["metrics"][k]["within_gene_acc"])]
    vs = [A["metrics"][k]["within_gene_acc"] for k in ks]
    ns = [A["metrics"][k]["n_pairs"] for k in ks]
    cols = [C[1] if k.startswith(("tx_", "n_ex", "read_")) else (C[2] if k == "mech" else C[0]) for k in ks]
    fig, ax = plt.subplots(figsize=(7.6, 0.42 * len(ks) + 1.5))
    yy = np.arange(len(ks))[::-1]
    ax.barh(yy, vs, height=0.62, color=cols, edgecolor=SURF, linewidth=2, zorder=2)
    ax.axvline(0.5, color=MUTED, lw=1.4, ls=(0, (3, 3)), zorder=4)
    for v, yv, n in zip(vs, yy, ns):
        ax.text(v + 0.008, yv, f"{v:.2f}  (n={n})", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(yy); ax.set_yticklabels([L[k] for k in ks])
    ax.set_xlim(0, 1.18); ax.set_xlabel("fraction of within-gene isoform pairs ordered correctly")
    ax.set_title(f"The 12 app genes: the mechanistic model can be run here\n"
                 f"{A['metrics']['mech']['n']} transcripts of {len({r['gene'] for r in A['rows']})} genes; "
                 "aqua = mechanistic, orange = trivial confound; dashed line = chance",
                 loc="left", fontsize=9.5)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    save(fig, "splice_ridge_app12")
    with open(os.path.join(FIG, "splice_ridge_app12.csv"), "w") as fh:
        fh.write("method,within_gene_acc,n_pairs,r2,pearson,spearman,n\n")
        for k in ks:
            m = A["metrics"][k]
            fh.write(f'{k},{m["within_gene_acc"]:.6g},{m["n_pairs"]},{m["r2"]:.6g},'
                     f'{m["pearson"]:.6g},{m["spearman"]:.6g},{m["n"]}\n')
