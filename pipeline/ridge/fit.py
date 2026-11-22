#!/usr/bin/env python
"""Fit interpretable ridge models (motifs.features) for MRL (Sample 2019 MPRA) and TE (Cao 2021 HEK293T),
compare with Optimus 5-Prime and the app's physics model on the same CV folds, write data/ridge.json + figures.

All comparison models are calibrated the same way: inside each CV fold, a 1-D linear map from the model's score to
the target is fit on the training folds and applied to the held-out fold (Optimus predicts MRL on the MPRA scale, the
physics model predicts initiations/min; neither is on the target's scale). r2 = 1 - SSE/SST on pooled held-out
predictions.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Ridge, RidgeCV
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.dirname(HERE)
ROOT = os.path.dirname(PIPE)
sys.path.insert(0, PIPE)
sys.path.insert(0, HERE)
import motifs  # noqa: E402
import feats  # noqa: E402

WORK = os.path.join(HERE, "work")
FIG = os.path.join(ROOT, "data", "figures")
OUT = os.path.join(ROOT, "data", "ridge.json")
ALPHAS = np.logspace(-3, 4, 29)
NBOOT = int(os.environ.get("NBOOT", "200"))
SEED = 0
# n_uAUG is dropped: it equals n_uAUG_strong + n_uAUG_adequate + n_uAUG_weak exactly (ridge would split its weight)
FEATS = [f for f in motifs.FEATURES if f != "n_uAUG"]


def r2(y, p):
    return float(1 - np.sum((y - p) ** 2) / np.sum((y - np.mean(y)) ** 2))


def usable(df):
    return [f for f in FEATS if f in df and df[f].notna().all() and df[f].std() > 0]


def folds(df, groups=None):
    if groups is None:
        return list(KFold(5, shuffle=True, random_state=SEED).split(df))
    return list(GroupKFold(5, shuffle=True, random_state=SEED).split(df, groups=groups))


def cv_ridge(X, y, fl):
    oof, alphas = np.zeros(len(y)), []
    for tr, te in fl:
        sc = StandardScaler().fit(X[tr])
        m = RidgeCV(alphas=ALPHAS).fit(sc.transform(X[tr]), y[tr])
        oof[te] = m.predict(sc.transform(X[te]))
        alphas.append(float(m.alpha_))
    return oof, alphas


def cv_calib(s, y, fl):
    oof = np.zeros(len(y))
    for tr, te in fl:
        m = LinearRegression().fit(s[tr, None], y[tr])
        oof[te] = m.predict(s[te, None])
    return oof


def boot_ci(X, y, alpha, groups=None, nboot=NBOOT):
    rng = np.random.default_rng(SEED)
    B = []
    if groups is not None:
        ug = np.unique(groups)
        idx_of = pd.Series(np.arange(len(groups))).groupby(groups).apply(np.array).to_dict()
    for _ in range(nboot):
        if groups is None:
            ii = rng.integers(0, len(y), len(y))
        else:
            ii = np.concatenate([idx_of[g] for g in rng.choice(ug, len(ug))])
        sc = StandardScaler().fit(X[ii])
        B.append(Ridge(alpha=alpha).fit(sc.transform(X[ii]), y[ii]).coef_)
    B = np.array(B)
    return np.percentile(B, 2.5, axis=0), np.percentile(B, 97.5, axis=0)


def within_gene(df, pred_cols, ycol="log2_te", min_diff=0.0):
    """All same-gene pairs with distinct 5' UTRs and |dy| > min_diff; fraction with sign(dpred) == sign(dy)."""
    res = {c: [0, 0, 0] for c in pred_cols}  # correct, wrong, tied
    npairs, ngenes = 0, 0
    for _, g in df.groupby("external_gene_id"):
        if len(g) < 2:
            continue
        v = g[[ycol, "utr"] + pred_cols].to_numpy()
        used = False
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                dy = v[i][0] - v[j][0]
                if v[i][1] == v[j][1] or abs(dy) <= min_diff:
                    continue
                npairs += 1
                used = True
                for k, c in enumerate(pred_cols):
                    dp = v[i][2 + k] - v[j][2 + k]
                    if abs(dp) < 1e-12:
                        res[c][2] += 1
                    else:
                        res[c][0 if np.sign(dp) == np.sign(dy) else 1] += 1
        ngenes += used
    out = {"n_pairs": npairs, "n_genes": ngenes, "min_abs_dlog2te": min_diff, "models": {}}
    for c, (a, b, t) in res.items():
        out["models"][c] = {"correct": a, "wrong": b, "tied": t,
                            "frac_correct_untied": a / (a + b) if a + b else None,
                            "frac_correct_ties_half": (a + 0.5 * t) / npairs if npairs else None}
    return out


def fit_dataset(df, y, groups, name):
    fs = usable(df)
    X = df[fs].to_numpy(float)
    fl = folds(df, groups)
    oof, al = cv_ridge(X, y, fl)
    res = {"n": int(len(y)), "features_used": fs, "fold_alphas": al, "r2_cv": r2(y, oof),
           "pearson_r_cv": float(np.corrcoef(y, oof)[0, 1])}
    preds = {"ridge": oof}
    preds["optimus"] = cv_calib(df["optimus"].to_numpy(float), y, fl)
    res["r2_optimus_same_split"] = r2(y, preds["optimus"])
    res["pearson_r_optimus"] = float(np.corrcoef(y, df["optimus"])[0, 1])
    if "phys_logperMin" in df:
        ph = df["phys_logperMin"].to_numpy(float)
        preds["physics"] = cv_calib(ph, y, fl)
        res["r2_physics_same_split"] = r2(y, preds["physics"])
        res["pearson_r_physics"] = float(np.corrcoef(y, ph)[0, 1])
        o2, _ = cv_ridge(np.column_stack([X, ph]), y, fl)
        preds["ridge+physics"] = o2
        res["r2_ridge_plus_physics"] = r2(y, o2)
    # ridge + optimus (does the interpretable set add anything on top of the CNN?)
    o3, _ = cv_ridge(np.column_stack([X, df["optimus"].to_numpy(float)]), y, fl)
    res["r2_ridge_plus_optimus"] = r2(y, o3)
    # final fit on all rows
    sc = StandardScaler().fit(X)
    alpha = float(np.median(al))
    m = Ridge(alpha=alpha).fit(sc.transform(X), y)
    lo, hi = boot_ci(X, y, alpha, groups)
    res["alpha"] = alpha
    res["intercept"] = float(m.intercept_)
    res["features"] = [{"name": f, "label": motifs.FEATURES[f], "coef": float(c), "ci": [float(a), float(b)]}
                       for f, c, a, b in zip(fs, m.coef_, lo, hi)]
    res["scaler"] = {"mean": dict(zip(fs, map(float, sc.mean_))), "std": dict(zip(fs, map(float, sc.scale_)))}
    res["coef_units"] = "change in target per +1 SD of the feature (standardized ridge coefficient)"
    res["ci_method"] = f"{NBOOT} bootstrap resamples ({'genes' if groups is not None else 'rows'}), fixed alpha"
    return res, preds


# ----------------------------------------------------------------------------- figures
def fig_forest(res, title, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    F = sorted(res["features"], key=lambda f: f["coef"])
    fig, ax = plt.subplots(figsize=(9.5, 0.32 * len(F) + 1.6))
    y = np.arange(len(F))
    c = ["#c0392b" if f["ci"][1] < 0 else "#2471a3" if f["ci"][0] > 0 else "#888888" for f in F]
    ax.hlines(y, [f["ci"][0] for f in F], [f["ci"][1] for f in F], color=c, lw=2)
    ax.scatter([f["coef"] for f in F], y, color=c, zorder=3, s=22)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y, [f["label"] for f in F], fontsize=8)
    ax.set_xlabel("ridge coefficient: target change per +1 SD of feature\n(95% bootstrap CI; grey = CI crosses 0)", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_scatter(y, preds, res, ylabel, title, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = [k for k in ["ridge", "optimus", "physics", "physics_fitted"] if k in preds]
    fig, axs = plt.subplots(1, len(names), figsize=(4 * len(names), 4), sharex=True, sharey=True)
    for ax, k in zip(np.atleast_1d(axs), names):
        ax.hexbin(preds[k], y, gridsize=60, bins="log", cmap="viridis", mincnt=1)
        lo, hi = np.percentile(np.concatenate([y, preds[k]]), [0.5, 99.5])
        ax.plot([lo, hi], [lo, hi], "w--", lw=1)
        ax.set_title(f"{k}  held-out r² = {r2(y, preds[k]):.3f}", fontsize=10)
        ax.set_xlabel(f"predicted {ylabel} (5-fold CV)")
    np.atleast_1d(axs)[0].set_ylabel(f"observed {ylabel}")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_within(wg_list, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(wg_list), figsize=(5.2 * len(wg_list), 4.4), sharey=True)
    for ax, wg in zip(np.atleast_1d(axs), wg_list):
        ks = list(wg["models"])
        v = [wg["models"][k]["frac_correct_untied"] or 0 for k in ks]
        n = [wg["models"][k]["correct"] + wg["models"][k]["wrong"] for k in ks]
        lab = [k.replace("pred_", "").replace("ctrl_", "control: ").replace("_", " ") for k in ks]
        bars = ax.bar(lab, v, color=["#2471a3", "#7d3c98", "#d68910", "#17a589", "#888888"][:len(ks)])
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
        for b, vv, nn in zip(bars, v, n):
            ax.text(b.get_x() + b.get_width() / 2, vv + 0.01, f"{vv:.2f}\n(n={nn})", ha="center", fontsize=8)
        ax.axhline(0.5, color="k", ls="--", lw=0.8)
        ax.set_ylim(0, 1)
        ax.set_title(f"pairs |Δlog2 TE| > {wg['min_abs_dlog2te']}: {wg['n_pairs']} pairs, {wg['n_genes']} genes",
                     fontsize=9)
    np.atleast_1d(axs)[0].set_ylabel("fraction of isoform pairs ordered correctly")
    fig.suptitle("Within-gene: which 5' UTR isoform has higher TE? (held-out genes, ties excluded)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------- app isoforms
def isoform_cds():
    """CDS (+3'UTR) for the app's isoforms from GENCODE v47 pc_transcripts (unversioned ENST match)."""
    import gzip
    want = {}
    d = json.load(open(os.path.join(ROOT, "data", "utr_data.json")))
    for g in d["genes"]:
        for i in g["isoforms"]:
            want[i["id"].split(".")[0]] = i
    cds = {}
    fa = os.path.join(HERE, "downloads", "gencode.v47.pc_transcripts.fa.gz")
    with gzip.open(fa, "rt") as fh:
        cur, seq = None, []
        for line in fh:
            if line[0] == ">":
                if cur:
                    cds[cur[0]] = (cur[1], "".join(seq))
                h = line[1:].strip().split("|")
                tid = h[0].split(".")[0]
                cur = (tid, h) if tid in want else None
                seq = []
            elif cur:
                seq.append(line.strip())
        if cur:
            cds[cur[0]] = (cur[1], "".join(seq))
    # build_data's cached Ensembl REST cDNA (mask_feature=1: lowercase 5' UTR, uppercase CDS+) is the same source
    # the app's isoforms came from -> preferred
    import glob
    cache = {}
    for f in glob.glob(os.path.join(PIPE, "cache", "*.json")):
        try:
            d = json.load(open(f))
        except Exception:  # noqa: BLE001
            continue
        if isinstance(d, list):
            for r in d:
                if isinstance(r, dict) and "seq" in r and "id" in r:
                    cache[r["id"].split(".")[0]] = r["seq"]
    out = {}
    for tid, iso in want.items():
        s = cache.get(tid)
        if s:
            k = len(s) - len(s.lstrip("acgtn"))
            u, c = s[:k].upper().replace("T", "U"), s[k:].upper().replace("T", "U")
            if u == iso["utr5"] and c:
                out[iso["id"]] = (iso, c, "ensembl_rest_cache")
                continue
        c = None
        if tid in cds:
            h, s = cds[tid]
            rng = {x.split(":")[0]: x.split(":")[1] for x in h if ":" in x}
            if "CDS" in rng:
                c = s[int(rng["CDS"].split("-")[0]) - 1:].replace("T", "U")
                u = s[:int(rng["CDS"].split("-")[0]) - 1].replace("T", "U")
                if u != iso["utr5"]:
                    c = None  # model differs from the app's isoform; fall back below
        out[iso["id"]] = (iso, c if c else iso["cds_context"][6:], "gencode_v47" if c else "cds_context_only")
    return out


def per_isoform(models):
    iso = isoform_cds()
    ids = list(iso)
    F = feats.many([iso[i][0]["utr5"] for i in ids], [iso[i][1] for i in ids], 4, chunksize=1)
    out = {}
    for i, f in zip(ids, F):
        rec = {"cds_source": iso[i][2]}
        for key, m in models.items():
            mu, sd = m["scaler"]["mean"], m["scaler"]["std"]
            contrib, pred = [], m["intercept"]
            for ft in m["features"]:
                n = ft["name"]
                v = f.get(n)
                if v is None:
                    v = mu[n]  # feature not computable -> mean (contribution 0)
                z = (v - mu[n]) / sd[n]
                contrib.append({"name": n, "label": ft["label"], "value": round(float(v), 4),
                                "contribution": round(ft["coef"] * z, 4)})
                pred += ft["coef"] * z
            rec[f"{key}_pred"] = round(float(pred), 4)
            rec[f"contrib_{key}"] = contrib
        rec["contrib"] = rec["contrib_mrl"]
        out[i] = rec
    return out


def main():
    os.makedirs(FIG, exist_ok=True)
    models, extra = {}, {}
    mrl = pd.read_csv(os.path.join(WORK, "mrl.csv.gz"))
    if os.path.exists(os.path.join(WORK, "mrl_physics.csv.gz")):
        mrl = mrl.merge(pd.read_csv(os.path.join(WORK, "mrl_physics.csv.gz")), on="utr", how="left")
    y = mrl["rl"].to_numpy(float)
    res, preds = fit_dataset(mrl, y, None, "mrl")
    pfo = os.path.join(WORK, "mrl_physics_fit_oof.csv.gz")
    if os.path.exists(pfo):
        preds["physics_fitted"] = pd.read_csv(pfo)["physics_fitted"].to_numpy(float)
        res["r2_physics_fitted_same_split"] = r2(y, preds["physics_fitted"])
    res.update({"dataset": "Sample et al. 2019 GEO GSE114002 / GSM4084997, varying-length (25-100 nt) eGFP "
                           "library, total_reads >= 10: all 'human' 5' UTR fragments + random-library test split "
                           "(top 100 by reads per length)",
                "target": "mean ribosome load (MRL, polysome profiling, HEK293T)", "cv": "5-fold (rows)",
                "cds_head": feats.EGFP_HEAD, "fold_method": feats.FOLD_METHOD,
                "leakage": "Optimus varying-length model was trained on the random-library training split only; "
                           "rows here are the human set and the random test split, so Optimus never saw them. "
                           "Optimus architecture/hyperparameters were however chosen on this library.",
                "applicability": "trained on 25-100 nt UTR fragments in an eGFP reporter; per_isoform mrl_pred for "
                                 "longer endogenous UTRs extrapolates (length / uAUG counts beyond the training range)",
                "physics_note": "physics = app/src/analytic.js expected() at app default params, log(initiations/min); "
                                "physics_fitted = 4 params refit per fold (see physics_fit)",
                "by_subset": {}})
    for s in mrl["subset"].unique():
        k = (mrl["subset"] == s).to_numpy()
        res["by_subset"][s] = {"n": int(k.sum()), "r2_ridge": r2(y[k], preds["ridge"][k]),
                               "r2_optimus": r2(y[k], preds["optimus"][k])}
        if "physics" in preds:
            res["by_subset"][s]["r2_physics"] = r2(y[k], preds["physics"][k])
    models["mrl"] = res
    fig_forest(res, f"MRL (Sample 2019 MPRA, n={res['n']}): ridge on interpretable 5' UTR features",
               os.path.join(FIG, "ridge_coef_mrl.png"))
    fig_scatter(y, preds, res, "MRL", "Sample 2019 MPRA, held-out (5-fold CV)", os.path.join(FIG, "pred_vs_obs_mrl.png"))
    pd.DataFrame({"observed": y, **{k: v for k, v in preds.items()}, "subset": mrl["subset"]}).to_csv(
        os.path.join(FIG, "pred_vs_obs_mrl.csv.gz"), index=False)

    if os.path.exists(os.path.join(WORK, "te.csv.gz")):
        te = pd.read_csv(os.path.join(WORK, "te.csv.gz"))
        if os.path.exists(os.path.join(WORK, "te_physics.csv.gz")):
            te = te.merge(pd.read_csv(os.path.join(WORK, "te_physics.csv.gz")), on="vid", how="left")
        y2 = te["log2_te"].to_numpy(float)
        res2, preds2 = fit_dataset(te, y2, te["external_gene_id"].to_numpy(), "te")
        res2.update({"dataset": "Cao et al. 2021 (Nat Commun) HEK293T transcript TE (Andreev 2015 Ribo-seq/RNA-seq), "
                                "5' UTR + CDS from GENCODE v19 pc_transcripts (exact versioned ID match)",
                     "target": "log2 TE (rpkm_riboseq / rpkm_rnaseq)", "cv": "5-fold grouped by gene",
                     "filters": "rnaseq_counts >= 50, riboseq_counts >= 10, te > 0, CDS starts ATG, 5' UTR >= 10 nt, "
                                "no sibling isoform with bit-identical TE (unresolved quantification)",
                     "n_genes": int(te["external_gene_id"].nunique()), "fold_method": feats.FOLD_METHOD,
                     "optimus_input": "last 100 nt of the 5' UTR",
                     "leakage": "none from Optimus (MPRA-trained). Ridge features/filters were not tuned on this set."})
        for c, p in preds2.items():
            te[f"pred_{c}"] = p
        pc = [f"pred_{c}" for c in preds2]
        # control, not a model: "longer CDS isoform has higher TE" (no 5' UTR information)
        te["ctrl_longer_CDS"] = te["cds_length"].astype(float)
        wc = pc + ["ctrl_longer_CDS"]
        res2["within_gene"] = [within_gene(te, wc, min_diff=0.0), within_gene(te, wc, min_diff=1.0)]
        res2["within_gene_note"] = ("pairs = same-gene isoforms with different 5' UTRs; out-of-fold predictions (gene held "
                                    "out); ctrl_longer_CDS = predicts the isoform with the longer CDS has higher TE "
                                    "(no 5' UTR information; flags isoform-quantification/length effects in the TE table)")
        models["te"] = res2
        fig_forest(res2, f"log2 TE (Cao 2021 HEK293T, n={res2['n']} transcripts, {res2['n_genes']} genes)",
                   os.path.join(FIG, "ridge_coef_te.png"))
        fig_scatter(y2, preds2, res2, "log2 TE", "Cao 2021 HEK293T TE, held-out genes (5-fold grouped CV)",
                    os.path.join(FIG, "pred_vs_obs_te.png"))
        fig_within(res2["within_gene"], os.path.join(FIG, "within_gene_te.png"))
        te[["vid", "external_gene_id", "log2_te"] + pc].to_csv(os.path.join(FIG, "pred_vs_obs_te.csv.gz"), index=False)
        extra["te"] = res2

    out = {"models": models}
    pf = os.path.join(WORK, "physics_fit.json")
    if os.path.exists(pf):
        out["physics_fit"] = json.load(open(pf))
    out["per_isoform"] = per_isoform(models)
    out["figures"] = sorted(f for f in os.listdir(FIG) if f.endswith(".png") and f.startswith(("ridge_coef_", "pred_vs_obs_", "within_gene_")))
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    for k, m in models.items():
        print(k, {x: m.get(x) for x in ["n", "r2_cv", "r2_optimus_same_split", "r2_physics_same_split", "r2_physics_fitted_same_split",
                                       "r2_ridge_plus_physics", "r2_ridge_plus_optimus"]})
        for f in sorted(m["features"], key=lambda f: -abs(f["coef"]))[:8]:
            print(f"   {f['name']:28s} {f['coef']:+.3f} [{f['ci'][0]:+.3f}, {f['ci'][1]:+.3f}]")
        if "within_gene" in m:
            for w in m["within_gene"]:
                print("   within", w["min_abs_dlog2te"], w["n_pairs"], w["n_genes"],
                      {c: round(v["frac_correct_untied"] or 0, 3) for c, v in w["models"].items()})


if __name__ == "__main__":
    main()
