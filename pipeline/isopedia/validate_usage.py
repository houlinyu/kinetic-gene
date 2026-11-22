#!/usr/bin/env python
"""Kinetic-model isoform fractions vs MEASURED Isopedia long-read isoform usage.

Inputs:  data/isoform_usage.json              (built by build_isoform_usage.py)
         data/figures/isopedia_pred.json      (node pipeline/isopedia/predict.mjs)
         data/figures/isopedia_join.csv       (app ENST -> Isopedia catalog isoform)
Outputs: data/figures/isopedia_validation.csv (one row per compared isoform)
         data/figures/isopedia_validation_summary.csv
         data/figures/isopedia_scatter.png, isopedia_bars.png

A comparison unit is one Isopedia catalog isoform (a splice pattern). If two app transcripts collapse
onto the same splice pattern, their predicted fractions are summed, because Isopedia cannot tell them
apart. Measured and predicted are each renormalised over the matched set within a gene, so the metric
is "of the mass the two resources share, is it distributed the same way".
"""
import json, os, sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
F = os.path.join(ROOT, "data", "figures")
MODELS = ["model", "promoter_only", "uniform"]
LABEL = {"model": "kinetic splicing model", "promoter_only": "CAGE promoter weights only (no splicing)",
         "uniform": "uniform over annotated isoforms"}


def main():
    usage = json.load(open(os.path.join(ROOT, "data", "isoform_usage.json")))
    pred = json.load(open(os.path.join(F, "isopedia_pred.json")))
    join = pd.read_csv(os.path.join(F, "isopedia_join.csv"))
    join = join[join.isoform_uid.notna()]

    rows = []
    for gene, gu in usage["genes"].items():
        if not gu.get("measured"):
            continue
        gj = join[join.gene == gene]
        P = pred[gene]
        for uid, sub in gj.groupby("isoform_uid"):
            bases = sorted(set(sub.app_base))
            meas = [gu["overall"][b] for b in bases if b in gu["overall"]]
            if not meas:
                continue
            r = dict(gene=gene, isoform_uid=uid, app_ids=";".join(sub.app_id),
                     app_names=";".join(sub.app_name.astype(str)), n_app=len(sub),
                     measured=float(np.mean(meas)))
            for m in MODELS:
                r["pred_" + m] = float(sum(P[m].get(i, 0.0) for i in sub.app_id))
            rows.append(r)
    d = pd.DataFrame(rows).sort_values(["gene", "measured"], ascending=[True, False])

    # renormalise within gene over the matched set
    for col in ["measured"] + ["pred_" + m for m in MODELS]:
        d[col + "_n"] = d.groupby("gene")[col].transform(lambda s: s / s.sum() if s.sum() > 0 else np.nan)
    d.to_csv(os.path.join(F, "isopedia_validation.csv"), index=False)

    def pairs(df, col):
        """within-gene pairwise ranking accuracy: does the model order the two isoforms correctly"""
        ok = tot = ties = 0
        for g, s in df.groupby("gene"):
            v = s[["measured", col]].values
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    if v[i][0] == v[j][0]:
                        continue
                    tot += 1
                    dp = v[i][1] - v[j][1]
                    if dp == 0:
                        ties += 1
                        ok += 0.5
                    elif (dp > 0) == (v[i][0] > v[j][0]):
                        ok += 1
        return ok, tot, ties

    summ = []
    for m in MODELS:
        c = "pred_" + m
        sub = d.dropna(subset=["measured_n", c + "_n"])
        pr = stats.pearsonr(sub["measured_n"], sub[c + "_n"])
        sp = stats.spearmanr(sub["measured_n"], sub[c + "_n"])
        ok, tot, ties = pairs(d, c)
        # raw (no renormalisation): predicted fractions are over app isoforms only
        prr = stats.pearsonr(d["measured"], d[c])
        summ.append(dict(predictor=m, label=LABEL[m], n_isoforms=len(sub), n_genes=d.gene.nunique(),
                         pearson_r=pr[0], pearson_p=pr[1], spearman_r=sp.statistic, spearman_p=sp.pvalue,
                         mae=float(np.abs(sub["measured_n"] - sub[c + "_n"]).mean()),
                         pairs_n=tot, pairs_correct=ok, pair_accuracy=(ok / tot if tot else np.nan),
                         pair_ties=ties,
                         pearson_r_unnormalised=prr[0]))
    S = pd.DataFrame(summ)
    S.to_csv(os.path.join(F, "isopedia_validation_summary.csv"), index=False)

    # per gene, for the model
    per = []
    for g, s in d.groupby("gene"):
        row = dict(gene=g, n_isoforms=len(s),
                   measured_share_of_gene=float(usage["genes"][g]["app_isoform_share_of_gene"]),
                   n_samples=usage["genes"][g]["n_samples"])
        for m in MODELS:
            c = "pred_" + m
            ok, tot, _ = pairs(s, c)
            row[m + "_pair_acc"] = ok / tot if tot else np.nan
            row[m + "_mae"] = float(np.abs(s["measured_n"] - s[c + "_n"]).mean())
            if len(s) >= 3 and s[c].std() > 0:
                row[m + "_pearson"] = float(stats.pearsonr(s["measured_n"], s[c + "_n"])[0])
            else:
                row[m + "_pearson"] = np.nan
            row[m + "_top1"] = bool(s.loc[s["measured"].idxmax(), c] == s[c].max() and s[c].max() > 0)
        per.append(row)
    PG = pd.DataFrame(per)
    PG.to_csv(os.path.join(F, "isopedia_per_gene.csv"), index=False)

    pd.set_option("display.width", 220)
    print("\n=== per compared isoform ===")
    print(d[["gene", "app_names", "n_app", "measured", "measured_n", "pred_model_n",
             "pred_promoter_only_n", "pred_uniform_n"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n=== overall ===")
    print(S.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n=== per gene ===")
    print(PG.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    return d, S, PG


if __name__ == "__main__":
    main()
