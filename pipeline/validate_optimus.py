#!/usr/bin/env python
"""Validate the numpy Optimus 5-Prime forward pass against GEO GSM4084997 (varying length 25-100).

Reproduces the split of modeling/Varying_length_25to100_modeling.ipynb:
  random set, total_reads >= 10; e_test = top-100-by-reads per length 25..100 (7600);
  e_train = the rest. Human set (total_reads >= 10), same top-100-per-length -> subhuman (7600).
Writes external/optimus/train_scaler.json (mean / population std of e_train.rl, as sklearn
StandardScaler) and prints Pearson r / r^2 (notebook reports r^2 = 0.838 random test, 0.778 human).
Note: pandas sort ties may make the selected sets differ slightly from the original notebook.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimus  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(HERE, "external", "geo", "GSM4084997_varying_length_25to100.csv.gz")


def top_per_len(d):
    parts = []
    for L in range(25, 101):
        t = d[d["len"] == L].sort_values("total_reads", ascending=False, kind="mergesort")
        parts.append(t.iloc[:100])
    return pd.concat(parts)


def main():
    df1 = pd.read_csv(GEO, index_col=0)
    rnd = df1[(df1["set"] == "random") & (df1["total_reads"] >= 10)]
    hum = df1[(df1["set"] == "human") & (df1["total_reads"] >= 10)]
    e_test = top_per_len(rnd)
    e_train = rnd.drop(index=e_test.index)
    sub_h = top_per_len(hum)
    print(f"random>=10: {len(rnd)}  e_test {len(e_test)}  e_train {len(e_train)}  human>=10 {len(hum)}  subhuman {len(sub_h)}")
    mean, std = float(e_train["rl"].mean()), float(e_train["rl"].std(ddof=0))
    with open(optimus.SCALER_JSON, "w") as fh:
        json.dump({"mean": mean, "std": std, "n": int(len(e_train)),
                   "note": "StandardScaler fit on e_train.rl (random set, total_reads>=10, minus e_test), GEO GSM4084997"}, fh)
    print(f"train scaler: mean {mean:.4f} std {std:.4f}")
    out = {}
    for name, d in [("random_test", e_test), ("human_test", sub_h)]:
        p = optimus.predict_mrl(list(d["utr"]))
        r = float(np.corrcoef(p, d["rl"])[0, 1])
        out[name] = {"n": len(d), "pearson_r": r, "r2": r * r}
        print(f"{name}: n={len(d)} pearson r={r:.4f} r^2={r*r:.4f}  pred mean {p.mean():.3f} obs mean {d['rl'].mean():.3f}")
    with open(os.path.join(HERE, "external", "optimus", "validation.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
