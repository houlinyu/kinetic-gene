#!/usr/bin/env python
"""Build feature tables for the ridge models.

  mrl : GEO GSM4084997 (Sample 2019 varying-length 25-100 nt library, eGFP reporter), total_reads >= 10
        (filter from Varying_length_25to100_modeling.ipynb). Rows kept = set 'human' (all) + the random-set test
        split e_test (top 100 by reads per length, as in the notebook / validate_optimus.py). Neither was used to
        train the Optimus varying-length model (trained on random e_train), so Optimus comparison is leakage-free.
  te  : Cao et al. 2021 HEK293T TE table (Andreev 2015 data) x GENCODE v19 pc_transcripts (versioned IDs match).
        Filters: rnaseq_counts >= 50, riboseq_counts >= 10, te > 0, CDS starts with ATG, 5' UTR >= 10 nt,
        not ambiguous (no sibling isoform of the same gene with a bit-identical TE).
Outputs pipeline/ridge/work/{mrl,te}.csv.gz with sequence, target, Optimus prediction and motifs.features().
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.dirname(HERE)
sys.path.insert(0, PIPE)
sys.path.insert(0, HERE)
import optimus  # noqa: E402
import feats  # noqa: E402

WORK = os.path.join(HERE, "work")
GEO = os.path.join(PIPE, "external", "geo", "GSM4084997_varying_length_25to100.csv.gz")
CAO_HEK = os.path.join(PIPE, "external", "df_counts_and_len.TE_sorted.HEK_Andrev2015.with_annot.txt")
GENCODE = os.path.join(HERE, "downloads", "gencode.v19.pc_transcripts.fa.gz")
NJ = int(os.environ.get("NJ", "8"))


def top_per_len(d):
    return pd.concat([d[d["len"] == L].sort_values("total_reads", ascending=False, kind="mergesort").iloc[:100]
                      for L in range(25, 101)])


def build_mrl():
    df = pd.read_csv(GEO, index_col=0)[["utr", "set", "total_reads", "rl", "len"]]
    df = df[df["total_reads"] >= 10]
    rnd = df[df["set"] == "random"]
    hum = df[df["set"] == "human"].copy()
    hum["subset"] = "human"
    et = top_per_len(rnd).copy()
    et["subset"] = "random_test"
    d = pd.concat([hum, et])
    d.index.name = "geo_row"
    d = d.reset_index()
    d = d.drop_duplicates("utr")
    print("mrl rows", len(d), d["subset"].value_counts().to_dict(), flush=True)
    d["optimus"] = optimus.predict_mrl(list(d["utr"]))
    F = feats.many(list(d["utr"]), [feats.EGFP_HEAD] * len(d), NJ)
    d = pd.concat([d.reset_index(drop=True), pd.DataFrame(F)], axis=1)
    d.to_csv(os.path.join(WORK, "mrl.csv.gz"), index=False)


def read_gencode():
    out = {}
    with gzip.open(GENCODE, "rt") as fh:
        hdr, seq = None, []
        for line in fh:
            if line[0] == ">":
                if hdr:
                    out[hdr[0]] = (hdr, "".join(seq))
                hdr, seq = line[1:].strip().split("|"), []
            else:
                seq.append(line.strip())
        if hdr:
            out[hdr[0]] = (hdr, "".join(seq))
    return out


def build_te():
    g = read_gencode()
    rows = []
    with open(CAO_HEK) as fh:
        header = fh.readline().split()
        for line in fh:
            p = line.split()
            r = dict(zip(header, p[1:]))
            r["vid"] = p[0]
            rows.append(r)
    t = pd.DataFrame(rows)
    for c in ["riboseq_counts", "rnaseq_counts", "te"]:
        t[c] = pd.to_numeric(t[c], errors="coerce")
    n0 = len(t)
    # ambiguity on the full table (before count filters), same rule as measured.py
    t["tekey"] = t["te"].map(lambda v: f"{v:.10g}" if v == v else "nan")
    dup = t[t["tekey"] != "nan"].groupby(["external_gene_id", "tekey"])["vid"].transform("size")
    t["ambiguous"] = False
    t.loc[dup.index, "ambiguous"] = dup.values > 1
    keep = (t["rnaseq_counts"] >= 50) & (t["riboseq_counts"] >= 10) & (t["te"] > 0) & ~t["ambiguous"]
    t = t[keep]
    utrs, cdss, ok = [], [], []
    for vid in t["vid"]:
        h = g.get(vid)
        u = c = None
        if h:
            hdr, s = h
            rng = {x.split(":")[0]: x.split(":")[1] for x in hdr if ":" in x}
            if "UTR5" in rng and "CDS" in rng:
                a, b = map(int, rng["UTR5"].split("-"))
                cs = int(rng["CDS"].split("-")[0])
                if a == 1 and b + 1 == cs:
                    u, c = s[:b], s[cs - 1:]
        good = u is not None and len(u) >= 10 and c[:3] == "ATG"
        utrs.append(u)
        cdss.append(c)
        ok.append(good)
    t = t.assign(utr=utrs, cds=cdss)[ok]
    print(f"te rows: table {n0} -> filtered+matched {len(t)} ({t['external_gene_id'].nunique()} genes)", flush=True)
    t["log2_te"] = np.log2(t["te"])
    t["optimus"] = optimus.predict_mrl([u[-100:] for u in t["utr"]])
    # cds passed whole (uORF stops searched through CDS+), fold uses its first 30 nt
    F = feats.many(list(t["utr"]), list(t["cds"]), NJ, chunksize=16)
    t = pd.concat([t.reset_index(drop=True), pd.DataFrame(F)], axis=1)
    t["cds_head"] = t["cds"].str[:60]
    t = t.drop(columns=["cds", "tekey"])
    t.to_csv(os.path.join(WORK, "te.csv.gz"), index=False)


if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    which = sys.argv[1:] or ["mrl", "te"]
    if "mrl" in which:
        build_mrl()
    if "te" in which:
        build_te()
