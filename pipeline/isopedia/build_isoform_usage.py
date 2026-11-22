#!/usr/bin/env python
"""Measured isoform usage for the 12 Kinetic Gene genes, from the Isopedia long-read atlas.

Reads (read-only) the user's own public-isoform project:
  <ISO>/table_wide2/context_table.parquet   usage per (sample, catalog isoform), 18,839-gene panel
  <ISO>/table_wide2/isoforms.parquet        catalog isoform -> locus, representative transcript, exons
  <ISO>/table_wide2/sample_meta.tsv         platform / project / tissue / depth per sample
  <ISO>/ref/gencode.v49.primary_assembly.annotation.gtf.gz   to map our ENSTs onto catalog isoforms

An Isopedia "isoform" is a splice pattern: GENCODE v49 transcripts that share an intron chain are one
isoform, so we join our app transcripts to the catalog by intron chain, not by transcript id (the app is
built on Ensembl 116 and many of its ENST ids do not exist in v49). Per gene and sample the catalog
isoforms carry usage fractions that sum to 1 together with one pooled "other:<gene>" row, which holds
every splice pattern that never reached 5% usage in 5% of samples.

Writes data/isoform_usage.json (+ app/public copy) and the CSVs behind the figures.
"""
import argparse, gzip, hashlib, json, os, re, sys
from collections import defaultdict
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ISO = "/data/mdl_data/hyu/projects/public_isoform/isopedia"
ATTR = re.compile(r'(\S+) "([^"]*)"')


def intron_chain(exons):
    ex = sorted(exons)
    if len(ex) < 2:
        return None
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))


def chain_key(chrom, strand, exons):
    """same hash the Isopedia context-table builder uses for a splice-pattern isoform"""
    ic = intron_chain(exons)
    if ic is None:
        return "mono:" + hashlib.sha1(f"{chrom}|{strand}|{tuple(sorted(exons))}".encode()).hexdigest()[:20]
    return hashlib.sha1(f"{chrom}|{strand}|{ic}".encode()).hexdigest()[:20]


def read_gtf_genes(path, genes):
    """gene_name in `genes` -> {transcript_id: (chrom, strand, [exons], gene_id, tx_name)}"""
    out = {}
    want = {f'gene_name "{g}"' for g in genes}
    ex = defaultdict(list)
    meta = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            if not any(w in f[8] for w in want):
                continue
            a = dict(ATTR.findall(f[8]))
            tid = a["transcript_id"]
            ex[tid].append((int(f[3]), int(f[4])))
            meta[tid] = (f[0], f[6], a.get("gene_id"), a.get("gene_name"), a.get("transcript_name"))
    for t, m in meta.items():
        out.setdefault(m[3], {})[t] = (m[0], m[1], ex[t], m[2], m[4])
    return out


def parse_exons_field(s):
    """'chr22:(+)[[s,e],...]' -> (chrom, strand, [(s,e)...])"""
    chrom, rest = s.split(":", 1)
    strand = rest[rest.index("(") + 1]
    arr = json.loads(rest[rest.index("["):])
    return chrom, strand, [tuple(x) for x in arr]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso", default=ISO)
    ap.add_argument("--table", default="table_wide2", help="which Isopedia context table build to use")
    ap.add_argument("--splice", default="app/public/splice_data.json")
    ap.add_argument("--out", default="data/isoform_usage.json")
    ap.add_argument("--csv-dir", default="data/figures")
    ap.add_argument("--min-tissue-samples", type=int, default=10)
    a = ap.parse_args()
    T = os.path.join(a.iso, a.table)

    app = json.load(open(a.splice))
    app_genes = {g["gene"]: g for g in app["genes"]}
    names = sorted(app_genes)
    print(f"app genes: {len(names)}  isoforms: {sum(len(g['isoforms']) for g in app['genes'])}")

    # ---- 1. measured usage, panel-filtered, for these gene names only
    ctx = pq.read_table(os.path.join(T, "context_table.parquet"),
                        columns=["model_id", "locus_id", "gene_name", "isoform_uid", "usage", "lineage"],
                        filters=[("gene_name", "in", names)]).to_pandas()
    for c in ("model_id", "locus_id", "gene_name", "isoform_uid", "lineage"):
        ctx[c] = ctx[c].astype(str)
    print(f"context rows: {len(ctx):,}  genes present: {ctx.gene_name.nunique()}  samples: {ctx.model_id.nunique()}")
    meta = pd.read_csv(os.path.join(T, "sample_meta.tsv"), sep="\t", dtype=str)
    meta = meta.set_index("model_id")

    # ---- 2. catalog isoform -> intron chain
    loci = set(ctx.locus_id)
    cat = pq.read_table(os.path.join(T, "isoforms.parquet")).to_pandas()
    cat = cat[cat.locus_id.astype(str).isin(loci)].copy()
    cat_chain = {}
    for r in cat.itertuples():
        chrom, strand, exons = parse_exons_field(r.exons)
        cat_chain[r.isoform_uid] = (chrom, strand, intron_chain(exons), chain_key(chrom, strand, exons), r.transcript_id)

    # ---- 3. our transcripts -> intron chain (via GENCODE v49 where the id exists; else via the
    #         app's own exon coordinates, so Ensembl-116-only ids still match structurally)
    gtf = read_gtf_genes(os.path.join(a.iso, "ref", "gencode.v49.primary_assembly.annotation.gtf.gz"), names)
    print("GENCODE v49 transcripts found per gene:", {g: len(v) for g, v in sorted(gtf.items())})

    def app_chain(g, iso):
        """intron chain of one app isoform from the app's genomic exon coordinates"""
        ex = []
        byid = {e["id"]: e for e in g["exons"]}
        for eid in iso["exon_ids"]:
            e = byid[eid]
            m = re.match(r"(chr[^:]+):(\d+)-(\d+)", e["g"])
            ex.append((int(m.group(2)), int(m.group(3))))
        return ("chr" + str(g["chrom"]), g["strand"], intron_chain(ex))

    # chain -> uid for this locus
    rows, join_log = [], []
    for gn in names:
        g = app_genes[gn]
        gctx = ctx[ctx.gene_name == gn]
        if gctx.empty:
            join_log.append(dict(gene=gn, status="gene absent from Isopedia panel"))
            continue
        locus = gctx.locus_id.iloc[0]
        uids = [u for u in gctx.isoform_uid.unique() if not u.startswith("other:")]
        chain2uid = {}
        for u in uids:
            chrom, strand, ic, ck, tx = cat_chain.get(u, (None, None, None, None, None))
            if ic is not None:
                chain2uid[ic] = u

        def fuzzy(ic, tol=10):
            """the builder merged chains whose boundaries agree within 10 bp, so retry loosely"""
            if ic is None:
                return None
            for c, u in chain2uid.items():
                if len(c) != len(ic):
                    continue
                if all(abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol for a, b in zip(c, ic)):
                    return u
            return None
        # v49 transcript ids of this gene, per chain
        v49 = gtf.get(gn, {})
        tx2chain = {t: intron_chain(v[2]) for t, v in v49.items()}
        for iso in g["isoforms"]:
            tid = iso["id"]
            base = tid.split(".")[0]
            ic = None
            how = None
            # prefer the v49 annotation of the same transcript id (version-insensitive)
            hit = [t for t in v49 if t.split(".")[0] == base]
            if hit:
                ic = tx2chain[hit[0]]
                how = "ENST present in GENCODE v49"
            else:
                ic = app_chain(g, iso)[2]
                how = "Ensembl-116-only id, matched by intron chain from app exon coordinates"
            uid = chain2uid.get(ic)
            if uid is None:
                uid = fuzzy(ic)
                if uid is not None:
                    how += " (+/-10 bp flank merge)"
            if uid is None:
                how += "; no catalog isoform with this intron chain -> pooled into 'other' (below the 5% catalog threshold)"
            rows.append(dict(gene=gn, locus_id=locus, app_id=tid, app_base=base,
                             app_name=iso.get("name"), isoform_uid=uid, join=how,
                             in_v49=bool(hit), n_introns=0 if ic is None else len(ic)))
        join_log.append(dict(gene=gn, status="measured", n_catalog=len(uids),
                             n_samples=int(gctx.model_id.nunique())))
    link = pd.DataFrame(rows)
    matched = link[link.isoform_uid.notna()]
    print(f"app isoforms joined to an Isopedia catalog isoform: {len(matched)} / {len(link)}")

    # ---- 4. aggregate usage: mean over samples (sums to 1 with 'other'), median as secondary
    def agg(sub):
        piv = sub.pivot_table(index="model_id", columns="isoform_uid", values="usage", aggfunc="mean")
        return piv.mean(axis=0), piv.median(axis=0), len(piv)

    genes_out, long_rows = {}, []
    for gn in names:
        gctx = ctx[ctx.gene_name == gn]
        if gctx.empty:
            genes_out[gn] = dict(measured=False,
                                 note="not in the Isopedia 18,839-gene panel: fewer than 20 supporting reads "
                                      "in >=2% of the 869 samples, so no usage fractions exist for it")
            continue
        mean, med, n = agg(gctx)
        lmap = {}
        for lin, sub in gctx.groupby("lineage"):
            if sub.model_id.nunique() < a.min_tissue_samples:
                continue
            m, _, nn = agg(sub)
            lmap[lin] = dict(n_samples=int(nn), usage=None, _m=m)
        sub = matched[matched.gene == gn]
        uid2app = defaultdict(list)
        for r in sub.itertuples():
            uid2app[r.isoform_uid].append(r.app_base)

        def by_enst(series):
            out = {}
            for uid, v in series.items():
                if uid.startswith("other:"):
                    continue
                for b in uid2app.get(uid, []):
                    out[b] = round(float(v), 6)
            return out

        other = float(mean.get(f"other:{gctx.locus_id.iloc[0]}", 0.0))
        covered = sum(float(mean[u]) for u in uid2app if u in mean.index)
        genes_out[gn] = dict(
            measured=True,
            locus_id=gctx.locus_id.iloc[0],
            n_samples=int(n),
            n_catalog_isoforms=int(gctx[~gctx.isoform_uid.str.startswith("other:")].isoform_uid.nunique()),
            n_app_isoforms=int((link.gene == gn).sum()),
            n_app_isoforms_matched=int(len(sub)),
            overall=by_enst(mean),
            overall_median=by_enst(med),
            pooled_other_fraction=round(other, 6),
            unmatched_catalog_fraction=round(float(1.0 - other - covered), 6),
            app_isoform_share_of_gene=round(covered, 6),
            tissues={k: dict(n_samples=v["n_samples"], usage=by_enst(v["_m"])) for k, v in lmap.items()},
            note=("usage fractions are means over samples; the denominator is all reads Isopedia assigned to "
                  "this gene, so our app isoforms' fractions do not sum to 1"),
        )
        for uid, v in mean.items():
            long_rows.append(dict(gene=gn, isoform_uid=uid, mean_usage=float(v),
                                  median_usage=float(med[uid]), n_samples=int(n),
                                  app_ids=";".join(uid2app.get(uid, [])),
                                  representative_tx=cat_chain.get(uid, (None,) * 5)[4]))

    plat = meta.loc[meta.index.intersection(ctx.model_id.unique())]
    out = dict(
        meta=dict(
            source="Isopedia (Zheng & Sedlazeck, bioRxiv 2026) — read-level index of 1,007 public human "
                   "long-read RNA-seq samples; queried and aggregated in the user's own public_isoform project "
                   f"at {a.iso} ({a.table}/)",
            dataset_description=f"{int(ctx.model_id.nunique())} samples passing depth >= 200,000 indexed reads, "
                                f"{(plat.platform == 'ONT').sum()} ONT and {(plat.platform == 'PB').sum()} PacBio, "
                                f"{plat.project.nunique()} BioProjects, {plat.std_tissue.nunique()} standardised tissues",
            how_quantified="`isopedia isoform -g gencode.v49.gtf` returns per-GENCODE-transcript, per-sample "
                           "supporting read counts (RC_FSM_JC = reads whose junction chain fully matches; "
                           "RC_EM_ISM = partial reads assigned by EM). Transcripts sharing an intron chain are "
                           "collapsed to one splice-pattern isoform taking max(RC_FSM_JC)+sum(RC_EM_ISM); usage = "
                           "isoform reads / gene reads in that sample. A (sample, gene) is measured at >= 20 gene "
                           "reads; the panel keeps genes measured in >= 2% of samples; catalog isoforms must reach "
                           "usage >= 0.05 in >= 5% of measured samples, everything else is pooled into 'other'. "
                           "Built by tools/isopedia_context_table.py in ~/IsoMap-RNAFM/isoform_usage_benchmark.",
            samples_or_tissues=sorted(set(ctx.lineage)),
            n=int(ctx.model_id.nunique()),
            denominator="all reads Isopedia assigned to the gene, i.e. every splice pattern it saw, including "
                        "isoforms this app does not model and the pooled 'other' bin. Per gene and sample, "
                        "catalog usage + other = 1; our app isoforms cover only part of that.",
            caveats=[
                "long-read quantification, not ground truth: EM assignment of partial (ISM) reads, 3'/5' "
                "degradation and library chemistry all shift usage",
                "annotation is GENCODE v49; the app is built on Ensembl 116, so app-only transcript ids are "
                "joined by intron chain, and transcripts that share an intron chain share one measured value",
                "heterogeneous public samples (tumour, cell line, tissue; ONT and PacBio mixed), not a matched "
                "panel, so 'overall' is a mean over a convenience sample",
                "isoforms below 5% usage in 95% of samples are invisible (pooled into 'other')",
                "mean is over samples, so a gene's measured fractions sum to 1 only when 'other' is included",
            ],
            citation_or_internal="Isopedia: Zheng & Sedlazeck, bioRxiv 2026 (MIT-licensed tool). Aggregation: "
                                 "internal (user's public_isoform project).",
            generated=pd.Timestamp.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            join=dict(key="unversioned Ensembl transcript id (strip '.N' from isoforms[].id)",
                      n_app_isoforms=int(len(link)), n_matched=int(len(matched))),
        ),
        genes=genes_out,
    )
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    os.makedirs(a.csv_dir, exist_ok=True)
    pd.DataFrame(long_rows).to_csv(os.path.join(a.csv_dir, "isopedia_usage_catalog.csv"), index=False)
    link.to_csv(os.path.join(a.csv_dir, "isopedia_join.csv"), index=False)
    print(json.dumps(join_log, indent=1))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
