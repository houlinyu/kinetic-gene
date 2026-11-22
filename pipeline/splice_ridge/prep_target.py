"""Aggregate Isopedia per-sample isoform usage into one usage fraction per transcript.

Source: <ISOPEDIA>/table_wide/context_table.parquet
  columns: model_id (sample), locus_id (ENSG), gene_name, isoform_uid, usage, lineage
  `usage` is that isoform's share of its locus in that sample.

Aggregation (read-weighted): per isoform, sum usage over all samples; then divide by the
locus total so the fractions of a locus sum to 1. Also keep n_pos = number of samples in
which the isoform has usage > 0 (used later as a "read support" confound baseline) and
n_samp = samples in which the locus is seen at all.
Output: work/usage.tsv
"""
import sys, collections, os
import pyarrow.parquet as pq

BUILD = os.environ.get("ISOPEDIA_BUILD", "table_wide2")   # table_wide2 = the 18,839-gene panel
B = f"/data/mdl_data/hyu/projects/public_isoform/isopedia/{BUILD}/"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work", "usage.tsv")

pf = pq.ParquetFile(B + "context_table.parquet")
usum = collections.Counter(); npos = collections.Counter()
locus = {}; gname = {}
locus_samples = collections.defaultdict(set)
cols = ["locus_id", "gene_name", "isoform_uid", "usage", "model_id"]
for i in range(pf.num_row_groups):
    t = pf.read_row_group(i, columns=cols)
    d = t.to_pydict()
    for lid, gn, iso, u, mid in zip(d["locus_id"], d["gene_name"], d["isoform_uid"], d["usage"], d["model_id"]):
        locus[iso] = lid; gname[iso] = gn
        if u and u > 0:
            usum[iso] += u; npos[iso] += 1
            locus_samples[lid].add(mid)
    print(f"  rowgroup {i+1}/{pf.num_row_groups} isoforms={len(locus)}", file=sys.stderr, flush=True)

gtot = collections.Counter()
for iso, s in usum.items():
    gtot[locus[iso]] += s

with open(OUT, "w") as fh:
    fh.write("isoform_uid\tlocus_id\tgene_name\tusage\tusage_sum\tn_pos\tn_locus_samples\n")
    for iso in locus:
        lid = locus[iso]; tot = gtot[lid]
        frac = (usum[iso] / tot) if tot > 0 else float("nan")
        fh.write(f"{iso}\t{lid}\t{gname[iso]}\t{frac:.6g}\t{usum[iso]:.6g}\t{npos[iso]}\t{len(locus_samples[lid])}\n")
print("wrote", OUT, len(locus), "isoforms", len(gtot), "loci", file=sys.stderr)
