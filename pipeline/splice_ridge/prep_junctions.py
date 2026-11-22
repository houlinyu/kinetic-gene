"""Dump every intron of every Isopedia isoform as (isoform, chrom, strand, exon_end, next_exon_start).
Coordinates are copied verbatim from table_wide/isoforms.parquet `exons`; prep_maxent.py decides
empirically which convention they are in (it tries both offsets and keeps the one that yields GT..AG).
Output: work/junctions.tsv
"""
import os, sys, ast
import pyarrow.parquet as pq
BUILD = os.environ.get("ISOPEDIA_BUILD", "table_wide2")   # table_wide2 = the 18,839-gene panel
B = f"/data/mdl_data/hyu/projects/public_isoform/isopedia/{BUILD}/"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work", "junctions.tsv")
t = pq.read_table(B + "isoforms.parquet", columns=["isoform_uid", "exons"]).to_pydict()
n = 0
with open(OUT, "w") as fh:
    fh.write("isoform_uid\tchrom\tstrand\texon_end\tnext_exon_start\n")
    for u, e in zip(t["isoform_uid"], t["exons"]):
        try:
            chrom, rest = e.split(":", 1)
            strand = rest[1]
            ex = sorted([(int(a), int(b)) for a, b in ast.literal_eval(rest[rest.index("["):])])
        except Exception:
            continue
        for i in range(len(ex) - 1):
            fh.write(f"{u}\t{chrom}\t{strand}\t{ex[i][1]}\t{ex[i+1][0]}\n"); n += 1
print("wrote", OUT, n, "introns", file=sys.stderr)
