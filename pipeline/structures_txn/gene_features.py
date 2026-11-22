"""Measure human exon / intron length distributions from the local GENCODE annotation.

Used instead of quoting a literature value: every number in txn_constants.json that
comes from here is reproducible with this script and the named annotation release.
One transcript per gene (Ensembl_canonical) so the statistic is not weighted by the
number of annotated isoforms.
"""
import gzip, json, sys, statistics as st
from collections import defaultdict

gtf = sys.argv[1]
exons = defaultdict(list)       # tx -> [(start,end)]
canon, pc = set(), set()
with gzip.open(gtf, "rt") as f:
    for line in f:
        if line[0] == "#":
            continue
        p = line.rstrip("\n").split("\t")
        if p[2] != "exon":
            continue
        a = p[8]
        if 'gene_type "protein_coding"' not in a or 'transcript_type "protein_coding"' not in a:
            continue
        tx = a.split('transcript_id "', 1)[1].split('"', 1)[0]
        if "Ensembl_canonical" in a:
            canon.add(tx)
        pc.add(tx)
        exons[tx].append((int(p[3]), int(p[4])))

ex_len, in_len, ntx = [], [], 0
for tx in canon:
    ivs = sorted(exons[tx])
    ntx += 1
    for s, e in ivs:
        ex_len.append(e - s + 1)
    for (s1, e1), (s2, e2) in zip(ivs, ivs[1:]):
        if s2 - e1 - 1 > 0:
            in_len.append(s2 - e1 - 1)

def stats(v):
    v = sorted(v)
    return {"n": len(v), "median_bp": st.median(v), "mean_bp": round(st.mean(v), 1),
            "p25_bp": v[len(v)//4], "p75_bp": v[3*len(v)//4], "min_bp": v[0], "max_bp": v[-1]}

out = {"annotation": gtf.split("/")[-1], "n_canonical_protein_coding_transcripts": ntx,
       "exons": stats(ex_len), "introns": stats(in_len),
       "method": "GENCODE exon records of protein_coding transcripts of protein_coding genes, "
                 "one Ensembl_canonical transcript per gene; introns = gaps between "
                 "consecutive exons of the same transcript"}
print(json.dumps(out, indent=1))
