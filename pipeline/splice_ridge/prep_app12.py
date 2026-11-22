"""Measured usage for the app's 12 genes, straight from the genome-wide Isopedia query table.

The curated table_wide catalogue keeps only a filtered subset of transcripts, so only 11 of the
app's annotated isoforms survive there. This instead reads every GENCODE v49 transcript of the
12 genes out of query/gencode_v49.isoform.tsv.gz (pre-grepped into app12_rows.tsv) and takes the
CPM_FSM_JC field (full-splice-match, junction-chain-assigned reads, CPM) of every sample.
usage(t) = sum_samples CPM(t) / sum_samples sum_t' CPM(t') over the gene.
Output: work/app12_usage.tsv   columns: transcript_id gene usage cpm n_pos n_samples exon_count length
"""
import sys, os, collections
SRC = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work", "app12_usage.tsv")
cpm = collections.Counter(); npos = collections.Counter()
gene = {}; meta = {}
nsamp = 0
for ln in open(SRC):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 14: continue
    tid, attrs, fmt = f[5], f[11], f[12]
    gn = dict(kv.split(":", 1) for kv in attrs.split(";") if ":" in kv).get("gene_name")
    if gn is None: continue
    keys = fmt.split(":")
    if "CPM_FSM_JC" not in keys: continue
    k = keys.index("CPM_FSM_JC")
    gene[tid] = gn; meta[tid] = (f[4], f[3])
    vals = f[13:]
    nsamp = max(nsamp, len(vals))
    for v in vals:
        p = v.split(":")
        if len(p) <= k: continue
        try: x = float(p[k])
        except ValueError: continue
        if x > 0: cpm[tid] += x; npos[tid] += 1
tot = collections.Counter()
for t, c in cpm.items(): tot[gene[t]] += c
with open(OUT, "w") as fh:
    fh.write("transcript_id\tgene\tusage\tcpm\tn_pos\tn_samples\texon_count\tlength\n")
    for t in gene:
        g = gene[t]; T = tot[g]
        fh.write(f"{t}\t{g}\t{(cpm[t]/T if T else 0):.6g}\t{cpm[t]:.6g}\t{npos[t]}\t{nsamp}\t{meta[t][0]}\t{meta[t][1]}\n")
print("wrote", OUT, len(gene), "transcripts,", len(tot), "genes,", nsamp, "samples", file=sys.stderr)
