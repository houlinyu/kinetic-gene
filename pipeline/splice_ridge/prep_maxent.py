"""Genuine MaxEntScan donor/acceptor scores for every Isopedia isoform.

Same scorer as pipeline/splice/build_splice.py: `maxentpy` (X. Dong), which ships the original
Yeo & Burge (2004) matrices -- score5 on the 9-mer -3..+6 around the donor, score3 on the 23-mer
-20..+3 around the acceptor. The reference examples are re-validated here, exactly as
build_splice.py does, and the run aborts if they do not reproduce.

Genome: /data/mdl_data/hyu/projects/asd_villages_sqtl/refs_nochrM/genome.fa (GRCh38, chr-named,
faidx-indexed, read-only, random access only -- the file is never loaded).
The exon-coordinate convention of table_wide/isoforms.parquet is decided empirically: both
candidate offsets are scored on a sample of introns and the one that yields canonical GT..AG is
kept; the chosen offset and its GT-AG rate are printed and stored in work/maxent_meta.json.
Output: work/maxent.tsv  (isoform_uid, mes_donor_min/mean, mes_acc_min/mean, n_scored)
"""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(HERE, "work")
FA = "/data/mdl_data/hyu/projects/asd_villages_sqtl/refs_nochrM/genome.fa"

class Faidx:
    def __init__(self, path):
        self.fh = open(path, "rb"); self.idx = {}
        for ln in open(path + ".fai"):
            n, L, off, lb, bb = ln.split("\t")[:5]
            self.idx[n] = (int(L), int(off), int(lb), int(bb))
    def seq(self, chrom, start, end):          # 1-based inclusive
        if chrom not in self.idx: return None
        L, off, lb, bb = self.idx[chrom]
        if start < 1 or end > L or end < start: return None
        s0 = off + (start - 1) // lb * bb + (start - 1) % lb
        e0 = off + (end - 1) // lb * bb + (end - 1) % lb
        self.fh.seek(s0)
        return self.fh.read(e0 - s0 + 1).decode().replace("\n", "").replace("\r", "")

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
rc = lambda s: s.translate(COMP)[::-1]

def windows(chrom, strand, ee, ns, off):
    """ee = exon end, ns = next exon start, both as stored; `off` shifts them to 1-based inclusive.
    Returns (donor 9-mer, acceptor 23-mer) on the transcribed strand."""
    E = ee + off[0]; S = ns + off[1]
    if strand == "+":
        d = fa.seq(chrom, E - 2, E + 6); a = fa.seq(chrom, S - 20, S + 2)
    else:
        # on the minus strand the gene runs right-to-left: S is the donor side, E the acceptor side
        d = fa.seq(chrom, S - 6, S + 2); a = fa.seq(chrom, E - 2, E + 20)
        d = rc(d) if d else None; a = rc(a) if a else None
    return d, a

fa = Faidx(FA)
rows = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(W, "junctions.tsv"))][1:]
print("introns", len(rows), file=sys.stderr)

# --- decide the coordinate convention
best = None
for off in [(0, 0), (0, 1), (1, 0), (1, 1)]:
    ok = tot = 0
    for r in rows[::997][:3000]:
        d, a = windows(r[1], r[2], int(r[3]), int(r[4]), off)
        if not d or not a or len(d) != 9 or len(a) != 23: continue
        tot += 1
        if d[3:5].upper() == "GT" and a[18:20].upper() == "AG": ok += 1
    rate = ok / tot if tot else 0
    print(f"  offset {off}: GT..AG {rate:.3f} on {tot}", file=sys.stderr)
    if best is None or rate > best[1]: best = (off, rate)
off, rate = best
if rate < 0.8:
    raise SystemExit(f"no exon-coordinate convention gives canonical GT..AG (best {rate:.3f}) - refusing "
                     f"to label these scores MaxEntScan")
print(f"using offset {off} (GT..AG {rate:.3f})", file=sys.stderr)

from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
# maxentpy reloads its matrices from disk on every call unless they are passed in (~18 ms vs
# ~10 us per score), so they are loaded once here. Same matrices, same scores.
M5, M3 = load_matrix5(), load_matrix3()
exp5, exp3 = 10.858313101356437, 2.8867730651152104
got5, got3 = maxent.score5("cagGTAAGT", matrix=M5), maxent.score3("ttccaaacgaacttttgtAGgga", matrix=M3)
print(f"  MaxEntScan validation: score5={got5!r} (expect {exp5!r}); score3={got3!r} (expect {exp3!r})",
      file=sys.stderr)
if abs(got5 - exp5) > 1e-9 or abs(got3 - exp3) > 1e-9:
    raise SystemExit("maxentpy does not reproduce its reference examples - refusing to label scores MaxEntScan")

cache5, cache3 = {}, {}
def s5(x):
    if x not in cache5:
        cache5[x] = maxent.score5(x, matrix=M5) if set(x.upper()) <= set("ACGT") else None
    return cache5[x]
def s3(x):
    if x not in cache3:
        cache3[x] = maxent.score3(x, matrix=M3) if set(x.upper()) <= set("ACGT") else None
    return cache3[x]

D = collections.defaultdict(list); A = collections.defaultdict(list)
bad = 0
for i, r in enumerate(rows):
    d, a = windows(r[1], r[2], int(r[3]), int(r[4]), off)
    if not d or not a or len(d) != 9 or len(a) != 23: bad += 1; continue
    x, z = s5(d.upper()), s3(a.upper())
    if x is not None: D[r[0]].append(x)
    if z is not None: A[r[0]].append(z)
    if i % 100000 == 0: print(f"  {i}/{len(rows)}", file=sys.stderr, flush=True)

out = os.path.join(W, "maxent.tsv")
with open(out, "w") as fh:
    fh.write("isoform_uid\tmes_donor_min\tmes_donor_mean\tmes_acc_min\tmes_acc_mean\tmes_n\n")
    for u in set(D) | set(A):
        d, a = D.get(u, []), A.get(u, [])
        if not d or not a: continue
        fh.write(f"{u}\t{min(d):.4f}\t{sum(d)/len(d):.4f}\t{min(a):.4f}\t{sum(a)/len(a):.4f}\t{len(d)}\n")
json.dump(dict(offset=list(off), gt_ag_rate=rate, unscored_windows=bad,
               scorer="maxentpy score5/score3 (Yeo & Burge 2004 matrices)", genome=FA),
          open(os.path.join(W, "maxent_meta.json"), "w"), indent=1)
print("wrote", out, "unscorable windows:", bad, file=sys.stderr)
