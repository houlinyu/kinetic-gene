#!/usr/bin/env python
"""Stage A of the TSS-weight pipeline: pull the FANTOM5 CAGE peaks that sit in our 12 loci.

Input (downloaded by fetch.sh, not regenerated here):
  hg38_fair+new_CAGE_peaks_phase1and2.bed.gz          hg38 coordinates of the robust CAGE peak set
  hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz   tag counts, peaks x 1829 CAGE libraries

The counts table is keyed by the peak's *hg19* name, so the BED file is the only place the hg38
coordinate lives; we join on that name. Output is a small JSON holding, for the peaks that fall
inside (or within --margin nt of) the gene loci in splice_data.json:
  { meta, samples: [name...], mapped: [library size...], peaks: [{id, chrom, start, end, strand,
    tss (1-based representative TSS = BED thickStart+1), counts: [...] }] }
Nothing is interpreted here; see weight_tss.py for the TSS assignment and the weights.
"""
import argparse
import gzip
import json
import os
import re
import sys
import urllib.parse

AP = argparse.ArgumentParser()
AP.add_argument('--splice-json', required=True)
AP.add_argument('--bed', required=True)
AP.add_argument('--counts', required=True)
AP.add_argument('--out', required=True)
AP.add_argument('--margin', type=int, default=2000,
                help='keep peaks this many nt outside the annotated locus too')
A = AP.parse_args()

with open(A.splice_json) as fh:
    data = json.load(fh)

loci = []
for g in data['genes']:
    chrom, span = g['locus'].split(':')
    lo, hi = (int(x) for x in span.split('-'))
    loci.append((g['gene'], 'chr' + str(chrom), lo - A.margin, hi + A.margin))
print(f'[extract] {len(loci)} loci', file=sys.stderr)


def in_loci(chrom, pos):
    for gene, c, lo, hi in loci:
        if c == chrom and lo <= pos <= hi:
            return gene
    return None


# ---- BED: peak name -> hg38 coordinate + representative TSS
keep = {}
with gzip.open(A.bed, 'rt') as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < 8:
            continue
        chrom, start, end, name, _score, strand, thick = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5], int(f[6])
        if in_loci(chrom, thick + 1) is None:
            continue
        keep[name] = {'id': name, 'chrom': chrom, 'start': start + 1, 'end': end,
                      'strand': strand, 'tss': thick + 1, 'gene_locus': in_loci(chrom, thick + 1)}
print(f'[extract] {len(keep)} CAGE peaks inside the loci', file=sys.stderr)

# ---- counts table: header, library sizes, and the rows we kept
SAMPLE = re.compile(r'^counts\.(.+)$')
samples, mapped, rows = None, None, {}
with gzip.open(A.counts, 'rt') as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[0] == '00Annotation':
            samples = [urllib.parse.unquote(SAMPLE.match(c).group(1)) for c in f[1:]]
            continue
        if f[0] == '01STAT:MAPPED':
            mapped = [int(x) for x in f[1:]]
            continue
        if f[0] in keep:
            rows[f[0]] = [int(x) for x in f[1:]]
            if len(rows) == len(keep):
                break
print(f'[extract] {len(samples)} CAGE libraries, {len(rows)} peak rows recovered', file=sys.stderr)

peaks = []
for name, rec in keep.items():
    if name not in rows:
        continue
    rec = dict(rec)
    rec['counts'] = rows[name]
    peaks.append(rec)
peaks.sort(key=lambda r: (r['chrom'], r['tss']))

os.makedirs(os.path.dirname(os.path.abspath(A.out)), exist_ok=True)
with open(A.out, 'w') as fh:
    json.dump({'meta': {'bed': os.path.basename(A.bed), 'counts': os.path.basename(A.counts),
                        'assembly': 'hg38', 'margin': A.margin},
               'samples': samples, 'mapped': mapped, 'peaks': peaks}, fh)
print(f'[extract] wrote {A.out}', file=sys.stderr)
