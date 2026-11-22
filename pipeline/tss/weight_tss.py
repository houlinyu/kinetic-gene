#!/usr/bin/env python
"""Stage B: turn the extracted FANTOM5 CAGE peaks into measured promoter-usage weights and merge
them into splice_data.json as `tss[i].weight` + `tss[i].weight_source`.

What a weight means
-------------------
`weight` is this promoter's share of the gene's CAGE tags, averaged over a 41-library human adult
tissue panel (the FANTOM5 "<tissue>, adult, poolN" libraries). CAGE counts 5' cap ends, so tags at
a TSS are a direct measurement of how often transcription starts there; the share of a gene's tags
is therefore its share of transcripts, i.e. exactly the promoter mixing weight the splicing model
needs. Weights are normalised to sum to 1 over the *matched* TSSs of each gene.

Matching
--------
splice_data.json positions are 0-based on a per-gene pre-mRNA axis in transcript orientation, so
gene.locus (1-based GRCh38 inclusive) converts them back to the genome:
    + strand: g = locus_lo + pos        - strand: g = locus_hi - pos
A CAGE peak is assigned to the annotated TSS it is nearest to, and only if the annotated TSS lies
inside the peak interval expanded by --window nt on each side (default 50) and the strands agree.
Assignment is nearest-TSS-wins, so one peak is never counted for two promoters; an annotated TSS
with no peak of its own keeps weight null (two annotated starts closer together than CAGE can
resolve leave the weaker one null rather than splitting a peak between them).

Nothing is estimated: a gene with no matched peak at all keeps every weight null and the model
falls back to uniform mixing with a flag.

Rerun: see pipeline/tss/README.md.
"""
import argparse
import datetime
import json
import os
import re
import statistics
import sys

PANEL_RE = re.compile(r'^(?P<tissue>[^.]+), adult, pool\d+\.')

SOURCE = {
    'key': 'fantom5_cage',
    'dataset': 'FANTOM5 phase 1+2 robust CAGE peak set, tag counts, hg38 reprocessed release '
               '(hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz + '
               'hg38_fair+new_CAGE_peaks_phase1and2.bed.gz)',
    'citation': 'Lizio M, Harshbarger J, Shimoji H, et al. Gateways to the FANTOM5 promoter level '
                'mammalian expression atlas. Genome Biology 2015;16:22. Data: FANTOM Consortium & '
                'the RIKEN PMI and CLST (DGT). A promoter-level mammalian expression atlas. '
                'Nature 2014;507:462-470.',
    'url': 'https://fantom.gsc.riken.jp/5/data/',
    'licence': 'CC BY 4.0 - "FANTOM5 and FANTOM6 data by RIKEN is licensed under a Creative '
               'Commons Attribution 4.0 International License" (https://fantom.gsc.riken.jp/5/data/; '
               'the data-file index declares "license":"http://creativecommons.org/licenses/by/4.0/").',
    'method': "CAGE (cap analysis of gene expression): 5'-cap tag counts at each robust CAGE peak, "
              'normalised to tags per million mapped reads with the library size given in the '
              "table's 01STAT:MAPPED row, then expressed as the peak's share of that gene's "
              'matched-peak signal.',
}

AP = argparse.ArgumentParser()
AP.add_argument('--cage', required=True, help='output of extract_cage.py')
AP.add_argument('--splice-json', action='append', required=True,
                help='splice_data.json to rewrite in place; repeat for each copy')
AP.add_argument('--window', type=int, default=50, help='match window in nt (each side)')
AP.add_argument('--report', required=True)
AP.add_argument('--dry-run', action='store_true')
A = AP.parse_args()

with open(A.cage) as fh:
    cage = json.load(fh)

panel = [(i, PANEL_RE.match(s).group('tissue')) for i, s in enumerate(cage['samples'])
         if PANEL_RE.match(s)]
if not panel:
    sys.exit('no "<tissue>, adult, poolN" libraries found in the CAGE table')
print(f'[weights] panel = {len(panel)} FANTOM5 adult tissue libraries', file=sys.stderr)
mapped = cage['mapped']

# peak -> tags per million, per panel library
for p in cage['peaks']:
    p['tpm'] = [p['counts'][i] / mapped[i] * 1e6 for i, _ in panel]

primary = A.splice_json[0]
with open(primary) as fh:
    data = json.load(fh)

report = {'generated': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
          'window_nt': A.window, 'panel': [t for _, t in panel], 'source': SOURCE, 'genes': []}

for g in data['genes']:
    chrom, span = g['locus'].split(':')
    lo, hi = (int(x) for x in span.split('-'))
    chrom = 'chr' + str(chrom)
    plus = g['strand'] == '+'

    def genomic(pos):
        return lo + pos if plus else hi - pos

    tss = g.get('tss') or []
    gpos = [genomic(t['pos']) for t in tss]
    cand = [p for p in cage['peaks'] if p['chrom'] == chrom and p['strand'] == g['strand']]

    # Assignment. A peak matches every annotated TSS inside its interval +/- window. When two or
    # more annotated starts fall in the same peak they are closer together than CAGE resolves (e.g.
    # TP53 T1/T2 are 3 nt apart), so the peak's signal is SPLIT EQUALLY between them -- a stated
    # modelling choice, flagged per TSS as `shared_peak`, not a measurement of either one alone.
    assigned = {i: [] for i in range(len(tss))}
    shared = {i: 0 for i in range(len(tss))}
    for p in cand:
        hits = [i for i, gp in enumerate(gpos)
                if p['start'] - A.window <= gp <= p['end'] + A.window]
        if not hits:
            continue
        for i in hits:
            assigned[i].append((p, 1.0 / len(hits)))
            if len(hits) > 1:
                shared[i] = max(shared[i], len(hits))

    # per-TSS mean TPM over the panel, and the per-tissue share (for the spread note)
    mean_tpm, per_tissue = {}, {}
    for i, peaks in assigned.items():
        if not peaks:
            continue
        vals = [sum(p['tpm'][k] * f for p, f in peaks) for k in range(len(panel))]
        per_tissue[i] = vals
        mean_tpm[i] = sum(vals) / len(vals)

    total = sum(mean_tpm.values())
    grec = {'gene': g['gene'], 'strand': g['strand'], 'locus': g['locus'],
            'n_tss': len(tss), 'n_matched': len(mean_tpm), 'tss': []}

    for i, t in enumerate(tss):
        peaks = assigned[i]
        if not peaks or total <= 0:
            t['weight'] = None
            t['weight_source'] = None
            why = ('no CAGE peak within ±%d nt on the %s strand in the panel' % (A.window, g['strand'])
                   if not peaks else 'gene has no CAGE signal in the panel')
            grec['tss'].append({'id': t['id'], 'pos': t['pos'], 'genomic': gpos[i],
                                'weight': None, 'reason': why})
            continue
        w = mean_tpm[i] / total
        # how much this promoter's share moves across the 41 tissues (weights are tissue-dependent)
        shares = []
        for k in range(len(panel)):
            tot_k = sum(per_tissue[j][k] for j in per_tissue)
            if tot_k > 0:
                shares.append(per_tissue[i][k] / tot_k)
        t['weight'] = round(w, 6)
        t['weight_source'] = {
            'dataset': SOURCE['dataset'],
            'citation': SOURCE['citation'],
            'url': SOURCE['url'],
            'licence': SOURCE['licence'],
            'method': SOURCE['method'],
            'tissue': 'mean over a %d-library FANTOM5 human adult tissue panel '
                      '("<tissue>, adult, poolN"); promoter usage is tissue-dependent, '
                      'per-tissue share of this promoter ranges %.3f-%.3f'
                      % (len(panel), min(shares) if shares else float('nan'),
                         max(shares) if shares else float('nan')),
            'cage_peaks': [p['id'] for p, _ in peaks],
            'offset_nt': int(min(abs(p['tss'] - gpos[i]) for p, _ in peaks)),
            'shared_peak': (None if not shared[i] else
                            '%d annotated starts fall in the same CAGE peak (closer than CAGE '
                            'resolves); the peak signal is split equally between them' % shared[i]),
            'mean_tpm': round(mean_tpm[i], 4),
            'window_nt': A.window,
        }
        grec['tss'].append({'id': t['id'], 'pos': t['pos'], 'genomic': gpos[i],
                            'weight': round(w, 6), 'mean_tpm': round(mean_tpm[i], 4),
                            'n_peaks': len(peaks), 'shared_with': shared[i] or None,
                            'offset_nt': t['weight_source']['offset_nt'],
                            'share_min': round(min(shares), 4) if shares else None,
                            'share_max': round(max(shares), 4) if shares else None,
                            'share_median': round(statistics.median(shares), 4) if shares else None})
    report['genes'].append(grec)

method_note = (
    'tss[i].weight is measured promoter usage, not an assumption: it is this TSS\'s share of the '
    'gene\'s FANTOM5 CAGE tags (5\'-cap counts at robust CAGE peaks), averaged over a %d-library '
    'human adult tissue panel ("<tissue>, adult, poolN"), and normalised to sum to 1 over the '
    'matched TSSs of the gene. A CAGE peak is assigned to the annotated TSS it is nearest to, and '
    'only when that TSS lies within ±%d nt of the peak interval on the same strand; annotated '
    'positions converted to GRCh38 from gene.locus. Where no peak matches, weight is null (never '
    'estimated) and app/src/splice/model_exact.js renormalises over the matched promoters, or '
    'falls back to uniform mixing with a flag if a gene has none. Promoter usage is '
    'tissue-dependent: the per-promoter share across the 41 tissues is given in '
    'tss[i].weight_source.tissue and in data/tss_weights.json. Pipeline: pipeline/tss/.'
    % (len(panel), A.window)
)

matched = sum(gr['n_matched'] for gr in report['genes'])
n_tss = sum(gr['n_tss'] for gr in report['genes'])
print(f'[weights] matched {matched}/{n_tss} annotated TSSs', file=sys.stderr)
for gr in report['genes']:
    print('  %-10s %d/%d  %s' % (gr['gene'], gr['n_matched'], gr['n_tss'],
                                 ' '.join('%s=%s' % (t['id'], 'null' if t['weight'] is None
                                                     else '%.3f' % t['weight'])
                                          for t in gr['tss'])), file=sys.stderr)

if A.dry_run:
    print('[weights] dry run, nothing written', file=sys.stderr)
    sys.exit(0)

with open(A.report, 'w') as fh:
    json.dump(report, fh, indent=1)
print(f'[weights] wrote {A.report}', file=sys.stderr)

# ---- merge: only meta.sources / meta.methods.tss_weight and the tss weight fields are touched
srcs = data['meta'].setdefault('sources', [])
if not any(s.get('key') == SOURCE['key'] for s in srcs):
    srcs.append({'key': SOURCE['key'], 'citation': SOURCE['citation'], 'url': SOURCE['url'],
                 'licence': SOURCE['licence']})
data['meta'].setdefault('methods', {})['tss_weight'] = method_note

weights = {g['gene']: {t['id']: (t.get('weight'), t.get('weight_source')) for t in (g.get('tss') or [])}
           for g in data['genes']}

for path in A.splice_json:
    with open(path) as fh:
        d = json.load(fh)
    d['meta'].setdefault('sources', [])
    if not any(s.get('key') == SOURCE['key'] for s in d['meta']['sources']):
        d['meta']['sources'].append({'key': SOURCE['key'], 'citation': SOURCE['citation'],
                                     'url': SOURCE['url'], 'licence': SOURCE['licence']})
    d['meta'].setdefault('methods', {})['tss_weight'] = method_note
    for g in d['genes']:
        w = weights.get(g['gene'], {})
        for t in (g.get('tss') or []):
            weight, wsrc = w.get(t['id'], (None, None))
            t['weight'] = weight
            t['weight_source'] = wsrc
    with open(path, 'w') as fh:
        json.dump(d, fh, separators=(',', ':'))
    print(f'[weights] merged into {path}', file=sys.stderr)
