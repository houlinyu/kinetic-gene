# `pipeline/tss/` — measured promoter usage (`tss[i].weight`)

Promoter choice decides which 5′ UTR a transcript carries, so it is the join between tab 2
("how the isoform is made") and tab 1 ("how well it is translated"). Until this pipeline existed
every `tss` entry in `splice_data.json` had `weight: null` and `app/src/splice/model_exact.js`
mixed a gene's promoters **uniformly** — an assumption that also capped any single-TSS isoform at
1/n (ATF4: 20%). These scripts replace that assumption with a measurement.

## What the weight is

`tss[i].weight` = this promoter's share of the gene's **FANTOM5 CAGE** tags, averaged over a
41-library human **adult tissue panel**, normalised to sum to 1 over the gene's matched TSSs.
CAGE sequences 5′ cap ends, so a tag at a TSS is a direct observation that transcription started
there; the share of a gene's tags is therefore the share of its transcripts.

- Dataset: FANTOM5 phase 1+2 **robust** CAGE peak set, hg38 reprocessed release —
  `hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz` (tag counts, 209 911 peaks × 1829
  libraries) joined to `hg38_fair+new_CAGE_peaks_phase1and2.bed.gz` (the counts table is keyed by
  the peak's *hg19* name, so the BED file is the only source of the hg38 coordinate).
- Citation: Lizio M *et al.* **Gateways to the FANTOM5 promoter level mammalian expression atlas.**
  *Genome Biology* 2015;16:22. Data: FANTOM Consortium and the RIKEN PMI and CLST (DGT),
  **A promoter-level mammalian expression atlas**, *Nature* 2014;507:462–470.
- <https://fantom.gsc.riken.jp/5/data/>
- **Licence: CC BY 4.0.** "FANTOM5 and FANTOM6 data by RIKEN is licensed under a Creative Commons
  Attribution 4.0 International License" (<https://fantom.gsc.riken.jp/5/data/>; the data-file
  index also declares `"license":"http://creativecommons.org/licenses/by/4.0/"`).
- Panel: the 41 libraries named `<tissue>, adult, poolN` (adipose … uterus), one library per
  tissue. Counts are converted to tags per million with each library's own `01STAT:MAPPED` total
  before averaging, so a big library does not dominate.
- **Promoter usage is tissue-dependent.** The weight shipped in the data is the panel mean;
  `tss[i].weight_source.tissue` and `data/tss_weights.json` carry the min–max of each promoter's
  share across the 41 tissues, which is wide for some genes (TP63: either promoter can be the only one used, share 0.000–1.000; CDKN2A T1 0.500–1.000, median 0.862; ESR1 T5 0.000–0.759).

## Matching, and what is deliberately left null

`splice_data.json` positions are 0-based on a per-gene pre-mRNA axis in transcript orientation, so
`gene.locus` (1-based GRCh38, inclusive) converts them back to the genome: `+` strand
`g = locus_lo + pos`, `−` strand `g = locus_hi − pos`.

A CAGE peak matches an annotated TSS when that TSS lies inside the **peak interval expanded by
±50 nt** on the same strand (`--window`, default 50). Where two or more annotated starts fall in
the same peak they are closer together than CAGE resolves (TP53 T1/T2 are 3 nt apart, DDIT3 T4/T5
1 nt), so the peak's signal is **split equally** between them; that is a stated modelling choice,
flagged per TSS as `weight_source.shared_peak`, not a measurement of either start alone.

A TSS with no peak of its own keeps **`weight: null`** — never an estimate. `model_exact.js` then
reports, in `result.tssWeighting`, either `mode: 'measured'` (weights used, unmatched promoters get
0 and the rest renormalised) or, if a gene has no weight at all, `mode: 'uniform (assumed)'`.

Match rate at ±50 nt: **36/57 annotated TSSs (63%)**, 12/12 genes with at least one match.
Weights are insensitive to the window: widening to ±100/±200 nt raises the match rate to 42/44 of
57 and moves the largest weights by <0.02 (ERBB2 T5 0.874 → 0.892 → 0.886).

## Rerun

```bash
# 1. download the FANTOM5 files (~121 MB, not committed)
bash pipeline/tss/fetch.sh ~/.cache/utrrace/fantom5

# 2. pull out the CAGE peaks that sit in the 12 loci (small intermediate, ~2 MB)
~/miniconda3/envs/utrrace/bin/python -I pipeline/tss/extract_cage.py \
  --splice-json data/splice_data.json \
  --bed    ~/.cache/utrrace/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz \
  --counts ~/.cache/utrrace/fantom5/hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz \
  --out    /tmp/cage_loci.json

# 3. weights -> both copies of splice_data.json + the provenance report
#    (--dry-run first: prints the per-gene match table and writes nothing)
~/miniconda3/envs/utrrace/bin/python -I pipeline/tss/weight_tss.py \
  --cage /tmp/cage_loci.json \
  --splice-json data/splice_data.json \
  --splice-json app/public/splice_data.json \
  --report data/tss_weights.json --window 50
```

`weight_tss.py` is a **merge** step, not a build step: it reads whatever `splice_data.json` already
contains and rewrites only `tss[i].weight`, `tss[i].weight_source`, `meta.sources[fantom5_cage]`
and `meta.methods.tss_weight`. Exons, introns, scores, isoforms, `measured` and `splice_findings`
are passed through untouched, so it is safe to re-run after `pipeline/splice/build_splice.py`
regenerates the file — and it **must** be re-run after it, or the weights are lost.
