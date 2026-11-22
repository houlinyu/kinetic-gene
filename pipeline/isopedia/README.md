# Measured isoform usage from Isopedia — and what it says about the splicing model

The ① Transcription page predicts **isoform fractions**. The validation that shipped first compared
predicted PSI on 17 alternative exons (ASCOT/GTEx, Pearson r = 0.16), which is the wrong unit: the
model's output is a distribution over transcripts, not an exon inclusion level. This directory
replaces it with an apples-to-apples test against **measured isoform usage** from Isopedia, a
read-level index of 1,007 public human long-read RNA-seq samples.

## What was read (read-only, nothing written inside it)

All under `/data/mdl_data/hyu/projects/public_isoform/isopedia` (the user's own public_isoform project):

| file | what it holds |
|---|---|
| `table_wide2/context_table.parquet` | 32.7 M rows: `usage` per (sample, catalog isoform) with `gene_name`, `locus_id`, `lineage`. The wide build: 18,839-gene panel, 105,837 catalog isoforms, 869 samples. |
| `table_wide2/isoforms.parquet` | catalog isoform → `locus_id`, exon coordinates, one representative GENCODE transcript |
| `table_wide2/sample_meta.tsv` | platform (547 ONT / 322 PacBio), project, tissue, depth per sample |
| `ref/gencode.v49.primary_assembly.annotation.gtf.gz` | exon structures, to place our ENSTs onto the catalog by intron chain |
| `build_wide2.log`, `build_v2.log` | the build parameters quoted below |
| `~/IsoMap-RNAFM/isoform_usage_benchmark/tools/isopedia_context_table.py` | the builder's docstring, which is the authoritative description of how usage is defined |

The narrow build (`table/`, 4,288-gene panel) covers only 5 of our 12 genes, so `table_wide2/`
(gene measured in ≥ 2 % of samples instead of ≥ 50 %) is used instead. Nothing was recomputed from
reads; everything here is an aggregation of an existing table.

## How Isopedia's usage is defined (quoted from the builder)

`isopedia isoform -g gencode.v49.gtf` returns, per GENCODE v49 transcript and per sample, supporting
read counts: `RC_FSM_JC` (reads whose junction chain fully matches) and `RC_EM_ISM` (partial reads
assigned by EM). Transcripts sharing an **intron chain** are collapsed into one *splice-pattern
isoform*, which takes `max(RC_FSM_JC) + sum(RC_EM_ISM)`. **usage = isoform reads / gene reads** in
that sample. A (sample, gene) is measured at ≥ 20 gene reads; samples need ≥ 200,000 indexed reads;
a catalog isoform must reach usage ≥ 0.05 in ≥ 5 % of measured samples, and everything else is pooled
into one `other:<gene>` row. **Per gene and sample, catalog usage + other = 1.**

## Rerun

```bash
cd /home/unix/yuhoulin/kinetic-gene
PY=~/miniconda3/envs/fungi/bin/python          # pandas, pyarrow, scipy, matplotlib
nice -n 19 $PY -I pipeline/isopedia/build_isoform_usage.py          # -> data/isoform_usage.json (~2 min)
nice -n 19 ~/miniconda3/envs/web/bin/node pipeline/isopedia/predict.mjs > data/figures/isopedia_pred.json
nice -n 19 $PY -I pipeline/isopedia/validate_usage.py               # -> validation CSVs, prints the tables
nice -n 19 $PY -I pipeline/isopedia/figures_usage.py                # -> isopedia_scatter.png, isopedia_bars.png
cp data/isoform_usage.json app/public/isoform_usage.json
```

`build_isoform_usage.py --table table` switches to the narrow 4,288-gene build;
`--min-tissue-samples N` sets how many samples a tissue needs to get its own breakdown (default 10).

## Joining Isopedia to the app

The app is built on **Ensembl 116**; Isopedia is annotated with **GENCODE v49**. 7 of the 57 transcript
ids in the 11 measured genes do not exist in v49 at all. So the join is **structural, not by id**: each app isoform's intron
chain (from its own genomic exon coordinates, or from its v49 entry when the id exists) is matched to
a catalog isoform's intron chain, exactly first, then allowing ±10 bp per boundary to mirror the
near-identical-chain merge the builder performs. `data/figures/isopedia_join.csv` records, per app
isoform, which route matched and which did not.

**Consequence:** two app transcripts that differ only in their ends can land on one catalog isoform
(PPP1R15A-201 and -206 do). Isopedia cannot separate them, so in the validation their predicted
fractions are summed and they share one measured value.

## Coverage

- 11 of 12 genes are measured. **ESR1 is absent** from the 18,839-gene panel entirely — fewer than
  20 supporting reads in ≥ 2 % of the 869 samples.
- **33 of 63 app isoforms** map onto a catalog isoform; 32 comparison units after collapsing.
  Of the 30 that do not: 6 are ESR1's (no data at all), and 24 are isoforms of measured genes whose
  splice pattern never reached 5 % usage in 5 % of samples, so Isopedia pooled them into `other` —
  the measurement says "rare", without giving a number. 17 of those 24 are annotated in v49, so the
  chain matching is not the reason they are missing; 7 are Ensembl-116-only ids whose chain is also
  not in the catalog.
- Per gene, the app's isoforms cover between 4 % (ERBB2, before the ±10 bp merge; 44 % after) and
  99 % (PPP1R15A) of the gene's reads — `app_isoform_share_of_gene` in the JSON.

## Result (the headline)

`data/figures/isopedia_validation_summary.csv`, 32 isoforms over 11 genes, both sides renormalised
within gene over the compared set:

| predictor | Pearson r | Spearman | MAE | within-gene pair accuracy |
|---|---|---|---|---|
| kinetic splicing model | 0.55 | 0.44 | 0.26 | **0.49** (17.5/36) |
| CAGE promoter weights only, no splicing | **0.78** | **0.71** | **0.19** | **0.68** (24.5/36) |
| uniform over annotated isoforms | 0.56 | 0.47 | 0.20 | 0.53 (34 of 36 pairs are ties → 0.5) |

**The splicing layer does not add anything; it subtracts.** Stripping the kinetics and keeping only
the measured CAGE promoter weights the model already consumes is better on every metric, and the
model's within-gene ranking is at the coin-flip baseline. Per gene
(`data/figures/isopedia_per_gene.csv`) the model inverts MAPT (predicts 2N4R dominant, measured says
0N3R), MDM2 (predicts the short MDM2-212/-211, measured says full-length MDM2-224/-202), CDKN2A and
TP63. n is small — 36 within-gene pairs — so the gap between 0.49 and 0.68 is suggestive, not
significant; what is solid is that the model is **not** above chance.

## Caveats

- **Quantification, not truth.** EM assignment of partial (ISM) reads, 5′/3′ degradation and library
  chemistry all move usage. Long-read usage estimates for genes with many similar isoforms (MDM2,
  TP53) are the least trustworthy.
- **Read-support thresholds hide the tail.** Anything below 5 % usage in 95 % of samples is invisible,
  so the comparison is biased toward the isoforms Isopedia could see.
- **Convenience sample.** 869 public samples: tumours, cell lines and tissues mixed, ONT and PacBio
  mixed, unbalanced across tissues (Blood 168, Brain 69, MAPT measured in only 28 samples). "Overall"
  is a mean over that sample, not over human tissue.
- **Annotation mismatch.** GENCODE v49 vs Ensembl 116; structural joining fixes the ids but a v49
  isoform our app does not model still absorbs real mass (ATF4's most-used splice pattern,
  `ENST00000960426`, is not one of the app's six).
- **Denominator.** Measured fractions are shares of *all* gene reads; the app's isoforms do not sum
  to 1. The headline metrics renormalise both sides within the compared set; `pearson_r_unnormalised`
  in the summary CSV is the un-renormalised version (model 0.46, promoter-only 0.65).
- MAPT, TP63, ERBB2 and ATF5 rest on 28–292 samples; ATF4, DDIT3, PPP1R15A, TP53 on 463–731.

## Files produced

- `data/isoform_usage.json`, copied to `app/public/isoform_usage.json`
- `data/figures/isopedia_scatter.png`, `isopedia_bars.png`
- `data/figures/isopedia_validation.csv` (per compared isoform), `isopedia_validation_summary.csv`,
  `isopedia_per_gene.csv`, `isopedia_usage_catalog.csv` (every catalog isoform of the 11 genes,
  whether or not the app models it), `isopedia_join.csv`, `isopedia_pred.json`

## JSON shape (for whoever wires the UI)

```
meta.{source, dataset_description, how_quantified, samples_or_tissues[], n, denominator,
      caveats[], citation_or_internal, generated, join:{key, n_app_isoforms, n_matched}}
genes.<GENE>.measured                       false for ESR1 (then only `note` is present)
genes.<GENE>.n_samples                      samples where this gene was measured
genes.<GENE>.n_catalog_isoforms             splice patterns Isopedia quantified for this gene
genes.<GENE>.n_app_isoforms / n_app_isoforms_matched
genes.<GENE>.overall.<ENST>                 mean usage over samples, 0-1, UNVERSIONED ENST
genes.<GENE>.overall_median.<ENST>          median over samples (robust; does not sum to 1)
genes.<GENE>.tissues.<Tissue>.{n_samples, usage.<ENST>}
genes.<GENE>.pooled_other_fraction          mass in the 'other' bin
genes.<GENE>.unmatched_catalog_fraction     mass in catalog isoforms the app does not model
genes.<GENE>.app_isoform_share_of_gene      mass the app's isoforms account for
```

Key is the **unversioned** ENST: join with `isoforms[].id.split('.')[0]`. A missing key means the
isoform is not in Isopedia's catalog — display it as "below detection (< 5 % usage)", not as 0.
