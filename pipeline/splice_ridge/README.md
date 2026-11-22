# `pipeline/splice_ridge/` — interpretable ridge on isoform usage

The ① Transcription counterpart of `pipeline/ridge/`: the same three-way comparison the
Translation page has (mechanistic model vs interpretable ridge vs measured data), but the
measured quantity is **a transcript's share of its gene** instead of reporter MRL.

Owns: `data/splice_ridge.json`, `app/public/splice_ridge.json`, `data/figures/splice_ridge_*`.
Reads (never writes) `/data/mdl_data/hyu/projects/public_isoform/isopedia/table_wide/`,
`data/splice_data.json`, `data/tss_weights.json`, `app/src/splice/model_exact.js`.

## Measured data

`table_wide2/context_table.parquet` from the user's Isopedia human long-read isoform atlas —
the same build `pipeline/isopedia/` uses, so the two agree: the 18,839-gene panel, 105,837
catalogue isoforms, 869 public PacBio/ONT samples, 32.7 M (isoform, sample) rows. Each row's
`usage` is that splice-pattern isoform's share of its gene's reads in that sample; transcripts
sharing an intron chain are collapsed, and whatever is left over sits in an `other:<gene>` row,
so catalogue usage + other = 1 per gene and sample (see `pipeline/isopedia/README.md`, which
quotes the builder).

Per transcript we sum `usage` over all samples and divide by the gene total **including the
`other` row**, so a gene's modelled fractions sum to 1 − (unmodelled remainder) rather than being
inflated to 1. `other` rows carry no features and never enter the fit. A locus is kept only if it
was seen in **≥ 20 samples** and if **≥ 2** of its transcripts survive the feature filters.

Set `ISOPEDIA_BUILD=table_wide` to run the whole thing against the narrower 4,288-gene build
instead. Because Isopedia is transcriptome-wide, the fit is **not** restricted to the app's 12 genes.

## Features (all interpretable, all per transcript)

Computed in `prep_features.py`; the human labels are the `FEATURES` dict at the top of
`fit.py` and are carried into `data/splice_ridge.json` under `models.usage.features[].label`.

| source | features |
|---|---|
| exon coordinates (`isoforms.parquet`) | `n_exons`, `log_tx_len`, `log_mean_exon`, `log_mean_intron`, `log_min_intron`, `log_max_intron`, `log_dist_dom_tss`, `is_dom_tss`, `is_longest`, `n_iso_locus_log`, and the alternative-event flags `cassette_skip`, `alt_5ss`, `alt_3ss`, `alt_first_exon` |
| `catalog_features.parquet` | `gc`, `utr5_len_log`, `orf_len_log`, `n_uatg`, `nmd_ptc`, `intron_retention`, `likely_5p_fragment`, `tss_gc`, `tss_cpg_island`, `tss_tata`, `tes_pas_canonical`, `frac_share_tss`, `frac_share_tes` |
| GRCh38 FASTA + `maxentpy` (`prep_maxent.py`) | `mes_donor_min`, `mes_donor_mean`, `mes_acc_min`, `mes_acc_mean` |
| `splice/isoform_splice_scores.parquet` (Pangolin) | `ss_p_min`, `ss_p_mean`, `ss_use_min`, `ss_use_mean` |

The alternative-event flags and `is_longest` are taken relative to the gene's **longest
transcript**, picked purely from exon structure and therefore without looking at usage.

### Splice-site strength: MaxEntScan (recomputed here) plus Pangolin

**MaxEntScan** is recomputed genome-wide by `prep_maxent.py` with the *same* scorer
`pipeline/splice/build_splice.py` uses — `maxentpy` (X. Dong), which ships the original
Yeo & Burge (2004) matrices: `score5` on the 9-mer −3..+6 around the donor, `score3` on the
23-mer −20..+3 around the acceptor. The two documented reference examples are re-validated at
the start of every run and the run aborts if they do not reproduce. `build_splice.py` pulls each
window from the Ensembl REST API one junction at a time, which does not scale to 730k introns,
so the windows are read here by random access out of a faidx-indexed GRCh38 FASTA
(`/data/mdl_data/hyu/projects/asd_villages_sqtl/refs_nochrM/genome.fa`, read-only, never loaded
into memory). The exon-coordinate convention of `isoforms.parquet` is *not* assumed: all four
candidate offsets are tried on a sample of introns and the one that yields canonical GT..AG is
kept, with the chosen offset and its GT..AG rate recorded in `work/maxent_meta.json` (offset
`(0,0)`, i.e. the stored coordinates are 1-based inclusive, 98.5 % GT..AG over 931k introns).
Per transcript this gives `mes_donor_min/mean` and `mes_acc_min/mean`.

One performance trap worth knowing: `maxent.score5/score3` reload their matrices from disk on
**every** call unless the matrix is passed in — 18 ms versus 10 µs, i.e. hours versus seconds over
931k introns. `prep_maxent.py` loads them once with `load_matrix5()/load_matrix3()` and passes them
in. Same matrices, identical scores (the reference examples still reproduce to the last digit).

**Pangolin** comes from `table_wide/splice/meta.json` = `{"model": "pangolin"}`: the shards give,
per transcript, the mean and minimum predicted splice-site probability and usage over the
junctions it uses (`ss_p_*`, `ss_use_*`). It is kept alongside MaxEntScan because it encodes
sequence context MaxEntScan does not; the two are reported separately, never merged.

## Fit and comparison

`RidgeCV` over 40 alphas in `1e-3 .. 1e4` on z-scored features. **Target = the raw usage
fraction, not its logit**, so that ridge, the mechanistic model and every baseline are on one
scale and one `r2` is comparable across them (a logit target would make the mechanistic model's
hard 0 predictions undefined). Cross-validation is `GroupKFold(5)` **grouped by gene**.
Coefficient CIs are the 2.5/97.5 percentiles of 200 bootstraps that resample whole genes.

Two compositional variants of the *same* held-out ridge predictions are reported beside it, since
the fractions of a gene must sum to 1: `ridge_renormalised` divides the held-out predictions by
their gene total, and `ridge_logit_softmax` refits the same features on `logit(usage)` under the
same folds and resolves the gene back to fractions with a within-gene softmax. Both are monotone
within a gene, so they are there to move `r2`, not the ranking score — `ridge_renormalised` still
differs from `ridge` by ~0.002 on the ranking only because clipping at 1e-6 ties predictions that
were distinct while negative.

Compared on the same held-out folds: ridge; ridge on the splice scores alone; `uniform`
(1/n isoforms); `longest_wins`; `dominant_tss_only` (gene mass split over the transcripts at the
dominant annotated TSS — the genome-wide stand-in for the CAGE baseline, since
`data/tss_weights.json` covers only the app's genes); and three deliberately trivial confounds
(transcript length, exon count, read support). The within-gene ranking score is over every pair
of transcripts of a gene whose measured usage differs, with 0.5 the no-information baseline.

The mechanistic model is `expected(gene, values())` from `app/src/splice/model_exact.js` after
`calibrate(data)`, at the app's default parameters, nothing fitted to Isopedia. It needs the
hand-built gene structure in `data/splice_data.json`, so it exists for **12 genes only**.
Note that `model_exact.js` imports `values` from `app/src/splice/params.js` — importing the
*translation* `app/src/params.js` instead silently yields NaN rates and all-zero fractions.

On those 12 genes the measured target is not the curated catalogue (which matches only 11 of the
app's isoforms) but `work/app12_usage.tsv`: every GENCODE v49 transcript of the 12 genes pulled
out of `query/gencode_v49.isoform.tsv.gz`, with usage = summed `CPM_FSM_JC` over the 1007
samples divided by the gene total. So that the ridge is a like-for-like comparison there,
`ridge_compact` is refit genome-wide on only the ten features that can also be derived from
`data/splice_data.json` (the four MaxEntScan aggregates, exon count, transcript length, isoform
count, longest-transcript flag, and the two dominant-TSS terms), **with the 12 genes' loci
removed from its training set**.

## What it says (run of 2026-11-22; `work/summary.txt` is the authoritative copy)

38,027 transcripts of 7,751 genes, 122,174 within-gene isoform pairs.

| method | r² | Pearson | Spearman | within-gene |
|---|---|---|---|---|
| ridge on logit, softmax per gene | −0.596 | 0.331 | 0.211 | **0.573** |
| ridge (held out, gene-grouped) | **0.209** | **0.457** | **0.363** | 0.561 |
| ridge, splice scores only | 0.014 | 0.120 | 0.137 | 0.538 |
| longest transcript wins | −2.280 | 0.312 | 0.211 | 0.518 |
| uniform fractions | −0.601 | 0.226 | 0.156 | 0.500 |
| dominant TSS only | −3.178 | 0.103 | 0.062 | 0.499 |
| CONFOUND transcript length | −349 | 0.072 | 0.084 | *0.548* |
| CONFOUND exon count | −7106 | −0.017 | −0.031 | 0.483 |
| CONFOUND read support | −4.8e5 | 0.029 | 0.128 | *0.780* |

**The confound check fires.** Plain transcript length orders 0.548 of the pairs correctly against
the ridge's 0.561: the whole 32-feature model is worth 1.3 points over "the longer transcript
wins". The paired head-to-head confirms the gap is real but small (ridge right / length wrong
22,052 vs the reverse 20,446, p = 6.7e-15). Read support reaches 0.780, but it is **not** an
independent confound — it is the number of samples the transcript is detected in, computed from
the same per-sample usage values as the target, so it is largely circular and is reported only to
show how much of the ranking is simply abundance.

**MaxEntScan splice-site strength barely moves isoform usage.** Its within-gene rank correlations
(`work/marginal.csv`) are +0.029 (mean donor), +0.004 (weakest donor), +0.006 (mean acceptor) and
−0.032 (weakest acceptor) — which is why `mes_acc_min` carries a *negative* ridge coefficient; that
sign is in the raw association, not an artefact of conditioning. Pangolin does better but still
small (+0.12 to +0.14). A ridge on the eight splice-score features alone reaches 0.538.

On the 12 app genes (56 transcripts, 124 pairs) the mechanistic model scores 0.653, the
like-for-like `ridge_compact` 0.661, the dominant-CAGE-promoter rule 0.629, MaxEntScan's weakest
donor 0.569 — and transcript length 0.649. At n = 124 none of these are separable; the app-gene
panel is a sanity check, not a benchmark.

## Rerun

```bash
cd /home/unix/yuhoulin/kinetic-gene
P=~/miniconda3/envs/fungi/bin/python          # has pyarrow (reads the Isopedia parquets)
Q=~/miniconda3/envs/utrrace_ml/bin/python     # has sklearn/scipy/matplotlib

R=~/miniconda3/envs/utrrace/bin/python                    # has maxentpy

nice -n 19 $P -I pipeline/splice_ridge/prep_target.py     # -> work/usage.tsv     (~6 min, 33M rows)
nice -n 19 $P -I pipeline/splice_ridge/prep_junctions.py  # -> work/junctions.tsv (~1 min)
nice -n 19 $R -I pipeline/splice_ridge/prep_maxent.py     # -> work/maxent.tsv    (~30 s)
nice -n 19 $P -I pipeline/splice_ridge/prep_features.py   # -> work/features.tsv  (~2 min)
PATH=~/miniconda3/envs/web/bin:$PATH nice -n 19 node \
    pipeline/splice_ridge/mech_eval.mjs > pipeline/splice_ridge/work/mech_pred.tsv

# measured usage of the 12 app genes, from the genome-wide query table (~10 min for the zcat)
I=/data/mdl_data/hyu/projects/public_isoform/isopedia
printf 'gene_name:%s\t' ATF4 ATF5 DDIT3 PPP1R15A BACE1 MDM2 ERBB2 ESR1 MAPT CDKN2A TP63 TP53 \
  | tr '\t' '\n' | sed 's/$/\t/' > /tmp/sr_pat.txt
nice -n 19 bash -c "zcat $I/query/gencode_v49.isoform.tsv.gz | grep -F -f /tmp/sr_pat.txt > /tmp/sr_app12.tsv"
nice -n 19 $Q -I pipeline/splice_ridge/prep_app12.py /tmp/sr_app12.tsv   # -> work/app12_usage.tsv

OMP_NUM_THREADS=8 nice -n 19 $Q -I pipeline/splice_ridge/fit.py   # -> work/{summary.txt,coefs.csv,...}
nice -n 19 $Q -I pipeline/splice_ridge/figures.py         # -> data/figures/splice_ridge_*
nice -n 19 $Q -I pipeline/splice_ridge/emit.py            # -> data/splice_ridge.json + app/public/
```

`work/` is intermediate and safe to delete; everything shipped is regenerated by the five
commands above.
