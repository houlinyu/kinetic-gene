"""Assemble data/splice_ridge.json from work/splice_ridge_core.json + work/app12.json, add
meta.methods / meta.caveats, and copy to app/public/."""
import os, json, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
W = os.path.join(HERE, "work")
core = json.load(open(os.path.join(W, "splice_ridge_core.json")))
app12 = json.load(open(os.path.join(W, "app12.json")))
core["models"]["usage"]["app12"] = dict(
    n=len(app12["rows"]), n_genes=len({r["gene"] for r in app12["rows"]}),
    note=("the app's mechanistic co-transcriptional splicing model can only be evaluated on the 12 genes "
          "for which data/splice_data.json carries a hand-built gene structure, so this block is restricted "
          "to them; MaxEntScan donor/acceptor scores here are the ones already in data/splice_data.json "
          "(maxentpy, Yeo & Burge 2004); the measured target there is every GENCODE v49 transcript "
          "of the 12 genes from the genome-wide Isopedia query table, and `ridge_compact` is the "
          "genome-wide ridge refit on only the features derivable from data/splice_data.json, with "
          "these 12 genes' loci removed from its training set"),
    methods=app12["metrics"], per_isoform=app12["rows"])
core["meta"] = dict(methods=[
 {"key": "target",
  "text": "Measured isoform usage from the Isopedia human long-read isoform atlas "
          "(/data/mdl_data/hyu/projects/public_isoform/isopedia, table_wide2 build - the same build "
          "pipeline/isopedia/ uses: an 18,839-gene panel, 105,837 catalogue isoforms, 869 public PacBio/ONT "
          "samples, 32.7M (isoform, sample) rows in table_wide2/context_table.parquet). Each row's usage is "
          "that splice-pattern isoform's share of its gene's reads in that sample; transcripts sharing an "
          "intron chain are collapsed and whatever is left over sits in an other:<gene> row. Per transcript "
          "we sum usage over all samples and divide by the gene total INCLUDING that other row, so a gene's "
          "modelled fractions sum to 1 minus the unmodelled remainder rather than being inflated to 1. Loci "
          "are kept only if seen in >=20 samples and if >=2 of their transcripts survive the filters."},
 {"key": "features",
  "text": "Per-transcript, all interpretable. Exon/intron counts and lengths, alternative-event flags "
          "(cassette skip / alt 5'ss / alt 3'ss / alt first exon, each relative to the gene's longest "
          "transcript as a structural reference chosen without looking at usage) and the distance from the "
          "gene's dominant annotated TSS are computed here from table_wide2/isoforms.parquet exon "
          "coordinates. GC, 5' UTR and ORF length, uAUG count, NMD/PTC flag, intron retention, TSS CpG/TATA "
          "and polyA-signal flags come from table_wide2/catalog_features.parquet."},
 {"key": "ss_score",
  "text": "MaxEntScan is recomputed genome-wide here with the same scorer pipeline/splice/build_splice.py "
          "uses - maxentpy (X. Dong), which ships the original Yeo & Burge (2004) matrices: score5 on the "
          "9-mer -3..+6 around the donor, score3 on the 23-mer -20..+3 around the acceptor. Both documented "
          "reference examples are re-validated at the start of every run and the run aborts if they do not "
          "reproduce. The windows are read by random access out of a faidx-indexed GRCh38 FASTA rather than "
          "from the Ensembl REST API, which does not scale to 730k introns; the exon-coordinate convention "
          "is chosen empirically as the offset that yields canonical GT..AG, and the chosen offset and its "
          "GT..AG rate are recorded in pipeline/splice_ridge/work/maxent_meta.json. Pangolin splice-site "
          "probabilities from table_wide2/splice (meta.json model=pangolin) are kept "
          "alongside as separate features, never merged with the MaxEntScan ones."},
 {"key": "fit",
  "text": "RidgeCV over 40 alphas in 1e-3..1e4 on z-scored features; target is the raw usage fraction "
          "(not its logit) so that ridge, the mechanistic model and every baseline live on one scale. "
          "Cross-validation is GroupKFold(5) grouped by gene, so no gene appears in both train and test. "
          "Coefficient CIs are 2.5/97.5 percentiles of 200 bootstraps that resample whole genes."},
 {"key": "mechanistic",
  "text": "expected(gene, values()) from app/src/splice/model_exact.js after calibrate(data), at the app's "
          "default parameters, run under node (pipeline/splice_ridge/mech_eval.mjs). Nothing in it is fitted "
          "to the Isopedia usage."},
 {"key": "within_gene",
  "text": "Every pair of transcripts inside a gene whose measured usage differs by more than 1e-6; a method "
          "scores 1 if it orders the pair correctly, 0.5 if it ties them. 0.5 is therefore the "
          "no-information baseline."}],
 caveats=[
  "Splice-site strength enters twice and the two are not interchangeable: MaxEntScan scores the 9-mer/23-mer "
  "windows alone, while Pangolin already encodes much of the wider sequence context.",
  "The genome-wide fit is restricted to transcripts with a Pangolin score, i.e. the Isopedia catalogue's "
  "spliced multi-exon isoforms; single-isoform loci carry no usage information and are dropped.",
  "Usage is pooled over a heterogeneous panel of public long-read samples (many tissues, two platforms). "
  "It is a panel-average usage, not a usage in any one tissue, and long-read 5'/3' completeness differs "
  "between platforms.",
  "The mechanistic-model comparison is on 12 genes only, which is far too few to separate the methods; "
  "treat its numbers as a sanity check, not a benchmark.",
  "'Read support' (number of samples in which a transcript is detected) is NOT an independent confound: it "
  "is computed from the same per-sample usage values as the target, so a transcript seen in more samples "
  "necessarily has a larger usage sum. It is reported to show how much of the within-gene signal is simply "
  "abundance, and must not be read as a model that beats the others. The non-circular confounds are "
  "transcript length and exon count.",
  "CAGE promoter weights (data/tss_weights.json, FANTOM5) exist only for the app's genes, so the "
  "genome-wide promoter baseline is the annotation-derived dominant-TSS share instead."])
out = os.path.join(ROOT, "data", "splice_ridge.json")
json.dump(core, open(out, "w"), indent=1)
shutil.copy(out, os.path.join(ROOT, "app", "public", "splice_ridge.json"))
print("wrote", out, os.path.getsize(out), "bytes -> app/public/")
