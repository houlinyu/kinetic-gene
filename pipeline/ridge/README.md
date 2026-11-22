# pipeline/ridge/ — interpretable ridge models of 5' UTR translation (RIDGE agent)

Writes `data/ridge.json` (contract: `data/SPEC_v2.md`) and figures in `data/figures/`.
Features come only from `pipeline/motifs.py` `features(utr, cds_head, fold)` (read-only import).

## Env (one-time)
```bash
~/miniconda3/bin/conda create -y -n utrrace_ml --clone utrrace
~/miniconda3/envs/utrrace_ml/bin/pip install scikit-learn pandas matplotlib
```

## Rerun (from repo root; ~15-25 min total, 8 cores)
```bash
cd ~/kinetic-gene/pipeline/ridge
PY="nice -n 19 $HOME/miniconda3/envs/utrrace_ml/bin/python"
mkdir -p downloads && curl -o downloads/gencode.v19.pc_transcripts.fa.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
curl -o downloads/gencode.v47.pc_transcripts.fa.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz
$PY prep.py mrl te                 # feature tables -> work/{mrl,te}.csv.gz   (NJ=8 workers)
NJ=4 $PY physics_prep.py mrl       # app physics model (node, app/src/analytic.js) at default params
NJ=4 $PY physics_prep.py te
$PY physics_fit.py                 # fit 4 physics params to MRL (5-fold CV) -> work/physics_fit.json
NBOOT=200 $PY fit.py               # CV, bootstrap CIs, figures, data/ridge.json
cd ../.. && nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/build_data.py --merge-only
```

## Data
- **MRL**: GEO GSM4084997 (Sample et al. 2019, varying-length 25-100 nt eGFP library; already in
  `pipeline/external/geo/`). Filter `total_reads >= 10` (as in `Varying_length_25to100_modeling.ipynb`). Rows = all
  `human` fragments + the notebook's random-library test split (top 100 by reads per length). Optimus
  (`Varying_length_25to100_model.hdf5`) was trained on the random training split only, so it never saw these rows.
  CDS head = eGFP start in the reporter, `ATGGGCGAATTAAGTAAGGGCGAGGAGCTGTTCACC` (Sample 2019 Online Methods, PMC7100133).
- **TE**: Cao et al. 2021 HEK293T table (`pipeline/external/df_counts_and_len.TE_sorted.HEK_Andrev2015.with_annot.txt`)
  x GENCODE v19 pc_transcripts (versioned IDs match 77.6k/79.1k rows). Filters: rnaseq_counts >= 50,
  riboseq_counts >= 10, te > 0, CDS starts ATG, 5' UTR >= 10 nt, drop transcripts whose TE is bit-identical to a
  sibling isoform (unresolved quantification). Target log2 TE. CV grouped by gene.
- **Fold** (`feats.py`): ViennaRNA on UTR + first 30 nt CDS; punp = RNAplfold-style (`pfl_fold_up`, W=150, L=100);
  mfe = UTR MFE with max bp span 100; cap_dG = MFE of first 50 nt.

## Models / evaluation (`fit.py`)
StandardScaler + RidgeCV (alphas 1e-3..1e4) inside each fold; pooled held-out r2. Optimus and physics scores are
mapped to the target by a 1-D linear fit on the training folds (same folds). Coefficients = final ridge on all rows
(alpha = median fold alpha), 95% CI from bootstrap (rows for MRL, genes for TE). `n_uAUG` dropped (exact sum of the
three context-specific counts). Within-gene test: all same-gene pairs with different 5' UTR sequences, out-of-fold
predictions (whole gene held out), fraction of pairs ordered correctly (ties excluded), also for |dlog2TE| > 1.
`per_isoform`: features of each app isoform (CDS from build_data's cached Ensembl REST cDNA, else GENCODE v47, else
the 4 nt of `cds_context`), contribution = coef x standardized value; `contrib` = MRL model, `contrib_te` = TE model.
