# pipeline/ — offline data build for UTR Race

Produces `data/utr_data.json` (contract: `data/SCHEMA.md`) and copies it to `app/public/utr_data.json`.

## Rerun

```bash
# env (one-time): ~/miniconda3/bin/conda create -y -n utrrace python=3.11
#                 ~/miniconda3/envs/utrrace/bin/pip install ViennaRNA requests numpy h5py pandas
cd ~/kinetic-gene
nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/build_data.py            # genes from pipeline/genes.txt
nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/build_data.py ATF4 MYC   # genes from CLI
#   --max-iso N    isoforms per gene (default 6)
#   --no-cache     ignore pipeline/cache/ (cached Ensembl REST responses) and refetch
#   --merge-only   reload data/utr_data.json, merge data/ridge.json (isoform.ridge, meta.ridge), write both copies
nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/test_motifs.py          # motif unit tests
```
To add genes, append symbols (one per line, `#` comments ok) to `pipeline/genes.txt`. Full 8-gene run: ~30 s
uncached; folding a ~2.8 kb UTR takes ~5 s. The script prints a per-gene table and exits non-zero if a sanity
check fails (layout2d / dotbracket / punp length == len, pair span <= 100, uORF/hairpin/motif positions in range and
motif seq matches, one kozak_main, motif cite in sources, start codon at cds_context[6:9]).

## What it does
- **Transcripts**: Ensembl REST (GRCh38) `/lookup/symbol` -> `/overlap/id/<gene>?feature=transcript` (biotype + tags;
  `gencode_basic` is written as `basic`) -> `POST /sequence/id?type=cdna;mask_feature=1` (5' UTR = leading lowercase).
  Kept: protein_coding, non-empty 5' UTR, no `*start_NF` tag, start codon = Ensembl-annotated start. Non-AUG annotated
  starts (NUG; e.g. MYC and VEGFA MANE use CUG) are kept and flagged with a tag `non_AUG_start:CUG`.
- **Dedup/selection**: identical UTR sequences collapsed to the best-tagged id (MANE_Select > Ensembl_canonical >
  basic). MANE_Select/canonical always kept, then greedy max-min 6-mer-Jaccard distance picks the most different UTRs
  up to `--max-iso`.
- **uORFs**: every AUG in the UTR scanned to its first in-frame stop (continuing into CDS/3'UTR). `end` is exclusive and
  includes the stop codon (may exceed `len` for overlapping ones); `aa` excludes the stop. `type`: `uORF` = stop inside
  the UTR; `overlapping` = stop past the main start, out of frame; `in-frame-extension` = in frame with the main start
  and no stop before it. If no stop exists in the cDNA, `end` = cDNA length and `aa` = null.
- **Kozak**: rule in `main_kozak.method` (strong = -3 purine AND +4 G; adequate = one; weak = neither). `score` = null.
- **Fold (v2, local)**: `punp` = RNAplfold unpaired probabilities (`RNA.pfl_fold_up`, W=150, L=100) on UTR + 30 nt
  CDS; `acc_start` = mean punp over -15..+14 around the start. Structure (`dotbracket`, `layout2d`, `hairpins`) = MFE
  with `max_bp_span=100` (an MEA from a span-limited partition function took ~55 s for 2.8 kb, so the span-limited MFE
  is used). `mfe`/`mfe_global` = unrestricted global MFE (kept); `mfe_local` = MFE of the local structure; `cap_dG` =
  MFE of the first 50 nt; `layout2d` = naview coordinates, centred and scaled to [-1, 1]. Full build ~60 s.
- **Motifs** (`pipeline/motifs.py`): `find_motifs(utr, cds_head)` -> `isoform.motifs`; `features(...)` = flat numeric
  features for the ridge agent (`FEATURES` = name -> label; `SOURCES` = cite keys, merged into `meta.sources`). Types
  and sources in `data/SCHEMA.md`; tests in `pipeline/test_motifs.py`.
- **Ridge**: if `data/ridge.json` exists its `per_isoform` entries become `isoform.ridge` and the rest `meta.ridge`.
- **Hairpins**: from each hairpin-closing pair, extend outward through stacks/bulges/interior loops until a multiloop
  or exterior-loop branch; keep stems with >= 3 bp; `dG` = `eval_structure` of that substring + substructure.
- **Blurbs / findings**: `pipeline/gene_notes.json` (blurbs written from the cited abstracts; sources merged into
  `meta.sources`; `findings` = literature pairwise comparisons mapped to Ensembl ids, checked to exist in the output).
- **measured** (`pipeline/measured.py`): Cao et al. 2021 TE tables in `pipeline/external/` (from
  github.com/zzz2010/5UTR_Optimizer, MIT), matched by unversioned ENST; `ambiguous` when a sibling isoform has the
  identical TE.
- **pred.optimus_mrl** (`pipeline/optimus.py`): numpy forward pass of Optimus 5-Prime `Varying_length_25to100_model.hdf5`
  (pjsample/human_5utr_modeling, GPL-3.0; weights + notebook in `pipeline/external/optimus/`). Input = whole UTR if
  <= 100 nt, else last 100 nt before the start codon; left-padded with N. Output un-scaled with the training-set MRL
  mean/sd. `python pipeline/validate_optimus.py` recomputes that scaler from GEO GSM4084997
  (`pipeline/external/geo/`) and reproduces the notebook's test r^2 (0.838 random / 0.778 human; ours 0.8375 / 0.7772).
  Needs `pip install h5py pandas` in `utrrace`.
- **tss_note**: Ensembl `lookup/id?expand=1` exon coordinates; heuristic described in `meta.methods.tss_note`.
- `utrlm_te` null.
