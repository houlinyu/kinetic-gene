# pipeline/splice/ — gene models + splice-site strengths + measured PSI for tab 2

Produces `data/splice_data.json` (contract: `data/SPEC_SPLICE.md`) and copies it to
`app/public/splice_data.json`. Owned by the SPLICE-DATA agent; nothing here touches tab 1.

## Rerun

```bash
cd ~/kinetic-gene
# one-time, into the existing utrrace env:
~/miniconda3/envs/utrrace/bin/pip install "git+https://github.com/kepbod/maxentpy.git"   # MaxEntScan
mkdir -p pipeline/splice/external && curl -sSL -o pipeline/splice/external/gtex_psi.csv.gz \
    https://snaptron.cs.jhu.edu/data/ascot/gtex_psi.csv.gz                               # ASCOT PSI (20 MB)

nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/splice/build_splice.py          # all genes (~40 s cached)
nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/splice/build_splice.py --quick   # first 2 genes only
nice -n 19 ~/miniconda3/envs/utrrace/bin/python pipeline/splice/build_splice.py MAPT TP53 # named genes
#   --no-cache   ignore pipeline/splice/cache/ (cached Ensembl REST responses) and refetch
```

The script prints a per-gene line (exon-type counts, introns, scored sites, TSSs, how many isoforms link to
tab 1, measured PSI values) and **exits non-zero if any validation fails** — nothing is written in that case.

To change the gene set, edit `pipeline/splice/genes_splice.json`:
`from_tab1: true` pulls that gene's isoform ids straight out of `data/utr_data.json` (which is what makes the
two tabs join), `extra` adds unversioned ENSTs, `event`/`cite`/`iso_names` are curated and must be backed by a
`cite` key that exists in `SOURCES` in `build_splice.py`.

## What it does

- **Gene models** — Ensembl REST (GRCh38), `/lookup/id/<ENST>?expand=1;utr=1`, cached in `pipeline/splice/cache/`
  exactly like `pipeline/build_data.py` (same generic User-Agent, same rate limit, md5-of-request filenames).
  The GENCODE GTF in `pipeline/discovery/dl/` was deliberately **not** used: it is v47 (= Ensembl 113) while tab 1
  is on Ensembl 116, so several tab-1 transcript ids are absent from it and the two tabs would not join.
- **Coordinates** — one pre-mRNA axis per gene, 0-based half-open, **transcript orientation**: 0 = the 5'-most
  annotated base of the locus. `gene.locus`/`gene.span` give the GRCh38 span, `exon.g` the exon's GRCh38
  coordinate (1-based inclusive). `intron.start` = first intron base = `donor.pos`; `intron.end` = first base of
  the next exon = `acceptor.pos`. Minus-strand genes are flipped, so everything reads 5'->3' left to right.
- **Exon classification** — derived by comparing the shipped transcripts (see `meta.methods.exon_type`), not a
  curated annotation: `first`/`last` beat `alt5`/`alt3` (shared 5' or 3' boundary with a different partner
  boundary) beat `cassette` (internal, absent from a transcript whose span covers it) and `constitutive`.
  Each exon also carries `transcripts` (which isoforms contain it) and `len`.
- **Documented exon names** — `EXON_ALIASES` in the script maps a GRCh38 coordinate to the exon number used in the
  primary literature (MAPT exons 2/3/10, ERBB2 exon 16, CDKN2A exons 1alpha/1beta, TP53 exon 4), each with its own
  `alias_cite`. Coordinate-keyed so it survives a change of transcript set; validation fails if an alias stops
  matching an exon of the expected class.
- **Splice-site strength** — MaxEntScan (Yeo & Burge 2004) via `maxentpy` 0.0.4, which ships the original
  maximum-entropy matrices. `score5` on the 9-mer -3..+6 (3 exonic + 6 intronic nt), `score3` on the 23-mer
  -20..+3 (20 intronic + 3 exonic nt), in bits. **Validation**: the build refuses to run unless `maxentpy`
  reproduces its documented reference examples to 1e-9 — `score5('cagGTAAGT') = 10.858313101356437` and
  `score3('ttccaaacgaacttttgtAGgga') = 2.8867730651152104`
  (https://github.com/kepbod/maxentpy). The Burge-lab web server (hollywood.mit.edu/burgelab/maxent) was
  unreachable from this machine, so those are the implementation's published examples, not a fresh run of the
  original Perl — `meta.methods.ss_score` says so. Independent sanity check on our own data: across all 150
  introns the median donor score is 8.4 and the median acceptor 9.2 bits, and 149/150 donors begin GT and
  150/150 acceptors end AG, which is what correct windows look like. The one exception is a real GC-AG intron in
  DDIT3 (donor `GAGGCGAGT`); the build prints it as a note rather than failing.
  Branch points are **not** scored: `intron.bp.pos` and `intron.bp.score` are always null.
- **Raw site sequences** — `intron.donor.seq` / `intron.acceptor.seq` are the exact scored windows, already in
  transcript orientation, so the UI can display them. Null (and `score` null) if the window runs off the fetched
  region or contains a non-ACGT base.
- **Measured PSI** — ASCOT GTEx cassette-exon table (Ling et al. 2020), joined to our exons by exact GRCh38
  coordinate + gene symbol; percentages divided by 100, `-1` (not measured) dropped. ASCOT's PSI file has no
  per-tissue sample count, so `n` is null (counts live in `ASCOT_accession_metadata.xlsx`, not parsed).
  VastDB (Tapial et al. 2017) whole-compendium PSI averages for MAPT exons 2 and 3 are hand-entered from the
  gene page with the exact values in a `note`. `measured.isoform_fraction` holds the one measured isoform
  fraction found (ERBB2 Delta16HER2, ~9% of HER2 mRNA, n = 46 breast carcinomas).
  Alternative-first-exon / alternative-promoter events (CDKN2A, TP63, TP53 delta133, MDM2 P1/P2) are in neither
  table, so their PSI is **null, not estimated**.
- **`splice_findings`** — top-level list of literature "which form is used more" pairs, read by the model's
  `allSpliceFindings()`. `higher`/`lower` are shipped isoform ids where the paper maps cleanly to a transcript,
  otherwise an event label; `higher_exon_id`/`lower_exon_id` point at the exon when the event is an exon.
  `value` is the paper's number or null.
- **Validation** (`validate()`, runs before anything is written) — exon positions ordered and inside `0..span`;
  exon types legal; every exon claimed by >= 1 transcript; TSS positions in range; every isoform's `tss_id`
  known, exon chain non-empty, non-overlapping and in order, and consistent with each exon's `transcripts`;
  every intron adjacent in some isoform's exon chain, with `donor.pos`/`acceptor.pos` equal to its bounds and
  9-/23-nt sequences; PSI values in 0..1; every `cite` key anywhere in the document present in `meta.sources`;
  every exon alias still matching an exon of its expected class; finding exon ids real.

## Caveats

- Exon types come from the transcripts we ship, so adding or removing a transcript can reclassify an exon.
  Loci where tab 1 picked fragmentary transcripts (ERBB2, MDM2) therefore show several `first` exons.
- `meta.release` is whatever Ensembl REST reports today; tab 1 was built on Ensembl 116.
- ASCOT and VastDB state no licence; cite Ling et al. 2020 / Tapial et al. 2017 respectively.
