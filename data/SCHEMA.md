# utr_data.json schema (contract between pipeline/ and app/)

Written by `pipeline/build_data.py`, read by `app/` from `app/public/utr_data.json`.
All sequences are RNA (A,C,G,U), 5'->3'. All positions are 0-based, half-open, relative to the 5' UTR start (cap = 0).

```jsonc
{
  "meta": { "generated": "ISO date", "gencode_or_ensembl_release": "...", "vienna_version": "...",
            "sources": [ { "key": "sample2019", "citation": "...", "url": "..." } ],
            "methods": { "optimus": "...", "measured": "...", "tss_note": "...", "fold": "...", "motifs": "..." },  // free-text method notes
            "ridge": { "models": {...}, "figures": ... } },   // data/ridge.json minus per_isoform (RIDGE agent, SPEC_v2); absent if no ridge.json
  "genes": [
    {
      "gene": "ATF4",
      "blurb": "one sentence, sourced (cite key) or empty",
      "cite": ["key", ...],
      "isoforms": [
        {
          "id": "ENST00000...",            // Ensembl transcript id (versioned)
          "name": "ATF4-201",
          "tags": ["MANE_Select", "basic"],
          "utr5": "AUG...",                 // full 5' UTR sequence (spliced)
          "len": 123,
          "cds_context": "GCCACCAUGG",       // -6..+4 around main AUG (A of AUG at index 6)
          "main_kozak": { "score": 0.0-1.0, "class": "strong|adequate|weak", "method": "..." },
          "uorfs": [ { "start": 10, "end": 70, "aa": 19, "kozak_class": "strong|adequate|weak",
                       "type": "uORF|overlapping|in-frame-extension" } ],
          "fold": { "method": "local folding: RNAplfold W=150 L=100 ...",   // see meta.methods.fold
                    "mfe": -23.4,            // unrestricted global MFE of the UTR (unchanged v1 field)
                    "mfe_global": -23.4,     // same value, explicit name
                    "mfe_local": -21.0,      // MFE of the span-restricted structure below (max bp span 100)
                    "dotbracket": "((((...))))...",   // LOCAL structure: MFE with max base-pair span 100 nt
                    "cap_dG": -5.1,          // MFE of first 50 nt
                    "layout2d": [[x,y], ...],          // naview coords of the local structure, centred, scaled to [-1,1]
                    "punp": [0.12, 0.98, ...],         // per-nt unpaired probability (RNAplfold W=150 L=100), 2 decimals, len = len
                    "acc_start": 0.51 },               // mean punp over the 30 nt around the main start (-15..+14, incl. CDS)
          "hairpins": [ { "start": 5, "end": 40, "stem_bp": 9, "dG": -12.3 } ],   // stem-loops from the LOCAL dotbracket, dG of that local structure
          "motifs": [ { "type": "uAUG|uORF_stop|nearcog|TOP|G4|TISU|kozak_main", "start": 24, "end": 27,   // 0-based half-open in the UTR
                        "seq": "AUG", "label": "uAUG (weak) -> uORF", "strength": "strong|adequate|weak|null",
                        "cite": "calvo2009", "note": "..." } ],   // pipeline/motifs.py find_motifs(); see meta.methods.motifs
          "ridge": { "mrl_pred": 0.0, "contrib": [ ... ] },   // data/ridge.json per_isoform entry (SPEC_v2); null/absent if none
          "gc": 0.62,
          "pred": { "optimus_mrl": 5.67,            // Optimus 5-Prime predicted mean ribosome load (MPRA scale), cite sample2019
                    "optimus_input": "full",        // "full" (len <= 100) | "last100" (100 nt upstream of the start codon)
                    "utrlm_te": null },             // null = not computed
          "measured": {                             // null if the transcript is in none of the Cao 2021 tables
            "metric": "TE (Ribo-seq/RNA-seq RPKM)",
            "values": { "HEK293T": 0.57, "PC3": 1.52, "muscle": null },   // null = absent / NaN in that table
            "cite": "cao2021",
            "ambiguous": false,                     // true: a sibling isoform has the identical TE (reads not resolved) -> gene-level
            "ambiguous_cells": [],                  // cell types where that happens
            "source_ids": { "HEK293T": "ENST...(their version)" },
            "note": "..." },
          "tss_note": {                             // null if exon structure unavailable
            "class": "reference|alt_first_exon|utr_splicing|length_only|alt_start_codon|other",
            "features": ["alt_first_exon", "alt_start_codon"],   // all criteria that apply
            "vs": "ENST...",                        // reference isoform (gene's first isoform: MANE/canonical)
            "tss": "12:68809002:+",                 // chrom:TSS(1-based GRCh38):strand
            "text": "alternative first exon / TSS vs MDM2-202" }
        }
      ]
    }
  ],
  "findings": [                                     // literature-reported pairwise isoform comparisons mapped to our ids
    { "id": "mdm2_p2", "title": "...", "cite": "brown1999", "gene": "MDM2",
      "better": "ENST... (isoform the paper says is translated better)", "worse": "ENST...",
      "evidence": "one sentence: how the paper's leader was matched to the Ensembl transcript" }
  ]
}
```

Rules: no invented values. Anything not computed or not sourced is `null` / omitted.

Field notes:
- `measured`: Cao et al. 2021 (Nat Commun 12:4138) transcript TE tables, matched by Ensembl transcript ID without
  version; their IDs are from an older Ensembl release, so a matching ID can still have a different model. Cells:
  HEK293T (Andreev 2015 data), PC3, muscle. `ambiguous` means identical TE for >= 2 isoforms of the gene in that table.
- `pred.optimus_mrl`: predicted for the UTR inside the MPRA reporter (eGFP), trained on 25-100 nt UTRs; isoforms sharing
  the last 100 nt get identical predictions. See `meta.methods.optimus` for the validation r.
- `motifs`: each type cites a verified source in `meta.sources` (uAUG/uORF_stop calvo2009, nearcog kearse2017 (only
  CUG/GUG/ACG/UUG in strong -3R/+4G context), TOP thoreen2012 TOP-like rule (>= 5 pyrimidines starting within 4 nt of
  the 5' end; Ensembl 5' ends are approximate), G4 huppert2005 (G3+N1-7 rule; kumari2007 = 5' UTR G4 represses
  translation), TISU elfakess2008 (AAS-AUG-GCGGC, UTR <= 30 nt), kozak_main kozak1986). `kozak_main` always exactly one;
  its `seq` and TISU `seq` extend into the CDS while `end` is clipped to `len`. PRTE and polyU are not reported (no
  verified definition). Numeric features for the ridge model: `motifs.features()` (labels in `motifs.FEATURES`).
- `fold` v2: structure is local (long-range pairs > 100 nt apart dropped), so `hairpins`/`layout2d` differ from v1.
- `tss_note` is a coordinate heuristic (see `meta.methods.tss_note`), not a curated annotation.
