# utr_data.json schema (contract between pipeline/ and app/)

Written by `pipeline/build_data.py`, read by `app/` from `app/public/utr_data.json`.
All sequences are RNA (A,C,G,U), 5'->3'. All positions are 0-based, half-open, relative to the 5' UTR start (cap = 0).

```jsonc
{
  "meta": { "generated": "ISO date", "gencode_or_ensembl_release": "...", "vienna_version": "...",
            "sources": [ { "key": "sample2019", "citation": "...", "url": "..." } ],
            "methods": { "optimus": "...", "measured": "...", "tss_note": "..." } },   // free-text method notes
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
          "fold": { "mfe": -23.4, "dotbracket": "((((...))))...", "cap_dG": -5.1,   // cap_dG = MFE of first 50 nt
                    "layout2d": [[x,y], ...] },                                    // one point per nt (RNAplot naview coords)
          "hairpins": [ { "start": 5, "end": 40, "stem_bp": 9, "dG": -12.3 } ],   // stem-loops from dotbracket, dG of that local structure
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
- `tss_note` is a coordinate heuristic (see `meta.methods.tss_note`), not a curated annotation.
