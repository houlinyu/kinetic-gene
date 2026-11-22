# Tab 2 spec: "How the isoform is made" (transcription + co-transcriptional splicing)

One app, two tabs: **Make it** (this spec) -> **Use it** (existing translation tab).
The join: tab 2 predicts the *relative abundance* of each isoform; tab 1 predicts *protein per mRNA* for each isoform.
Protein output per gene = sum over isoforms of (abundance x protein-per-mRNA). Both tabs share the gene picker.

## Model rules (one sentence each; same style as app/src/sim.js)
- RNA Pol II initiates at a promoter, pauses near the TSS, then elongates at a settable rate.
- A splice site can only be used once Pol II has transcribed it (co-transcriptional).
- Each 5' and 3' splice site has a strength score; commitment to a site happens at a rate that scales with its strength.
- Competing sites race: a weak upstream site wins only if it commits before a stronger downstream site appears.
- Intron removal takes time; slow Pol II gives weak/alternative sites more time (kinetic coupling).
Every parameter is sourced or labelled "assumed", same rule as app/src/params.js.

## Ownership
- SPLICE-DATA agent: `pipeline/splice/` + `data/splice_data.json` + `app/public/splice_data.json`.
- SPLICE-MODEL agent: `app/src/splice/model.js` (+ `model_exact.js` if an expected-value version is feasible) and `app/src/splice/params.js`.
- SPLICE-UI agent: `app/src/splice/view.js` and tab wiring in `app/src/main.js` / `index.html` / `style.css`.
- EVIDENCE: extend `data/evidence.json` with a `splice_genes` / `splice_findings` section.

## `splice_data.json` shape (0-based, genomic-order along the transcript)
```jsonc
{ "meta": { "release": "...", "sources": [{key, citation, url}], "methods": {"ss_score": "...", "psi": "..."} },
  "genes": [ { "gene": "...", "strand": "+",
    "tss": [ { "id": "...", "pos": 0, "support": "...", "cite": "..." } ],
    "exons": [ { "id": "E1", "start": 0, "end": 120, "type": "constitutive|cassette|alt5|alt3|first|last" } ],
    "introns": [ { "id": "I1", "start": 120, "end": 1500,
        "donor": { "pos": 120, "seq": "...", "score": 0.0, "scale": "MaxEntScan 5'" },
        "acceptor": { "pos": 1500, "seq": "...", "score": 0.0, "scale": "MaxEntScan 3'" },
        "bp": { "pos": null, "score": null } } ],
    "isoforms": [ { "id": "ENST...", "name": "...", "tss_id": "...", "exon_ids": ["E1","E3"], "utr5_id": "ENST... (link to tab 1)" } ],
    "measured": { "psi": [ { "event": "E2 inclusion", "value": 0.0, "tissue": "...", "n": 0, "cite": "..." } ],
                  "isoform_fraction": [ { "id": "ENST...", "value": 0.0, "tissue": "...", "cite": "..." } ] } } ] }
```
Rules: no invented values; anything unmeasured is null. Splice-site scores must come from a real, cited scoring method.
