# pipeline/targets — oncogene focus + therapeutic levers

Owner: ONCOTARGET-DATA agent. Output: `data/targets.json`, copied to `app/public/targets.json`.
Contract: `data/SPEC_TARGET.md`. Evidence bar: `data/evidence.json` tiers (A/B, same definitions).
Regenerate with `python3 pipeline/targets/build_targets.py` (no network, no deps — the file is hand
curated and the script only serialises it, so the JSON is never silently rewritten by a scraper).

## How this was curated (2026-11-22)

1. Every candidate in `SPEC_TARGET.md` (and the lead's candidate list) was turned into a PubMed query
   and resolved with NCBI E-utilities: `esearch` (relevance sort) then `efetch rettype=abstract
   retmode=text`. The curator read the fetched text; nothing was taken from memory.
2. A row ships only if the paper the curator read **experimentally tested that gene's own 5' leader or
   its own splicing event**. Transcriptome-wide-only evidence ships as tier B and says so.
3. Each row carries a verbatim quote of **≤ 15 words** from the text actually fetched, `quote_source`
   (abstract unless stated), the PubMed URL, what was measured, and in which cells/species.
4. All quotes were re-checked mechanically against the fetched text after the file was written
   (unicode-normalised substring match + word count). 50/50 quotes match; 0 exceed 15 words.
5. Anything that could not be verified was dropped and recorded in `targets.json → rejected` with the
   reason. Nothing is marked "assumed" and shipped anyway.

`param` fields are **our** mapping of a drug class onto a model knob (`app/src/params.js`,
`app/src/splice/params.js`). No cited paper mentions our parameters; each lever says so in its caveat.
No effect sizes are precomputed — per SPEC_TARGET.md the app scores live from the exact models.

## What shipped

- **14 cancer genes** (12 oncogene, 1 tumour suppressor + MDM4): MDM2, ERBB2, PIM1, ODC1, MYC, KRAS,
  MET, BCL2L1, MCL1, AR, BRAF, CCND1, CDKN2A, MDM4. Nine are isoform-resolved (MDM2, MET, BCL2L1,
  MCL1, AR, BRAF, CCND1, CDKN2A, MDM4); ERBB2, PIM1, ODC1, MYC, KRAS are leader-level only.
- **7 levers**: eIF4A inhibition (`helicase`↓, `hpFree`↓), ISR/eIF2 (`tc`↓ under stress, ISRIB back up),
  mTOR (`loadRate`↓), splice-switching ASO (`site_score`↓ at the masked site), SF3B1 modulation
  (`commitScale`↓, `removalTime`↑), CLK/SRPK (`scoreSlope`↓), Pol II elongation (`elongRate`↓,
  `pauseEscape`↓). Each lever carries machine-readable `param` + `also_params` (real keys from
  `app/src/params.js` / `app/src/splice/params.js`), a single `direction` (`up`/`down`) so the UI can
  derive the extreme, `model` (`translation`/`splice`), and the prose nuance in `param_note` /
  `direction_note` — never inside `param` itself, which `app/src/target/panel.js` tokenises.
- **3 reference genes** (GAPDH, ACTB, HBB) justified by data, not assertion — see below.
- **`differential_evidence`**: per lever, what preferential sensitivity was actually shown, in which
  system, the conflicts between papers, and an explicit `what_the_app_may_claim`.

## Reference set: why these three, and why not the usual suspects

The comparator must be the *insensitive* class, demonstrated, not assumed.

- **GAPDH and ACTB** come from Reiter 2008 (PMID 18638456): in rat liver, AMPK activation repressed the
  TOP mRNAs rpS6, rpS8 and eEF1α while "neither global protein synthesis nor translation of mRNAs
  encoding GAPDH or beta-actin was changed". That is a measured dissociation on the same axis our
  `loadRate` lever represents. It also settles a common assumption: β-actin did **not** behave as a TOP
  mRNA here, so ACTB is usable as a reference rather than excluded.
- **HBB** is the classic efficiently translated standard: it is the control message the ODC 5'UTR
  constructs were compared against (Van Steeg 1991). This is comparator *usage*, not a measurement of
  the HBB leader — the weakest of the three justifications, and flagged as such in the JSON.
- **Ribosomal-protein and eEF1A genes were deliberately excluded**: they are the TOP class and are the
  *most* mTOR-sensitive transcripts (Thoreen 2012, Reiter 2008). They belong in the sensitive arm as
  positive controls (recorded under `differential_evidence.mtor`), never in the bulk reference set.
- Reference leaders still need the pipeline to confirm their computed features (length, uAUG, TOP, G4)
  once the genes are built — the literature here supports *insensitivity*, not "unstructured".

## What was rejected and why (full list in `targets.json → rejected`)

| rejected | reason |
|---|---|
| EGFRvIII | genomic deletion of exons 2–7, not a regulated splicing choice; neither model represents it, and no primary paper was read in time |
| KRAS4A/4B splicing | no readable primary paper found this session; KRAS ships on its structured-leader/eIF4A evidence only |
| start-codon-fidelity lever (`kozakStrong/Adeq/Weak`) | no agent with a verified published effect on start-codon selection; SPEC_TARGET.md already marked it "verify". Dropped, not invented |
| BCL2, CTNNB1 | mechanism unverified in this session |
| VEGFA, TP53 | verified mechanism is IRES/cap-independent (`evidence.json`), which the scanning model does not implement |
| MYC P0/P1/P2 leader claims | published leaders are not mapped to Ensembl transcripts and MANE MYC starts at CUG (`evidence.json`); MYC ships gene-level, tier B |
| "a drug shifts the MCL1/BCL2L1 ratio" | the primary papers establish the isoforms by overexpression; no drug shifted the endogenous ratio |
| a phase number for zotatifin | the papers read say "development candidate", "first-in-class" and report patient biopsies; no phase stated, so none claimed |
| DDIT3 as a second ISR source | the PMID resolved to an unrelated paper; ISR differential ships on Vattem 2004 alone |
| GAPDH leader-feature citation | the search returned an unrelated trypanosome paper; GAPDH ships on the insensitivity data only |

## Honesty rules this file must keep

- The app's output is **"predicted differential sensitivity (model)"** — never efficacy, never a patient claim.
- Repeat the failure modes: no validation against isoform-level TE, no IRES/RBP/stalling mechanisms,
  kinetics assumed, and (new) **the mTOR lever must not be ranked by leader structure** — Thoreen 2012
  found no evidence for that, only for TOP motifs.
- The eIF4A lever is a structure-unwinding caricature: Iwasaki 2016 showed rocaglates are
  gain-of-function clamps on polypurine tracts and that 5'UTR structure is "only a minor determinant".
  Any therapeutic-window display driven by our structure score must say this.

## Genes not yet in `pipeline/genes.txt` (the lead must add them; this agent did not edit that file)

Leader model: **PIM1, ODC1, KRAS, GAPDH, ACTB, HBB**.
Splice model (`pipeline/splice/`): **MET, BCL2L1, MCL1, AR, BRAF, CCND1, MDM4, CDKN2A**.
Already present: MDM2, ERBB2, MYC. Same list is machine-readable in `targets.json → pipeline_requests`.
