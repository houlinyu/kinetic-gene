# Kinetic Gene — one physics from promoter to protein

Built at Sundai Hack 150 ("Building 3D worlds", MIT, 2026-11-22). Live: https://kinetic-gene.vercel.app

Both halves of gene expression are the same physics: a molecular machine moves 5′→3′ along a template, and decisions get made in real time as it passes. A splice site and a start codon are both races against that machine.

- **① Transcription** — RNA Pol II elongates, pausing and speed set the timing, and the spliceosome commits to sites only once they have been transcribed. Output: which isoforms exist, and in what proportion.
- **② Translation** — 43S complexes load at the cap and scan each real 5′ UTR (folded with ViennaRNA), stalling at hairpins, initiating at uORFs, reinitiating or being lost. Output: protein per mRNA per isoform.
- **Therapeutic levers** on both pages, split by where the drug acts: splice-switching ASOs, SF3B1 and CLK/SRPK inhibitors and elongation compounds on transcription; eIF4A inhibitors and ISR modulation on translation.

Every example carries a primary citation, and 3D views are built from real cryo-EM structures at true scale.

## What is validated, and what is not

This is a mechanistic model, stated honestly in the app itself:

- **Reproduced** (one parameter set, no per-gene tuning): removing ATF4's CDS-overlapping uORF2 raises basal output ~33× (Vattem & Wek measured ≈35×); mutating the ERBB2 uORF de-represses 4.4× (Child 1999: ~5×); ATF4 rises ×3.5 under simulated eIF2α-P while the median mRNA falls ×0.42; slowing Pol II shifts cassette-exon inclusion in the direction de la Mata 2003 reported.
- **5 of 6** literature direction checks pass at our chosen parameter point; **CHOP fails** and the failure is shown, because the model lacks the mechanism Palam 2011 invokes.
- **Not validated**: predicted splicing does not match measured PSI on alternative exons (17 paired GTEx events: Pearson r = 0.16, barely better than predicting the mean; MDM2 exon 10 is sign-inverted). The translation model ranks isoforms at chance against endogenous TE.
- **Robustness is narrow**: sampling all 16 exposed parameters over their declared ranges, the ATF4 induction holds in 3.5%, and all five findings simultaneously in 0.70%.
- **Dynamic range is unphysical** (an unbounded hairpin stall), so output ratios are displayed as rank order, capped at ×100.
- The therapeutic "window" is largely leader length: against length-matched controls it falls from ~360× to ~5×.

A by-product worth noting: in the Cao 2021 isoform-level TE table, the isoform with more Ribo-seq reads wins 80% of within-gene comparisons although read counts cancel in the ratio — so that benchmark largely measures quantification, not biology.

## Run

```bash
cd app && npm i && npm run dev        # Vite + three.js
```

Rebuild data (Python ≥3.11): see `pipeline/README.md`. Third-party inputs are re-downloaded, not redistributed.

## Data and credits

- Transcripts/exons: Ensembl REST, GENCODE (GRCh38). Folding: ViennaRNA (Lorenz et al. 2011). Splice sites: MaxEntScan (Yeo & Burge 2004), validated against the published reference scores.
- Measured: GTEx PSI via ASCOT (Ling et al. 2020); translation efficiency from Cao et al. 2021 via https://github.com/zzz2010/5UTR_Optimizer (MIT); Optimus 5-Prime (Sample et al. 2019, GPL-3.0) re-implemented and validated at r² = 0.84 — weights not redistributed.
- Structures: 48S scanning complex PDB 8PJ1 (Petrychenko et al. 2025); 48S+eIF4F 8OZ0 (Brito Querido et al. 2024); Pol II elongation complex 6GMH (Vos et al. 2018, mammalian); human spliceosome C* 5XJC (Zhang et al. 2017); U1 snRNP 4PJO (Kondo et al. 2015).
- Full source list in the app and in `data/*.json` → `meta.sources`.
- Starter kit: the Sundai 3D-worlds workshop (https://github.com/hunterh37/sundai-3dworlds).

License: MIT (code). Data files keep their sources' terms.
