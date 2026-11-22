# UTR Race: one physics for 5′ UTR translation

Built at Sundai Hack 150 ("Building 3D worlds", MIT, 2026-11-22).

Alternative 5′ UTR isoforms of one gene race side by side. Ribosomes (43S scanning complexes) load at the cap, scan through real Ensembl 5′ UTRs folded with ViennaRNA, get trapped by uORFs, stall at hairpins, and reinitiate after short uORFs. Under stress (low eIF2·TC) the race changes. ATF4 goes **up** while most mRNAs go down.

## More than a teaching toy

Most of what we know about 5′ UTR translation comes from separate single-gene reporter studies. Here every published finding must be reproduced by **one shared physical model with one parameter set**. Nothing is tuned per gene. This is similar to how AlphaFold's structures are held to physical constraints.

- **Scorecard:** published findings written as direction checks, computed with an exact expected-value version of the model (`app/src/analytic.js`, `app/src/findings.js`).
- **Feasibility map:** the region of parameter space where *all* findings hold at once. For example, ATF4 induction, its uORF1/uORF2 mutants and global repression jointly fix the TC re-acquisition per scanned nucleotide (λ = tcRate / scanSpeed). Scanning speed by itself cannot be identified from these data.
- **Failures are informative:** CHOP (DDIT3) induction is *not* reproduced at any setting. Its known mechanism, ribosome stalling inside the uORF, is not in the model.
- **Three layers per isoform:** the mechanism (this model), a statistical prediction (Optimus 5-Prime CNN, Sample et al. 2019, re-implemented in numpy and validated at r² = 0.84 on the paper's test set), and measured translation efficiency where available (Cao et al. 2021, HEK293T/PC3).

Honest limits: kinetic rates are inferred, not measured. Checks are on the direction of each effect, not its size. Reporter findings are not always endogenous behaviour. The MFE structure is a single prediction.

## Run

```bash
cd app && npm i && npm run dev        # Vite + three.js
```

Rebuild the data (Python ≥ 3.11, `pip install ViennaRNA requests numpy h5py`); see `pipeline/README.md`. Third-party inputs are re-downloaded into `pipeline/external/` and are not redistributed here.

## Data and credits

- Transcripts: Ensembl REST (GRCh38). Folding: ViennaRNA (Lorenz et al. 2011).
- Optimus 5-Prime: Sample et al., Nat Biotechnol 2019 (https://github.com/pjsample/human_5utr_modeling, GPL-3.0). Weights are not included. Only derived predictions ship in `utr_data.json`.
- Measured TE: Cao et al., Nat Commun 2021, via https://github.com/zzz2010/5UTR_Optimizer (MIT).
- Findings: Vattem & Wek 2004 (ATF4); MDM2 P1/P2 leaders; full source list in `utr_data.json` → `meta.sources` and in the app.
- Starter kit and prompts: the Sundai 3D-worlds workshop (https://github.com/hunterh37/sundai-3dworlds).

License: MIT (code). Data files keep their sources' terms.
