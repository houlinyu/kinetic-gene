# v2 spec: local structure, motifs, ridge fit (contract between agents)

Owners (do not edit files you don't own):
- STRUCT agent: `pipeline/motifs.py` (new, write FIRST, within ~15 min), `pipeline/build_data.py`, `data/SCHEMA.md`.
- RIDGE agent: `pipeline/ridge/` (new dir), `data/ridge.json`, `data/figures/`. Imports `pipeline/motifs.py` (read-only).
- UI agent: `app/src/main.js`, `world.js`, `style.css`, `index.html`, `scorecard.js`, new UI files.
- Lead (me): `app/src/sim.js`, `analytic.js`, `findings.js`, `params.js`, this spec.

## `pipeline/motifs.py` API (pure functions, RNA or DNA input, positions 0-based within the 5' UTR)
```python
find_motifs(utr: str, cds_head: str) -> list[dict]   # cds_head = first >=30 nt of CDS incl. start codon
features(utr: str, cds_head: str, fold: dict | None = None) -> dict[str, float]  # flat numeric features for ridge
```
Motif dict: `{type, start, end, seq, label, strength (str|None), cite (key), note}`. Types (each needs a cite key with a real source in `meta.sources`; drop a type if no source verifies it):
- `uAUG` (upstream AUG; strength = Kozak class strong/adequate/weak; -3 R / +4 G rule, same as main_kozak)
- `uORF_stop` (stop codon ending a uORF)
- `nearcog` (upstream CUG/GUG/ACG/AUU/AUA... restrict to CUG, GUG, ACG, UUG; strength = context class; only report strong context)
- `TOP` (5' terminal oligopyrimidine: C at position 0 followed by >= 4 pyrimidines; note: depends on exact TSS, Ensembl 5' ends are approximate)
- `TISU` (translation initiator of short UTRs: SAASAUGGCGGC around the main AUG, only for UTR <= 30 nt; report on the main start)
- `G4` (putative G-quadruplex: (G{3,}N{1,7}){3,}G{3,})
- `PRTE` (pyrimidine-rich translational element, if a sourced definition exists)
- `polyU` (U{6,})
- `kozak_main` (always one entry: the -6..+4 context of the main start with its class)

## Fold fields (per isoform), STRUCT agent
- `fold.method`: e.g. "RNAplfold W=150 L=100 + MEA" (local folding; long-range pairs dropped).
- `fold.dotbracket`, `fold.layout2d` (one [x,y] per nt, normalized to [-1,1], from the local structure), `fold.mfe_global` (old global MFE kept for reference), `fold.cap_dG`.
- `fold.punp`: per-nt unpaired probability (array of floats, 2 decimals), `fold.acc_start`: mean unpaired prob of 30 nt around the main start.
- `hairpins` re-derived from the local structure (same shape as before).

## Ridge outputs, RIDGE agent
`data/ridge.json`:
```jsonc
{ "models": {
    "mrl": { "dataset": "Sample 2019 GSE114002 (human library + 25-100 nt random)", "n": 0, "cv": "5-fold", "r2_cv": 0.0, "r2_optimus_same_split": 0.0,
             "features": [ { "name": "n_uAUG_strong", "label": "upstream AUG, strong context", "coef": 0.0, "ci": [0, 0] } ], "intercept": 0.0, "scaler": { "mean": {}, "std": {} } },
    "te":  { "dataset": "Cao 2021 HEK293T TE (endogenous)", "cv": "5-fold grouped by gene", ... } },
  "physics_fit": { "dataset": "...", "params": {}, "r2_cv": 0.0, "note": "" },   // optional, if time
  "per_isoform": { "ENST...": { "mrl_pred": 0.0, "te_pred": 0.0, "contrib": [ { "name": "...", "label": "...", "value": 0.0, "contribution": 0.0 } ] } } }
```
Figures (PNG + the data behind them) in `data/figures/`: coefficient forest plot per model; predicted vs observed scatter (held-out) for ridge vs Optimus vs physics.
STRUCT's `build_data.py` merges `per_isoform` into `isoform.ridge` if `data/ridge.json` exists, and copies `models` (minus per_isoform) into `meta.ridge`.
