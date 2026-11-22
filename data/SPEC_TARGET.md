# Spec: oncogene focus + therapeutic levers

Framing: the model's parameters ARE the targets of real drugs. Turning a knob = applying a drug class.
The useful output is a THERAPEUTIC WINDOW: which transcripts lose the most output when a knob is turned,
relative to a reference set (housekeeping / bulk). That is hypothesis generation, never a claim of efficacy.

## Levers (each needs a real drug/tool compound + a primary source; drop any you cannot verify)
| lever | model parameter(s) | example agent | what it should do |
|---|---|---|---|
| eIF4A inhibition | `helicase` down | silvestrol / rocaglates / zotatifin (eFT226) | hurts mRNAs with structured 5' leaders most |
| ISR / eIF2 modulation | `tc` (and `loadTCexp`) | ISRIB, or stress inducers | changes uORF-containing mRNAs; can go either way |
| mTOR / cap-dependence | `loadRate` | rapalogs, mTOR kinase inhibitors | hits TOP-motif mRNAs most |
| start-codon fidelity | `kozakWeak`/`kozakAdeq` | (tool compounds; verify) | changes leaky scanning |
| splice-site choice | splice model: site strength / commitment | splice-switching ASOs (nusinersen as proof of concept), SF3B1 inhibitors, CLK/SRPK inhibitors | shifts isoform fractions |
| elongation rate | splice model: elongation | DRB / flavopiridol (tool compounds) | kinetic coupling -> isoform shift |

## Required honesty
- Every lever row shows its agent, mechanism and citation; "model parameter" is stated explicitly.
- Output is labelled "predicted differential sensitivity (model)", with the known failure modes repeated:
  no validation against isoform-level TE (chance within genes), no IRES/RBP/stalling mechanisms, kinetics assumed.
- Oncogene examples must pass the same evidence bar as tab 1 (data/evidence.json tiers; a primary paper that
  experimentally tested that gene's leader or splicing event).

## Output files
- `data/targets.json`: { levers:[{id,label,param,direction,agent,mechanism,cite,url}],
    genes:[{gene, role:"oncogene|tumour suppressor|reference", disease, cite, url, why_targetable}],
    reference_set:[gene ids used as the "bulk" comparator and why] }
- Scores are computed live in the app from the existing exact models; do not precompute effect sizes.
