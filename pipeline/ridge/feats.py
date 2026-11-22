"""Fold + feature extraction shared by all ridge datasets (uses pipeline/motifs.py, read-only).

Fold (ViennaRNA, default params, 37 C), on RNA = utr + cds_head[:30]:
  punp   : RNAplfold-style unpaired probability per nt (RNA.pfl_fold_up, u=1, W=150, L=100)
  mfe    : MFE of the 5' UTR alone with max base-pair span 100 nt (local; long-range pairs dropped)
  cap_dG : MFE of the first 50 nt of the UTR
Then motifs.features(utr, cds_head, fold) -> dict.
"""
import os
import sys
from multiprocessing import Pool

import RNA

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import motifs  # noqa: E402

# eGFP CDS head in the Sample 2019 reporter (Online Methods, PMC7100133: "atgggcgaattaagtaagggcgaggagctgttcacc")
EGFP_HEAD = "AUGGGCGAAUUAAGUAAGGGCGAGGAGCUGUUCACC"
FOLD_METHOD = ("ViennaRNA %s: punp = RNAplfold-style unpaired prob (pfl_fold_up u=1 W=150 L=100) on UTR+first 30 nt "
               "CDS; mfe = UTR MFE with max bp span 100; cap_dG = MFE of first 50 nt" % RNA.__version__)


def fold(utr, cds_head):
    utr = utr.upper().replace("T", "U")
    seq = utr + cds_head.upper().replace("T", "U")[:30]
    up = RNA.pfl_fold_up(seq, 1, 150, 100)
    punp = [float(up[i][1]) for i in range(1, len(seq) + 1)]
    md = RNA.md()
    md.max_bp_span = 100
    mfe = RNA.fold_compound(utr, md).mfe()[1] if utr else 0.0
    cap = RNA.fold(utr[:50])[1] if utr else 0.0
    return {"mfe": mfe, "cap_dG": cap, "punp": punp}


def one(args):
    utr, cds = args
    try:
        return motifs.features(utr, cds, fold(utr, cds))
    except Exception as e:  # noqa: BLE001
        return {"_error": str(e)}


def many(utrs, cdss, n_jobs=8, chunksize=64):
    with Pool(n_jobs) as p:
        return p.map(one, list(zip(utrs, cdss)), chunksize=chunksize)
