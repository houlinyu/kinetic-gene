"""Run the app's deterministic scanning model (app/src/analytic.js expected(iso, params)) on many UTRs via node.

Isoform objects are built with the same rules as build_data.py (find_uorfs, find_hairpins, kozak_class, imported
read-only); hairpins come from a local MFE fold (max bp span 100, same as feats.fold), not the global MFE.
"""
import json
import os
import subprocess
import sys
from multiprocessing import Pool

import RNA

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import build_data as bd  # noqa: E402

NODE = os.path.expanduser("~/miniconda3/envs/web/bin/node")


def iso(args):
    utr, cds = args
    utr, cds = utr.upper().replace("T", "U"), cds.upper().replace("T", "U")
    md = RNA.md()
    md.max_bp_span = 100
    db, _ = RNA.fold_compound(utr, md).mfe()
    return {"len": len(utr), "tags": [], "main_kozak": {"class": bd.kozak_class(utr + cds, len(utr))},
            "uorfs": bd.find_uorfs(utr, cds), "hairpins": bd.find_hairpins(db, utr)}


def build(utrs, cdss, n_jobs=8):
    with Pool(n_jobs) as p:
        return p.map(iso, list(zip(utrs, cdss)), chunksize=32)


def run(isos, params=None):
    inp = "\n".join(json.dumps(i) for i in isos) + "\n"
    cmd = [NODE, os.path.join(HERE, "physics_eval.mjs")] + ([json.dumps(params)] if params else [])
    out = subprocess.run(cmd, input=inp, capture_output=True, text=True, check=True).stdout
    return [tuple(float(x) for x in l.split("\t")) for l in out.strip().split("\n")]
