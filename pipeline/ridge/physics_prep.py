#!/usr/bin/env python
"""Physics-model score per row: app/src/analytic.js expected(iso, defaults) -> pMain, flux, log(perMin).
usage: physics_prep.py mrl|te   (reads work/<ds>.csv.gz, writes work/<ds>_physics.csv.gz)"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import feats  # noqa: E402
import physics  # noqa: E402

ds = sys.argv[1]
nj = int(os.environ.get("NJ", "4"))
d = pd.read_csv(os.path.join(HERE, "work", f"{ds}.csv.gz"))
if ds == "mrl":
    key, cds = "utr", [feats.EGFP_HEAD] * len(d)
else:
    # cds_head column = first 60 nt of the CDS; enough for physics uORF typing only if the stop is inside it,
    # so re-read full CDS from GENCODE v19
    import prep
    g = prep.read_gencode()
    cds = []
    for vid in d["vid"]:
        hdr, s = g[vid]
        rng = {x.split(":")[0]: x.split(":")[1] for x in hdr if ":" in x}
        cds.append(s[int(rng["CDS"].split("-")[0]) - 1:])
    key = "vid"
isos = physics.build(list(d["utr"]), cds, nj)
r = np.array(physics.run(isos))
pd.DataFrame({key: d[key], "phys_pMain": r[:, 0], "phys_flux": r[:, 1],
              "phys_logperMin": np.log(r[:, 2] + 1e-6)}).to_csv(os.path.join(HERE, "work", f"{ds}_physics.csv.gz"), index=False)
import gzip, json  # noqa: E401
with gzip.open(os.path.join(HERE, "work", f"{ds}_isos.jsonl.gz"), "wt") as fh:
    fh.write("\n".join(json.dumps(i) for i in isos) + "\n")
print("done", ds, len(d))
