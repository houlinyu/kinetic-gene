#!/bin/bash
# Discovery chain (nice 19, 8 procs). Inputs: pipeline/ridge/downloads/gencode.v47.pc_transcripts.fa.gz,
# pipeline/discovery/dl/gencode.v47.basic.annotation.gtf.gz (curl -O from ftp.ebi.ac.uk/.../release_47/).
set -e
cd "$(dirname "$0")/../.."
N=~/miniconda3/envs/web/bin/node; PY=~/miniconda3/envs/utrrace/bin/python
nice -n 19 $PY pipeline/discovery/prep.py                    # ~7 min: groups + 2 fold methods
nice -n 19 $N pipeline/discovery/feasible.mjs                # feasible parameter grid (findings.score)
nice -n 19 $N --max-old-space-size=8000 pipeline/discovery/evaluate.mjs   # analytic model, all points
nice -n 19 $PY pipeline/discovery/calls.py > pipeline/discovery/work/calls.log   # data/discovery.json
nice -n 19 ~/miniconda3/envs/fungi/bin/python -I pipeline/discovery/figs.py      # data/figures/discovery_*
