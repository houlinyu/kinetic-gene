#!/usr/bin/env python
"""Fit 4 physics-model parameters (app/src/analytic.js) to MRL (Sample 2019 rows in work/mrl.csv.gz).
Score = log(perMin) from expected(iso, params); objective = r2 of a linear map score -> MRL on the training folds.
Nelder-Mead from the app defaults, params clipped to the slider ranges in app/src/params.js. 5-fold CV (same
KFold(5, shuffle, seed 0) as fit.py) -> held-out r2; then a final fit on all rows. Writes work/physics_fit.json."""
import gzip
import json
import os
import subprocess

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.model_selection import KFold

HERE = os.path.dirname(os.path.abspath(__file__))
NODE = os.path.expanduser("~/miniconda3/envs/web/bin/node")
NAMES = ["kozakWeak", "kozakAdeq", "reinitMax", "hpKT"]
DEF = np.array([0.25, 0.6, 0.8, 4.0])            # app/src/params.js values at time of fit
LO, HI = np.array([0, 0, 0, 1.0]), np.array([1, 1, 1, 10.0])
INP = gzip.open(os.path.join(HERE, "work", "mrl_isos.jsonl.gz"), "rt").read()
y = pd.read_csv(os.path.join(HERE, "work", "mrl.csv.gz"), usecols=["rl"])["rl"].to_numpy(float)


def score(x):
    p = dict(zip(NAMES, map(float, np.clip(x, LO, HI))))
    out = subprocess.run([NODE, os.path.join(HERE, "physics_eval.mjs"), json.dumps(p)], input=INP,
                         capture_output=True, text=True, check=True).stdout
    return np.log(np.array([float(l.split("\t")[2]) for l in out.strip().split("\n")]) + 1e-6)


def r2lin(s, yy, s_tr=None, y_tr=None):
    s_tr = s if s_tr is None else s_tr
    y_tr = yy if y_tr is None else y_tr
    b, a = np.polyfit(s_tr, y_tr, 1)
    p = a + b * s
    return 1 - np.sum((yy - p) ** 2) / np.sum((yy - yy.mean()) ** 2)


def fit(idx):
    f = lambda x: -r2lin(score(x)[idx], y[idx])  # noqa: E731
    r = minimize(f, DEF, method="Nelder-Mead", options={"maxfev": 80, "xatol": 1e-3, "fatol": 1e-4})
    return np.clip(r.x, LO, HI)


oof = np.zeros(len(y))
s0 = score(DEF)
for tr, te in KFold(5, shuffle=True, random_state=0).split(y):
    x = fit(tr)
    s = score(x)
    b, a = np.polyfit(s[tr], y[tr], 1)
    oof[te] = a + b * s[te]
    print("fold", dict(zip(NAMES, np.round(x, 3))), flush=True)
r2cv = float(1 - np.sum((y - oof) ** 2) / np.sum((y - y.mean()) ** 2))
xf = fit(np.arange(len(y)))
res = {"dataset": "Sample 2019 MRL rows used for models.mrl (same rows, same 5 folds)",
       "params_fitted": NAMES, "params": dict(zip(NAMES, map(float, xf))), "defaults": dict(zip(NAMES, DEF.tolist())),
       "r2_cv": r2cv, "r2_cv_default_params": None,
       "r2_in_sample_default_params": float(r2lin(s0, y)), "r2_in_sample_fitted": float(r2lin(score(xf), y)),
       "note": "score = log(initiations/min) from analytic.js expected(); linear map to MRL fit on training folds; "
               "Nelder-Mead (<=80 evals) from app defaults, params clipped to the app slider ranges; hairpins from "
               "a local MFE fold (max bp span 100); all other params at app defaults."}
oof_def = np.zeros(len(y))
for tr, te in KFold(5, shuffle=True, random_state=0).split(y):
    b, a = np.polyfit(s0[tr], y[tr], 1)
    oof_def[te] = a + b * s0[te]
res["r2_cv_default_params"] = float(1 - np.sum((y - oof_def) ** 2) / np.sum((y - y.mean()) ** 2))
pd.DataFrame({"physics_fitted": oof}).to_csv(os.path.join(HERE, "work", "mrl_physics_fit_oof.csv.gz"), index=False)
json.dump(res, open(os.path.join(HERE, "work", "physics_fit.json"), "w"), indent=1)
print(res)
