"""Numpy re-implementation of Optimus 5-Prime (Sample et al. 2019, Nat Biotech 37:803),
variable-length 25-100 nt model (modeling/saved_models/Varying_length_25to100_model.hdf5,
https://github.com/pjsample/human_5utr_modeling, GPL-3.0).

Architecture / encoding taken from modeling/Varying_length_25to100_modeling.ipynb and the HDF5
model_config (Keras 2.2.4, TF backend):
  input  : 100 x 4 one-hot, order A,C,G,T; 'N' = all zeros; U -> T.
           UTR left-padded with N to 100 nt (notebook: utr100 = 75*'N' + utr, then [-100:]),
           i.e. the 3' end of the UTR (next to the start codon) is always at position 99.
  layers : Conv1D(120, k=8, 'same', relu) x3 -> Flatten -> Dense(40, relu) -> Dense(1, linear)
           (dropout layers are identity at inference).
           Keras/TF 'same' with even kernel 8 pads 3 before / 4 after.
  output : standard-scaled MRL. The notebook trains on StandardScaler(e_train.rl) and, in
           test_data(), inverse-transforms with a scaler fit on the *evaluated* set's observed rl.
           For new sequences (no observed MRL) we inverse-transform with the training-set scaler
           (mean/std of e_train.rl, recomputed from GEO GSM4084997 with the notebook's filters);
           see TRAIN_RL_MEAN / TRAIN_RL_STD and validate_optimus.py.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WEIGHTS = os.path.join(HERE, "external", "optimus", "Varying_length_25to100_model.hdf5")
SCALER_JSON = os.path.join(HERE, "external", "optimus", "train_scaler.json")

_W = None


def _load():
    global _W
    if _W is None:
        import h5py
        with h5py.File(WEIGHTS, "r") as f:
            mw = f["model_weights"]
            names = [n.decode() if isinstance(n, bytes) else n for n in mw.attrs["layer_names"]]
            layers = []
            for n in names:
                g = mw[n]
                wn = [w.decode() if isinstance(w, bytes) else w for w in g.attrs["weight_names"]]
                if wn:
                    d = {w.split("/")[-1].split(":")[0]: np.array(g[w], dtype=np.float64) for w in wn}
                    layers.append((n, d["kernel"], d["bias"]))
        assert [l[0].split("_")[0] for l in layers] == ["conv1d"] * 3 + ["dense"] * 2, layers
        _W = layers
    return _W


def one_hot(seqs, L=100):
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    X = np.zeros((len(seqs), L, 4))
    for i, s in enumerate(seqs):
        s = s.upper().replace("U", "T")
        s = ("N" * L + s)[-L:]
        for j, c in enumerate(s):
            k = idx.get(c)
            if k is not None:
                X[i, j, k] = 1.0
    return X


def _conv_same_relu(X, K, b):
    # X (n, L, cin), K (k, cin, cout); TF 'SAME': pad_before = (k-1)//2, pad_after = k-1-pad_before
    k = K.shape[0]
    pb = (k - 1) // 2
    pa = k - 1 - pb
    n, L, _ = X.shape
    Xp = np.pad(X, ((0, 0), (pb, pa), (0, 0)))
    out = np.zeros((n, L, K.shape[2]))
    for t in range(k):
        out += Xp[:, t:t + L, :] @ K[t]
    return np.maximum(out + b, 0.0)


def predict_scaled(seqs, batch=2048):
    W = _load()
    res = []
    for s in range(0, len(seqs), batch):
        X = one_hot(seqs[s:s + batch])
        for _, K, b in W[:3]:
            X = _conv_same_relu(X, K, b)
        h = X.reshape(X.shape[0], -1)
        h = np.maximum(h @ W[3][1] + W[3][2], 0.0)
        res.append((h @ W[4][1] + W[4][2]).reshape(-1))
    return np.concatenate(res) if res else np.zeros(0)


def train_scaler():
    with open(SCALER_JSON) as fh:
        d = json.load(fh)
    return d["mean"], d["std"]


def predict_mrl(seqs):
    m, s = train_scaler()
    return predict_scaled(seqs) * s + m
