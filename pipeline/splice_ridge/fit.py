"""Interpretable ridge on isoform usage, with mechanistic + baseline comparisons.

Target  : a transcript's measured share of its gene (Isopedia long-read panel).
CV      : GroupKFold(5) grouped by gene (locus_id) -- no gene is ever in train and test.
Compared: ridge | uniform | longest-transcript-wins | dominant-TSS-only | trivial confounds
          and, on the 12 app genes only, the app's mechanistic co-transcriptional model.
"""
import os, sys, json, math, csv, collections
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold
from scipy import stats as sps

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
W = os.path.join(HERE, os.environ.get("SR_WORK", "work"))
os.makedirs(W, exist_ok=True)
FIG = os.path.join(ROOT, "data", "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(20261122)
MIN_SAMPLES = 20   # locus must be seen in >=20 Isopedia samples for its usage to be used

FEATURES = {
    "mes_donor_min":   "Weakest donor used (MaxEntScan 5' bits)",
    "mes_donor_mean":  "Mean donor strength of the junctions used (MaxEntScan 5')",
    "mes_acc_min":     "Weakest acceptor used (MaxEntScan 3' bits)",
    "mes_acc_mean":    "Mean acceptor strength of the junctions used (MaxEntScan 3')",
    "ss_p_min":        "Weakest splice site used (Pangolin site probability)",
    "ss_p_mean":       "Mean splice-site strength of the junctions used (Pangolin)",
    "ss_use_min":      "Weakest predicted site usage of the junctions used (Pangolin)",
    "ss_use_mean":     "Mean predicted site usage of the junctions used (Pangolin)",
    "n_exons":         "Number of exons",
    "log_tx_len":      "Transcript length (log10 nt)",
    "log_mean_exon":   "Mean exon length (log10 nt)",
    "log_mean_intron": "Mean intron length (log10 nt)",
    "log_min_intron":  "Shortest intron (log10 nt)",
    "log_max_intron":  "Longest intron (log10 nt)",
    "utr5_len_log":    "5' UTR length (log10 nt)",
    "log_dist_dom_tss":"Distance from the gene's dominant TSS (log10 nt)",
    "is_dom_tss":      "Starts at the gene's dominant TSS",
    "cassette_skip":   "Skips a cassette exon of the reference transcript",
    "alt_5ss":         "Uses an alternative 5' splice site",
    "alt_3ss":         "Uses an alternative 3' splice site",
    "alt_first_exon":  "Alternative first exon",
    "is_longest":      "Is the gene's longest transcript (structural reference)",
    "gc":              "GC content",
    "orf_len_log":     "ORF length (log10 nt)",
    "n_uatg":          "Upstream AUGs in the 5' UTR",
    "nmd_ptc":         "Premature stop codon (NMD candidate)",
    "intron_retention":"Retains an intron",
    "likely_5p_fragment":"Flagged as a likely 5' fragment",
    "tss_cpg_island":  "TSS in a CpG island",
    "tss_tata":        "TATA box at the TSS",
    "tss_gc":          "GC content around the TSS",
    "tes_pas_canonical":"Canonical polyA signal at the 3' end",
    "frac_share_tss":  "Fraction of the gene's transcripts sharing this TSS",
    "frac_share_tes":  "Fraction of the gene's transcripts sharing this 3' end",
    "n_iso_locus_log": "Number of annotated isoforms of the gene (log10)",
}

def read_tsv(p):
    with open(p) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return list(r)

def fnum(x, d=float("nan")):
    try:
        v = float(x)
        return v if math.isfinite(v) else d
    except Exception:
        return d

# ---------------------------------------------------------------- load + merge
SRC = os.path.join(HERE, "work")
import time; T0=time.time()
def tick(m): print(f"[{time.time()-T0:7.1f}s] {m}", file=sys.stderr, flush=True)
feat = read_tsv(os.path.join(SRC, "features.tsv"))
usage = {r["isoform_uid"]: r for r in read_tsv(os.path.join(SRC, "usage.tsv"))}
MES = {}
_pm = os.path.join(SRC, "maxent.tsv")
if os.path.exists(_pm):
    MES = {r["isoform_uid"]: r for r in read_tsv(_pm)}
tick("loaded tables")
rows = []
for r in feat:
    u = usage.get(r["isoform_uid"])
    if not u: continue
    y = fnum(u["usage"])
    if not math.isfinite(y): continue
    if int(u["n_locus_samples"]) < MIN_SAMPLES: continue
    if not int(r["has_pangolin"]): continue
    d = dict(r)
    d["y"] = y
    d["n_pos"] = int(u["n_pos"])
    d["gene_name"] = u["gene_name"]
    d["utr5_len_log"] = math.log10(1 + fnum(r["utr5_len"], 0))
    d["orf_len_log"] = math.log10(1 + fnum(r["orf_len_nt"], 0))
    d["n_iso_locus_log"] = math.log10(float(r["n_iso_locus"]))
    m = MES.get(r["isoform_uid"])
    if MES and not m: continue       # keep one consistent feature matrix
    for k in ("mes_donor_min", "mes_donor_mean", "mes_acc_min", "mes_acc_mean"):
        d[k] = fnum(m[k]) if m else float("nan")
    rows.append(d)
# only loci where >=2 isoforms survived the filters (a usage share needs competitors)
cnt = collections.Counter(r["locus_id"] for r in rows)
rows = [r for r in rows if cnt[r["locus_id"]] >= 2]
print(f"rows={len(rows)} genes={len({r['locus_id'] for r in rows})}", file=sys.stderr)

names = list(FEATURES)
X = np.array([[fnum(r[n], 0.0) for n in names] for r in rows], float)
X = np.where(np.isfinite(X), X, 0.0)
y = np.array([r["y"] for r in rows], float)
groups = np.array([r["locus_id"] for r in rows])
gidx = collections.defaultdict(list)
for i, g in enumerate(groups): gidx[g].append(i)

ALPHAS = np.logspace(-3, 4, 15)

def zfit(Xtr, ytr):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1.0
    m = RidgeCV(alphas=ALPHAS).fit((Xtr - mu) / sd, ytr)
    return mu, sd, m

# ---------------------------------------------------------------- grouped CV
NF = 5
gkf = GroupKFold(n_splits=NF)
oof = np.full(len(y), np.nan)
alphas_used = []
tick("start CV")
for tr, te in gkf.split(X, y, groups):
    mu, sd, m = zfit(X[tr], y[tr]); tick("fold")
    oof[te] = m.predict((X[te] - mu) / sd)
    alphas_used.append(float(m.alpha_))

def metrics(pred, truth):
    ok = np.isfinite(pred) & np.isfinite(truth)
    p, t = pred[ok], truth[ok]
    if len(p) < 3 or np.std(p) == 0:
        return dict(r2=float("nan"), pearson=float("nan"), spearman=float("nan"), n=int(len(p)))
    ss = 1 - np.sum((t - p) ** 2) / np.sum((t - t.mean()) ** 2)
    return dict(r2=float(ss), pearson=float(sps.pearsonr(p, t)[0]),
                spearman=float(sps.spearmanr(p, t)[0]), n=int(len(p)))

# ---------------------------------------------------------------- baselines
niso = np.array([len(gidx[g]) for g in groups], float)
pred_uniform = 1.0 / niso
pred_longest = np.array([float(r["is_longest"]) for r in rows])
# dominant-TSS-only: split the gene's mass equally over the transcripts at its dominant TSS
dom_n = collections.Counter()
for r in rows:
    if int(r["is_dom_tss"]): dom_n[r["locus_id"]] += 1
pred_domtss = np.array([ (1.0 / dom_n[r["locus_id"]]) if int(r["is_dom_tss"]) and dom_n[r["locus_id"]] else 0.0
                         for r in rows])
conf = {
    "confound_tx_length":  np.array([fnum(r["log_tx_len"]) for r in rows]),
    "confound_n_exons":    np.array([fnum(r["n_exons"]) for r in rows]),
    "confound_read_support": np.array([float(r["n_pos"]) for r in rows]),
}
# splice-only ridge (the four Pangolin features alone), same folds
def renorm(p):
    """project a prediction onto the simplex of each gene (fractions of a gene must sum to 1)"""
    q = np.clip(np.asarray(p, float), 1e-6, None); out = np.empty_like(q)
    for g, idx in gidx.items():
        t = q[idx].sum(); out[idx] = q[idx] / t if t > 0 else 1.0 / len(idx)
    return out

# ridge on logit(usage), resolved back to fractions by a within-gene softmax: the compositional
# version of the same model, same folds, same features
ylog = np.log(np.clip(y, 1e-3, 1 - 1e-3) / (1 - np.clip(y, 1e-3, 1 - 1e-3)))
oof_lg = np.full(len(y), np.nan)
for tr, te in gkf.split(X, y, groups):
    mu_, sd_, m_ = zfit(X[tr], ylog[tr]); oof_lg[te] = m_.predict((X[te] - mu_) / sd_)
soft = np.empty_like(oof_lg)
for g, idx in gidx.items():
    e = np.exp(oof_lg[idx] - oof_lg[idx].max()); soft[idx] = e / e.sum()

ss_names = ["ss_p_min", "ss_p_mean", "ss_use_min", "ss_use_mean",
            "mes_donor_min", "mes_donor_mean", "mes_acc_min", "mes_acc_mean"]
Xs = X[:, [names.index(n) for n in ss_names]]
oof_ss = np.full(len(y), np.nan)
for tr, te in gkf.split(Xs, y, groups):
    mu, sd, m = zfit(Xs[tr], y[tr]); oof_ss[te] = m.predict((Xs[te] - mu) / sd)

PRED = collections.OrderedDict([
    ("ridge", oof), ("ridge_renormalised", renorm(oof)), ("ridge_logit_softmax", soft),
    ("ridge_splice_only", oof_ss),
    ("uniform", pred_uniform), ("longest_wins", pred_longest), ("dominant_tss_only", pred_domtss),
] + list(conf.items()))
tick("baselines done")
MET = {k: metrics(v, y) for k, v in PRED.items()}
tick("metrics done")

# ---------------------------------------------------------------- within-gene ranking
def within_gene_acc(pred):
    """All isoform pairs inside a gene with a real usage difference: does `pred` order them right?
    Ties in pred count as half (chance), so 0.5 is the no-information baseline."""
    ok = n = 0.0
    for g, idx in gidx.items():
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                if abs(y[i] - y[j]) < 1e-6: continue
                if not (np.isfinite(pred[i]) and np.isfinite(pred[j])): continue
                n += 1
                d = (pred[i] - pred[j]) * (y[i] - y[j])
                ok += 1.0 if d > 0 else (0.5 if pred[i] == pred[j] else 0.0)
    return (ok / n if n else float("nan")), int(n)

tick("start within-gene")
def head_to_head(pa, pb):
    """paired (McNemar-style) comparison on the same within-gene pairs: how often A is right and
    B wrong, vs the reverse, with the normal-approximation p-value of the difference"""
    a_b = b_a = 0
    for g, idx in gidx.items():
        for u in range(len(idx)):
            for v in range(u + 1, len(idx)):
                i, j = idx[u], idx[v]
                if abs(y[i] - y[j]) < 1e-6: continue
                sa = (pa[i] - pa[j]) * (y[i] - y[j]); sb = (pb[i] - pb[j]) * (y[i] - y[j])
                ra, rb = sa > 0, sb > 0
                if ra and not rb: a_b += 1
                elif rb and not ra: b_a += 1
    n = a_b + b_a
    z = (a_b - b_a) / math.sqrt(n) if n else float("nan")
    return dict(a_right_b_wrong=a_b, b_right_a_wrong=b_a,
                p=float(2 * sps.norm.sf(abs(z))) if n else float("nan"))

WG = {}
for k, v in PRED.items():
    a, n = within_gene_acc(np.asarray(v, float))
    WG[k] = dict(acc=a, n_pairs=n, se=(math.sqrt(0.25 / n) if n else float("nan")))

# Marginal (and within-gene) rank correlation of every feature with usage, so a ridge
# coefficient whose sign is surprising can be checked against the raw association.
MARG = {}
yc_all = y.copy()
for g_, idx in gidx.items(): yc_all[idx] -= y[idx].mean()
for k in names:
    v = X[:, names.index(k)]
    vc = v.copy()
    for g_, idx in gidx.items(): vc[idx] -= v[idx].mean()
    MARG[k] = dict(label=FEATURES[k],
                   marginal_spearman=float(sps.spearmanr(v, y)[0]),
                   within_gene_spearman=float(sps.spearmanr(vc, yc_all)[0]))
with open(os.path.join(W, "marginal.csv"), "w") as fh:
    fh.write("name,label,marginal_spearman,within_gene_spearman\n")
    for k, d in sorted(MARG.items(), key=lambda kv: -abs(kv[1]["within_gene_spearman"])):
        fh.write(f'{k},"{d["label"]}",{d["marginal_spearman"]:.6g},{d["within_gene_spearman"]:.6g}\n')
tick("marginals done")

H2H = {f"ridge_vs_{k}": head_to_head(np.asarray(oof, float), np.asarray(PRED[k], float))
       for k in ("confound_tx_length", "confound_n_exons", "longest_wins", "ridge_splice_only")}
tick("head-to-head done")

# ---------------------------------------------------------------- final coefs + bootstrap over genes
tick("within-gene done")
mu, sd, mfull = zfit(X, y)
coef = mfull.coef_.copy()
glist = list(gidx)
B = int(os.environ.get("SR_BOOT", "200"))
boot = np.zeros((B, len(names)))
for b in range(B):
    gs = RNG.choice(len(glist), len(glist), replace=True)
    idx = np.concatenate([gidx[glist[k]] for k in gs])
    try:
        m2, s2, mb = zfit(X[idx], y[idx]); boot[b] = mb.coef_
    except Exception:
        boot[b] = np.nan
    if b % 10 == 0: tick(f"boot {b}")
lo = np.nanpercentile(boot, 2.5, 0); hi = np.nanpercentile(boot, 97.5, 0)

order = np.argsort(-np.abs(coef))
feat_out = [dict(name=names[i], label=FEATURES[names[i]], coef=float(coef[i]),
                 ci=[float(lo[i]), float(hi[i])]) for i in order]

# ---------------------------------------------------------------- app-gene (12) track
tick("bootstrap done")
app = json.load(open(os.path.join(ROOT, "data", "splice_data.json")))
mech = {}
for ln in open(os.path.join(SRC, "mech_pred.tsv")).read().splitlines()[1:]:
    g, iso, p = ln.split("\t"); mech[iso] = float(p)
tid2row = {}
for i, r in enumerate(rows): tid2row[r["transcript_id"]] = i
usage_by_tid = {}
for u in usage.values(): pass
tid_usage = {}
feat_by_tid = {}
tid_locus = {}
for r in feat:
    uu = usage.get(r["isoform_uid"])
    if uu and math.isfinite(fnum(uu["usage"])):
        tid_usage[r["transcript_id"]] = fnum(uu["usage"])
        tid_usage[r["transcript_id"].split(".")[0]] = fnum(uu["usage"])   # version-insensitive join
        feat_by_tid[r["transcript_id"].split(".")[0]] = r
        tid_locus[r["transcript_id"].split(".")[0]] = r["locus_id"]

# a ridge prediction for the app's transcripts that never saw those genes in training
app_loci = set()
for g in app["genes"]:
    for iso in g["isoforms"]:
        lid = tid_locus.get(iso["id"].split(".")[0])
        if lid: app_loci.add(lid)
keep_i = [i for i, r in enumerate(rows) if r["locus_id"] not in app_loci]
mu_a, sd_a, m_a = zfit(X[keep_i], y[keep_i])
def app_ridge(tid):
    r = feat_by_tid.get(tid.split(".")[0])
    if not r: return float("nan")
    d = dict(r)
    d["utr5_len_log"] = math.log10(1 + fnum(r["utr5_len"], 0))
    d["orf_len_log"] = math.log10(1 + fnum(r["orf_len_nt"], 0))
    d["n_iso_locus_log"] = math.log10(float(r["n_iso_locus"]))
    v = np.array([[fnum(d.get(n), 0.0) for n in names]], float)
    v = np.where(np.isfinite(v), v, 0.0)
    return float(m_a.predict((v - mu_a) / sd_a)[0])

# preferred app-gene target: every GENCODE v49 transcript of the 12 genes, straight from the
# genome-wide Isopedia query table (work/app12_usage.tsv from prep_app12.py). The curated
# catalogue keeps only a filtered subset, so it alone matches ~11 of the app's isoforms.
app_usage, app_support = {}, {}
_p12 = os.path.join(SRC, "app12_usage.tsv")
if os.path.exists(_p12):
    for r in read_tsv(_p12):
        for k in (r["transcript_id"], r["transcript_id"].split(".")[0]):
            app_usage[k] = fnum(r["usage"]); app_support[k] = float(r["n_pos"])
    tick(f"app12 usage: {len(app_usage)//2} transcripts")

# A like-for-like ridge for the 12 app genes: refit genome-wide on ONLY the features that can
# also be computed from data/splice_data.json, with the app's own genes excluded from training.
COMPACT = ["mes_donor_min", "mes_donor_mean", "mes_acc_min", "mes_acc_mean",
           "n_exons", "log_tx_len", "n_iso_locus_log", "is_longest", "log_dist_dom_tss", "is_dom_tss"]
ci = [names.index(n) for n in COMPACT if n in names]


def compact_feats(g, iso, exlen):
    """the COMPACT features, computed from data/splice_data.json the same way prep_features.py
    computes them from the Isopedia exon coordinates"""
    intr = [I for I in g["introns"] if iso["id"] in I["transcripts"]]
    dn = [I["donor"]["score"] for I in intr if I["donor"]["score"] is not None]
    ac = [I["acceptor"]["score"] for I in intr if I["acceptor"]["score"] is not None]
    lens = {i["id"]: sum(exlen.get(e, 0) for e in i["exon_ids"]) for i in g["isoforms"]}
    longest = max(lens, key=lambda k: lens[k])
    ex = sorted([(e["start"], e["end"]) for e in g["exons"] if e["id"] in iso["exon_ids"]])
    starts = collections.Counter()
    for i2 in g["isoforms"]:
        e2 = sorted([(e["start"], e["end"]) for e in g["exons"] if e["id"] in i2["exon_ids"]])
        if e2: starts[e2[0][0]] += 1
    dom = starts.most_common(1)[0][0] if starts else (ex[0][0] if ex else 0)
    return dict(_mes_donor_min=(min(dn) if dn else float("nan")),
                _mes_donor_mean=(float(np.mean(dn)) if dn else float("nan")),
                _mes_acc_min=(min(ac) if ac else float("nan")),
                _mes_acc_mean=(float(np.mean(ac)) if ac else float("nan")),
                _n_exons=float(len(iso["exon_ids"])),
                _log_tx_len=math.log10(max(lens[iso["id"]], 1)),
                _n_iso_locus_log=math.log10(len(g["isoforms"])),
                _is_longest=float(iso["id"] == longest),
                _log_dist_dom_tss=math.log10(1 + abs((ex[0][0] if ex else 0) - dom)),
                _is_dom_tss=float(bool(ex) and ex[0][0] == dom))

app_rows = []
for g in app["genes"]:
    tssw = {t["id"]: t.get("weight") for t in g["tss"]}
    mx = max([v for v in tssw.values() if v is not None] or [0])
    exlen = {e["id"]: e["len"] for e in g["exons"]}
    for iso in g["isoforms"]:
        tid = iso["id"]; bare = tid.split(".")[0]
        yv = app_usage.get(tid, app_usage.get(bare, tid_usage.get(tid, tid_usage.get(bare))))
        if yv is None or not math.isfinite(yv): continue
        intr = [I for I in g["introns"] if tid in I["transcripts"]]
        dn = [I["donor"]["score"] for I in intr if I["donor"]["score"] is not None]
        ac = [I["acceptor"]["score"] for I in intr if I["acceptor"]["score"] is not None]
        w = tssw.get(iso.get("tss_id"))
        app_rows.append(dict(gene=g["gene"], tid=tid, y=yv,
            read_support=app_support.get(tid, app_support.get(bare, float("nan"))),
            mech=mech.get(tid, float("nan")),
            mes_donor_min=(min(dn) if dn else float("nan")), mes_donor_mean=(float(np.mean(dn)) if dn else float("nan")),
            mes_acc_min=(min(ac) if ac else float("nan")), mes_acc_mean=(float(np.mean(ac)) if ac else float("nan")),
            n_exons=len(iso["exon_ids"]), tx_len=sum(exlen.get(e, 0) for e in iso["exon_ids"]),
            cage=(w if w is not None else float("nan")),
            is_dom_cage=(1.0 if (w is not None and mx > 0 and abs(w - mx) < 1e-12) else 0.0),
            ridge=app_ridge(tid),
            **compact_feats(g, iso, exlen)))
mu_c, sd_c, m_c = zfit(X[np.ix_(keep_i, ci)], y[keep_i])
Vc = np.array([[r["_" + n] for n in COMPACT] for r in app_rows], float)
Vc = np.where(np.isfinite(Vc), Vc, 0.0)
pc = m_c.predict((Vc - mu_c) / sd_c) if len(app_rows) else np.array([])
for i, r in enumerate(app_rows): r["ridge_compact"] = float(pc[i])

agid = collections.defaultdict(list)
for i, r in enumerate(app_rows): agid[r["gene"]].append(i)
for gname_, idx in agid.items():
    for i in idx: app_rows[i]["uniform"] = 1.0 / len(idx)
ay = np.array([r["y"] for r in app_rows])

def app_wg(key):
    ok = n = 0.0
    for gname_, idx in agid.items():
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                pi, pj = app_rows[i][key], app_rows[j][key]
                if abs(ay[i] - ay[j]) < 1e-6 or not (math.isfinite(pi) and math.isfinite(pj)): continue
                n += 1
                d = (pi - pj) * (ay[i] - ay[j])
                ok += 1.0 if d > 0 else (0.5 if pi == pj else 0.0)
    return (ok / n if n else float("nan")), int(n)

APP = {}
for key in ["mech", "ridge", "mes_donor_min", "mes_acc_min", "mes_donor_mean", "mes_acc_mean", "cage", "is_dom_cage", "tx_len", "n_exons", "read_support", "uniform", "ridge_compact"]:
    arr = np.array([r[key] for r in app_rows], float)
    a, n = app_wg(key)
    APP[key] = dict(**metrics(arr, ay), within_gene_acc=a, n_pairs=n)

json.dump(dict(rows=[{k: (None if isinstance(v, float) and not math.isfinite(v) else v) for k, v in r.items()}
                     for r in app_rows], metrics=APP),
          open(os.path.join(W, "app12.json"), "w"), indent=1)

# ---------------------------------------------------------------- outputs
# per_isoform is shipped to the browser, so it carries the 12 app genes in full plus an evenly
# spaced sample of the rest; the complete table is work/preds.tsv.
APP_GENES = {g["gene"] for g in app["genes"]}
CAP = 1200
others = [i for i, r in enumerate(rows) if r["gene_name"] not in APP_GENES]
step = max(1, len(others) // CAP)
keep_rows = set(i for i, r in enumerate(rows) if r["gene_name"] in APP_GENES) | set(others[::step])
per_iso = {}
Z = (X - mu) / sd
for i, r in enumerate(rows):
    if i not in keep_rows: continue
    contrib = sorted([dict(name=names[k], label=FEATURES[names[k]], value=float(X[i, k]),
                           contribution=float(Z[i, k] * coef[k])) for k in range(len(names))],
                     key=lambda d: -abs(d["contribution"]))[:8]
    per_iso[r["transcript_id"]] = dict(pred=float(oof[i]), measured=float(y[i]),
                                       gene=r["gene_name"], contrib=contrib)

out = dict(models=dict(usage=dict(
    dataset=("Isopedia human long-read isoform atlas (table_wide2 build, 869 PacBio/ONT samples); "
             "per-transcript usage = that transcript's read-weighted share of its gene, aggregated "
             "over all samples in which the gene was detected"),
    n=len(rows), n_genes=len({r["locus_id"] for r in rows}),
    cv=f"GroupKFold({NF}) grouped by gene; alphas chosen per fold by RidgeCV {sorted(set(alphas_used))}",
    target="transcript usage fraction (raw, not logit) so that ridge, the mechanistic model and the baselines are all on one scale",
    r2_cv=MET["ridge"]["r2"], pearson=MET["ridge"]["pearson"], spearman=MET["ridge"]["spearman"],
    within_gene_acc=WG["ridge"]["acc"], within_gene_pairs=WG["ridge"]["n_pairs"],
    features=feat_out, head_to_head=H2H, marginal=MARG,
    baselines={k: dict(**MET[k], within_gene_acc=WG[k]["acc"]) for k in PRED if k != "ridge"},
    app12=APP)), per_isoform=per_iso,
    meta=dict(methods=[], caveats=[]))
json.dump(out, open(os.path.join(W, "splice_ridge_core.json"), "w"), indent=1)

with open(os.path.join(W, "summary.txt"), "w") as fh:
    fh.write(f"n={len(rows)} isoforms, {len({r['locus_id'] for r in rows})} genes, pairs={WG['ridge']['n_pairs']}\n")
    fh.write(f"{'method':24s} {'r2':>8s} {'pearson':>8s} {'spearman':>9s} {'within-gene':>12s}\n")
    for k in PRED:
        m = MET[k]
        fh.write(f"{k:24s} {m['r2']:8.3f} {m['pearson']:8.3f} {m['spearman']:9.3f} {WG[k]['acc']:12.3f}\n")
    fh.write("\nhead-to-head (same pairs):\n")
    for k, v in H2H.items():
        fh.write(f"  {k:36s} ridge right/other wrong {v['a_right_b_wrong']:6d} | "
                 f"other right/ridge wrong {v['b_right_a_wrong']:6d} | p={v['p']:.2g}\n")
    fh.write("\n--- 12 app genes (mechanistic model available) ---\n")
    fh.write(f"n={len(app_rows)} isoforms, {len(agid)} genes\n")
    fh.write(f"{'method':18s} {'r2':>8s} {'spearman':>9s} {'within-gene':>12s} {'pairs':>6s}\n")
    for k, m in APP.items():
        fh.write(f"{k:18s} {m['r2']:8.3f} {m['spearman']:9.3f} {m['within_gene_acc']:12.3f} {m['n_pairs']:6d}\n")
print(open(os.path.join(W, "summary.txt")).read())
np.save(os.path.join(W, "oof.npy"), oof)
with open(os.path.join(W, "preds.tsv"), "w") as fh:
    fh.write("transcript_id\tlocus_id\tgene_name\tmeasured\t" + "\t".join(PRED) + "\n")
    for i, r in enumerate(rows):
        fh.write(f"{r['transcript_id']}\t{r['locus_id']}\t{r['gene_name']}\t{y[i]:.6g}\t"
                 + "\t".join(f"{np.asarray(v,float)[i]:.6g}" for v in PRED.values()) + "\n")
np.save(os.path.join(W, "boot.npy"), boot)
with open(os.path.join(W, "coefs.csv"), "w") as fh:
    fh.write("name,label,coef,ci_lo,ci_hi\n")
    for d in feat_out:
        fh.write(f'{d["name"]},"{d["label"]}",{d["coef"]:.6g},{d["ci"][0]:.6g},{d["ci"][1]:.6g}\n')
with open(os.path.join(W, "within_gene.csv"), "w") as fh:
    fh.write("method,accuracy,n_pairs,se\n")
    for k in PRED: fh.write(f"{k},{WG[k]['acc']:.6g},{WG[k]['n_pairs']},{WG[k]['se']:.6g}\n")
print("done", file=sys.stderr)
