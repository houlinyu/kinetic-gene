#!/usr/bin/env python
"""Build data/utr_data.json (see data/SCHEMA.md) for a list of human genes.

Usage:
    python pipeline/build_data.py                 # genes from pipeline/genes.txt
    python pipeline/build_data.py ATF4 MYC        # genes from the CLI
    python pipeline/build_data.py --max-iso 6 --no-cache ATF4

Data source: Ensembl REST API (GRCh38). Folding: ViennaRNA (RNAplfold W=150 L=100 + span-limited MFE; global MFE kept). --merge-only merges data/ridge.json.
Ensembl responses are cached in pipeline/cache/ so reruns are fast and offline-capable.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys
import time

import requests
import RNA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measured  # noqa: E402  (Cao et al. 2021 TE tables)
import optimus  # noqa: E402  (numpy Optimus 5-Prime, Sample et al. 2019)
import motifs  # noqa: E402  (sequence motifs / features, data/SPEC_v2.md)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT_DATA = os.path.join(ROOT, "data", "utr_data.json")
OUT_APP = os.path.join(ROOT, "app", "public", "utr_data.json")
CACHE = os.path.join(HERE, "cache")

SERVER = "https://rest.ensembl.org"
HEADERS = {"User-Agent": "utr-race-pipeline/0.1 (python-requests)",
           "Content-Type": "application/json", "Accept": "application/json"}
USE_CACHE = True

KOZAK_METHOD = ("rule-based: strong = purine (A/G) at -3 AND G at +4; adequate = exactly one of "
                "these; weak = neither (positions relative to first base of the start codon = +1; "
                "main start = Ensembl-annotated start codon, non-AUG flagged in tags; a -3 position "
                "upstream of the cap counts as not purine). score not computed (null).")

# --------------------------------------------------------------------------- Ensembl REST
_last_call = [0.0]


def _get(path, params=None, method="GET", body=None):
    key = hashlib.md5(json.dumps([method, path, params, body], sort_keys=True).encode()).hexdigest()
    cpath = os.path.join(CACHE, key + ".json")
    if USE_CACHE and os.path.exists(cpath):
        with open(cpath) as fh:
            return json.load(fh)
    for attempt in range(6):
        wait = 0.1 - (time.time() - _last_call[0])  # <= ~10 req/s (limit is 15)
        if wait > 0:
            time.sleep(wait)
        _last_call[0] = time.time()
        if method == "GET":
            r = requests.get(SERVER + path, params=params, headers=HEADERS, timeout=60)
        else:
            r = requests.post(SERVER + path, params=params, headers=HEADERS, json=body, timeout=120)
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 2)) + 0.5)
            continue
        if r.status_code >= 500:
            time.sleep(2 * (attempt + 1))
            continue
        r.raise_for_status()
        data = r.json()
        os.makedirs(CACHE, exist_ok=True)
        with open(cpath, "w") as fh:
            json.dump(data, fh)
        return data
    raise RuntimeError(f"Ensembl request failed repeatedly: {path}")


def ensembl_release():
    rel = _get("/info/data")["releases"]
    asm = _get("/info/assembly/homo_sapiens")["assembly_name"]
    return f"Ensembl {rel[0]} ({asm})"


def fetch_gene(symbol):
    """Return list of dicts: id, name, tags, utr5 (RNA), cds (RNA from AUG to end of cDNA)."""
    g = _get(f"/lookup/symbol/homo_sapiens/{symbol}", {"expand": 0})
    gid = g["id"]
    feats = _get(f"/overlap/id/{gid}", {"feature": "transcript"})
    tx = {}
    for f in feats:
        if f.get("Parent") != gid or f.get("biotype") != "protein_coding":
            continue
        tags = ["basic" if t == "gencode_basic" else t for t in (f.get("tag") or [])]
        tx[f["transcript_id"]] = {"id": f"{f['transcript_id']}.{f['version']}",
                                  "name": f.get("external_name"), "tags": tags}
    # exon structure + genomic CDS start (for tss_note)
    full = _get(f"/lookup/id/{gid}", {"expand": 1})
    for tr in full.get("Transcript", []):
        if tr["id"] not in tx or not tr.get("Translation"):
            continue
        strand = tr["strand"]
        exons = sorted(((e["start"], e["end"]) for e in tr["Exon"]), key=lambda e: e[0], reverse=strand < 0)
        cds0 = tr["Translation"]["start"] if strand > 0 else tr["Translation"]["end"]
        tx[tr["id"]].update(strand=strand, exons=exons, cds_gstart=cds0, chrom=tr["seq_region_name"])
    ids = sorted(tx)
    seqs = {}
    for i in range(0, len(ids), 50):
        res = _get("/sequence/id", {"type": "cdna", "mask_feature": 1}, method="POST",
                   body={"ids": ids[i:i + 50]})
        for rec in res:
            seqs[rec["query"]] = rec["seq"]
    out = []
    for tid in ids:
        s = seqs.get(tid)
        if not s:
            continue
        n5 = 0
        while n5 < len(s) and s[n5].islower():
            n5 += 1
        utr = s[:n5].upper().replace("T", "U")
        rest = s[n5:].upper().replace("T", "U")
        t = tx[tid]
        # complete 5' UTR: non-empty, CDS starts with AUG, no start_NF tags
        if n5 == 0 or len(rest) < 4:
            continue
        if rest[:3] != "AUG":
            # Ensembl annotates some starts as non-AUG (e.g. MYC/VEGFA MANE use CUG).
            # Keep them (annotated start codon is used as the main start) but flag it.
            if rest[1:3] != "UG":
                continue
            t["tags"] = t["tags"] + [f"non_AUG_start:{rest[:3]}"]
        if any("start_NF" in x for x in t["tags"]):
            continue
        if set(utr) - set("ACGU"):
            continue
        t.update(utr5=utr, cds=rest)
        if "exons" in t:
            t["utr_segs"] = utr_segments(t)
            if sum(hi - lo + 1 for lo, hi in t["utr_segs"]) != len(utr):
                t["utr_segs"] = None  # structure inconsistent with sequence; tss_note -> null
        out.append(t)
    return gid, out


# --------------------------------------------------------------------------- Optimus 5-Prime
NOTES = {"genes": {}, "sources": {}}


def optimus_method():
    m, sd = optimus.train_scaler()
    v = json.load(open(os.path.join(HERE, "external", "optimus", "validation.json")))
    return ("Optimus 5-Prime variable-length model (Sample et al. 2019; trained on 25-100 nt random 5' UTRs in "
            "an MPRA, HEK293T). Input: whole UTR if len <= 100 (optimus_input='full'), else the 100 nt "
            "immediately upstream of the start codon ('last100'); U->T, left-padded with N (zeros) to 100. "
            f"Output inverse-scaled with the training-set MRL mean {m:.4f} / sd {sd:.4f} (recomputed from GEO "
            "GSM4084997 with the notebook's split). Validation of this numpy forward pass on the notebook's "
            f"held-out sets: random test Pearson r = {v['random_test']['pearson_r']:.3f} (n={v['random_test']['n']}), "
            f"human test r = {v['human_test']['pearson_r']:.3f} (n={v['human_test']['n']}). Predicted MRL is for "
            "the UTR in the MPRA reporter context (eGFP CDS), not the endogenous mRNA; UTRs < 25 nt are outside "
            "the training length range.")


# --------------------------------------------------------------------------- TSS / first-exon heuristic
TSS_METHOD = ("Heuristic vs the gene's reference isoform (first chosen: MANE_Select/canonical), from Ensembl "
              "exon coordinates (lookup/id?expand=1). 5' UTR genomic segments = exons clipped at the annotated "
              "CDS start. class: reference; alt_first_exon = first exons do not overlap genomically; "
              "utr_splicing = first exons overlap but the set of 5' UTR introns inside the region covered by "
              "both UTRs differs; length_only = same start codon and one UTR sequence is a 3'-suffix of the "
              "other (TSS shifted within the same exon structure); alt_start_codon = otherwise, the annotated "
              "start codon differs; other = none of these. features lists every criterion that applies.")


def utr_segments(t):
    """Genomic (lo, hi) inclusive intervals of the 5' UTR in transcript order."""
    segs, c, strand = [], t["cds_gstart"], t["strand"]
    for lo, hi in t["exons"]:
        if strand > 0:
            if lo >= c:
                break
            segs.append((lo, min(hi, c - 1)))
            if hi >= c:
                break
        else:
            if hi <= c:
                break
            segs.append((max(lo, c + 1), hi))
            if lo <= c:
                break
    return segs


def _introns(segs, strand):
    out = set()
    for a, b in zip(segs, segs[1:]):
        out.add((a[1], b[0]) if strand > 0 else (b[1], a[0]))
    return out


def tss_note(t, ref):
    if not t.get("utr_segs") or not ref.get("utr_segs"):
        return None
    st = t["strand"]
    tss = t["utr_segs"][0][1] if st < 0 else t["utr_segs"][0][0]
    base = {"vs": ref["id"], "tss": f"{t['chrom']}:{tss}:{'+' if st > 0 else '-'}"}
    if t is ref:
        return dict(base, **{"class": "reference", "features": [],
                             "text": "reference isoform for comparison (MANE Select / canonical)"})
    fe_t, fe_r = t["exons"][0], ref["exons"][0]
    first_overlap = fe_t[0] <= fe_r[1] and fe_r[0] <= fe_t[1]
    same_start = t["cds_gstart"] == ref["cds_gstart"]
    span = lambda s: (min(x[0] for x in s), max(x[1] for x in s))
    a, b = span(t["utr_segs"]), span(ref["utr_segs"])
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    it = {x for x in _introns(t["utr_segs"], st) if lo <= x[0] and x[1] <= hi}
    ir = {x for x in _introns(ref["utr_segs"], st) if lo <= x[0] and x[1] <= hi}
    suffix = t["utr5"].endswith(ref["utr5"]) or ref["utr5"].endswith(t["utr5"])
    feats = []
    if not first_overlap:
        feats.append("alt_first_exon")
    if first_overlap and it != ir:
        feats.append("utr_splicing")
    if not same_start:
        feats.append("alt_start_codon")
    if same_start and suffix:
        feats.append("length_only")
    rname = ref.get("name") or ref["id"]
    if "alt_first_exon" in feats:
        cls, txt = "alt_first_exon", f"alternative first exon / TSS vs {rname}"
    elif "utr_splicing" in feats:
        cls, txt = "utr_splicing", f"same first exon as {rname}, different splicing within the 5' UTR"
    elif "length_only" in feats:
        d = len(t["utr5"]) - len(ref["utr5"])
        cls, txt = "length_only", (f"same exon structure as {rname}; TSS shifted, UTR "
                                   f"{'longer' if d > 0 else 'shorter'} by {abs(d)} nt")
    elif "alt_start_codon" in feats:
        cls, txt = "alt_start_codon", f"annotated start codon differs from {rname}"
    else:
        cls, txt = "other", f"differs from {rname} (TSS within same first exon, no UTR intron difference)"
    if not same_start and cls != "alt_start_codon":
        txt += "; annotated start codon also differs"
    return dict(base, **{"class": cls, "features": feats, "text": txt})


# --------------------------------------------------------------------------- selection
def tag_score(tags):
    return (3 * ("MANE_Select" in tags) + 2 * ("Ensembl_canonical" in tags)
            + 1 * ("MANE_Plus_Clinical" in tags) + 1 * ("basic" in tags)
            + 0.5 * ("gencode_primary" in tags))


def kmers(s, k=6):
    return {s[i:i + k] for i in range(max(1, len(s) - k + 1))}


def utr_distance(a, b):
    ka, kb = kmers(a), kmers(b)
    jac = len(ka & kb) / max(1, len(ka | kb))
    return 1.0 - jac


def select_isoforms(txs, max_iso):
    # dedupe identical UTRs -> keep best-tagged id (tie: lowest id)
    best = {}
    for t in txs:
        cur = best.get(t["utr5"])
        if cur is None or (-tag_score(t["tags"]), t["id"]) < (-tag_score(cur["tags"]), cur["id"]):
            best[t["utr5"]] = t
    pool = sorted(best.values(), key=lambda t: (-tag_score(t["tags"]), t["id"]))
    chosen = []
    for t in pool:  # always keep MANE_Select / canonical
        if ("MANE_Select" in t["tags"] or "Ensembl_canonical" in t["tags"]) and len(chosen) < max_iso:
            chosen.append(t)
    if not chosen and pool:
        chosen.append(pool[0])
    rest = [t for t in pool if t not in chosen]
    while rest and len(chosen) < max_iso:  # greedy max-min UTR distance (6-mer Jaccard)
        nxt = max(rest, key=lambda t: (round(min(utr_distance(t["utr5"], c["utr5"]) for c in chosen), 6),
                                       tag_score(t["tags"])))
        chosen.append(nxt)
        rest.remove(nxt)
    return chosen, len(txs), len(best)


# --------------------------------------------------------------------------- features
def kozak_class(full, a):
    """full = UTR+CDS RNA; a = index of A of AUG."""
    m3 = full[a - 3] if a - 3 >= 0 else "N"
    p4 = full[a + 3] if a + 3 < len(full) else "N"
    n = (m3 in "AG") + (p4 == "G")
    return ["weak", "adequate", "strong"][n]


def find_uorfs(utr, cds):
    full = utr + cds
    L = len(utr)
    res = []
    for i in range(L - 2):
        if utr[i:i + 3] != "AUG":
            continue
        stop = None
        for j in range(i + 3, len(full) - 2, 3):
            if full[j:j + 3] in ("UAA", "UAG", "UGA"):
                stop = j
                break
        inframe = (L - i) % 3 == 0
        if stop is None:
            end, aa = len(full), None
            typ = "in-frame-extension" if inframe else "overlapping"
        else:
            end = stop + 3
            aa = (stop - i) // 3
            if inframe and stop >= L:
                typ = "in-frame-extension"
            elif end <= L:
                typ = "uORF"
            else:
                typ = "overlapping"
        res.append({"start": i, "end": end, "aa": aa, "kozak_class": kozak_class(full, i), "type": typ})
    return res


def find_hairpins(db, seq, min_bp=3):
    pt = [-1] * len(db)
    st = []
    for i, c in enumerate(db):
        if c == "(":
            st.append(i)
        elif c == ")":
            j = st.pop()
            pt[i], pt[j] = j, i
    out = []
    for k, l in enumerate(pt):
        if l <= k:
            continue
        if any(pt[x] != -1 for x in range(k + 1, l)):
            continue  # not a hairpin-closing pair
        i, j = k, l
        while True:  # extend outward through stacks / bulges / interior loops
            a = i - 1
            while a >= 0 and pt[a] == -1:
                a -= 1
            b = j + 1
            while b < len(pt) and pt[b] == -1:
                b += 1
            if a >= 0 and b < len(pt) and pt[a] == b:
                i, j = a, b
            else:
                break
        sub = db[i:j + 1]
        bp = sub.count("(")
        if bp < min_bp:
            continue
        dG = RNA.fold_compound(seq[i:j + 1]).eval_structure(sub)
        out.append({"start": i, "end": j + 1, "stem_bp": bp, "dG": round(dG, 2)})
    return out



def pairs_of(db):
    st, out = [], []
    for i, c in enumerate(db):
        if c == "(":
            st.append(i)
        elif c == ")":
            out.append((st.pop(), i))
    return out


def layout(db):
    RNA.cvar.rna_plot_type = RNA.PLOT_TYPE_NAVIEW
    c = RNA.get_xy_coordinates(db)
    xy = [(c.get(i).X, c.get(i).Y) for i in range(len(db))]
    cx = sum(p[0] for p in xy) / len(xy)
    cy = sum(p[1] for p in xy) / len(xy)
    m = max(max(abs(p[0] - cx), abs(p[1] - cy)) for p in xy) or 1.0
    return [[round((x - cx) / m, 4), round((y - cy) / m, 4)] for x, y in xy]


PLF_W, PLF_L = 150, 100
FOLD_METHOD = (f"local folding: per-nt unpaired probabilities from RNAplfold (RNA.pfl_fold_up, W={PLF_W}, "
               f"L={PLF_L}, ulength 1) on the UTR + first 30 nt of CDS (punp reported for the UTR only); structure "
               f"(dotbracket, layout2d, hairpins) = MFE with maximal base-pair span {PLF_L} nt (RNAfold --maxBPspan "
               f"{PLF_L}), so long-range pairs are dropped. mfe / mfe_global = unrestricted global MFE (old field, "
               "kept for reference); mfe_local = MFE of the span-restricted structure. ViennaRNA default energy "
               "parameters, 37 C.")
ACC_FLANK = 15


def local_fold(utr, cds):
    md = RNA.md()
    md.max_bp_span = PLF_L
    db, mfe_local = RNA.fold_compound(utr, md).mfe()
    seq = utr + cds[:30]
    up = RNA.pfl_fold_up(seq, 1, PLF_W, PLF_L)  # 1-based; up[i][1] = P(i unpaired)
    pu = [up[i + 1][1] for i in range(len(seq))]
    L = len(utr)
    win = pu[max(0, L - ACC_FLANK):L + ACC_FLANK]
    return db, mfe_local, [round(x, 2) for x in pu[:L]], sum(win) / len(win)


def analyse(t):
    utr, cds = t["utr5"], t["cds"]
    L = len(utr)
    t0 = time.time()
    db_g, mfe = RNA.fold_compound(utr).mfe()
    db, mfe_local, punp, acc = local_fold(utr, cds)
    fold_s = time.time() - t0
    cap_db, cap_mfe = RNA.fold_compound(utr[:50]).mfe()
    ctx = ("N" * max(0, 6 - L) + utr[-6:]) + cds[:4]
    return {
        "id": t["id"],
        "name": t["name"],
        "tags": t["tags"],
        "utr5": utr,
        "len": L,
        "cds_context": ctx,
        "main_kozak": {"score": None, "class": kozak_class(utr + cds, L), "method": KOZAK_METHOD},
        "uorfs": find_uorfs(utr, cds),
        "fold": {"method": FOLD_METHOD, "mfe": round(mfe, 2), "mfe_global": round(mfe, 2),
                 "mfe_local": round(mfe_local, 2), "dotbracket": db, "cap_dG": round(cap_mfe, 2),
                 "layout2d": layout(db), "punp": punp, "acc_start": round(acc, 3)},
        "hairpins": find_hairpins(db, utr),
        "gc": round(sum(utr.count(x) for x in "GC") / L, 4),
        "pred": {"optimus_mrl": None, "optimus_input": None, "utrlm_te": None},
        "measured": measured.lookup(t["id"]),
        "tss_note": None,
        "motifs": motifs.find_motifs(utr, cds),
    }, fold_s


MOTIF_METHOD = ("pipeline/motifs.py find_motifs(); positions 0-based half-open in the 5' UTR. uAUG (strength = Kozak "
                "class, -3 purine / +4 G; cite calvo2009) with uORF_stop for uORFs ending in the UTR; nearcog = CUG/GUG/"
                "ACG/UUG in strong context only (kearse2017); TOP = >= 5 pyrimidines starting within 4 nt of the 5' end "
                "(TOP-like rule of thoreen2012; Ensembl 5' ends are approximate); G4 = (G3+N1-7){3,}G3+ "
                "(huppert2005 folding rule; kumari2007 for 5' UTR translational repression), not structure-tested; "
                "TISU = AAS-AUG-GCGGC at the main start for UTR <= 30 nt (elfakess2008); kozak_main = -6..+4 context "
                "(kozak1986). Dropped (no definition verified in time): PRTE, polyU.")


def merge_ridge(data, path=None):
    """Merge data/ridge.json (RIDGE agent) into isoform.ridge + meta.ridge. Returns problems list."""
    path = path or os.path.join(ROOT, "data", "ridge.json")
    if not os.path.exists(path):
        return []
    rj = json.load(open(path))
    per = rj.get("per_isoform", {}) or {}
    probs, n = [], 0
    for g in data["genes"]:
        for iso in g["isoforms"]:
            r = per.get(iso["id"]) or per.get(iso["id"].split(".")[0])
            iso["ridge"] = r
            n += r is not None
    data["meta"]["ridge"] = {k: v for k, v in rj.items() if k != "per_isoform"}
    keys = {i["id"] for g in data["genes"] for i in g["isoforms"]} | {i["id"].split(".")[0] for g in data["genes"] for i in g["isoforms"]}
    missing = [k for k in per if k not in keys]
    if missing:
        probs.append(f"ridge.json: {len(missing)} per_isoform ids not in data (e.g. {missing[:3]})")
    print(f"merged ridge.json: {n} isoforms with ridge, models {list(rj.get('models', {}))}")
    return probs


def write_out(data):
    os.makedirs(os.path.dirname(OUT_DATA), exist_ok=True)
    os.makedirs(os.path.dirname(OUT_APP), exist_ok=True)
    with open(OUT_DATA, "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    shutil.copyfile(OUT_DATA, OUT_APP)
    print(f"\nwrote {OUT_DATA} ({os.path.getsize(OUT_DATA)/1e6:.2f} MB) and {OUT_APP}")


# --------------------------------------------------------------------------- main
def main():
    global USE_CACHE
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="*")
    ap.add_argument("--genes-file", default=os.path.join(HERE, "genes.txt"))
    ap.add_argument("--max-iso", type=int, default=6)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--merge-only", action="store_true",
                    help="reload data/utr_data.json, merge data/ridge.json, write both copies")
    args = ap.parse_args()
    if args.merge_only:
        data = json.load(open(OUT_DATA))
        probs = merge_ridge(data)
        write_out(data)
        if probs:
            print("RIDGE MERGE PROBLEMS:\n  " + "\n  ".join(probs))
            sys.exit(1)
        return
    USE_CACHE = not args.no_cache
    genes = args.genes or [l.split("#")[0].strip() for l in open(args.genes_file)]
    genes = [g for g in genes if g]

    global NOTES
    with open(os.path.join(HERE, "gene_notes.json")) as fh:
        NOTES = json.load(fh)
    out_genes, failed, problems = [], [], []
    print(f"{'gene':10s} {'iso':>3s} {'pc_tx':>5s} {'uniqUTR':>7s}  isoform details (name len #uORF[upstream/overlap/ext] MFE capdG fold_s)")
    for g in genes:
        try:
            gid, txs = fetch_gene(g)
            chosen, n_tx, n_uniq = select_isoforms(txs, args.max_iso)
            if not chosen:
                raise RuntimeError("no protein-coding transcript with complete 5' UTR")
            isos = []
            details = []
            opt_in = [t["utr5"] if len(t["utr5"]) <= 100 else t["utr5"][-100:] for t in chosen]
            opt = optimus.predict_mrl(opt_in)
            for t, mrl in zip(chosen, opt):
                rec, fs = analyse(t)
                rec["pred"]["optimus_mrl"] = round(float(mrl), 3)
                rec["pred"]["optimus_input"] = "full" if rec["len"] <= 100 else "last100"
                rec["tss_note"] = tss_note(t, chosen[0])
                isos.append(rec)
                L = rec["len"]
                # sanity checks
                if len(rec["fold"]["layout2d"]) != L:
                    problems.append(f"{g} {rec['id']}: layout2d len {len(rec['fold']['layout2d'])} != {L}")
                if len(rec["fold"]["dotbracket"]) != L:
                    problems.append(f"{g} {rec['id']}: dotbracket len mismatch")
                if len(rec["cds_context"]) != 10 or rec["cds_context"][6:9] != t["cds"][:3]:
                    problems.append(f"{g} {rec['id']}: bad cds_context {rec['cds_context']}")
                for u in rec["uorfs"]:
                    if not (0 <= u["start"] < L and u["end"] > u["start"]) or rec["utr5"][u["start"]:u["start"] + 3] != "AUG":
                        problems.append(f"{g} {rec['id']}: uORF out of range {u}")
                if len(rec["fold"]["punp"]) != L or not all(0 <= x <= 1 for x in rec["fold"]["punp"]):
                    problems.append(f"{g} {rec['id']}: punp len/range bad")
                if not (0 <= rec["fold"]["acc_start"] <= 1):
                    problems.append(f"{g} {rec['id']}: acc_start out of range")
                if any(abs(i - j) > PLF_L for i, j in pairs_of(rec["fold"]["dotbracket"])):
                    problems.append(f"{g} {rec['id']}: base pair span > {PLF_L}")
                if sum(m["type"] == "kozak_main" for m in rec["motifs"]) != 1:
                    problems.append(f"{g} {rec['id']}: kozak_main count != 1")
                for m in rec["motifs"]:
                    if not (0 <= m["start"] < m["end"] <= L) or m["cite"] not in motifs.SOURCES:
                        problems.append(f"{g} {rec['id']}: motif out of range / uncited {m}")
                    elif m["type"] in ("uAUG", "uORF_stop", "nearcog", "G4", "TOP") and rec["utr5"][m["start"]:m["end"]] != m["seq"]:
                        problems.append(f"{g} {rec['id']}: motif seq mismatch {m}")
                for h in rec["hairpins"]:
                    if not (0 <= h["start"] < h["end"] <= L):
                        problems.append(f"{g} {rec['id']}: hairpin out of range {h}")
                ty = [u["type"] for u in rec["uorfs"]]
                details.append(f"{rec['name']}{'*' if 'MANE_Select' in rec['tags'] else ''} "
                               f"{L}nt {len(ty)}[{ty.count('uORF')}/{ty.count('overlapping')}/"
                               f"{ty.count('in-frame-extension')}] {rec['fold']['mfe']}/{rec['fold']['mfe_local']} "
                               f"{rec['fold']['cap_dG']} hp{len(rec['hairpins'])} m{len(rec['motifs'])}"
                               + (f" fold={fs:.1f}s" if L > 1000 else ""))
            gn = NOTES["genes"].get(g, {})
            out_genes.append({"gene": g, "blurb": gn.get("blurb", ""), "cite": gn.get("cite", []), "isoforms": isos})
            print(f"{g:10s} {len(isos):3d} {n_tx:5d} {n_uniq:7d}  " + " | ".join(details))
        except Exception as e:  # keep going for other genes
            failed.append((g, repr(e)))
            print(f"{g:10s} FAILED: {e!r}")

    data = {
        "meta": {
            "generated": datetime.date.today().isoformat(),
            "gencode_or_ensembl_release": ensembl_release(),
            "vienna_version": RNA.__version__,
            "sources": [
                {"key": "ensembl", "citation": "Ensembl REST API (rest.ensembl.org), Homo sapiens GRCh38",
                 "url": "https://rest.ensembl.org"},
                {"key": "viennarna", "citation": "Lorenz R. et al. ViennaRNA Package 2.0. Algorithms Mol Biol 6:26 (2011)",
                 "url": "https://doi.org/10.1186/1748-7188-6-26"},
                {"key": "sample2019", "citation": "Sample PJ et al. Human 5' UTR design and variant effect prediction from a "
                 "massively parallel translation assay. Nat Biotechnol 37:803-809 (2019). Optimus 5-Prime "
                 "Varying_length_25to100_model.hdf5 (github.com/pjsample/human_5utr_modeling, GPL-3.0), forward pass "
                 "re-implemented in numpy (pipeline/optimus.py)",
                 "url": "https://www.nature.com/articles/s41587-019-0164-5"},
                {"key": "cao2021", "citation": "Cao J et al. High-throughput 5' UTR engineering for enhanced protein "
                 "production in non-viral gene therapies. Nat Commun 12:4138 (2021). TE tables from "
                 "github.com/zzz2010/5UTR_Optimizer (MIT)",
                 "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8260622/"},
            ] + [dict(key=k, **v) for k, v in NOTES["sources"].items()]
              + [dict(key=k, **v) for k, v in motifs.SOURCES.items() if k not in NOTES["sources"]],
            "methods": {
                "optimus": optimus_method(),
                "measured": measured.NOTE,
                "tss_note": TSS_METHOD,
                "fold": FOLD_METHOD,
                "motifs": MOTIF_METHOD,
            },
        },
        "genes": out_genes,
        "findings": NOTES.get("findings", []),
    }
    ids = {i["id"]: g["gene"] for g in out_genes for i in g["isoforms"]}
    keys = {x["key"] for x in data["meta"]["sources"]}
    for f in data["findings"]:
        for side in ("better", "worse"):
            if ids.get(f[side]) != f["gene"]:
                problems.append(f"finding {f['id']}: {side} {f[side]} not among {f['gene']} isoforms")
        if f["cite"] not in keys:
            problems.append(f"finding {f['id']}: cite {f['cite']} not in meta.sources")
    problems += merge_ridge(data)
    write_out(data)
    print("* = MANE_Select; uORF counts [fully-upstream/overlapping-CDS/in-frame-extension]")
    if failed:
        print("FAILED genes:", failed)
    if problems:
        print("SANITY PROBLEMS:\n  " + "\n  ".join(problems))
        sys.exit(1)
    print("sanity checks passed (layout2d/dotbracket/punp length == len, pair span <= 100, motif/uORF/hairpin positions in range, start codon at cds_context[6:9])")


if __name__ == "__main__":
    main()
