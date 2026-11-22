"""Discovery calls from model_out.json (exact analytic model, many feasible parameter points x 2 fold methods).

A call is kept only if it holds at EVERY parameter point (default + all feasible grid points) AND for BOTH
fold methods ("robust"). Call types, per same-CDS isoform pair (A, B) of a gene:
  leader_effect : calm output differs >= 2-fold (|log2| >= 1) in the same direction everywhere.
  stress_flip   : calm A < B (by >= 1.25x) but under stress (tc=0.15) A > B (by >= 1.25x) everywhere:
                  the ISR is predicted to switch which leader dominates protein output.
  isr_escape    : A's stress/calm ratio >= 1 (not repressed or induced) while B's <= 0.6, everywhere.
Plus a model-validation check against Cao 2021 measured TE (within-gene, unambiguous values only).
Writes data/discovery.json and data/figures/discovery_*.png (+ .json with the plotted numbers).
"""
import itertools
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = os.path.join(ROOT, "pipeline/discovery/work")
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
import measured  # noqa: E402  read-only import (Cao 2021 TE tables)

FOLDS = ("global", "local")
EPS = 1e-9


def arr(r, cond):
    return np.array([r[cond][f] for f in FOLDS])  # shape (2, npoints)


def abl(r, key):
    return np.array([r["abl"][f][key] for f in FOLDS])


QC = QS = None  # 25th percentile of calm / stress output over all isoforms, per fold x point (set in main)


def holds(kind, ca, cb, sa, sb):
    lc = np.log2(np.maximum(ca, EPS) / np.maximum(cb, EPS))
    ls = np.log2(np.maximum(sa, EPS) / np.maximum(sb, EPS))
    if kind == "leader_effect":
        return bool((lc >= 1).all() and (ca >= QC).all())
    if kind == "stress_flip":
        return bool((lc <= -math.log2(1.25)).all() and (ls >= math.log2(1.25)).all() and (sa >= QS).all())
    ra, rb = sa / np.maximum(ca, EPS), sb / np.maximum(cb, EPS)
    return bool((ra >= 1).all() and (rb <= 0.6).all() and (sa >= QS).all() and (ca > EPS).all())


def mechanism(kind, ra, rb):
    """Which model ingredient the call needs: re-test with all upstream AUGs disabled / all hairpins melted."""
    noU = holds(kind, abl(ra, "noU_calm"), abl(rb, "noU_calm"), abl(ra, "noU_stress"), abl(rb, "noU_stress"))
    noH = holds(kind, abl(ra, "noH_calm"), abl(rb, "noH_calm"), abl(ra, "noH_stress"), abl(rb, "noH_stress"))
    if noU and noH:
        both = holds(kind, abl(ra, "noUH_calm"), abl(rb, "noUH_calm"), abl(ra, "noUH_stress"), abl(rb, "noUH_stress"))
        return "other (length / Kozak / cap-proximal AUG)" if both else "redundant (uORF or structure each sufficient)"
    if not noU and noH:
        return "uORF/uAUG"
    if noU and not noH:
        return "hairpin/structure"
    if not noU and not noH:
        return "uORF+structure"


def summarize_iso(i):
    u = i["uorfs"]
    return {
        "id": i["id"], "name": i["name"], "len": i["len"], "tags": [t for t in i["tags"] if t in ("MANE_Select", "Ensembl_canonical", "basic")],
        "tsl": i["tsl"], "main_kozak": i["main_kozak"]["class"],
        "n_uORF": sum(x["type"] == "uORF" for x in u),
        "n_overlapping": sum(x["type"] == "overlapping" for x in u),
        "n_inframe_ext": sum(x["type"] == "in-frame-extension" for x in u),
        "uorfs": [{k: x[k] for k in ("start", "end", "aa", "kozak_class", "type")} for x in u][:12],
        "min_hairpin_dG_global": min([h["dG"] for h in i["fold"]["global"]["hairpins"]], default=None),
        "cap_hairpin_dG_global": min([h["dG"] for h in i["fold"]["global"]["hairpins"] if h["start"] < 15], default=None),
        "mfe_global": i["fold"]["global"]["mfe"], "gc": i["gc"],
    }


def hc(i):
    return "MANE_Select" in i["tags"] or i["tsl"] in ("1", "2")


def main():
    iso_doc = json.load(open(os.path.join(W, "isoforms.json")))
    mo = json.load(open(os.path.join(W, "model_out.json")))
    feas = json.load(open(os.path.join(W, "feasible.json")))
    M = mo["iso"]
    global QC, QS
    QC = np.percentile(np.array([arr(r, "calm") for r in M.values()]), 25, axis=0)
    QS = np.percentile(np.array([arr(r, "stress") for r in M.values()]), 25, axis=0)
    npts = len(mo["points"])
    calls = {"leader_effect": [], "stress_flip": [], "isr_escape": []}
    n_pairs = n_pairs_hc = 0
    all_ratio_default = []
    for g in iso_doc["groups"]:
        for i in g["isoforms"]:
            r = M[i["id"]]
            c, s = r["calm"]["global"][0], r["stress"]["global"][0]
            if c > EPS:
                all_ratio_default.append(s / c)
        for A, B in itertools.permutations(g["isoforms"], 2):
            ca, cb = arr(M[A["id"]], "calm"), arr(M[B["id"]], "calm")
            sa, sb = arr(M[A["id"]], "stress"), arr(M[B["id"]], "stress")
            if A["id"] < B["id"]:
                n_pairs += 1
                n_pairs_hc += hc(A) and hc(B)
            if (ca <= EPS).any() or (cb <= EPS).any():
                continue
            lc = np.log2(ca / cb)
            ls = np.log2(np.maximum(sa, EPS) / np.maximum(sb, EPS))
            ra, rb = sa / ca, sb / cb
            base = {"gene": g["gene"], "A": summarize_iso(A), "B": summarize_iso(B),
                    "both_high_conf": bool(hc(A) and hc(B)),
                    "default": {"calm_A": round(float(ca[0, 0]), 4), "calm_B": round(float(cb[0, 0]), 4),
                                "stress_A": round(float(sa[0, 0]), 4), "stress_B": round(float(sb[0, 0]), 4),
                                "log2_calm_A_over_B": round(float(lc[0, 0]), 3),
                                "log2_stress_A_over_B": round(float(ls[0, 0]), 3),
                                "stress_over_calm_A": round(float(ra[0, 0]), 3), "stress_over_calm_B": round(float(rb[0, 0]), 3)},
                    "range_over_points_and_folds": {"log2_calm": [round(float(lc.min()), 3), round(float(lc.max()), 3)],
                                                    "log2_stress": [round(float(ls.min()), 3), round(float(ls.max()), 3)],
                                                    "ratio_A": [round(float(ra.min()), 3), round(float(ra.max()), 3)],
                                                    "ratio_B": [round(float(rb.min()), 3), round(float(rb.max()), 3)]}}
            # each ordered pair; A = the isoform the call is about (better / escaping / winning under stress)
            mA, mB = M[A["id"]], M[B["id"]]
            if holds("leader_effect", ca, cb, sa, sb):
                calls["leader_effect"].append({**base, "score": float(lc.min()), "mechanism": mechanism("leader_effect", mA, mB)})
            if holds("stress_flip", ca, cb, sa, sb):
                calls["stress_flip"].append({**base, "score": float(min(-lc.max(), ls.min())), "mechanism": mechanism("stress_flip", mA, mB)})
            if holds("isr_escape", ca, cb, sa, sb):
                calls["isr_escape"].append({**base, "score": float(ra.min() / rb.max()), "mechanism": mechanism("isr_escape", mA, mB)})
    induced = []
    for g in iso_doc["groups"]:
        for i in g["isoforms"]:
            r = M[i["id"]]
            c, st = arr(r, "calm"), arr(r, "stress")
            if (c > EPS).all() and (st / c >= 1.5).all() and (st >= QS).all():
                sib = [x for x in g["isoforms"] if x is not i]
                induced.append({"gene": g["gene"], "A": summarize_iso(i), "high_conf": hc(i),
                                "stress_over_calm_default": round(float(st[0, 0] / c[0, 0]), 3),
                                "stress_over_calm_range": [round(float((st / c).min()), 3), round(float((st / c).max()), 3)],
                                "calm_default": round(float(c[0, 0]), 4),
                                "siblings": [{"name": x["name"], "len": x["len"],
                                              "stress_over_calm_default": round(M[x["id"]]["stress"]["global"][0] / max(M[x["id"]]["calm"]["global"][0], EPS), 3),
                                              "calm_default": round(M[x["id"]]["calm"]["global"][0], 4)} for x in sib],
                                "score": float((st / c).min())})
    calls["induced_isoform"] = induced
    for k in calls:
        calls[k].sort(key=lambda x: (-x.get("both_high_conf", x.get("high_conf", False)), -x["score"]))

    # ---- validation vs Cao 2021 measured TE (within same-CDS groups; unambiguous only)
    val = {}
    for cell in ("HEK293T", "PC3", "muscle"):
        agree = tot = 0
        pts = []
        for g in iso_doc["groups"]:
            ms = []
            for i in g["isoforms"]:
                m = measured.lookup(i["id"])
                if not m or m["values"].get(cell) is None or cell in m.get("ambiguous_cells", []) or m["values"][cell] <= 0:
                    continue
                ms.append((i, m["values"][cell]))
            for (a, ta), (b, tb) in itertools.combinations(ms, 2):
                pa, pb = M[a["id"]]["calm"]["global"][0], M[b["id"]]["calm"]["global"][0]
                if pa <= EPS or pb <= EPS or ta == tb:
                    continue
                lp, lm = math.log2(pa / pb), math.log2(ta / tb)
                if abs(lp) < 0.25:
                    continue
                tot += 1
                agree += (lp > 0) == (lm > 0)
                pts.append([round(lp, 3), round(lm, 3)])
        rho = None
        if len(pts) >= 5:
            x, y = np.array(pts).T
            rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
            rho = float(np.corrcoef(rx, ry)[0, 1])
        # binomial two-sided p for sign agreement
        from math import comb
        p = None
        if tot:
            k = max(agree, tot - agree)
            p = min(1.0, 2 * sum(comb(tot, j) for j in range(k, tot + 1)) / 2 ** tot)
        val[cell] = {"pairs": tot, "sign_agree": agree, "frac": round(agree / tot, 3) if tot else None,
                     "binom_p_two_sided": p, "spearman_log2": round(rho, 3) if rho is not None else None, "points": pts}

    # positive controls / known ISR-regulated genes present in the same-CDS set: per-isoform ranges
    CTRL = {"ATF4": "vattem2004 (fitted finding, MANE only)", "ATF5": "ISR-induced via uORFs (Zhou et al. 2008 JBC; not used in fitting)",
            "DDIT3": "induced in vivo via uORF stalling, a mechanism the model lacks (palam2011)",
            "SLC35A4": "ISR-resistant uORF gene (Andreev et al. 2015 eLife)", "IFRD1": "ISR-resistant (Andreev et al. 2015 eLife)",
            "UCP2": "ISR-resistant (Andreev et al. 2015 eLife)", "PPP1R15B": "ISR-resistant (Andreev et al. 2015 eLife)"}
    controls = []
    for g in iso_doc["groups"]:
        if g["gene"] not in CTRL:
            continue
        rows = []
        for i in g["isoforms"]:
            r = M[i["id"]]
            c, st = arr(r, "calm"), arr(r, "stress")
            ratio = st / np.maximum(c, EPS)
            rows.append({"name": i["name"], "id": i["id"], "len": i["len"], "tsl": i["tsl"], "mane": "MANE_Select" in i["tags"],
                         "calm_range": [round(float(c.min()), 3), round(float(c.max()), 3)],
                         "stress_over_calm_range": [round(float(ratio.min()), 3), round(float(ratio.max()), 3)],
                         "stress_above_floor_everywhere": bool((st >= QS).all())})
        controls.append({"gene": g["gene"], "note": CTRL[g["gene"]], "isoforms": rows})
    ratios = np.array(all_ratio_default)
    doc = {
        "meta": {
            "what": "Model-based predictions for same-protein isoforms that differ only in the 5' leader (GENCODE v47).",
            "model": "app/src/analytic.js expected() (deterministic expected value of the sim.js scanning model), perMin output",
            "stress": "tc = 0.15 (findings.js STRESS_TC); calm tc = 1",
            "params_default": {k: v for k, v in feas["default"].items() if k not in ("passed", "total", "failed")},
            "robustness": f"every call holds at all {npts} parameter points (default + every grid point where all "
                          f"{feas['default']['total']} mechanism findings pass; grid scanSpeed x tcRate x loadTCexp) "
                          "and for both fold methods (global MFE hairpins; local MFE max_bp_span=100 hairpins)",
            "points": mo["points"],
            "feasible_grid": {"n_grid": len(feas["grid"]), "n_feasible": len(feas["feasible"]),
                              "lambda_range": [min(q["lambda"] for q in feas["feasible"]), max(q["lambda"] for q in feas["feasible"])],
                              "loadTCexp_range": [min(q["loadTCexp"] for q in feas["feasible"]), max(q["loadTCexp"] for q in feas["feasible"])]},
            "input": "GENCODE v47 pc_transcripts.fa (UTR5/CDS ranges) + gencode.v47.basic.annotation.gtf (tags, TSL)",
            "filters": "protein_coding, tag basic, no mRNA_start_NF/cds_start_NF/readthrough, not PAR_Y, UTR5 from nt 1, CDS starts AUG, "
                       f"UTR length {iso_doc['filters']['minlen']}..{iso_doc['filters']['maxlen']} nt; grouped by gene + identical CDS nucleotide "
                       "sequence; identical UTRs collapsed; groups with >= 2 distinct UTRs",
            "output_floor": "the isoform a call is about must reach the 25th percentile of all isoforms' output (calm for leader_effect; stress for stress_flip, isr_escape, induced_isoform) at every point, so ratios of near-zero outputs are not called",
            "high_conf": "both isoforms MANE_Select or transcript_support_level 1-2 (field both_high_conf); lists sorted high-conf first",
            "stats": {**iso_doc["stats"], "pairs": n_pairs, "pairs_both_high_conf": n_pairs_hc,
                      "median_stress_over_calm_default": round(float(np.median(ratios)), 3),
                      "frac_isoforms_ratio_ge1_default": round(float((ratios >= 1).mean()), 4)},
            "counts": {k: len(v) for k, v in calls.items()},
            "counts_high_conf": {k: sum(x.get("both_high_conf", x.get("high_conf", False)) for x in v) for k, v in calls.items()},
            "mechanism_note": "mechanism = which ablation (all upstream AUGs disabled in both isoforms; all hairpins melted) makes the call fail at some point",
            "mechanism_counts": {k: {m: sum(x.get("mechanism") == m for x in v) for m in sorted({x.get("mechanism") for x in v if x.get("mechanism")})} for k, v in calls.items()},
            "caveats": [
                "Predictions from a scanning model with assumed (mostly unsourced) parameters; robustness across the feasible region is not validation.",
                "GENCODE 5' ends are approximate TSSs; real leaders vary by a few to many nt.",
                "The model has no IRES, RBP, m6A, uORF peptide-stalling (e.g. CHOP) or non-AUG uORF starts; calls depending on those are out of scope.",
                "Same-CDS isoforms share their whole CDS, so Ribo-seq CDS reads cannot be assigned between them: the Cao 2021 TE validation uses only values not flagged ambiguous and is weak evidence either way.",
            ],
        },
        "validation_cao2021": val,
        "controls": controls,
        "calls": {k: v[:200] for k, v in calls.items()},
    }
    out = os.path.join(ROOT, "data/discovery.json")
    json.dump(doc, open(out, "w"), indent=1)
    print(json.dumps(doc["meta"]["stats"]), json.dumps(doc["meta"]["counts"]), json.dumps(doc["meta"]["counts_high_conf"]))
    print(json.dumps(doc["meta"]["mechanism_counts"]))
    for x in calls["induced_isoform"][:40]:
        print("  IND", x["gene"], x["A"]["name"], x["A"]["len"], "hc", x["high_conf"], "ratio", x["stress_over_calm_range"], "calm", x["calm_default"],
              [(s_["name"], s_["stress_over_calm_default"], s_["calm_default"]) for s_ in x["siblings"]][:4])
    for c, v in val.items():
        print(c, {k: v[k] for k in v if k != "points"})
    for k, v in calls.items():
        print("==", k)
        if k == "induced_isoform":
            continue
        for x in v[:15]:
            d = x["default"]
            print(f"  {x['gene']:10s} {x['A']['name']:14s} vs {x['B']['name']:14s} hc={x['both_high_conf']} "
                  f"calm {d['calm_A']:.3g}/{d['calm_B']:.3g} stress {d['stress_A']:.3g}/{d['stress_B']:.3g} "
                  f"mech={x.get('mechanism')} lenA={x['A']['len']} lenB={x['B']['len']} uORF/ovl A={x['A']['n_uORF']}/{x['A']['n_overlapping']} B={x['B']['n_uORF']}/{x['B']['n_overlapping']}")


if __name__ == "__main__":
    main()
