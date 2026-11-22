#!/usr/bin/env python
"""Reload the published struct assets and sanity-check them.

Run:  ~/miniconda3/envs/utrrace/bin/python -I pipeline/structures/validate.py
Exit code 0 = all checks pass.  Prints one line per check; FAIL lines are errors.
"""
from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "app/public/struct")

REQUIRED_PARTS = {"18S_rRNA", "40S_proteins", "mRNA", "tRNAi"}
fails: list[str] = []


def ck(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


def check_struct(path: str, primary: bool) -> None:
    name = os.path.basename(path)
    print(f"\n--- {name} ({os.path.getsize(path)/1e6:.2f} MB) ---")
    d = json.load(open(path))
    m, g = d["meta"], d["geometry"]

    ck(d.get("format") == "utr-race/struct-v1", f"{name}: format tag")
    for k in ("pdb_id", "title", "method", "organism", "citation", "doi", "pubmed",
              "url", "license", "parts_mapping", "dimensions_nm", "n_mrna_nt_resolved"):
        ck(k in m and m[k] not in (None, "", {}), f"{name}: meta.{k} present")
    ck(len(m["pdb_id"]) == 4, f"{name}: pdb_id is a 4-character code ({m['pdb_id']})")
    ck(m["license"] == "PDB data are freely available", f"{name}: license string")
    ck(isinstance(m["resolution_A"], (int, float)) and 1.0 < m["resolution_A"] < 10.0,
       f"{name}: resolution {m['resolution_A']} A in a plausible cryo-EM range")
    R = m["rotation_matrix_row_major"]
    ck(len(R) == 9, f"{name}: rotation matrix has 9 elements")
    det = (R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6])
           + R[2] * (R[3] * R[7] - R[4] * R[6]))
    ck(abs(det - 1.0) < 1e-3, f"{name}: rotation is a proper rotation (det={det:.6f})")
    for i in range(3):
        nrm = math.sqrt(sum(R[3 * i + j] ** 2 for j in range(3)))
        ck(abs(nrm - 1.0) < 1e-3, f"{name}: rotation row {i} is unit length ({nrm:.6f})")

    # traces
    tr = g["traces"]
    ck(len(tr) > 10, f"{name}: {len(tr)} trace polylines")
    npts = 0
    bad_len = bad_rng = bad_short = 0
    for t in tr:
        if len(t["xyz"]) != 3 * t["n"]:
            bad_len += 1
        if t["n"] < 2:
            bad_short += 1
        npts += t["n"]
        for v in t["xyz"]:
            if not (-400.0 < v < 400.0):
                bad_rng += 1
                break
    ck(bad_len == 0, f"{name}: every trace xyz length == 3*n ({bad_len} bad)")
    ck(bad_short == 0, f"{name}: every trace has >= 2 points ({bad_short} bad)")
    ck(bad_rng == 0, f"{name}: all coordinates within +/-400 A ({bad_rng} traces out of range)")
    print(f"      {npts} trace points total")

    parts = {t["part"] for t in tr}
    ck(REQUIRED_PARTS <= parts, f"{name}: required parts present "
                                f"(missing: {sorted(REQUIRED_PARTS - parts)})")
    declared = set(m["parts_mapping"]["parts"])
    ck(parts <= declared, f"{name}: every trace part is declared in parts_mapping")

    # consecutive points inside a polyline must be close (it is a backbone trace)
    worst = 0.0
    for t in tr:
        x = t["xyz"]
        for i in range(0, len(x) - 3, 3):
            dd = math.dist(x[i:i + 3], x[i + 3:i + 6])
            worst = max(worst, dd)
    ck(worst < 45.0, f"{name}: largest step inside a polyline {worst:.1f} A (< 45)")

    # surface spheres
    s = g["surface_spheres"]
    ck(len(s["xyzr"]) == 4 * s["count"], f"{name}: sphere array length == 4*count")
    radii = s["xyzr"][3::4]
    ck(all(2.0 < r < 40.0 for r in radii),
       f"{name}: sphere radii in 2-40 A (min {min(radii):.1f}, max {max(radii):.1f})")
    ck(s["count"] > 100, f"{name}: {s['count']} surface spheres")

    # mRNA channel
    c = g["mrna_channel"]
    ck(c is not None, f"{name}: mrna_channel present")
    n = c["n_nt_resolved"]
    ck(n == m["n_mrna_nt_resolved"], f"{name}: n_mrna_nt_resolved agrees ({n})")
    ck(len(c["res_numbers"]) == n, f"{name}: res_numbers length == n ({n})")
    ck(len(c["path_xyz_5to3_A"]) == 3 * n, f"{name}: mRNA path length == 3*n")
    ck(len(c["sequence_resolved_5to3"]) == n, f"{name}: mRNA sequence length == n")
    ck(set(c["sequence_resolved_5to3"]) <= set("ACGU"),
       f"{name}: mRNA sequence is RNA letters only")
    ck(c["res_numbers"] == sorted(c["res_numbers"]),
       f"{name}: mRNA residue numbers increase 5'->3'")
    ck(3.0 < c["mean_PP_spacing_A"] < 9.0,
       f"{name}: mean P-P spacing {c['mean_PP_spacing_A']} A is plausible for ssRNA")
    for end in ("end_5prime_exit_channel", "end_3prime_entry_channel"):
        ck(len(c[end]["xyz_A"]) == 3, f"{name}: {end} has an xyz")
    # the mRNA path should run roughly along +x by construction
    xs = c["path_xyz_5to3_A"][0::3]
    ck(xs[-1] > xs[0], f"{name}: mRNA 3' end is at larger x than the 5' end "
                       f"({xs[0]:.0f} -> {xs[-1]:.0f} A)")

    # dimensions
    dm = m["dimensions_nm"]
    ck(10.0 < dm["40S_max_diameter_nm"] < 40.0,
       f"{name}: 40S max diameter {dm['40S_max_diameter_nm']} nm is plausible")
    ck(all(5.0 < v < 40.0 for v in dm["40S_bbox_nm"]),
       f"{name}: 40S bbox {dm['40S_bbox_nm']} nm")
    ck(dm["40S_max_diameter_nm"] >= max(dm["40S_principal_axis_extents_nm"]) - 0.5,
       f"{name}: max diameter >= largest principal extent")

    # flat point-cloud compat view consumed by app/src/struct3d.js
    ck(d.get("pdb_id") == m["pdb_id"], f"{name}: top-level pdb_id alias")
    ck(d.get("citation") == m["citation"], f"{name}: top-level citation alias")
    ck(d.get("dimensions_nm") == m["dimensions_nm"], f"{name}: top-level dimensions_nm alias")
    cp = d["parts"]
    ck(isinstance(cp, list) and len(cp) >= 4, f"{name}: {len(cp)} compat parts")
    tot = 0
    for p in cp:
        ck(p["n"] == len(p["points"]), f"{name}: compat part {p['name']} n == len(points)")
        want = 4 if p["kind"] == "spheres" else 3
        ck(all(len(q) == want for q in p["points"]),
           f"{name}: compat part {p['name']} points are length-{want}")
        ck(all(all(isinstance(v, (int, float)) and -400 < v < 400 for v in q[:3])
               for q in p["points"]),
           f"{name}: compat part {p['name']} coordinates finite and in range")
        tot += p["n"]
    ck(tot >= 500, f"{name}: {tot} compat point-cloud points")
    ck({p["name"] for p in cp} <= declared,
       f"{name}: every compat part name is declared in parts_mapping")
    ck(len(d["channel_path"]) == n, f"{name}: channel_path length == n mRNA nt")
    ck(all(len(q) == 3 for q in d["channel_path"]), f"{name}: channel_path points are xyz")
    flat = [v for q in d["channel_path"] for v in q]
    ck(flat == c["path_xyz_5to3_A"], f"{name}: channel_path agrees with geometry.mrna_channel")

    if primary:
        ck(m["pdb_id"] == "8PJ1", "primary asset is PDB 8PJ1")
        ck(n == 47, f"primary asset has 47 modelled mRNA nt (got {n})")
        ck(c["entry_edge_offset_from_p_site_nt"]["value"] == 16,
           "primary asset: entry edge is +16 nt from the geometric P site "
           "(published toeprint leading edge +16..+18)")
        ck({"eIF1", "eIF1A", "eIF2", "eIF3"} <= parts,
           "primary asset carries eIF1, eIF1A, eIF2 and eIF3")


def check_constants() -> None:
    p = os.path.join(OUT, "constants.json")
    print(f"\n--- constants.json ({os.path.getsize(p)/1e3:.1f} kB) ---")
    d = json.load(open(p))
    cs = d["constants"]
    ck(len(cs) >= 8, f"{len(cs)} constants")
    names = set()
    for c in cs:
        for k in ("name", "value", "unit", "source"):
            ck(k in c, f"constant {c.get('name','?')}: has {k}")
        ck("citation" in c["source"] and "url" in c["source"],
           f"constant {c['name']}: source has citation + url")
        ck(c["source"]["url"].startswith("https://"),
           f"constant {c['name']}: url is https")
        ck(c["name"] not in names, f"constant {c['name']}: name is unique")
        names.add(c["name"])
        v = c["value"]
        ck(isinstance(v, (int, float)) or (isinstance(v, list) and len(v) == 2),
           f"constant {c['name']}: value is a number or a [lo, hi] range")
    # cross-check the measured constants against the primary asset
    a = json.load(open(os.path.join(OUT, "43s.json")))
    by = {c["name"]: c for c in cs}
    ck(by["mrna_nt_modelled_in_48S"]["value"] == a["meta"]["n_mrna_nt_resolved"],
       "constants: mrna_nt_modelled_in_48S matches 43s.json")
    ck(abs(by["s40_max_dimension"]["value"]
           - a["meta"]["dimensions_nm"]["40S_max_diameter_nm"]) < 0.1,
       "constants: s40_max_dimension matches 43s.json")
    ck(abs(by["mrna_pp_spacing_in_40S"]["value"]
           - a["geometry"]["mrna_channel"]["mean_PP_spacing_A"]) < 0.01,
       "constants: mrna_pp_spacing matches 43s.json")

    # app/src/struct3d.js picks its six scale slots by regex over name+label, first
    # match wins. Guard the ordering so a later edit cannot silently mis-wire the UI.
    import re
    pats = {"rise": r"rise",
            "bpturn": r"(bp|base.?pairs?).*(per|_)turn|per_turn",
            "diam": r"helix.*diam|diam.*helix|a.?form.*diam",
            "ribo": r"40s.*(max|dim|size|diam)|(max|dim).*40s",
            "chan": r"channel", "foot": r"footprint"}
    want = {"rise": "a_form_rna_rise_per_bp", "bpturn": "a_form_rna_bp_per_turn",
            "diam": "a_form_rna_helix_diameter", "ribo": "s40_max_dimension",
            "chan": "mrna_nt_in_40S_channel", "foot": "riboseq_footprint_80S_yeast"}
    for slot, pat in pats.items():
        rx = re.compile(pat, re.I)
        first = next((c["name"] for c in cs
                      if isinstance(c["value"], (int, float))
                      and rx.search(f"{c['name']} {c.get('label','')}")), None)
        ck(first == want[slot],
           f"constants: struct3d.js slot '{slot}' resolves to {want[slot]} (got {first})")


def check_readme() -> None:
    p = os.path.join(OUT, "README.json")
    print(f"\n--- README.json ({os.path.getsize(p)/1e3:.1f} kB) ---")
    d = json.load(open(p))
    ids = set()
    for s in d["structures"]:
        for k in ("pdb_id", "title", "method", "resolution_A", "organism",
                  "citation", "doi", "url", "asset", "why_chosen"):
            ck(k in s and s[k], f"README {s.get('pdb_id','?')}: has {k}")
        ck(os.path.exists(os.path.join(OUT, s["asset"])),
           f"README {s['pdb_id']}: asset {s['asset']} exists")
        a = json.load(open(os.path.join(OUT, s["asset"])))
        ck(a["meta"]["pdb_id"] == s["pdb_id"],
           f"README {s['pdb_id']}: asset meta.pdb_id agrees")
        ck(abs(a["meta"]["resolution_A"] - s["resolution_A"]) < 0.05,
           f"README {s['pdb_id']}: asset resolution agrees")
        ids.add(s["pdb_id"])
    ck("license" in d and d["license"].startswith("PDB data are freely available"),
       "README: license statement")
    ck("8PJ1" in ids, "README documents the primary structure 8PJ1")


def main() -> int:
    check_struct(os.path.join(OUT, "43s.json"), primary=True)
    sec = os.path.join(OUT, "48s_eif4f.json")
    if os.path.exists(sec):
        check_struct(sec, primary=False)
    check_constants()
    check_readme()
    print()
    if fails:
        print(f"{len(fails)} CHECK(S) FAILED:")
        for f in fails:
            print("  -", f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
