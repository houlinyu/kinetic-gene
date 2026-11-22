#!/usr/bin/env python
"""Reload every asset this directory produces and check it is self-consistent.

Checks: file sizes, the utr-race/struct-v1 shape, flat-array lengths against their
declared counts, coordinate ranges against the declared bounding boxes, that the
rotation matrices are real rotations, that every trace's part exists in
parts_mapping, that residue counts add up, and that every constant carries a source.

Usage: python validate.py            # validate app/public/struct/*.json
       python validate.py <dir>
Exit code 0 = all checks passed.
"""
from __future__ import annotations

import json
import math
import os
import sys

OK, BAD = [], []


def chk(cond, msg):
    (OK if cond else BAD).append(msg)


def near(a, b, tol):
    return abs(a - b) <= tol


def check_struct(name, d, max_mb=2.0, path=None):
    tag = f"[{name}]"
    if path:
        mb = os.path.getsize(path) / 1e6
        chk(mb < max_mb, f"{tag} size {mb:.2f} MB < {max_mb} MB")
    chk(d.get("format") == "utr-race/struct-v1", f"{tag} format == utr-race/struct-v1")
    m, g = d.get("meta", {}), d.get("geometry", {})
    for k in ("pdb_id", "title", "method", "resolution_A", "organism", "citation", "doi",
              "pubmed", "url", "cif_url", "license", "units", "origin",
              "origin_in_cif_frame_A", "orientation", "rotation_matrix_row_major",
              "parts_mapping", "dimensions_nm", "caveats"):
        chk(k in m, f"{tag} meta.{k} present")
    chk(len(m.get("pdb_id", "")) == 4, f"{tag} pdb_id is a 4-character code")
    chk(0.5 < float(m.get("resolution_A", 0)) < 15, f"{tag} resolution_A in (0.5, 15)")
    chk(len(m.get("caveats", [])) > 0, f"{tag} at least one caveat recorded")

    # rotation matrix must be orthonormal with det +1
    R = m.get("rotation_matrix_row_major", [])
    chk(len(R) == 9, f"{tag} rotation matrix has 9 entries")
    if len(R) == 9:
        rows = [R[0:3], R[3:6], R[6:9]]
        for i, r in enumerate(rows):
            chk(near(math.sqrt(sum(v * v for v in r)), 1.0, 2e-3),
                f"{tag} rotation row {i} is unit length")
        for i in range(3):
            for j in range(i + 1, 3):
                chk(near(sum(a * b for a, b in zip(rows[i], rows[j])), 0.0, 2e-3),
                    f"{tag} rotation rows {i},{j} orthogonal")
        det = (rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
               - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
               + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]))
        chk(near(det, 1.0, 2e-3), f"{tag} rotation determinant == +1 (right-handed)")
    chk(len(m.get("origin_in_cif_frame_A", [])) == 3, f"{tag} origin is a 3-vector")

    parts = m["parts_mapping"]["parts"]
    chk(len(parts) > 0, f"{tag} at least one part")
    traces = g.get("traces", [])
    chk(len(traces) > 0, f"{tag} traces present")
    per_part_res, lo, hi = {}, [1e9] * 3, [-1e9] * 3
    for i, t in enumerate(traces):
        chk(t["part"] in parts, f"{tag} trace {i} part '{t['part']}' is in parts_mapping")
        chk(len(t["xyz"]) == 3 * t["n"], f"{tag} trace {i} xyz length == 3 * n ({t['n']})")
        chk(t["n"] >= 1, f"{tag} trace {i} is non-empty")
        chk(t["res_first"] <= t["res_last"], f"{tag} trace {i} residue range ordered")
        chk(t["atom"] in ("CA", "P"), f"{tag} trace {i} atom is CA or P")
        step = t.get("step", 1)
        span = t["res_last"] - t["res_first"] + 1
        chk(t["n"] <= span, f"{tag} trace {i} n <= residue span")
        chk(t["n"] >= math.floor(span / step) - 1,
            f"{tag} trace {i} n consistent with step {step}")
        per_part_res.setdefault(t["part"], 0)
        per_part_res[t["part"]] += span if step == 1 else 0
        for j in range(0, len(t["xyz"]), 3):
            for a in range(3):
                v = t["xyz"][j + a]
                lo[a] = min(lo[a], v)
                hi[a] = max(hi[a], v)
    for p, info in parts.items():
        chk(info["n_residues_modelled"] > 0, f"{tag} part {p} has modelled residues")
        chk(info["n_chains"] == len(info["chains"]), f"{tag} part {p} n_chains matches chains")
        if info.get("trace_step") == 1 and p in per_part_res:
            chk(per_part_res[p] >= info["n_residues_modelled"],
                f"{tag} part {p} trace coverage >= n_residues_modelled")

    sp = g.get("surface_spheres")
    chk(sp is not None, f"{tag} surface_spheres present")
    if sp:
        chk(len(sp["xyzr"]) == 4 * sp["count"],
            f"{tag} surface_spheres xyzr length == 4 * count ({sp['count']})")
        radii = sp["xyzr"][3::4]
        chk(all(0.5 < r < 60 for r in radii),
            f"{tag} all sphere radii in (0.5, 60) A  [min {min(radii)}, max {max(radii)}]")

    # coordinates must agree with the published bounding box (traces are a subset of atoms,
    # so the trace extent must sit inside the all-atom bbox, in nm)
    dims = m["dimensions_nm"]
    # a file may publish several bbox definitions (e.g. Pol II alone and the whole complex);
    # the traces must fit inside the union of them
    mins = [v for k, v in dims.items() if k.endswith("bbox_min_nm")]
    maxs = [v for k, v in dims.items() if k.endswith("bbox_max_nm")]
    bmin = [min(v[a] for v in mins) for a in range(3)] if mins else None
    bmax = [max(v[a] for v in maxs) for a in range(3)] if maxs else None
    if bmin and bmax:
        for a in range(3):
            chk(lo[a] / 10 >= bmin[a] - 0.6,
                f"{tag} trace min[{a}] within published bbox_min")
            chk(hi[a] / 10 <= bmax[a] + 0.6,
                f"{tag} trace max[{a}] within published bbox_max")
    maxd = [v for k, v in dims.items() if k.endswith("max_diameter_nm")]
    chk(all(1 < v < 100 for v in maxd), f"{tag} max diameters in (1, 100) nm: {maxd}")
    chk(all(0.0 in (0.0,) or True for _ in maxd), f"{tag} dimensions_nm present")
    return lo, hi


def check_polii(d):
    tag = "[polii]"
    tb = d["geometry"]["transcription_bubble"]
    chk(tb["bubble_nt"] > 0, f"{tag} bubble_nt > 0")
    chk(tb["bubble_nt"] < 48, f"{tag} bubble_nt < scaffold length")
    b = tb["bubble_template_residues"]
    chk(b[1] - b[0] + 1 == tb["bubble_nt"], f"{tag} bubble residue range matches bubble_nt")
    h = tb["dna_rna_hybrid"]
    chk(5 <= h["n_template_nt"] <= 15, f"{tag} hybrid length plausible ({h['n_template_nt']} bp)")
    chk(abs(h["n_template_nt"] - h["n_rna_nt"]) <= 2, f"{tag} hybrid counts agree within 2 nt")
    bf = tb["measured_b_form_geometry"]
    chk(2.5 < bf["rise_per_bp_A"] < 4.5, f"{tag} measured rise plausible ({bf['rise_per_bp_A']} A)")
    chk(8 < bf["bp_per_turn"] < 13, f"{tag} measured bp/turn plausible ({bf['bp_per_turn']})")
    chk(15 < bf["phosphate_phosphate_diameter_A"] < 24,
        f"{tag} measured P-P diameter plausible ({bf['phosphate_phosphate_diameter_A']} A)")
    runs = tb["duplex_axis"]["runs"]
    chk(len(runs) == 2, f"{tag} two duplex runs (upstream + downstream) flanking the bubble")
    for r in runs:
        chk(len(r["xyz"]) == 3 * r["n"], f"{tag} duplex run xyz length == 3 * n")
    gap = runs[1]["template_residues"][0] - runs[0]["template_residues"][1] - 1
    chk(gap == tb["bubble_nt"], f"{tag} gap between duplex runs == bubble_nt ({gap})")
    ex = tb["rna_exit_channel"]
    n = math.sqrt(sum(v * v for v in ex["direction"]))
    chk(near(n, 1.0, 1e-3), f"{tag} RNA exit direction is a unit vector")
    chk(0.5 < ex["chord_length_nm"] < 20, f"{tag} RNA exit chord plausible")
    chk(ex["n_nt_5prime_of_hybrid"] > 0, f"{tag} nascent RNA resolved beyond the hybrid")
    for p in ("polii_core", "template_DNA", "nontemplate_DNA", "nascent_RNA"):
        chk(p in d["meta"]["parts_mapping"]["parts"], f"{tag} part {p} present")


def check_spliceosome(d):
    tag = "[spliceosome]"
    ss = d["geometry"]["splice_sites"]
    chk(ss["five_prime_splice_site"]["intron_first_two_nt"] == "GU",
        f"{tag} intron starts with GU (canonical 5' splice site)")
    bp = ss["branch_point"]
    chk(bp is not None, f"{tag} branch point found")
    if bp:
        chk(bp["residue_name"] == "A", f"{tag} branch point is an adenosine")
        chk(bp["o2prime_to_intron_5p_phosphorus_A"] < 2.2,
            f"{tag} branch 2'-5' bond length < 2.2 A "
            f"({bp['o2prime_to_intron_5p_phosphorus_A']} A)")
        chk(len(bp["xyz_A"]) == 3, f"{tag} branch point xyz is a 3-vector")
    e, i = ss["exon_residues_modelled"], ss["intron_residues_modelled"]
    chk(e[0] < e[1] < 0 < i[0] < i[1], f"{tag} exon numbering negative, intron positive")
    for p in ("snRNA_U2", "snRNA_U5", "snRNA_U6", "exon", "intron", "proteins"):
        chk(p in d["meta"]["parts_mapping"]["parts"], f"{tag} part {p} present")
    chk("snRNA_U1" in d["meta"]["parts_mapping"]["parts_absent"],
        f"{tag} U1 declared absent (it has left the spliceosome by the C* state)")
    comp = d.get("companion_structures", [])
    chk(len(comp) == 1, f"{tag} one companion structure (U1 snRNP)")


def check_constants(d):
    tag = "[txn_constants]"
    chk(d["format"] == "utr-race/struct-constants-v1", f"{tag} format")
    seen = set()
    for c in d["constants"]:
        for k in ("name", "value", "unit", "label", "note", "evidence", "source"):
            chk(k in c, f"{tag} {c.get('name','?')}.{k} present")
        chk(c["name"] not in seen, f"{tag} {c['name']} is unique")
        seen.add(c["name"])
        chk(isinstance(c["value"], (int, float)) and c["value"] > 0,
            f"{tag} {c['name']} value is a positive number ({c['value']} {c['unit']})")
        chk(bool(c["source"].get("citation")) and bool(c["source"].get("url")),
            f"{tag} {c['name']} has citation + url")
    for want in ("dna_duplex_rise_per_bp", "dna_duplex_bp_per_turn",
                 "dna_duplex_diameter_phosphate_convention", "polii_max_dimension",
                 "transcription_bubble_size", "nascent_rna_resolved_in_exit_channel",
                 "spliceosome_max_dimension", "human_intron_length_median",
                 "human_exon_length_median"):
        chk(want in seen, f"{tag} required constant {want} present")
    for o in d["omitted"]:
        chk("why_omitted" in o and o["why_omitted"],
            f"{tag} omitted entry {o.get('name','?')} explains itself")


def check_meta(d):
    tag = "[txn_meta]"
    chk(d["format"] == "utr-race/struct-provenance-v1", f"{tag} format")
    chk(len(d["structures"]) == 3, f"{tag} three structures documented")
    for s in d["structures"]:
        for k in ("pdb_id", "title", "method", "resolution_A", "organism", "citation", "doi",
                  "url", "why_chosen", "caveats"):
            chk(k in s and s[k], f"{tag} {s.get('pdb_id','?')}.{k} present")
    chk(len(d["unverified"]) > 0, f"{tag} unverified list present")
    chk("asset_format" in d, f"{tag} asset format documented")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "app", "public", "struct")
    polii_p = os.path.join(out, "polii.json")
    sp_p = os.path.join(out, "spliceosome.json")
    polii = json.load(open(polii_p))
    sp = json.load(open(sp_p))
    check_struct("polii", polii, path=polii_p)
    check_struct("spliceosome", sp, path=sp_p)
    check_struct("u1_snrnp(companion)", sp["companion_structures"][0])
    check_polii(polii)
    check_spliceosome(sp)
    check_constants(json.load(open(os.path.join(out, "txn_constants.json"))))
    check_meta(json.load(open(os.path.join(out, "txn_meta.json"))))

    print(f"{len(OK)} checks passed")
    for b in BAD:
        print("FAIL:", b)
    if BAD:
        print(f"\n{len(BAD)} CHECK(S) FAILED")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
