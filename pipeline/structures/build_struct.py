#!/usr/bin/env python
"""Derive a lightweight web asset from a real cryo-EM mmCIF of a human 48S
translation pre-initiation complex.

Input : pipeline/structures/dl/<id>.cif   (downloaded from RCSB, untrusted dir)
Output: app/public/struct/43s.json

Nothing here invents coordinates or numbers: every value is either read from
the deposited model or computed from its coordinates. Provenance (PDB id,
citation, DOI/PubMed) is read out of the mmCIF itself where possible.

Run:  ~/miniconda3/envs/utrrace/bin/python -I pipeline/structures/build_struct.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import date

import gemmi
import numpy as np

# ---------------------------------------------------------------- part mapping
# Matched against the mmCIF _entity.pdbx_description, in order; first hit wins.
# (substring, case-insensitive)  ->  part name
PART_RULES: list[tuple[str, str]] = [
    ("18s rrna", "18S_rRNA"),
    ("initiator met-trna", "tRNAi"),
    ("trnaimet", "tRNAi"),
    ("trnai", "tRNAi"),
    ("mrna", "mRNA"),
    ("initiation factor 1a", "eIF1A"),
    ("initiation factor 1,", "eIF1"),
    ("initiation factor 2 subunit", "eIF2"),
    ("initiation factor 3 subunit", "eIF3"),
    ("initiation factor 4a", "eIF4A"),
    ("initiation factor 4 gamma", "eIF4G"),
    ("initiation factor 4e", "eIF4E"),
    ("initiation factor 5b", "eIF5B"),
    ("initiation factor 5", "eIF5"),
    ("40s ribosomal protein", "40S_proteins"),
    ("small ribosomal subunit protein", "40S_proteins"),  # newer naming in some entries
    ("receptor of activated protein c kinase 1", "40S_proteins"),  # RACK1, 40S-bound
    ("60s ribosomal protein l41", "40S_proteins"),  # eL41, modelled with the 40S here
    ("ubiquitin", "40S_proteins"),  # the ubiquitin moiety of eS31/RPS27a
]
# exact-description overrides (eIF1 has a description without a trailing comma)
EXACT_RULES = {
    "eukaryotic translation initiation factor 1": "eIF1",
}

# the "40S" for dimension/centring purposes
CORE_40S = ("18S_rRNA", "40S_proteins")

# decimation step (residues) per part; 1 = keep every residue
STEP = {
    "18S_rRNA": 1,
    "mRNA": 1,
    "tRNAi": 1,
    "eIF1": 1,
    "eIF1A": 1,
    "40S_proteins": 2,
    "eIF2": 2,
    "eIF3": 3,
    "eIF4A": 2,
    "eIF4G": 2,
    "eIF4E": 1,
    "eIF5": 1,
    "eIF5B": 2,
    "other": 3,
}
# backbone atom used for the trace
TRACE_ATOM = {"rna": "P", "protein": "CA"}
# split a polyline when consecutive kept trace atoms are further apart than this
GAP_A = {"rna": 14.0, "protein": 12.0}
SURFACE_WINDOW = 10  # residues per surface sphere
VDW_PAD_A = 3.0  # padding added to the RMS atom spread of a surface window


# Some mmCIF citation blocks are stale (deposited before the paper was paginated).
# These replacements were each checked against PubMed esummary and Crossref on
# 2026-11-22; the key is the PDB id and `pubmed` is the record they were read from.
CITATION_FIXUP = {
    "8PJ1": {"journal": "Nat Struct Mol Biol", "volume": "32", "page_first": "62",
             "page_last": "72", "year": "2025", "pubmed": "39289545",
             "verified": "PubMed esummary 39289545 and Crossref doi "
                         "10.1038/s41594-024-01378-4 (2026-11-22); the mmCIF citation "
                         "block still carries year 2024 with no volume or pages",
             "correction": "An Author Correction was published: Nat Struct Mol Biol "
                           "33:724 (2026), doi:10.1038/s41594-026-01777-9. Not read "
                           "here; it does not affect the deposited coordinates."},
    "8OZ0": {"journal": "Nat Struct Mol Biol", "volume": "31", "page_first": "455",
             "page_last": "464", "year": "2024", "pubmed": "38287194",
             "verified": "PubMed esummary 38287194 (2026-11-22)"},
}


def classify(desc: str) -> str:
    d = desc.strip().strip("'\"").lower()
    if d in EXACT_RULES:
        return EXACT_RULES[d]
    for sub, part in PART_RULES:
        if sub in d:
            return part
    # "initiation factor 1A, X-chromosomal" is caught above; plain eIF1 may lack a comma
    if d.startswith("eukaryotic translation initiation factor 1"):
        return "eIF1A" if "1a" in d else "eIF1"
    return "other"


def read_entity_desc(path: str) -> tuple[dict[str, str], dict[str, str], dict]:
    """entity_id -> description, auth_chain -> entity_id, plus header info."""
    blk = gemmi.cif.read(path).sole_block()
    ent = {}
    for r in blk.find("_entity.", ["id", "pdbx_description"]):
        ent[r.str(0)] = r.str(1)
    ch2ent = {}
    for r in blk.find("_atom_site.", ["auth_asym_id", "label_entity_id"]):
        ch2ent.setdefault(r.str(0), r.str(1))

    def one(tag, default=None):
        v = blk.find_value(tag)
        if v is None:
            return default
        v = gemmi.cif.as_string(v).strip()
        return None if v in ("?", ".", "") else v

    hdr = {
        "pdb_id": (one("_entry.id") or "").upper(),
        "title": one("_struct.title"),
        "method": one("_exptl.method"),
        "resolution_A": one("_em_3d_reconstruction.resolution"),
        "resolution_method": one("_em_3d_reconstruction.resolution_method"),
        "deposition_date": one("_pdbx_database_status.recvd_initial_deposition_date"),
        "citation": {
            "title": one("_citation.title"),
            "journal": one("_citation.journal_abbrev"),
            "volume": one("_citation.journal_volume"),
            "page_first": one("_citation.page_first"),
            "year": one("_citation.year"),
            "doi": one("_citation.pdbx_database_id_DOI"),
            "pubmed": one("_citation.pdbx_database_id_PubMed"),
        },
    }
    authors = [r.str(0) for r in blk.find("_citation_author.", ["name"])]
    hdr["citation"]["authors"] = authors
    return ent, ch2ent, hdr


def chain_kind(ch) -> str:
    pt = str(ch.get_polymer().check_polymer_type())
    return "rna" if "Rna" in pt or "Dna" in pt else "protein"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    cif = a.cif or os.path.join(root, "pipeline/structures/dl/8pj1.cif")
    out = a.out or os.path.join(root, "app/public/struct/43s.json")

    ent_desc, ch2ent, hdr = read_entity_desc(cif)
    st = gemmi.read_structure(cif)
    st.setup_entities()
    model = st[0]

    # ---------------------------------------------------- gather per-chain data
    chains = []  # dict: name, part, desc, kind, resnums, trace (N,3), all_atoms (M,3), atom_res_idx
    for ch in model:
        eid = ch2ent.get(ch.name)
        desc = (ent_desc.get(eid, "") or "").strip("'\"")
        part = classify(desc)
        kind = chain_kind(ch)
        want = TRACE_ATOM[kind]
        resnums, trace, allxyz, allres = [], [], [], []
        for ri, res in enumerate(ch):
            at = res.find_atom(want, "*")
            if at is None:  # fall back to first atom so nothing is silently dropped
                at = res[0] if len(res) else None
            if at is None:
                continue
            resnums.append(res.seqid.num)
            trace.append([at.pos.x, at.pos.y, at.pos.z])
            for atom in res:
                if atom.element == gemmi.Element("H"):
                    continue
                allxyz.append([atom.pos.x, atom.pos.y, atom.pos.z])
                allres.append(len(resnums) - 1)
        if not trace:
            continue
        chains.append(
            dict(name=ch.name, part=part, desc=desc, kind=kind,
                 resnums=np.array(resnums), trace=np.array(trace, float),
                 allxyz=np.array(allxyz, float), allres=np.array(allres))
        )

    by_part: dict[str, list] = {}
    for c in chains:
        by_part.setdefault(c["part"], []).append(c)
    missing = [p for p in ("18S_rRNA", "40S_proteins", "mRNA", "tRNAi",
                           "eIF1", "eIF1A", "eIF2", "eIF3") if p not in by_part]

    # ------------------------------------------------- origin: 40S centre of mass
    core = [c for c in chains if c["part"] in CORE_40S]
    core_atoms = np.vstack([c["allxyz"] for c in core])
    # mass-free centroid of all modelled non-H atoms of 18S rRNA + 40S proteins
    origin = core_atoms.mean(axis=0)

    # ------------------------------------------------------------- orientation
    mrna = by_part.get("mRNA", [])
    rot = np.eye(3)
    orient_note = "none (deposited frame, only translated to the 40S centroid)"
    if mrna:
        m = mrna[0]
        P = m["trace"] - origin
        # channel axis = 1st principal component of the mRNA phosphate path,
        # signed so that 5'->3' is +x
        Pc = P - P.mean(axis=0)
        u, s, vt = np.linalg.svd(Pc, full_matrices=False)
        ax = vt[0]
        if np.dot(P[-1] - P[0], ax) < 0:
            ax = -ax
        xh = ax / np.linalg.norm(ax)
        # +z towards the initiator tRNA (the intersubunit/40S-interface side)
        if by_part.get("tRNAi"):
            ref = np.vstack([c["trace"] for c in by_part["tRNAi"]]).mean(axis=0) - origin
        else:
            ref = np.array([0.0, 0.0, 1.0])
        zh = ref - np.dot(ref, xh) * xh
        zh /= np.linalg.norm(zh)
        yh = np.cross(zh, xh)
        rot = np.vstack([xh, yh, zh])  # rows = new basis in old coords; p' = rot @ p
        if np.linalg.det(rot) < 0:
            rot[1] *= -1
        orient_note = (
            "+x = 1st principal component of the mRNA phosphate path (signed 5'->3'); "
            "+z = component of (tRNAi centroid - origin) orthogonal to +x; +y = z x x. "
            "p_out = R @ (p_cif - origin)."
        )

    def xf(p: np.ndarray) -> np.ndarray:
        return (rot @ (p - origin).T).T

    for c in chains:
        c["trace_t"] = xf(c["trace"])
        c["allxyz_t"] = xf(c["allxyz"])

    # ------------------------------------------------------------------ traces
    traces = []
    for c in sorted(chains, key=lambda c: (c["part"], c["name"])):
        step = STEP.get(c["part"], 3)
        keep = np.arange(0, len(c["trace_t"]), step)
        if keep[-1] != len(c["trace_t"]) - 1:
            keep = np.append(keep, len(c["trace_t"]) - 1)
        pts = c["trace_t"][keep]
        nums = c["resnums"][keep]
        gap = GAP_A[c["kind"]] * max(1, step)
        # split into continuous polylines at chain breaks
        cuts = [0]
        for i in range(1, len(pts)):
            if np.linalg.norm(pts[i] - pts[i - 1]) > gap:
                cuts.append(i)
        cuts.append(len(pts))
        for j in range(len(cuts) - 1):
            s, e = cuts[j], cuts[j + 1]
            if e - s < 2:
                continue
            traces.append({
                "chain": c["name"], "part": c["part"], "label": c["desc"],
                "atom": TRACE_ATOM[c["kind"]], "step": step,
                "res_first": int(nums[s]), "res_last": int(nums[e - 1]),
                "n": int(e - s),
                "xyz": [round(float(v), 1) for v in pts[s:e].ravel()],
            })

    # ----------------------------------------------------------- surface blobs
    # Windows must not straddle a chain break, or a sphere would bridge a missing
    # loop and balloon. So split each chain into spatially continuous runs first.
    spheres = []
    for c in chains:
        if c["part"] not in CORE_40S:
            continue
        tr_ = c["trace_t"]
        runs, start = [], 0
        for i in range(1, len(tr_)):
            if np.linalg.norm(tr_[i] - tr_[i - 1]) > GAP_A[c["kind"]]:
                runs.append((start, i))
                start = i
        runs.append((start, len(tr_)))
        for r0, r1 in runs:
            for s in range(r0, r1, SURFACE_WINDOW):
                e = min(s + SURFACE_WINDOW, r1)
                sel = (c["allres"] >= s) & (c["allres"] < e)
                pts = c["allxyz_t"][sel]
                if len(pts) < 4:
                    continue
                cen = pts.mean(axis=0)
                r = float(np.sqrt(((pts - cen) ** 2).sum(axis=1).mean())) + VDW_PAD_A
                spheres.append((cen, r, c["part"]))
    rmax = max(r for _, r, _ in spheres)
    assert rmax < 40.0, f"surface sphere radius {rmax:.1f} A is too large"
    sph_xyzr = [round(float(v), 1) for cen, r, _ in spheres for v in (*cen, r)]

    # -------------------------------------------------------------- dimensions
    core_t = np.vstack([c["allxyz_t"] for c in core])
    bb_min, bb_max = core_t.min(axis=0), core_t.max(axis=0)
    # max diameter: exact over the decimated backbone trace of the 40S core
    core_tr = np.vstack([c["trace_t"] for c in core])
    d2 = 0.0
    pair = (0, 0)
    for i in range(0, len(core_tr), 1):
        dd = ((core_tr[i + 1:] - core_tr[i]) ** 2).sum(axis=1)
        if len(dd) and dd.max() > d2:
            d2 = float(dd.max())
            pair = (i, i + 1 + int(dd.argmax()))
    max_diam_A = float(np.sqrt(d2))

    dims = {
        "frame": "oriented frame described in meta.orientation; nm",
        "40S_definition": "18S rRNA + 40S ribosomal proteins (incl. RACK1, eL41, "
                          "the ubiquitin moiety of eS31) as modelled in this entry",
        "40S_bbox_nm": [round(float(v) / 10, 2) for v in (bb_max - bb_min)],
        "40S_bbox_min_nm": [round(float(v) / 10, 2) for v in bb_min],
        "40S_bbox_max_nm": [round(float(v) / 10, 2) for v in bb_max],
        "40S_bbox_note": "axis-aligned over all modelled non-H atoms of the 40S core",
        "40S_max_diameter_nm": round(max_diam_A / 10, 2),
        "40S_max_diameter_note": "largest distance between any two 40S backbone trace "
                                 "atoms (18S P, protein CA); all-atom extent is ~0.5 nm larger",
        "40S_radius_of_gyration_nm": round(
            float(np.sqrt(((core_t - core_t.mean(axis=0)) ** 2).sum(axis=1).mean())) / 10, 2),
    }
    # tighter, orientation-independent extents: along the 40S's own principal axes.
    # Computed over the SAME point set as the max diameter (backbone trace atoms) so
    # the two numbers are directly comparable.
    cc = core_tr - core_tr.mean(axis=0)
    _, _, pvt = np.linalg.svd(cc, full_matrices=False)
    proj = cc @ pvt.T
    dims["40S_principal_axis_extents_nm"] = [
        round(float(v) / 10, 2) for v in (proj.max(axis=0) - proj.min(axis=0))]
    dims["40S_principal_axis_note"] = ("extents along the 40S core's own principal axes, "
                                       "over the backbone trace atoms (same point set as "
                                       "40S_max_diameter_nm); orientation-independent, so "
                                       "the axis-aligned all-atom bbox above is larger")

    # whole-complex extent too (48S)
    all_t = np.vstack([c["allxyz_t"] for c in chains])
    dims["48S_bbox_nm"] = [round(float(v) / 10, 2) for v in (all_t.max(axis=0) - all_t.min(axis=0))]

    # ------------------------------------------------------------- mRNA channel
    mrna_block = None
    if mrna:
        m = mrna[0]
        pts = m["trace_t"]
        nums = [int(v) for v in m["resnums"]]
        seq = []
        for res in model[m["name"]]:
            seq.append(res.name[-1] if len(res.name) > 1 else res.name)
        seglen = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        # geometric P site: mRNA residue closest to the tRNAi anticodon (res 34-36)
        psite = None
        if by_part.get("tRNAi"):
            t = by_part["tRNAi"][0]
            ac = t["trace_t"][(t["resnums"] >= 34) & (t["resnums"] <= 36)]
            if len(ac):
                dmat = np.linalg.norm(pts[:, None, :] - ac[None, :, :], axis=2).min(axis=1)
                k = int(dmat.argmin())
                psite = {"mrna_res": nums[k], "nt": seq[k],
                         "dist_to_tRNAi_anticodon_A": round(float(dmat[k]), 1),
                         "note": "mRNA residue whose P atom is closest to a tRNAi "
                                 "anticodon (res 34-36) P atom; in this open scanning "
                                 "state the codon-anticodon duplex is NOT formed, so this "
                                 "is a geometric proxy only"}
        # what lines each mRNA end (independent support for the entry/exit labels)
        def near(p, cutoff=20.0, k=6):
            hits = []
            for c in chains:
                if c["part"] == "mRNA":
                    continue
                d = np.linalg.norm(c["allxyz_t"] - p, axis=1).min()
                if d < cutoff:
                    hits.append((round(float(d), 1), c["desc"]))
            return [f"{n} ({d} A)" for d, n in sorted(hits)[:k]]

        mrna_block = {
            "chain": m["name"],
            "n_nt_resolved": len(nums),
            "res_numbering": "author numbering of this entry's mRNA construct "
                             "(entity sequence is 255 nt; 0 is skipped)",
            "res_first": nums[0], "res_last": nums[-1],
            "res_numbers": nums,
            "sequence_resolved_5to3": "".join(seq),
            "path_xyz_5to3_A": [round(float(v), 1) for v in pts.ravel()],
            "contour_length_nm": round(float(seglen.sum()) / 10, 2),
            "mean_PP_spacing_A": round(float(seglen.mean()), 2),
            "end_5prime_exit_channel": {
                "res": nums[0], "nt": seq[0],
                "xyz_A": [round(float(v), 1) for v in pts[0]],
                "neighbours_40S_proteins": near(pts[0]),
            },
            "end_3prime_entry_channel": {
                "res": nums[-1], "nt": seq[-1],
                "xyz_A": [round(float(v), 1) for v in pts[-1]],
                "neighbours_40S_proteins": near(pts[-1]),
            },
            "direction_note": "The 40S translocates 5'->3', so mRNA is drawn IN at its "
                              "downstream (3') end through the entry channel and pushed "
                              "OUT at its upstream (5') end through the exit channel.",
            "p_site_geometric": psite,
        }
        # where is the AUG of the construct, and how far is the entry-channel edge
        # from the geometric P site? (cross-checks the published toeprint leading edge)
        s = "".join(seq)
        k = s.find("AUG")
        if k >= 0:
            mrna_block["start_codon_in_resolved_path"] = {
                "res_A_of_AUG": nums[k], "resolved_nt_of_codon": s[k:k + 3],
                "note": "first AUG inside the resolved stretch, author numbering",
            }
        elif s.endswith("AU"):
            mrna_block["start_codon_in_resolved_path"] = {
                "res_A_of_AUG": nums[-2],
                "resolved_nt_of_codon": "AU (G not modelled)",
                "note": "the construct's ACCACC-AUG-G start codon; only its A and U are "
                        "modelled, at the 3' (entry-channel) edge of the resolved path",
            }
        if psite:
            mrna_block["entry_edge_offset_from_p_site_nt"] = {
                "value": nums[-1] - psite["mrna_res"],
                "note": "nt from the geometric P-site residue to the 3'-most modelled "
                        "nucleotide (the leading/entry edge of the 40S). Compare the "
                        "published toeprint leading edge of +16 to +18 nt (Pisareva & "
                        "Pisarev 2016, Nucleic Acids Res 44:4252, PMID 27067542).",
            }
        dims["mRNA_span_in_channel_nt"] = len(nums)
        dims["mRNA_end_to_end_nm"] = round(float(np.linalg.norm(pts[-1] - pts[0])) / 10, 2)

    # ------------------------------------------------------------------- output
    c = hdr["citation"]
    fix = CITATION_FIXUP.get(hdr["pdb_id"])
    if fix:
        c.update({k: v for k, v in fix.items()
                  if k in ("journal", "volume", "page_first", "page_last", "year", "pubmed")})
        c["fixup_note"] = fix["verified"]
        if "correction" in fix:
            c["correction"] = fix["correction"]
    auth = c.get("authors") or []
    auth_s = ", ".join(auth[:6]) + (" et al." if len(auth) > 6 else "")
    pages = c.get("page_first") or ""
    if pages and c.get("page_last"):
        pages += "-" + c["page_last"]
    cit = (f"{auth_s} {c.get('title','').rstrip('.')}. {c.get('journal')} "
           f"{c.get('volume') or ''}{':' + pages if pages else ''} "
           f"({c.get('year')}).").replace("  ", " ")
    if not (c.get("volume") and pages):
        print("  WARNING incomplete citation (no volume/pages); add a CITATION_FIXUP entry")

    parts_summary = {}
    for p, cs in sorted(by_part.items()):
        parts_summary[p] = {
            "chains": sorted(x["name"] for x in cs),
            "n_chains": len(cs),
            "n_residues_modelled": int(sum(len(x["resnums"]) for x in cs)),
            "trace_step": STEP.get(p, 3),
            "entity_descriptions": sorted({x["desc"] for x in cs}),
        }

    doc = {
        "format": "utr-race/struct-v1",
        "meta": {
            "pdb_id": hdr["pdb_id"],
            "title": hdr["title"],
            "resolution_A": float(hdr["resolution_A"]) if hdr["resolution_A"] else None,
            "resolution_method": hdr["resolution_method"],
            "method": hdr["method"],
            "organism": "Homo sapiens",
            "organism_note": "all polymer entities of this entry are annotated Homo sapiens "
                             "in RCSB (reconstituted from human components); the mRNA is a "
                             "synthetic model mRNA construct, not a cellular transcript",
            "state": hdr["title"],
            "citation": cit,
            "citation_fields": c,
            "doi": c.get("doi"),
            "pubmed": c.get("pubmed"),
            "url": f"https://www.rcsb.org/structure/{hdr['pdb_id']}",
            "cif_url": f"https://files.rcsb.org/download/{hdr['pdb_id']}.cif.gz",
            "license": "PDB data are freely available",
            "generated": date.today().isoformat(),
            "generated_by": "pipeline/structures/build_struct.py (gemmi "
                            + gemmi.__version__ + ")",
            "units": "all coordinates in Angstroms; all 'dimensions_nm' values in nanometres",
            "origin": "centroid of all modelled non-H atoms of 18S rRNA + 40S ribosomal "
                      "proteins (unweighted, i.e. geometric centroid, not mass-weighted)",
            "origin_in_cif_frame_A": [round(float(v), 2) for v in origin],
            "orientation": orient_note,
            "rotation_matrix_row_major": [round(float(v), 6) for v in rot.ravel()],
            "parts_mapping": {
                "rule": "substring match on mmCIF _entity.pdbx_description, first hit wins",
                "rules": [{"match": s, "part": p} for s, p in PART_RULES],
                "parts": parts_summary,
                "parts_absent": missing,
            },
            "dimensions_nm": dims,
            "n_mrna_nt_resolved": mrna_block["n_nt_resolved"] if mrna_block else 0,
            "caveats": [
                "Coordinates are decimated backbone traces, not full atomic detail.",
                "The mRNA is the synthetic construct used in vitro by the depositors.",
                "In this entry (48S-1, open codon-scanning state) the start codon is NOT "
                "in the P site: author residues 18/19 are the A and U of the ACCACC-AUG-G "
                "start codon and sit at the 3' (entry-channel) edge of the resolved path.",
                "Residue numbering is the depositors' construct numbering, not "
                "AUG-relative; 0 is skipped.",
            ],
        },
        "geometry": {
            "traces": traces,
            "traces_note": "one entry per continuous modelled stretch; xyz is a flat "
                           "[x0,y0,z0,x1,...] polyline in Angstroms, 0.1 A precision",
            "surface_spheres": {
                "note": "one sphere per up to %d consecutive residues of 18S rRNA + 40S "
                        "proteins; centre = mean of the window's non-H atoms, radius = "
                        "RMS atom spread + %.1f A padding. Draw as an InstancedMesh of "
                        "low-poly spheres for a cheap blob." % (SURFACE_WINDOW, VDW_PAD_A),
                "stride": "flat [x,y,z,r] per sphere, Angstroms",
                "count": len(spheres),
                "xyzr": sph_xyzr,
            },
            "mrna_channel": mrna_block,
        },
    }

    # ------------------------------------------------- flat point-cloud compat view
    # app/src/struct3d.js reads a generic shape: top-level `parts` (name + points,
    # each point [x,y,z] or [x,y,z,r]), a channel path as a plain point array, and
    # top-level pdb_id / citation / dimensions_nm. Provided here so the asset drops
    # straight into that renderer. It is a REDUNDANT VIEW of `geometry` above:
    # the 40S core comes from surface_spheres (so it carries radii), everything else
    # from the backbone traces. Same frame, same Angstrom units.
    compat_parts = []
    for p in CORE_40S:
        pl = [[round(float(v), 1) for v in (*cen, r)]
              for cen, r, part in spheres if part == p]
        if pl:
            compat_parts.append({"name": p, "kind": "spheres",
                                 "n": len(pl), "points": pl})
    COMPAT_STEP = {"mRNA": 1, "tRNAi": 1, "eIF1": 1, "eIF1A": 1}
    for p in sorted(by_part):
        if p in CORE_40S:
            continue
        st_ = COMPAT_STEP.get(p, 2)
        pl = []
        for c in by_part[p]:
            for q in c["trace_t"][::st_]:
                pl.append([round(float(v), 1) for v in q])
        if pl:
            compat_parts.append({"name": p, "kind": "trace_points",
                                 "n": len(pl), "points": pl})
    doc["parts"] = compat_parts
    doc["parts_note"] = ("flat point-cloud view of `geometry`, for generic renderers; "
                         "40S core points are [x,y,z,radius] surface spheres, the rest "
                         "are [x,y,z] backbone trace atoms. Angstroms, same frame.")
    if mrna_block:
        doc["channel_path"] = [[round(float(v), 1) for v in q] for q in mrna[0]["trace_t"]]
        doc["channel_path_note"] = ("mRNA phosphate path 5'->3', Angstroms; same data as "
                                    "geometry.mrna_channel.path_xyz_5to3_A")
    doc["pdb_id"] = doc["meta"]["pdb_id"]
    doc["citation"] = doc["meta"]["citation"]
    doc["dimensions_nm"] = dims  # alias of meta.dimensions_nm

    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(doc, f, separators=(",", ":"))
    print(f"wrote {out}  {os.path.getsize(out)/1e6:.2f} MB")
    print(f"  traces={len(traces)}  trace_points="
          f"{sum(t['n'] for t in traces)}  spheres={len(spheres)}")
    print(f"  parts: " + ", ".join(f"{p}:{v['n_chains']}ch" for p, v in parts_summary.items()))
    print(f"  40S max diameter {dims['40S_max_diameter_nm']} nm  bbox {dims['40S_bbox_nm']} nm")
    if mrna_block:
        print(f"  mRNA {mrna_block['n_nt_resolved']} nt resolved "
              f"({mrna_block['res_first']}..{mrna_block['res_last']}), "
              f"contour {mrna_block['contour_length_nm']} nm")
        print(f"  P-site proxy: {mrna_block['p_site_geometric']}")
        print(f"  5' end near: {mrna_block['end_5prime_exit_channel']['neighbours_40S_proteins']}")
        print(f"  3' end near: {mrna_block['end_3prime_entry_channel']['neighbours_40S_proteins']}")
    if missing:
        print("  WARNING parts absent:", missing)
    return 0


if __name__ == "__main__":
    sys.exit(main())
