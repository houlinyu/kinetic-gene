#!/usr/bin/env python
"""Build lightweight 3D assets for the transcription+splicing tab from real PDB entries.

Inputs : pipeline/structures_txn/dl/*.cif  (see fetch.sh)
Outputs: app/public/struct/polii.json, spliceosome.json, txn_meta.json

Format follows `utr-race/struct-v1`, the convention established by
pipeline/structures/build_struct.py for app/public/struct/43s.json: top-level
{format, meta, geometry}; meta carries the PDB provenance, the parts mapping, the
origin/orientation (with an explicit rotation matrix) and dimensions_nm; geometry
carries per-stretch backbone `traces`, a `surface_spheres` point cloud
(flat [x,y,z,r] quadruples) and one domain-specific path block.

Everything written out is either (a) a coordinate taken straight from the deposited
mmCIF, or (b) a quantity measured from those coordinates by this script. Nothing is
invented; every derivation is implemented here so it can be audited.
"""
from __future__ import annotations

import json
import os
import numpy as np
import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DL = os.path.join(HERE, "dl")
OUT = os.path.join(ROOT, "app", "public", "struct")

FORMAT = "utr-race/struct-v1"
GENERATED = "2026-11-22"
GENERATED_BY = f"pipeline/structures_txn/build_txn_structures.py (gemmi {gemmi.__version__})"
LICENSE = "PDB data are freely available"
UNITS = "all coordinates in Angstroms; all 'dimensions_nm' values in nanometres"
WINDOW = 10          # residues per surface sphere (matches 43s.json)
SPHERE_PAD = 3.0     # A padding on the RMS atom spread (matches 43s.json)

# --------------------------------------------------------------- provenance

META = {
    "6GMH": {
        "pdb_id": "6GMH",
        "title": "Structure of activated transcription complex Pol II-DSIF-PAF-SPT6",
        "resolution_A": 3.1,
        "resolution_method": "FSC 0.143 CUT-OFF",
        "method": "ELECTRON MICROSCOPY",
        "organism": "Sus scrofa RNA polymerase II; Homo sapiens DSIF (SPT4/SPT5), PAF1 "
                    "complex and SPT6; synthetic nucleic-acid scaffold",
        "organism_note": "MAMMALIAN, NOT HUMAN Pol II. The paper title states 'an activated "
                         "elongation complex of Sus scrofa Pol II and Homo sapiens DSIF, PAF "
                         "and SPT6'. Pig Pol II is the standard source for mammalian Pol II "
                         "elongation complexes; the elongation factors here are human. Label "
                         "the enzyme 'mammalian Pol II' in the UI, not 'human Pol II'.",
        "state": "activated transcription elongation complex (post pause-release), with a "
                 "complete 48-bp transcription scaffold and nascent RNA",
        "citation": "Vos SM, Farnung L, Boehning M, Wigge C, Linden A, Urlaub H, Cramer P. "
                    "Structure of activated transcription complex Pol II-DSIF-PAF-SPT6. "
                    "Nature 560:607-612 (2018).",
        "citation_fields": {
            "title": "Structure of activated transcription complex Pol II-DSIF-PAF-SPT6.",
            "journal": "Nature",
            "volume": "560",
            "page_first": "607",
            "year": "2018",
            "doi": "10.1038/s41586-018-0440-4",
            "pubmed": "30135578",
            "authors": ["Vos, S.M.", "Farnung, L.", "Boehning, M.", "Wigge, C.",
                        "Linden, A.", "Urlaub, H.", "Cramer, P."],
        },
        "doi": "10.1038/s41586-018-0440-4",
        "pubmed": "30135578",
        "url": "https://www.rcsb.org/structure/6GMH",
        "cif_url": "https://files.rcsb.org/download/6GMH.cif.gz",
    },
    "5XJC": {
        "pdb_id": "5XJC",
        "title": "Cryo-EM structure of the human spliceosome just prior to exon ligation at "
                 "3.6 angstrom (C* complex)",
        "resolution_A": 3.6,
        "resolution_method": "FSC 0.143 CUT-OFF",
        "method": "ELECTRON MICROSCOPY",
        "organism": "Homo sapiens (snRNAs and proteins); pre-mRNA substrate derived from the "
                    "Human adenovirus 2 major late leader",
        "organism_note": "HUMAN spliceosome - no yeast substitution was needed for this asset.",
        "state": "C* complex: catalytic step-II-activated spliceosome, after branching "
                 "(5' exon cleaved, intron lariat formed) and before exon ligation",
        "citation": "Zhang X, Yan C, Hang J, Finci LI, Lei J, Shi Y. An Atomic Structure of "
                    "the Human Spliceosome. Cell 169:918-929 (2017).",
        "citation_fields": {
            "title": "An Atomic Structure of the Human Spliceosome.",
            "journal": "Cell",
            "volume": "169",
            "page_first": "918",
            "year": "2017",
            "doi": "10.1016/j.cell.2017.04.033",
            "pubmed": "28502770",
            "authors": ["Zhang, X.", "Yan, C.", "Hang, J.", "Finci, L.I.", "Lei, J.", "Shi, Y."],
        },
        "doi": "10.1016/j.cell.2017.04.033",
        "pubmed": "28502770",
        "url": "https://www.rcsb.org/structure/5XJC",
        "cif_url": "https://files.rcsb.org/download/5XJC.cif.gz",
    },
    "4PJO": {
        "pdb_id": "4PJO",
        "title": "Minimal U1 snRNP (human U1 snRNP bound to a 5' splice site RNA)",
        "resolution_A": 3.3,
        "resolution_method": "X-ray, highest-resolution shell reported by depositors",
        "method": "X-RAY DIFFRACTION",
        "organism": "Homo sapiens",
        "organism_note": "Engineered minimal U1 snRNP: the U1 snRNA four-helix junction is "
                         "replaced by a kissing loop and stem-loop 4 is shortened. The bound "
                         "10-mer 5'-AGGUAAGUCC-3' is a 5' splice site.",
        "state": "U1 snRNP with a 5' splice site RNA base-paired to the 5' end of U1 snRNA "
                 "(the earliest step of splice-site recognition)",
        "citation": "Kondo Y, Oubridge C, van Roon AM, Nagai K. Crystal structure of human U1 "
                    "snRNP, a small nuclear ribonucleoprotein particle, reveals the mechanism "
                    "of 5' splice site recognition. eLife 4:e04986 (2015).",
        "citation_fields": {
            "title": "Crystal structure of human U1 snRNP, a small nuclear ribonucleoprotein "
                     "particle, reveals the mechanism of 5' splice site recognition.",
            "journal": "eLife",
            "volume": "4",
            "page_first": "e04986",
            "year": "2015",
            "doi": "10.7554/eLife.04986",
            "pubmed": "25555158",
            "authors": ["Kondo, Y.", "Oubridge, C.", "van Roon, A.M.", "Nagai, K."],
        },
        "doi": "10.7554/eLife.04986",
        "pubmed": "25555158",
        "url": "https://www.rcsb.org/structure/4PJO",
        "cif_url": "https://files.rcsb.org/download/4PJO.cif.gz",
    },
}

# --------------------------------------------------------------- helpers


def load(pdb_id):
    st = gemmi.read_structure(os.path.join(DL, f"{pdb_id.lower()}.cif"))
    st.setup_entities()
    return st


def is_protein(ch):
    return str(ch.get_polymer().check_polymer_type()).endswith("PeptideL")


def atom_xyz(a):
    return np.array([a.pos.x, a.pos.y, a.pos.z])


def rep_atoms(ch, kind, lo=None, hi=None):
    """Ordered (resnum, resname, xyz) of the representative backbone atom."""
    names = ["CA"] if kind == "protein" else ["P", "C1'"]
    out = []
    for r in ch:
        if lo is not None and not (lo <= r.seqid.num <= hi):
            continue
        for n in names:
            a = r.find_atom(n, "*")
            if a is not None:
                out.append((r.seqid.num, r.name, atom_xyz(a)))
                break
    return out


def all_atoms(ch, lo=None, hi=None):
    """(resnum, [xyz...]) for every modelled non-H atom, used for surface spheres."""
    out = []
    for r in ch:
        if lo is not None and not (lo <= r.seqid.num <= hi):
            continue
        pts = [atom_xyz(a) for a in r if a.element != gemmi.Element("H")]
        if pts:
            out.append((r.seqid.num, np.asarray(pts)))
    return out


def flat(xyz, nd=1):
    out = []
    for p in xyz:
        out.extend(round(float(v), nd) for v in p)
    return out


def runs(pts, max_gap=1):
    """Split an ordered residue list into continuous modelled stretches."""
    segs, cur = [], [pts[0]]
    for a, b in zip(pts, pts[1:]):
        if b[0] - a[0] <= max_gap:
            cur.append(b)
        else:
            segs.append(cur)
            cur = [b]
    segs.append(cur)
    return segs


def decimate(seq, step):
    if step <= 1:
        return seq
    keep = seq[::step]
    if keep[-1] is not seq[-1]:
        keep = keep + [seq[-1]]
    return keep


def trace_entries(chain_name, part, label, pts, kind, step):
    atom = "CA" if kind == "protein" else "P"
    out = []
    for seg in runs(pts):
        d = decimate(seg, step)
        out.append({
            "chain": chain_name, "part": part, "label": label, "atom": atom, "step": step,
            "res_first": seg[0][0], "res_last": seg[-1][0], "n": len(d),
            "xyz": flat([p[2] for p in d]),
        })
    return out


def spheres_for(residues):
    """One sphere per <=WINDOW consecutive residues; centre = mean atom, r = RMS + pad."""
    out = []
    for i in range(0, len(residues), WINDOW):
        win = residues[i:i + WINDOW]
        pts = np.vstack([a for _, a in win])
        c = pts.mean(0)
        r = float(np.sqrt(((pts - c) ** 2).sum(1).mean())) + SPHERE_PAD
        out.append((c, r))
    return out


def extent(points):
    a = np.asarray(points)
    lo, hi = a.min(0), a.max(0)
    s = a[:: max(1, len(a) // 3500)]
    d = np.sqrt(((s[:, None, :] - s[None, :, :]) ** 2).sum(-1)).max()
    pc = a - a.mean(0)
    _, _, vt = np.linalg.svd(pc, full_matrices=False)
    pr = pc @ vt.T
    return {
        "bbox_nm": [round(float(v) / 10, 2) for v in (hi - lo)],
        "bbox_min_nm": [round(float(v) / 10, 2) for v in lo],
        "bbox_max_nm": [round(float(v) / 10, 2) for v in hi],
        "max_diameter_nm": round(float(d) / 10, 2),
        "radius_of_gyration_nm": round(float(np.sqrt((pc ** 2).sum(1).mean())) / 10, 2),
        "principal_axis_extents_nm": [round(float(v) / 10, 2)
                                      for v in (pr.max(0) - pr.min(0))],
    }


def frame(x_dir, z_hint):
    """Right-handed orthonormal frame: +x = x_dir, +z ~ z_hint, +y = z x x."""
    x = x_dir / np.linalg.norm(x_dir)
    z = z_hint - x * float(z_hint @ x)
    z /= np.linalg.norm(z)
    y = np.cross(z, x)
    return np.vstack([x, y, z])          # rows; p_out = R @ (p_cif - origin)


def apply_frame(pts, origin, R):
    return [(R @ (p - origin)) for p in pts]


# --------------------------------------------------------------- Pol II

POLII = {"A": "RPB1", "B": "RPB2", "C": "RPB3", "D": "RPB4", "E": "RPABC1/RPB5",
         "F": "RPABC2/RPB6", "G": "RPB7", "H": "RPABC3/RPB8", "I": "RPB9",
         "J": "RPABC5/RPB10", "K": "RPB11", "L": "RPABC4/RPB12"}
FACTORS = {"Y": "SPT4 (DSIF)", "Z": "SPT5 (DSIF)", "M": "SPT6", "V": "PAF1", "U": "LEO1",
           "Q": "CTR9", "W": "WDR61", "X": "CDC73"}
CLAMP = (1, 346)
PAIR_CUT = 12.0     # C1'-C1' cutoff: a Watson-Crick pair is ~10.5 A, a broken one is >12


def build_polii():
    st = load("6GMH")
    m = st[0]
    ch = {c.name: c for c in m}

    def c1(name, lo=None, hi=None):
        out = []
        for r in ch[name]:
            if lo is not None and not (lo <= r.seqid.num <= hi):
                continue
            a = r.find_atom("C1'", "*")
            if a is not None:
                out.append((r.seqid.num, r.name, atom_xyz(a)))
        return out

    tmpl, nont, rna = (rep_atoms(ch[k], "na") for k in ("T", "N", "P"))
    tc1, nc1, rc1 = c1("T"), c1("N"), c1("P")

    # -- base-pair register between the two DNA strands, determined from geometry.
    # The scaffold is antiparallel, so template i pairs with non-template (k - i) for a
    # single offset k. Find the k that maximises the number of C1'-C1' contacts.
    nmap = {n: x for n, _, x in nc1}
    best = None
    for k in range(2, 100):
        hits = [(n, nmap[k - n], np.linalg.norm(x - nmap[k - n]))
                for n, _, x in tc1 if (k - n) in nmap]
        good = [h for h in hits if h[2] <= PAIR_CUT]
        if best is None or len(good) > best[1]:
            best = (k, len(good))
    reg_k = best[0]
    pairs = [(n, reg_k - n, float(np.linalg.norm(x - nmap[reg_k - n])), (x + nmap[reg_k - n]) / 2)
             for n, _, x in tc1 if (reg_k - n) in nmap]
    paired = sorted([p for p in pairs if p[2] <= PAIR_CUT])
    unpaired = sorted(n for n, _, d, _ in pairs if d > PAIR_CUT)
    nont_nums = {p[0] for p in nont}
    nont_missing = [i for i in range(min(nont_nums), max(nont_nums) + 1) if i not in nont_nums]
    bubble = sorted(set(unpaired) | {reg_k - i for i in nont_missing if 1 <= reg_k - i <= 48})

    # -- duplex axis: contiguous runs of base-pair midpoints
    axis_runs, cur = [], [paired[0]]
    for a, b in zip(paired, paired[1:]):
        if b[0] - a[0] == 1:
            cur.append(b)
        else:
            axis_runs.append(cur)
            cur = [b]
    axis_runs.append(cur)
    axis_runs = [r for r in axis_runs if len(r) > 1]
    longest = max(axis_runs, key=len)

    # -- measured duplex geometry on the longest duplex run.
    # The C1'-C1' midpoints precess around the true helix axis, so a straight
    # point-to-point distance over-estimates the rise. Fit the axis by PCA over the run
    # and project, which removes the precession.
    mids = np.vstack([p[3] for p in longest])
    ax = np.linalg.svd(mids - mids.mean(0), full_matrices=False)[2][0]
    if (mids[-1] - mids[0]) @ ax < 0:
        ax = -ax
    proj = mids @ ax
    rise_A = float((proj[-1] - proj[0]) / (len(longest) - 1))
    rises = [float(b - a) for a, b in zip(proj, proj[1:])]
    pmap_t = {n: x for n, _, x in rep_atoms(ch["T"], "na")}
    pmap_n = {n: x for n, _, x in rep_atoms(ch["N"], "na")}
    pp = [float(np.linalg.norm(pmap_t[n] - pmap_n[reg_k - n]))
          for n, _, _, _ in longest if n in pmap_t and (reg_k - n) in pmap_n]
    # helical twist about the fitted axis, from successive C1'->C1' vectors
    tmap = {n: x for n, _, x in tc1}
    tw = []
    for a, b in zip(longest, longest[1:]):
        va = nmap[reg_k - a[0]] - tmap[a[0]]
        vb = nmap[reg_k - b[0]] - tmap[b[0]]
        va = va - ax * (va @ ax)
        vb = vb - ax * (vb @ ax)
        c = float(np.clip((va @ vb) / (np.linalg.norm(va) * np.linalg.norm(vb)), -1, 1))
        tw.append(np.degrees(np.arccos(c)))

    # -- DNA:RNA hybrid
    rr = np.asarray([p[2] for p in rc1])
    tt = np.asarray([p[2] for p in tc1])
    hyb_t = sorted(n for n, _, x in tc1 if np.linalg.norm(rr - x, axis=1).min() <= PAIR_CUT)
    hyb_r = sorted(n for n, _, x in rc1 if np.linalg.norm(tt - x, axis=1).min() <= PAIR_CUT)

    # -- origin and oriented frame
    core_pts = [p for cid in POLII if cid in ch for p in rep_atoms(ch[cid], "protein")]
    fac_pts = [p for cid in FACTORS if cid in ch for p in rep_atoms(ch[cid], "protein")]
    origin = np.asarray([p[2] for p in core_pts]).mean(0)
    # +x = downstream -> upstream DNA travel, from the axis midpoint of the two duplex runs
    down = min(axis_runs, key=lambda r: r[0][0])
    up = max(axis_runs, key=lambda r: r[0][0])
    x_dir = up[-1][3] - down[0][3]
    rna_sorted = sorted(rna, key=lambda p: p[0])
    exit_from = next(p for p in rna_sorted if p[0] == min(hyb_r))
    exit_to = rna_sorted[0]
    R = frame(x_dir, exit_to[2] - exit_from[2])

    def to_frame(pts):
        return [(n, nm, R @ (x - origin)) for n, nm, x in pts]

    def to_frame_atoms(res):
        return [(n, np.asarray([R @ (a - origin) for a in arr])) for n, arr in res]

    traces, parts, sph_res = [], {}, []
    spec = [
        ("polii_core", "RNA polymerase II (12 subunits)", list(POLII), "protein", 2, True),
        ("clamp", "RPB1 clamp region (approximate)", ["A"], "protein", 2, False),
        ("elongation_factors", "DSIF, PAF1 complex and SPT6", list(FACTORS), "protein", 2, True),
        ("template_DNA", "template DNA strand", ["T"], "na", 1, True),
        ("nontemplate_DNA", "non-template DNA strand", ["N"], "na", 1, True),
        ("nascent_RNA", "nascent RNA transcript", ["P"], "na", 1, True),
    ]
    for part, label, cids, kind, step, in_surface in spec:
        lo, hi = (CLAMP if part == "clamp" else (None, None))
        chains, nres, descs = [], 0, []
        for cid in cids:
            if cid not in ch:
                continue
            pts = to_frame(rep_atoms(ch[cid], kind, lo, hi))
            if not pts:
                continue
            chains.append(cid)
            nres += len(pts)
            descs.append(POLII.get(cid) or FACTORS.get(cid) or
                         {"T": "Template DNA", "N": "Non-template DNA", "P": "RNA"}[cid])
            traces += trace_entries(cid, part, label, pts, kind, step)
            if in_surface:
                for seg in runs(to_frame_atoms(all_atoms(ch[cid], lo, hi))):
                    sph_res.append(seg)
        parts[part] = {"chains": chains, "n_chains": len(chains),
                       "n_residues_modelled": nres, "trace_step": step,
                       "entity_descriptions": descs}
    parts["clamp"]["residue_range"] = list(CLAMP)
    parts["clamp"]["warning"] = (
        "APPROXIMATE. The residue range is the clamp definition for the homologous "
        "S. cerevisiae Rpb1 (Cramer, Bushnell & Kornberg, Science 292:1863, 2001) applied "
        "to this entry's RPB1 numbering; it was not independently validated against this "
        "model. Use as a highlight only, never as a stated domain boundary.")

    spheres = [s for seg in sph_res for s in spheres_for(seg)]
    xyzr = []
    for c, r in spheres:
        xyzr.extend([round(float(c[0]), 1), round(float(c[1]), 1), round(float(c[2]), 1),
                     round(float(r), 1)])

    core_f = [p[2] for p in to_frame(core_pts)]
    all_f = core_f + [p[2] for p in to_frame(fac_pts + tmpl + nont + rna)]
    na_f = [p[2] for p in to_frame(tmpl + nont + rna)]

    e_core, e_all = extent(core_f), extent(all_f)
    dims = {
        "frame": "oriented frame described in meta.orientation; nm",
        "polii_definition": "the 12 Pol II subunits RPB1-RPB12 as modelled in this entry",
        **{f"polii_{k}": v for k, v in e_core.items() if k != "frame"},
        "polii_max_diameter_note": "largest distance between any two Pol II backbone trace "
                                   "atoms (protein CA); all-atom extent is ~0.5 nm larger",
        **{f"complex_{k}": v for k, v in e_all.items() if k != "frame"},
        "complex_definition": "Pol II + DSIF + PAF1 complex + SPT6 + the 48-bp scaffold",
        "scaffold_end_to_end_nm": round(float(np.linalg.norm(
            np.asarray(na_f).max(0) - np.asarray(na_f).min(0))) / 10, 2),
    }

    exit_f = {n: x for n, _, x in to_frame(rna_sorted)}
    v = exit_f[exit_to[0]] - exit_f[exit_from[0]]
    axis_f = {}
    for n, _, _, mid in paired:
        axis_f[n] = R @ (mid - origin)

    doc = {
        "format": FORMAT,
        "meta": {
            **META["6GMH"],
            "license": LICENSE,
            "generated": GENERATED,
            "generated_by": GENERATED_BY,
            "units": UNITS,
            "origin": "centroid of the backbone trace atoms (CA) of the 12 Pol II subunits",
            "origin_in_cif_frame_A": [round(float(x), 2) for x in origin],
            "orientation": "+x = direction of DNA travel through the enzyme, from the midpoint "
                           "of the downstream duplex to the midpoint of the upstream duplex "
                           "(Pol II moves along the template towards +x in this frame); "
                           "+z = component of the RNA exit direction orthogonal to +x; "
                           "+y = z x x. p_out = R @ (p_cif - origin).",
            "rotation_matrix_row_major": [round(float(v), 6) for v in R.reshape(-1)],
            "parts_mapping": {
                "rule": "explicit auth_asym_id -> part assignment (this entry has a small, "
                        "fully annotated chain set; see pipeline/structures_txn/README.md)",
                "parts": parts,
                "parts_absent": [],
            },
            "dimensions_nm": dims,
            "n_scaffold_nt_modelled": {"template_DNA": len(tmpl), "nontemplate_DNA": len(nont),
                                       "nascent_RNA": len(rna)},
            "caveats": [
                "Coordinates are decimated backbone traces, not full atomic detail.",
                "The enzyme is Sus scrofa (pig) Pol II, not human; the elongation factors are "
                "human. Do not label the polymerase 'human' in the UI.",
                "The nucleic-acid scaffold is a synthetic construct used in vitro by the "
                "depositors, not a cellular gene.",
                "Part of the single-stranded non-template strand in the bubble "
                f"(residues {nont_missing[0]}-{nont_missing[-1]}) is disordered and absent from "
                "the deposited model, so a bubble drawn from these coordinates alone has a gap.",
                "Only 21 nt of the 46-nt nascent RNA entity are modelled; the RNA emerging "
                "beyond the exit channel is not resolved.",
                "The 'clamp' part is an approximate, literature-derived residue range (see "
                "meta.parts_mapping.parts.clamp.warning).",
            ],
        },
        "geometry": {
            "traces": traces,
            "traces_note": "one entry per continuous modelled stretch; xyz is a flat "
                           "[x0,y0,z0,x1,...] polyline in Angstroms, 0.1 A precision",
            "surface_spheres": {
                "note": f"one sphere per up to {WINDOW} consecutive residues of Pol II, the "
                        "elongation factors and the nucleic-acid scaffold; centre = mean of the "
                        f"window's non-H atoms, radius = RMS atom spread + {SPHERE_PAD} A "
                        "padding. Draw as an InstancedMesh of low-poly spheres.",
                "stride": "flat [x,y,z,r] per sphere, Angstroms",
                "count": len(spheres),
                "xyzr": xyzr,
            },
            "transcription_bubble": {
                "template_chain": "T",
                "nontemplate_chain": "N",
                "rna_chain": "P",
                "base_pair_register": {
                    "rule": f"template residue i pairs with non-template residue {reg_k} - i",
                    "derivation": "the single antiparallel offset that maximises the number of "
                                  f"template-C1' to non-template-C1' contacts <= {PAIR_CUT} A",
                    "n_pairs_within_cutoff": len(paired),
                },
                "bubble_nt": len(bubble),
                "bubble_template_residues": [bubble[0], bubble[-1]],
                "bubble_derivation": "template nucleotides with no non-template partner within "
                                     f"{PAIR_CUT} A, plus the template positions whose "
                                     "non-template partners are not modelled at all; measured "
                                     "from the deposited coordinates of PDB 6GMH",
                "bubble_caveat": "This is the bubble AS MODELLED in 6GMH, i.e. a measurement of "
                                 "this structure. It is not a value quoted from the paper, and "
                                 "the disordered non-template segment means the true "
                                 "single-stranded extent could differ by a nucleotide or two.",
                "nontemplate_not_modelled": {"residues": [nont_missing[0], nont_missing[-1]],
                                             "n_nt": len(nont_missing)},
                "dna_rna_hybrid": {
                    "n_template_nt": len(hyb_t), "n_rna_nt": len(hyb_r),
                    "template_residues": [hyb_t[0], hyb_t[-1]],
                    "rna_residues": [hyb_r[0], hyb_r[-1]],
                    "derivation": f"nucleotides whose C1' lies within {PAIR_CUT} A of a C1' of "
                                  "the other strand; a distance criterion, not a base-pair "
                                  "annotation, so the two counts can differ by one",
                },
                "duplex_axis": {
                    "note": "midpoints of C1'(template)-C1'(non-template) pairs, ordered by "
                            "template residue, split into contiguous duplex runs; the gap "
                            "between runs is the bubble",
                    "runs": [{"template_residues": [r[0][0], r[-1][0]], "n": len(r),
                              "xyz": flat([axis_f[p[0]] for p in r])} for r in axis_runs],
                },
                "measured_b_form_geometry": {
                    "run_used": {"template_residues": [longest[0][0], longest[-1][0]],
                                 "n_bp": len(longest)},
                    "rise_per_bp_A": round(rise_A, 2),
                    "rise_per_bp_sd_A": round(float(np.std(rises)), 2),
                    "rise_per_bp_note": "axial rise along the PCA-fitted duplex axis of this run; the duplex is bound and bent by the enzyme, so it is not a canonical B-DNA value",
                    "twist_per_bp_deg": round(float(np.mean(tw)), 1),
                    "bp_per_turn": round(360.0 / float(np.mean(tw)), 1),
                    "phosphate_phosphate_diameter_A": round(float(np.mean(pp)), 1),
                    "phosphate_phosphate_diameter_sd_A": round(float(np.std(pp)), 1),
                    "derivation": "axis = PCA of the C1'-C1' base-pair midpoints of this run; "
                                  "rise = axial span / (n_bp - 1) along that axis (projecting "
                                  "removes the precession of the midpoints about the axis); "
                                  "twist = angle between successive C1'->C1' vectors projected "
                                  "perpendicular to that axis; diameter = P-to-P distance "
                                  "across a base pair (phosphate convention, so it is smaller "
                                  "than the van-der-Waals width)",
                },
                "rna_exit_channel": {
                    "start_res": exit_from[0], "end_res": exit_to[0],
                    "start_xyz_A": [round(float(x), 1) for x in exit_f[exit_from[0]]],
                    "end_xyz_A": [round(float(x), 1) for x in exit_f[exit_to[0]]],
                    "direction": [round(float(x), 4) for x in v / np.linalg.norm(v)],
                    "chord_length_nm": round(float(np.linalg.norm(v)) / 10, 2),
                    "n_nt_5prime_of_hybrid": int(min(hyb_r) - min(n for n, _, _ in rna_sorted)),
                    "derivation": "from the 5'-most nucleotide of the DNA:RNA hybrid to the "
                                  "5'-most modelled nucleotide of the transcript; the nt count "
                                  "is what is RESOLVED threading the exit channel in this "
                                  "model, not the full nuclease-protected length",
                },
                "direction_note": "Pol II moves along the template DNA towards +x in this "
                                  "frame; DNA enters as the downstream duplex at -x, is opened "
                                  "in the bubble, and re-anneals as the upstream duplex at +x. "
                                  "The nascent RNA leaves sideways through the exit channel.",
            },
        },
    }
    return doc


# --------------------------------------------------------------- spliceosome

SNRNA = {"H": ("snRNA_U2", "U2 snRNA"), "B": ("snRNA_U5", "U5 snRNA"),
         "F": ("snRNA_U6", "U6 snRNA")}


def build_spliceosome():
    st = load("5XJC")
    m = st[0]
    ch = {c.name: c for c in m}
    ent_desc = {}
    doc_cif = gemmi.cif.read(os.path.join(DL, "5xjc.cif")).sole_block()
    e2d = {r[0]: r[1].strip("'\";\n ") for r in doc_cif.find("_entity.", ["id", "pdbx_description"])}
    for r in doc_cif.find("_atom_site.", ["auth_asym_id", "label_entity_id"]):
        ent_desc.setdefault(r[0], e2d.get(r[1], "?"))

    g_all = rep_atoms(ch["G"], "na")
    exon = [p for p in g_all if p[0] < 0]
    intron = [p for p in g_all if p[0] > 0]

    # branch point: the adenosine whose O2' is covalently bonded to the 5'-terminal
    # phosphate of intron residue +1 (the 2'-5' linkage of the lariat). Measured.
    g1 = [r for r in ch["G"] if r.seqid.num == 1][0]
    pa = g1.find_atom("P", "*")
    bp = None
    for r in ch["G"]:
        o2 = r.find_atom("O2'", "*")
        if o2 is not None and pa is not None and o2.pos.dist(pa.pos) < 2.2:
            bp = (r, o2.pos.dist(pa.pos))
            break

    origin = np.asarray([p[2] for ccc in m if True for p in rep_atoms(
        ccc, "protein" if is_protein(ccc) else "na")]).mean(0)
    # +x = 5'->3' principal axis of the pre-mRNA path; +z towards the branch point
    pm = np.asarray([p[2] for p in g_all])
    _, _, vt = np.linalg.svd(pm - pm.mean(0), full_matrices=False)
    x_dir = vt[0]
    if (pm[-1] - pm[0]) @ x_dir < 0:
        x_dir = -x_dir
    bp_xyz = atom_xyz(bp[0].find_atom("C1'", "*")) if bp else pm.mean(0)
    R = frame(x_dir, bp_xyz - origin)

    def tf(pts):
        return [(n, nm, R @ (x - origin)) for n, nm, x in pts]

    def tfa(res):
        return [(n, np.asarray([R @ (a - origin) for a in arr])) for n, arr in res]

    traces, parts, sph_res, all_f = [], {}, [], []

    def add(part, label, cid, kind, step, lo=None, hi=None, surface=True):
        pts = tf(rep_atoms(ch[cid], kind, lo, hi))
        if not pts:
            return 0
        traces.extend(trace_entries(cid, part, label, pts, kind, step))
        if surface:
            for seg in runs(tfa(all_atoms(ch[cid], lo, hi))):
                sph_res.append(seg)
        all_f.extend(p[2] for p in pts)
        return len(pts)

    for cid, (part, label) in SNRNA.items():
        n = add(part, label, cid, "na", 1)
        parts[part] = {"chains": [cid], "n_chains": 1, "n_residues_modelled": n,
                       "trace_step": 1, "entity_descriptions": [ent_desc.get(cid, label)]}
    n = add("exon", "5' exon", "G", "na", 1, lo=-10_000, hi=-1)
    parts["exon"] = {"chains": ["G"], "n_chains": 1, "n_residues_modelled": n, "trace_step": 1,
                     "entity_descriptions": ["pre-mRNA (5' exon portion)"],
                     "residue_range": [exon[0][0], exon[-1][0]],
                     "note": "negative author numbering = 5' exon; in the C* state it is "
                             "already cleaved from the intron at the 5' splice site"}
    n = add("intron", "intron lariat", "G", "na", 1, lo=1, hi=10_000)
    parts["intron"] = {"chains": ["G"], "n_chains": 1, "n_residues_modelled": n, "trace_step": 1,
                       "entity_descriptions": ["pre-mRNA (intron portion)"],
                       "residue_range": [intron[0][0], intron[-1][0]],
                       "note": "author numbering +1 at the 5' splice site; branched (lariat) "
                               "through a 2'-5' bond to the branch-point adenosine"}
    prot_chains, nprot, descs = [], 0, []
    for c in m:
        if is_protein(c):
            k = add("proteins", "spliceosomal proteins", c.name, "protein", 3)
            if k:
                prot_chains.append(c.name)
                nprot += k
                d = ent_desc.get(c.name, "?")
                if d not in descs:
                    descs.append(d)
    parts["proteins"] = {"chains": sorted(prot_chains), "n_chains": len(prot_chains),
                         "n_residues_modelled": nprot, "trace_step": 3,
                         "entity_descriptions": descs}

    spheres = [s for seg in sph_res for s in spheres_for(seg)]
    xyzr = []
    for c, r in spheres:
        xyzr.extend([round(float(c[0]), 1), round(float(c[1]), 1), round(float(c[2]), 1),
                     round(float(r), 1)])

    e_all = extent(all_f)
    rna_f = [p[2] for p in tf(g_all)]
    dims = {
        "frame": "oriented frame described in meta.orientation; nm",
        "spliceosome_definition": "all modelled chains of the C* complex in this entry "
                                  "(including the exon junction complex and AQR)",
        **{f"spliceosome_{k}": v for k, v in e_all.items() if k != "frame"},
        "spliceosome_max_diameter_note": "largest distance between any two backbone trace atoms "
                                         "(protein CA, RNA P); all-atom extent is ~0.5 nm larger",
        "pre_mRNA_end_to_end_nm": round(float(np.linalg.norm(
            np.asarray(rna_f).max(0) - np.asarray(rna_f).min(0))) / 10, 2),
    }

    doc = {
        "format": FORMAT,
        "meta": {
            **META["5XJC"],
            "license": LICENSE,
            "generated": GENERATED,
            "generated_by": GENERATED_BY,
            "units": UNITS,
            "origin": "centroid of all modelled backbone trace atoms (protein CA, nucleic P)",
            "origin_in_cif_frame_A": [round(float(x), 2) for x in origin],
            "orientation": "+x = 1st principal component of the pre-mRNA phosphate path "
                           "(signed 5'->3'); +z = component of (branch-point adenosine - "
                           "origin) orthogonal to +x; +y = z x x. p_out = R @ (p_cif - origin).",
            "rotation_matrix_row_major": [round(float(v), 6) for v in R.reshape(-1)],
            "parts_mapping": {
                "rule": "explicit auth_asym_id -> part assignment; the pre-mRNA chain G is "
                        "split by author residue number (negative = 5' exon, positive = intron)",
                "parts": parts,
                "parts_absent": ["snRNA_U1", "snRNA_U4", "three_prime_exon"],
            },
            "dimensions_nm": dims,
            "state_detail": {
                "name": "C* complex",
                "species_label": "HUMAN (not yeast)",
                "intron_visible": True,
                "branch_point_visible": True,
                "three_prime_exon_resolved": False,
                "note": "U1 and U4 snRNA have left the spliceosome by this stage, so they are "
                        "absent by biology, not by disorder. The 3' exon of the 275-nt pre-mRNA "
                        "entity is not modelled (the highest modelled pre-mRNA residue is +150, "
                        "still intron).",
            },
            "caveats": [
                "Coordinates are decimated backbone traces, not full atomic detail.",
                "Protein traces use every 3rd CA to keep the file small; do not measure "
                "distances along protein traces.",
                "The pre-mRNA substrate is an adenovirus-derived model substrate used in vitro, "
                "not a human transcript.",
                "Only 84 of the 275 pre-mRNA nucleotides are modelled (18 exon + 66 intron).",
                "This is one state of a cycle. Do not animate it as if the whole splicing "
                "reaction were captured in this one structure.",
            ],
        },
        "geometry": {
            "traces": traces,
            "traces_note": "one entry per continuous modelled stretch; xyz is a flat "
                           "[x0,y0,z0,x1,...] polyline in Angstroms, 0.1 A precision",
            "surface_spheres": {
                "note": f"one sphere per up to {WINDOW} consecutive residues of every modelled "
                        "chain; centre = mean of the window's non-H atoms, radius = RMS atom "
                        f"spread + {SPHERE_PAD} A padding. Draw as an InstancedMesh.",
                "stride": "flat [x,y,z,r] per sphere, Angstroms",
                "count": len(spheres),
                "xyzr": xyzr,
            },
            "splice_sites": {
                "pre_mRNA_chain": "G",
                "res_numbering": "author numbering of this entry's pre-mRNA construct: "
                                 "negative = 5' exon, +1 = first intron nucleotide",
                "exon_residues_modelled": [exon[0][0], exon[-1][0]],
                "intron_residues_modelled": [intron[0][0], intron[-1][0]],
                "five_prime_splice_site": {
                    "intron_first_two_nt": "".join(
                        r.name for r in ch["G"] if r.seqid.num in (1, 2)),
                    "derivation": "intron residues +1 and +2 in PDB 5XJC",
                    "note": "the canonical GU dinucleotide at the 5' splice site",
                },
                "branch_point": None if bp is None else {
                    "residue_number": bp[0].seqid.num,
                    "residue_name": bp[0].name,
                    "xyz_A": [round(float(x), 1) for x in
                              R @ (atom_xyz(bp[0].find_atom("C1'", "*")) - origin)],
                    "o2prime_to_intron_5p_phosphorus_A": round(float(bp[1]), 2),
                    "derivation": "the only nucleotide whose O2' lies within covalent bonding "
                                  "distance of the 5'-terminal phosphorus of intron residue +1 "
                                  "in PDB 5XJC, i.e. the 2'-5' branch of the lariat; measured "
                                  "from the deposited coordinates",
                },
                "direction_note": "In this frame the pre-mRNA runs 5'->3' towards +x and the "
                                  "branch point lies towards +z. The 5' exon 3'-end and the "
                                  "intron 5'-end are no longer joined: branching has happened.",
            },
        },
    }
    return doc


# --------------------------------------------------------------- U1 snRNP


def build_u1():
    st = load("4PJO")
    m = st[0]
    ch = {c.name: c for c in m}
    doc_cif = gemmi.cif.read(os.path.join(DL, "4pjo.cif")).sole_block()
    e2d = {r[0]: r[1].strip("'\";\n ") for r in doc_cif.find("_entity.", ["id", "pdbx_description"])}
    c2e, ent_desc = {}, {}
    for r in doc_cif.find("_atom_site.", ["auth_asym_id", "label_entity_id"]):
        c2e.setdefault(r[0], r[1])
    for cid, eid in c2e.items():
        ent_desc[cid] = e2d.get(eid, "?")

    # 4PJO holds four copies in the asymmetric unit. Keep exactly one: for each polymer
    # entity, the single chain closest to the 5' splice site RNA chain X.
    ref = np.asarray([p[2] for p in rep_atoms(ch["X"], "na")])
    best = {}
    for c in m:
        kind = "protein" if is_protein(c) else "na"
        pts = rep_atoms(c, kind)
        if not pts or c.name == "X":
            continue
        a = np.asarray([p[2] for p in pts])
        d = float(np.linalg.norm(a[:, None, :] - ref[None, :, :], axis=-1).min())
        eid = c2e.get(c.name, c.name)
        if eid not in best or d < best[eid][0]:
            best[eid] = (d, c.name, kind)
    keep = [(n, k) for d, n, k in best.values() if d <= 45.0]

    pts_all = [p for n, k in keep for p in rep_atoms(ch[n], k)] + rep_atoms(ch["X"], "na")
    origin = np.asarray([p[2] for p in pts_all]).mean(0)
    ss = rep_atoms(ch["X"], "na")
    x_dir = ss[-1][2] - ss[0][2]
    R = frame(x_dir, np.asarray([p[2] for p in pts_all]).mean(0) - ss[0][2] + np.array([0, 0, 1e-6]))

    def tf(pts):
        return [(n, nm, R @ (x - origin)) for n, nm, x in pts]

    def tfa(res):
        return [(n, np.asarray([R @ (a - origin) for a in arr])) for n, arr in res]

    traces, parts, sph_res, all_f = [], {}, [], []

    def add(part, label, cid, kind, step):
        pts = tf(rep_atoms(ch[cid], kind))
        traces.extend(trace_entries(cid, part, label, pts, kind, step))
        for seg in runs(tfa(all_atoms(ch[cid]))):
            sph_res.append(seg)
        all_f.extend(p[2] for p in pts)
        return len(pts)

    groups = {
        "u1_proteins": ("U1-70K, U1-C and the Sm ring", "protein", 2),
        "snRNA_U1": ("U1 snRNA (engineered minimal construct)", "na", 1),
    }
    for part, (label, kind, step) in groups.items():
        cids = [n for n, k in keep if k == kind]
        nres = sum(add(part, label, cid, kind, step) for cid in cids)
        parts[part] = {"chains": sorted(cids), "n_chains": len(cids),
                       "n_residues_modelled": nres, "trace_step": step,
                       "entity_descriptions": sorted({ent_desc.get(c, "?") for c in cids})}
    n = add("five_prime_splice_site_RNA", "5' splice site RNA", "X", "na", 1)
    parts["five_prime_splice_site_RNA"] = {
        "chains": ["X"], "n_chains": 1, "n_residues_modelled": n, "trace_step": 1,
        "entity_descriptions": [ent_desc.get("X", "?")],
        "sequence": "".join(r.name for r in ch["X"]),
        "note": "10-mer 5' splice site RNA base-paired to the 5' end of U1 snRNA"}

    spheres = [s for seg in sph_res for s in spheres_for(seg)]
    xyzr = []
    for c, r in spheres:
        xyzr.extend([round(float(c[0]), 1), round(float(c[1]), 1), round(float(c[2]), 1),
                     round(float(r), 1)])
    e_all = extent(all_f)

    return {
        "format": FORMAT,
        "meta": {
            **META["4PJO"],
            "license": LICENSE,
            "generated": GENERATED,
            "generated_by": GENERATED_BY,
            "units": UNITS,
            "origin": "centroid of the backbone trace atoms of the selected copy",
            "origin_in_cif_frame_A": [round(float(x), 2) for x in origin],
            "orientation": "+x = 5'->3' direction of the bound 5' splice site RNA (chain X); "
                           "+z = component of (particle centroid - 5'SS 5' end) orthogonal to "
                           "+x; +y = z x x. p_out = R @ (p_cif - origin).",
            "rotation_matrix_row_major": [round(float(v), 6) for v in R.reshape(-1)],
            "parts_mapping": {
                "rule": "4PJO contains four copies of the particle in the asymmetric unit; for "
                        "each polymer entity the single chain closest to 5'SS RNA chain X was "
                        "kept (cutoff 45 A), so exactly one copy is exported",
                "parts": parts,
                "parts_absent": [],
            },
            "dimensions_nm": {
                "frame": "oriented frame described in meta.orientation; nm",
                "u1_definition": "one copy of the minimal U1 snRNP plus the bound 5'SS RNA",
                **{f"u1_{k}": v for k, v in e_all.items() if k != "frame"},
            },
            "caveats": [
                "ENGINEERED CONSTRUCT: the U1 snRNA four-helix junction is replaced by a "
                "kissing loop and stem-loop 4 is shortened, so this is smaller than a native "
                "U1 snRNP. Do not quote its size as the size of U1 snRNP.",
                "U1-A and the SL2 region are not part of this minimal particle.",
                "X-ray structure, so this is a crystal copy, not a particle in solution.",
                "Shown as the step BEFORE the C* complex in spliceosome.json; the two "
                "structures are in different coordinate frames and must not be superposed.",
            ],
            "role": "5' splice site recognition by base pairing with the 5' end of U1 snRNA",
        },
        "geometry": {
            "traces": traces,
            "traces_note": "one entry per continuous modelled stretch; xyz is a flat "
                           "[x0,y0,z0,x1,...] polyline in Angstroms, 0.1 A precision",
            "surface_spheres": {
                "note": f"one sphere per up to {WINDOW} consecutive residues; centre = mean of "
                        f"the window's non-H atoms, radius = RMS atom spread + {SPHERE_PAD} A.",
                "stride": "flat [x,y,z,r] per sphere, Angstroms",
                "count": len(spheres),
                "xyzr": xyzr,
            },
        },
    }


# --------------------------------------------------------------- main


def write(name, doc):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    with open(p, "w") as f:
        json.dump(doc, f, separators=(",", ":"))
    print(f"wrote {p}  {os.path.getsize(p)/1024:.0f} KB")
    return doc


if __name__ == "__main__":
    polii = write("polii.json", build_polii())
    sp = build_spliceosome()
    sp["companion_structures"] = [build_u1()]
    write("spliceosome.json", sp)
    g = polii["geometry"]["transcription_bubble"]
    print("register k =", g["base_pair_register"]["rule"])
    print("bubble nt =", g["bubble_nt"], g["bubble_template_residues"])
    print("hybrid =", g["dna_rna_hybrid"]["n_template_nt"], "/", g["dna_rna_hybrid"]["n_rna_nt"])
    print("B-form measured:", json.dumps(g["measured_b_form_geometry"], indent=1)[:600])
    print("exit:", {k: v for k, v in g["rna_exit_channel"].items() if k != "derivation"})
    print("polii dims:", {k: v for k, v in polii["meta"]["dimensions_nm"].items()
                          if k.startswith("polii_") and "note" not in k and "def" not in k})
    print("spliceo dims:", {k: v for k, v in sp["meta"]["dimensions_nm"].items()
                            if "max_diameter" in k or "bbox_nm" in k})
    print("branch point:", sp["geometry"]["splice_sites"]["branch_point"])
    print("u1 dims:", sp["companion_structures"][0]["meta"]["dimensions_nm"])
