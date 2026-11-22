#!/usr/bin/env python
"""Write app/public/struct/txn_constants.json and txn_meta.json.

Every constant is pulled out of the already-built assets (so it cannot drift from the
geometry the app draws) or out of the GENCODE measurement in gene_features.py.
Nothing is a remembered textbook value: anything that would need a literature quote we
could not read is listed in `omitted` instead of being guessed.

Run AFTER build_txn_structures.py and gene_features.py.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "app", "public", "struct")
GTF = os.path.join(ROOT, "pipeline", "discovery", "dl", "gencode.v47.basic.annotation.gtf.gz")
GENERATED = "2026-11-22"

PDB_SRC = {
    "6GMH": {
        "citation": "Vos SM, Farnung L, Boehning M, Wigge C, Linden A, Urlaub H, Cramer P. "
                    "Structure of activated transcription complex Pol II-DSIF-PAF-SPT6. "
                    "Nature 560:607-612 (2018). doi:10.1038/s41586-018-0440-4, PMID 30135578. "
                    "Deposited coordinates: PDB 6GMH.",
        "url": "https://www.rcsb.org/structure/6GMH",
    },
    "5XJC": {
        "citation": "Zhang X, Yan C, Hang J, Finci LI, Lei J, Shi Y. An Atomic Structure of the "
                    "Human Spliceosome. Cell 169:918-929 (2017). "
                    "doi:10.1016/j.cell.2017.04.033, PMID 28502770. "
                    "Deposited coordinates: PDB 5XJC.",
        "url": "https://www.rcsb.org/structure/5XJC",
    },
    "4PJO": {
        "citation": "Kondo Y, Oubridge C, van Roon AM, Nagai K. Crystal structure of human U1 "
                    "snRNP, a small nuclear ribonucleoprotein particle, reveals the mechanism "
                    "of 5' splice site recognition. eLife 4:e04986 (2015). "
                    "doi:10.7554/eLife.04986, PMID 25555158. "
                    "Deposited coordinates: PDB 4PJO.",
        "url": "https://www.rcsb.org/structure/4PJO",
    },
    "GENCODE": {
        "citation": "GENCODE release 47, basic gene annotation for GRCh38 "
                    "(gencode.v47.basic.annotation.gtf.gz), measured by "
                    "pipeline/structures_txn/gene_features.py.",
        "url": "https://www.gencodegenes.org/human/release_47.html",
    },
}


def c(name, value, unit, label, note, evidence, src, **extra):
    d = {"name": name, "value": value, "unit": unit, "label": label, "note": note,
         "evidence": evidence, "source": PDB_SRC[src]}
    d.update(extra)
    return d


def main():
    polii = json.load(open(os.path.join(OUT, "polii.json")))
    sp = json.load(open(os.path.join(OUT, "spliceosome.json")))
    u1 = sp["companion_structures"][0]
    tb = polii["geometry"]["transcription_bubble"]
    bf = tb["measured_b_form_geometry"]
    pd = polii["meta"]["dimensions_nm"]
    sd = sp["meta"]["dimensions_nm"]

    gf = json.loads(subprocess.check_output(
        [sys.executable, "-I", os.path.join(HERE, "gene_features.py"), GTF]).decode())

    run = bf["run_used"]["template_residues"]
    nbp = bf["run_used"]["n_bp"]
    dna_where = (f"downstream DNA duplex of PDB 6GMH, template residues {run[0]}-{run[1]} "
                 f"({nbp} bp)")

    constants = [
        c("dna_duplex_rise_per_bp", bf["rise_per_bp_A"], "A",
          "DNA axial rise per base pair (measured)",
          "Measured in a real double-stranded DNA duplex bound to Pol II, not a canonical "
          "B-DNA literature value. The duplex is held and slightly bent by the enzyme, so "
          f"treat this as 'the rise in this structure' (sd {bf['rise_per_bp_sd_A']} A).",
          f"Measured from the {dna_where}: {bf['derivation']}",
          "6GMH", sd_A=bf["rise_per_bp_sd_A"], n_bp=nbp,
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("dna_duplex_bp_per_turn", bf["bp_per_turn"], "bp/turn",
          "DNA base pairs per helical turn (measured)",
          f"From the measured mean twist of {bf['twist_per_bp_deg']} deg/bp "
          "(360 / twist). Same caveat as the rise: this is this duplex, in this complex.",
          f"Measured from the {dna_where}: {bf['derivation']}",
          "6GMH", twist_per_bp_deg=bf["twist_per_bp_deg"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("dna_duplex_diameter_phosphate_convention", bf["phosphate_phosphate_diameter_A"], "A",
          "DNA duplex diameter, phosphate-to-phosphate convention (measured)",
          "CONVENTION-DEPENDENT. This is the P-to-P distance across a base pair, which is the "
          "backbone-to-backbone width. The commonly quoted ~20 A for B-DNA is a "
          "van-der-Waals-inclusive convention and is therefore larger; we could not read a "
          "citable source for that number, so it is not given here. Use this one for a tube "
          "radius and add the vdW shell if you want a solid-looking helix.",
          f"Measured from the {dna_where}: mean P(template)-P(non-template) distance across "
          f"each base pair (sd {bf['phosphate_phosphate_diameter_sd_A']} A)",
          "6GMH", sd_A=bf["phosphate_phosphate_diameter_sd_A"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("polii_max_dimension", pd["polii_max_diameter_nm"], "nm",
          "RNA polymerase II maximum dimension",
          "The 12-subunit Pol II core alone (RPB1-RPB12), elongation factors excluded. "
          "Backbone trace atoms only, so the all-atom extent is ~0.5 nm larger. "
          f"Axis-aligned bounding box in the asset's own frame: {pd['polii_bbox_nm']} nm.",
          "Largest distance between any two Pol II CA atoms modelled in PDB 6GMH",
          "6GMH", bbox_nm=pd["polii_bbox_nm"],
          radius_of_gyration_nm=pd["polii_radius_of_gyration_nm"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("polii_elongation_complex_max_dimension", pd["complex_max_diameter_nm"], "nm",
          "Activated Pol II elongation complex maximum dimension",
          "Pol II + DSIF + PAF1 complex + SPT6 + the 48-bp nucleic-acid scaffold, i.e. the "
          "whole particle in this entry. Use this, not the Pol II number, when the view shows "
          "the elongation factors.",
          "Largest distance between any two backbone trace atoms modelled in PDB 6GMH",
          "6GMH", bbox_nm=pd["complex_bbox_nm"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("transcription_bubble_size", tb["bubble_nt"], "nt",
          "Transcription bubble: unpaired DNA nucleotides",
          "MEASURED FROM THE MODEL, not quoted from a paper. "
          f"Template residues {tb['bubble_template_residues'][0]}-"
          f"{tb['bubble_template_residues'][1]} have no non-template partner. "
          f"{tb['nontemplate_not_modelled']['n_nt']} nt of the non-template strand "
          "(residues "
          f"{tb['nontemplate_not_modelled']['residues'][0]}-"
          f"{tb['nontemplate_not_modelled']['residues'][1]}) are disordered and absent from the "
          "deposited model, so the true single-stranded extent could differ by a nucleotide "
          "or two. Do not present this as a published bubble size.",
          tb["bubble_derivation"], "6GMH",
          template_residues=tb["bubble_template_residues"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("dna_rna_hybrid_length", tb["dna_rna_hybrid"]["n_template_nt"], "bp",
          "DNA:RNA hybrid in the Pol II active site",
          "Template nucleotides paired with the nascent transcript. A distance criterion "
          "(C1'-C1' <= 12 A), not a base-pair annotation, so it can be off by one.",
          tb["dna_rna_hybrid"]["derivation"], "6GMH",
          template_residues=tb["dna_rna_hybrid"]["template_residues"],
          rna_residues=tb["dna_rna_hybrid"]["rna_residues"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("nascent_rna_resolved_in_exit_channel",
          tb["rna_exit_channel"]["n_nt_5prime_of_hybrid"], "nt",
          "Nascent RNA resolved beyond the hybrid, in the exit channel",
          "This is what is MODELLED threading the exit channel in this structure. It is NOT "
          "the nuclease-protected length of nascent RNA, and it is a lower bound: only 21 of "
          "the 46 nt of the RNA entity are modelled at all.",
          tb["rna_exit_channel"]["derivation"], "6GMH",
          chord_length_nm=tb["rna_exit_channel"]["chord_length_nm"],
          exit_direction_in_asset_frame=tb["rna_exit_channel"]["direction"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("spliceosome_max_dimension", sd["spliceosome_max_diameter_nm"], "nm",
          "Spliceosome (human C* complex) maximum dimension",
          "All modelled chains of the C* complex, including the exon junction complex and AQR. "
          f"Axis-aligned bounding box in the asset's own frame: {sd['spliceosome_bbox_nm']} nm. "
          "Backbone trace atoms only, so the all-atom extent is ~0.5 nm larger.",
          "Largest distance between any two backbone trace atoms modelled in PDB 5XJC",
          "5XJC", bbox_nm=sd["spliceosome_bbox_nm"],
          radius_of_gyration_nm=sd["spliceosome_radius_of_gyration_nm"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("u1_snrnp_max_dimension", u1["meta"]["dimensions_nm"]["u1_max_diameter_nm"], "nm",
          "Minimal U1 snRNP + 5' splice site, maximum dimension",
          "ENGINEERED CONSTRUCT: the U1 snRNA four-helix junction is replaced by a kissing loop "
          "and stem-loop 4 is shortened, and U1-A is absent, so this is smaller than a native "
          "U1 snRNP. Do not quote it as the size of U1 snRNP.",
          "Largest distance between any two backbone trace atoms of the exported copy of "
          "PDB 4PJO",
          "4PJO", bbox_nm=u1["meta"]["dimensions_nm"]["u1_bbox_nm"],
          measured_by="pipeline/structures_txn/build_txn_structures.py"),
        c("human_intron_length_median", gf["introns"]["median_bp"], "bp",
          "Median human intron length",
          "Measured genome-wide from the GENCODE annotation, one Ensembl_canonical transcript "
          f"per protein-coding gene ({gf['n_canonical_protein_coding_transcripts']} transcripts, "
          f"{gf['introns']['n']} introns). The distribution is extremely skewed: mean "
          f"{gf['introns']['mean_bp']} bp, IQR {gf['introns']['p25_bp']}-"
          f"{gf['introns']['p75_bp']} bp, max {gf['introns']['max_bp']} bp. Use the median for "
          "a 'typical' intron and show the IQR if the UI implies a range.",
          gf["method"], "GENCODE",
          mean_bp=gf["introns"]["mean_bp"], p25_bp=gf["introns"]["p25_bp"],
          p75_bp=gf["introns"]["p75_bp"], max_bp=gf["introns"]["max_bp"],
          n=gf["introns"]["n"], annotation=gf["annotation"],
          measured_by="pipeline/structures_txn/gene_features.py"),
        c("human_exon_length_median", gf["exons"]["median_bp"], "bp",
          "Median human exon length",
          "Same annotation and transcript set as the intron constant "
          f"({gf['exons']['n']} exons). Mean {gf['exons']['mean_bp']} bp, IQR "
          f"{gf['exons']['p25_bp']}-{gf['exons']['p75_bp']} bp; the mean is inflated by long "
          "terminal exons carrying 3' UTRs, so the median is the number to show.",
          gf["method"], "GENCODE",
          mean_bp=gf["exons"]["mean_bp"], p25_bp=gf["exons"]["p25_bp"],
          p75_bp=gf["exons"]["p75_bp"], max_bp=gf["exons"]["max_bp"],
          n=gf["exons"]["n"], annotation=gf["annotation"],
          measured_by="pipeline/structures_txn/gene_features.py"),
    ]

    omitted = [
        {"name": "b_form_dna_rise_per_bp_canonical", "wanted": "3.4 A (the textbook value)",
         "why_omitted": "We could not read a citable source for it. Europe PMC's fullTextXML "
                        "endpoint returned HTTP 500 for every article we tried, and the "
                        "publisher pages for the classic references (Wing et al., Nature 1980; "
                        "Wang, PNAS 1979) returned 403 or a reCAPTCHA interstitial. "
                        "Use `dna_duplex_rise_per_bp`, which is measured from PDB 6GMH, and "
                        "label it as measured."},
        {"name": "b_form_dna_bp_per_turn_canonical", "wanted": "10.5 bp/turn in solution",
         "why_omitted": "Same reason. `dna_duplex_bp_per_turn` (10.0, measured from 6GMH) is "
                        "the enzyme-bound duplex, not the solution value; do not relabel it."},
        {"name": "b_form_dna_diameter_vdw", "wanted": "~20 A van-der-Waals width",
         "why_omitted": "Convention-dependent and unsourced. "
                        "`dna_duplex_diameter_phosphate_convention` is given instead."},
        {"name": "transcription_bubble_size_published", "wanted": "the bubble size in nt as "
                                                                  "stated in a primary paper",
         "why_omitted": "Barnes et al. (Mol Cell 2015) and Kettenberger et al. (Mol Cell 2004) "
                        "were both behind 403s. `transcription_bubble_size` is measured from "
                        "PDB 6GMH instead and is labelled as such."},
        {"name": "nascent_rna_nuclease_protected_length", "wanted": "the length of nascent RNA "
                                                                    "protected by Pol II",
         "why_omitted": "No source we could read. The asset reports only what is modelled in "
                        "the exit channel (`nascent_rna_resolved_in_exit_channel`), which is a "
                        "different quantity."},
        {"name": "pol_ii_elongation_rate", "wanted": "nt/s in vivo",
         "why_omitted": "Not needed for scaling the 3D view and not verified here. Do not put "
                        "a speed in the UI from this file."},
    ]

    doc = {
        "format": "utr-race/struct-constants-v1",
        "note": "Physical constants for the transcription + splicing view. Every entry is "
                "MEASURED - from the deposited coordinates of the cited PDB entry, or from the "
                "cited GENCODE release - by the script named in `measured_by`, and can be "
                "regenerated. Nothing here is a remembered textbook value: numbers that would "
                "have needed a literature quote we could not actually read are in `omitted` "
                "with the reason. Companion to constants.json (the translation-side constants "
                "from pipeline/structures/); the two files do not overlap.",
        "generated": GENERATED,
        "generated_by": "pipeline/structures_txn/build_txn_constants.py",
        "constants": constants,
        "omitted": omitted,
    }
    p = os.path.join(OUT, "txn_constants.json")
    json.dump(doc, open(p, "w"), indent=1)
    print(f"wrote {p}  {os.path.getsize(p)/1024:.0f} KB  ({len(constants)} constants, "
          f"{len(omitted)} omitted)")

    # ---------------------------------------------------------- txn_meta.json
    def prov(role, asset, d, why, contents):
        m = d["meta"]
        return {
            "role": role, "asset": asset, "pdb_id": m["pdb_id"], "title": m["title"],
            "method": m["method"], "resolution_A": m["resolution_A"],
            "resolution_method": m["resolution_method"],
            "organism": m["organism"], "organism_note": m["organism_note"],
            "state": m["state"], "contents": contents, "why_chosen": why,
            "citation": m["citation"], "doi": m["doi"], "pubmed": m["pubmed"],
            "url": m["url"], "cif_url": m["cif_url"],
            "paper_url": "https://doi.org/" + m["doi"],
            "caveats": m["caveats"],
        }

    meta = {
        "format": "utr-race/struct-provenance-v1",
        "generated": GENERATED,
        "generated_by": "pipeline/structures_txn/build_txn_constants.py",
        "what": "Provenance for the experimental structures behind the transcription + splicing "
                "view (app/public/struct/polii.json and spliceosome.json). All geometry is "
                "derived from these deposited models; no coordinates were invented or modelled "
                "by us. Show this when the 3D view displays a polymerase or a spliceosome.",
        "license": "PDB data are freely available (wwPDB; no restrictions on use or "
                   "redistribution).",
        "license_url": "https://www.rcsb.org/pages/usage-policy",
        "attribution_line_short": "Mammalian Pol II elongation complex, PDB 6GMH (Vos et al., "
                                  "Nature 2018); human spliceosome C* complex, PDB 5XJC (Zhang "
                                  "et al., Cell 2017); human U1 snRNP with a 5' splice site, "
                                  "PDB 4PJO (Kondo et al., eLife 2015) - PDB data are freely "
                                  "available",
        "asset_format": {
            "name": "utr-race/struct-v1",
            "note": "Same format as app/public/struct/43s.json, produced by a different "
                    "pipeline (pipeline/structures_txn/). Read a file as "
                    "{format, meta, geometry}.",
            "meta": "pdb_id, title, method, resolution_A, resolution_method, organism, "
                    "organism_note, state, citation(+citation_fields), doi, pubmed, url, "
                    "cif_url, license, generated, generated_by, units, origin, "
                    "origin_in_cif_frame_A, orientation, rotation_matrix_row_major, "
                    "parts_mapping{rule, parts{<part>:{chains, n_chains, n_residues_modelled, "
                    "trace_step, entity_descriptions}}, parts_absent}, dimensions_nm, caveats",
            "geometry": "traces[] (one per continuous modelled stretch: chain, part, label, "
                        "atom, step, res_first, res_last, n, xyz flat polyline in Angstroms at "
                        "0.1 A precision), traces_note, surface_spheres{note, stride, count, "
                        "xyzr flat [x,y,z,r] quadruples in Angstroms}, plus one "
                        "structure-specific path block: `transcription_bubble` in polii.json, "
                        "`splice_sites` in spliceosome.json",
            "units": "Angstroms everywhere except keys ending in _nm",
            "orientation_convention": "p_out = R @ (p_cif - origin): translate by the asset's "
                                      "origin, then rotate by rotation_matrix_row_major (3x3, "
                                      "row-major). Translation and rotation only - no scaling, "
                                      "so distances in the files are true physical distances.",
            "extra_key": "spliceosome.json carries `companion_structures`: [ <a full "
                         "utr-race/struct-v1 document for PDB 4PJO> ]. It is a SEPARATE "
                         "structure in its own frame and must not be superposed on the C* "
                         "complex.",
        },
        "parts": {
            "polii.json": ["polii_core", "clamp", "elongation_factors", "template_DNA",
                           "nontemplate_DNA", "nascent_RNA"],
            "spliceosome.json": ["snRNA_U2", "snRNA_U5", "snRNA_U6", "exon", "intron",
                                 "proteins"],
            "spliceosome.json#companion_structures[0]": ["u1_proteins", "snRNA_U1",
                                                         "five_prime_splice_site_RNA"],
        },
        "structures": [
            prov("primary", "polii.json", polii,
                 "A mammalian Pol II ELONGATION complex with a complete 48-bp scaffold: "
                 "template DNA (48 nt modelled), non-template DNA (37 nt) and nascent RNA "
                 "(21 nt) are all resolved, so the transcription bubble, the DNA entry and "
                 "exit duplexes and the RNA exit channel are real geometry rather than a "
                 "drawing. At 3.1 A it is the highest-resolution mammalian elongation complex "
                 "with a scaffold this long, and it is the activated (post pause-release) "
                 "state, which is the state a transcribing polymerase is actually in.",
                 "Pol II (12 subunits), DSIF (SPT4/SPT5), PAF1 complex (PAF1, LEO1, CTR9, "
                 "WDR61, CDC73), SPT6, template DNA, non-template DNA, nascent RNA"),
            prov("primary", "spliceosome.json", sp,
                 "A HUMAN spliceosome (no yeast substitution needed) in the C* state: after "
                 "branching, before exon ligation. The 5' exon, the intron lariat and the "
                 "branch-point adenosine are all modelled, and the 2'-5' branch can be "
                 "verified directly in the coordinates (O2' of A144 is 1.63 A from the "
                 "5'-terminal phosphorus of intron nucleotide +1), so the branch point in the "
                 "asset is measured, not assumed. U2, U5 and U6 snRNA are all present.",
                 "PRP8, BRR2/SNRNP200, AQR, SNRNP200, RBM22, SLU7, PRP17, PRP22/DHX8, the "
                 "exon junction complex, Sm cores, 46 protein chains in total, plus U2, U5 and "
                 "U6 snRNA and the pre-mRNA (18 nt of 5' exon + 66 nt of intron modelled)"),
            prov("companion", "spliceosome.json#companion_structures[0]", u1,
                 "The only step the C* complex cannot show is how a 5' splice site is first "
                 "recognised. 4PJO is a human U1 snRNP crystallised with a 10-mer 5' splice "
                 "site RNA (5'-AGGUAAGUCC-3') base-paired to the 5' end of U1 snRNA, which is "
                 "exactly that step. Included as a companion, in its own frame.",
                 "U1-70K, U1-C, the seven Sm proteins, U1 snRNA (engineered minimal "
                 "construct) and the bound 5' splice site RNA"),
        ],
        "rejected_alternatives": [
            {"pdb_id": "9EHZ", "title": "Cryo-EM structure of Human RNA polymerase II "
                                        "Elongation Complex in an Intermediate Translocation "
                                        "State",
             "resolution_A": 2.6, "organism": "Homo sapiens",
             "citation": "Florez Ariza AJ, Lue NZ, et al. Nat Struct Mol Biol (2025).",
             "url": "https://www.rcsb.org/structure/9EHZ",
             "why_rejected": "Genuinely HUMAN Pol II and higher resolution (2.6 A), but the "
                             "scaffold is short and poorly resolved: only 22 nt of template "
                             "DNA, 14 nt of non-template DNA and 11 nt of RNA are modelled, so "
                             "there is almost no upstream or downstream duplex to draw. If the "
                             "UI must say 'human polymerase', swap to this entry and accept a "
                             "much shorter DNA."},
            {"pdb_id": "5FLM", "title": "Structure of transcribing mammalian RNA polymerase II",
             "resolution_A": 3.4, "organism": "Sus scrofa",
             "citation": "Bernecky C, Herzog F, Baumeister W, Plitzko JM, Cramer P. Nature "
                         "(2016).",
             "url": "https://www.rcsb.org/structure/5FLM",
             "why_rejected": "Pol II alone with a 39-bp scaffold - cleaner, but lower "
                             "resolution and less RNA modelled (14 nt) than 6GMH. A good "
                             "fallback if the elongation factors clutter the view; 6GMH's "
                             "`elongation_factors` part can simply be hidden instead."},
            {"pdb_id": "8XSO", "title": "RNA polymerase II elongation complex transcribing "
                                        "genomic DNA extracted from human nuclei",
             "resolution_A": 2.7, "organism": "Homo sapiens",
             "citation": "Kujirai T, Kato J, et al. Nat Commun (2025).",
             "url": "https://www.rcsb.org/structure/8XSO",
             "why_rejected": "Human and high resolution, but a 37-nt scaffold with only 11 nt "
                             "of RNA."},
            {"pdb_id": "6ICZ", "title": "Cryo-EM structure of a human post-catalytic "
                                        "spliceosome (P complex) at 3.0 angstrom",
             "resolution_A": 3.0, "organism": "Homo sapiens",
             "citation": "Zhang X, Zhan X, et al. Cell Res (2019).",
             "url": "https://www.rcsb.org/structure/6ICZ",
             "why_rejected": "Higher resolution than 5XJC, but it is the P (post-catalytic) "
                             "state - the exons are already ligated, so it cannot show "
                             "splicing about to happen."},
            {"pdb_id": "5MQF", "title": "Cryo-EM structure of a human spliceosome activated "
                                        "for step 2 of splicing (C* complex)",
             "resolution_A": 5.9, "organism": "Homo sapiens",
             "citation": "Bertram K, Agafonov DE, et al. Nature (2017).",
             "url": "https://www.rcsb.org/structure/5MQF",
             "why_rejected": "Same state as 5XJC but at 5.9 A."},
        ],
        "unverified": [
            "The `clamp` part of polii.json uses RPB1 residues 1-346, the clamp definition for "
            "the homologous S. cerevisiae Rpb1 (Cramer, Bushnell & Kornberg, Science 292:1863, "
            "2001) applied to this entry's numbering. We did not validate that mapping against "
            "6GMH, and we could not read that paper's text. Treat it as a highlight, not a "
            "domain boundary, and do not print the residue numbers in the UI.",
            "No published transcription-bubble size, DNA:RNA hybrid length, canonical B-DNA "
            "rise/twist/diameter, or nuclease-protected nascent-RNA length is quoted anywhere "
            "in these files: every publisher page and the Europe PMC full-text endpoint failed "
            "(403 / reCAPTCHA / HTTP 500). See txn_constants.json `omitted`.",
            "Resolutions and titles were read from the RCSB REST API "
            "(data.rcsb.org/rest/v1/core/entry/<id>); the organism note for 6GMH comes from "
            "the abstract of the primary paper, read via the Europe PMC search API.",
        ],
        "rerun": "see pipeline/structures_txn/README.md",
    }
    p = os.path.join(OUT, "txn_meta.json")
    json.dump(meta, open(p, "w"), indent=1)
    print(f"wrote {p}  {os.path.getsize(p)/1024:.0f} KB")


if __name__ == "__main__":
    main()
