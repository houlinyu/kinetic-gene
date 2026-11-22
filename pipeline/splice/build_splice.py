#!/usr/bin/env python
"""Build data/splice_data.json (contract: data/SPEC_SPLICE.md) and copy it to app/public/.

Gene models: Ensembl REST (GRCh38), same server/caching style as pipeline/build_data.py.
Splice-site strength: MaxEntScan (Yeo & Burge 2004) via the `maxentpy` implementation.
Measured splicing: ASCOT GTEx cassette-exon PSI table (Ling et al. 2020) + VastDB (Tapial et al. 2017).

Usage:
    python pipeline/splice/build_splice.py              # all genes in genes_splice.json
    python pipeline/splice/build_splice.py --quick      # only the first 2 genes (schema-correct stub)
    python pipeline/splice/build_splice.py --no-cache   # refetch from Ensembl
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import shutil
import sys
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT_DATA = os.path.join(ROOT, "data", "splice_data.json")
OUT_APP = os.path.join(ROOT, "app", "public", "splice_data.json")
UTR_DATA = os.path.join(ROOT, "data", "utr_data.json")
CACHE = os.path.join(HERE, "cache")
GENES_CFG = os.path.join(HERE, "genes_splice.json")
ASCOT_PSI = os.path.join(HERE, "external", "gtex_psi.csv.gz")

SERVER = "https://rest.ensembl.org"
HEADERS = {"User-Agent": "utr-race-pipeline/0.1 (python-requests)",
           "Content-Type": "application/json", "Accept": "application/json"}
USE_CACHE = True

# --------------------------------------------------------------------------- sources
SOURCES = [
    {"key": "ensembl", "citation": "Ensembl 116 (GRCh38.p14) gene/transcript annotation, Ensembl REST API.",
     "url": "https://rest.ensembl.org"},
    {"key": "yeo2004", "citation": "Yeo G, Burge CB. Maximum entropy modeling of short sequence motifs with "
     "applications to RNA splicing signals. Journal of Computational Biology 2004;11:377-394.",
     "url": "https://doi.org/10.1089/1066527041410418"},
    {"key": "maxentpy", "citation": "maxentpy (X. Dong), a Python implementation of MaxEntScan score5ss/score3ss "
     "shipping the original Yeo & Burge maximum-entropy matrices; v0.0.4.",
     "url": "https://github.com/kepbod/maxentpy"},
    {"key": "ling2020", "citation": "Ling JP, Wilks C, Charles R, Leavey PJ, Ghosh D, Jiang L, Santiago CP, Pang B, "
     "Venkataraman A, Clark BS, Nellore A, Langmead B, Blackshaw S. ASCOT identifies key regulators of neuronal "
     "subtype-specific splicing. Nature Communications 2020;11:137.",
     "url": "https://doi.org/10.1038/s41467-019-14020-5"},
    {"key": "tapial2017", "citation": "Tapial J, Ha KCH, Sterne-Weiler T, et al. An atlas of alternative splicing "
     "profiles and functional associations reveals new regulatory programs and genes that simultaneously express "
     "multiple major isoforms. Genome Research 2017;27:1759-1768 (VastDB).",
     "url": "https://vastdb.crg.eu"},
    {"key": "goedert1989", "citation": "Goedert M, Spillantini MG, Jakes R, Rutherford D, Crowther RA. Multiple "
     "isoforms of human microtubule-associated protein tau: sequences and localization in neurofibrillary tangles "
     "of Alzheimer's disease. Neuron 1989;3:519-526.", "url": "https://doi.org/10.1016/0896-6273(89)90210-9"},
    {"key": "arikan2002", "citation": "Arikan MC, Memmott J, Broderick JA, Lafyatis R, Screaton G, Stamm S, "
     "Andreadis A. Modulation of the membrane-binding projection domain of tau protein: splicing regulation of "
     "exon 3. Brain Res Mol Brain Res 2002;101:109-121.", "url": "https://doi.org/10.1016/S0169-328X(02)00178-X"},
    {"key": "andreadis1992", "citation": "Andreadis A, Brown WM, Kosik KS. Structure and novel exons of the human "
     "tau gene. Biochemistry 1992;31:10626-10633.", "url": "https://doi.org/10.1021/bi00158a027"},
    {"key": "quelle1995", "citation": "Quelle DE, Zindy F, Ashmun RA, Sherr CJ. Alternative reading frames of the "
     "INK4a tumor suppressor gene encode two unrelated proteins capable of inducing cell cycle arrest. Cell "
     "1995;83:993-1000 (mouse Ink4a).", "url": "https://doi.org/10.1016/0092-8674(95)90214-7"},
    {"key": "stott1998", "citation": "Stott FJ, Bates S, James MC, et al. The alternative product from the human "
     "CDKN2A locus, p14(ARF), participates in a regulatory feedback loop with p53 and MDM2. EMBO J "
     "1998;17:5001-5014.", "url": "https://doi.org/10.1093/emboj/17.17.5001"},
    {"key": "yang1998", "citation": "Yang A, Kaghad M, Wang Y, et al. p63, a p53 homolog at 3q27-29, encodes "
     "multiple products with transactivating, death-inducing, and dominant-negative activities. Mol Cell "
     "1998;2:305-316.", "url": "https://doi.org/10.1016/S1097-2765(00)80275-0"},
    {"key": "bourdon2005", "citation": "Bourdon JC, Fernandes K, Murray-Zmijewski F, Liu G, Diot A, Xirodimas DP, "
     "Saville MK, Lane DP. p53 isoforms can regulate p53 transcriptional activity. Genes Dev 2005;19:2122-2137.",
     "url": "https://doi.org/10.1101/gad.1339905"},
    {"key": "kwong1998", "citation": "Kwong KY, Hung MC. A novel splice variant of HER2 with increased "
     "transformation activity. Mol Carcinog 1998;23:62-68.",
     "url": "https://doi.org/10.1002/(SICI)1098-2744(199810)23:2<62::AID-MC2>3.0.CO;2-O"},
    {"key": "castiglioni2006", "citation": "Castiglioni F, Tagliabue E, Campiglio M, Pupa SM, Balsari A, Menard S. "
     "Role of exon-16-deleted HER2 in breast carcinomas. Endocr Relat Cancer 2006;13:221-232.",
     "url": "https://doi.org/10.1677/erc.1.01047"},
    {"key": "caffrey2006", "citation": "Caffrey TM, Joachim C, Paracchini S, Esiri MM, Wade-Martins R. "
     "Haplotype-specific expression of exon 10 at the human MAPT locus. Hum Mol Genet 2006;15:3529-3537.",
     "url": "https://doi.org/10.1093/hmg/ddl429"},
    {"key": "landers1997", "citation": "Landers JE, Cassel SL, George DL. Translational enhancement of mdm2 "
     "oncogene expression in human tumor cells containing a stabilized wild-type p53 protein. Cancer Res "
     "1997;57:3562-3568.", "url": "https://pubmed.ncbi.nlm.nih.gov/9270029/"},
    {"key": "tanahashi2001", "citation": "Tanahashi H, Tabira T. Three novel alternatively spliced isoforms of the "
     "human beta-site amyloid precursor protein cleaving enzyme (BACE) and their effect on amyloid beta-peptide "
     "production. Neurosci Lett 2001;307:9-12.", "url": "https://doi.org/10.1016/S0304-3940(01)01912-7"},
]

SS_METHOD = ("MaxEntScan maximum-entropy splice-site models (Yeo & Burge 2004, yeo2004), computed with the "
             "maxentpy implementation v0.0.4 (maxentpy), which ships the original matrices. score5 on the 9-mer "
             "-3..+6 around the donor (3 exonic + 6 intronic nt, GT at +1+2); score3 on the 23-mer -20..+3 around "
             "the acceptor (20 intronic + 3 exonic nt, AG at -2-1). Units: bits (log2 odds vs a 1st-order "
             "background), higher = stronger; typical real sites ~0-12 (5') and ~0-14 (3'). Validation: maxentpy "
             "reproduces the documented reference examples exactly - score5('cagGTAAGT') = 10.858313101356437 and "
             "score3('ttccaaacgaacttttgtAGgga') = 2.8867730651152104 (maxentpy README, "
             "https://github.com/kepbod/maxentpy); the Burge-lab web server (hollywood.mit.edu/burgelab/maxent) was "
             "unreachable from this machine, so the reference values are the implementation's published examples "
             "rather than a fresh run of the original Perl. Sites whose window runs off the fetched region, or "
             "contains a non-ACGT base, get score null. Branch points are not scored (bp.pos/bp.score always null).")

PSI_METHOD = ("Measured PSI comes from the ASCOT GTEx cassette-exon table (ling2020; "
              "snaptron.cs.jhu.edu/data/ascot/gtex_psi.csv.gz), joined to our exons by exact GRCh38 coordinate "
              "(chrN:start-end, 1-based inclusive) and gene symbol; ASCOT percentages are divided by 100 so value "
              "is a fraction 0-1, and -1 (not measured) is dropped. ASCOT's PSI file carries no per-tissue sample "
              "count, so n is null (counts are in ASCOT_accession_metadata.xlsx, not parsed here). VastDB "
              "(tapial2017) PSI averages are whole-compendium averages, not one tissue, and are recorded with "
              "tissue 'VastDB human compendium (average)' and n null. Alternative-first-exon / alternative-promoter "
              "events (CDKN2A, TP63, TP53 delta133, MDM2 P1/P2) are not cassette exons and are in neither table: "
              "their psi and isoform_fraction are null, not estimated. No value here is computed or inferred by us.")

EXON_TYPE_METHOD = ("Exon types are derived by comparing the transcripts shipped for that gene (not a curated "
                    "annotation): 'first'/'last' = the 5'-most / 3'-most exon of at least one transcript; "
                    "'alt5' = shares its 5' boundary with another exon that has a different 3' (donor) boundary; "
                    "'alt3' = shares its 3' boundary with another exon that has a different 5' (acceptor) boundary; "
                    "'cassette' = internal exon absent from at least one transcript whose span covers it; "
                    "'constitutive' = internal exon present in every transcript whose span covers it. "
                    "first/last win over the others, then alt5/alt3, then cassette.")

# Documented exon names, keyed by GRCh38 coordinate so they do not depend on the shipped transcript set.
# Each is checked against the built model by validate() (coordinate must exist, type must match).
EXON_ALIASES = {
    "MAPT": {
        "chr17:45971859-45971945": {"alias": "exon 2", "cite": "andreadis1992", "expect": "cassette",
                                    "note": "87 nt = the 29-aa N-terminal insert (goedert1989)"},
        "chr17:45974385-45974471": {"alias": "exon 3", "cite": "andreadis1992", "expect": "cassette",
                                    "note": "87 nt = the second 29-aa N-terminal insert; never included without "
                                            "exon 2 (arikan2002)"},
        "chr17:46010310-46010402": {"alias": "exon 10", "cite": "ensembl", "expect": "cassette",
                                    "note": "93 nt microtubule-binding-repeat cassette: present in the three 4R "
                                            "transcripts shipped here and absent from the 3R one, which is what "
                                            "defines 3R vs 4R tau"}},
    "ERBB2": {
        "chr17:39719787-39719834": {"alias": "exon 16", "cite": "kwong1998", "expect": "cassette",
                                    "note": "48 nt, in frame, the 16th exon of MANE ENST00000269571; skipping it "
                                            "gives Delta16HER2 (kwong1998)"}},
    "CDKN2A": {
        "chr9:21974678-21974857": {"alias": "exon 1alpha", "cite": "stott1998", "expect": "first",
                                   "note": "first exon of the p16INK4a transcript"},
        "chr9:21994139-21994392": {"alias": "exon 1beta", "cite": "stott1998", "expect": "first",
                                   "note": "first exon of the p14ARF transcript, ~19.3 kb further 5'"}},
    "TP53": {
        "chr17:7675994-7676272": {"alias": "exon 4", "cite": "bourdon2005", "expect": "cassette",
                                  "note": "the internal promoter that makes the delta133p53 variants lies in the "
                                          "intron immediately 3' of this exon (bourdon2005)"}},
}

# --------------------------------------------------------------------------- Ensembl REST
_last_call = [0.0]


def _get(path, params=None):
    key = hashlib.md5(json.dumps(["GET", path, params], sort_keys=True).encode()).hexdigest()
    cpath = os.path.join(CACHE, key + ".json")
    if USE_CACHE and os.path.exists(cpath):
        with open(cpath) as fh:
            return json.load(fh)
    for attempt in range(6):
        wait = 0.12 - (time.time() - _last_call[0])
        if wait > 0:
            time.sleep(wait)
        _last_call[0] = time.time()
        r = requests.get(SERVER + path, params=params, headers=HEADERS, timeout=120)
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
    return f"Ensembl {rel[0]} (GRCh38)"


# --------------------------------------------------------------------------- sequence helpers
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def maxent_scores():
    from maxentpy import maxent
    return maxent.score5, maxent.score3


def validate_maxent():
    """Reproduce the implementation's documented reference examples before using it."""
    score5, score3 = maxent_scores()
    exp5, exp3 = 10.858313101356437, 2.8867730651152104
    got5, got3 = score5("cagGTAAGT"), score3("ttccaaacgaacttttgtAGgga")
    ok = abs(got5 - exp5) < 1e-9 and abs(got3 - exp3) < 1e-9
    print(f"  MaxEntScan validation: score5('cagGTAAGT')={got5!r} (expect {exp5!r}); "
          f"score3('ttccaaacgaacttttgtAGgga')={got3!r} (expect {exp3!r}) -> {'OK' if ok else 'MISMATCH'}")
    if not ok:
        raise SystemExit("MaxEntScan implementation does not reproduce its reference examples - refusing to "
                         "label scores 'MaxEntScan'. Fall back to a PWM and change meta.methods.ss_score.")
    return score5, score3


# --------------------------------------------------------------------------- gene model
def transcript_model(enst):
    """Ensembl exon list for one transcript -> dict."""
    t = _get(f"/lookup/id/{enst}", {"expand": 1, "utr": 1})
    exons = sorted(({"start": e["start"], "end": e["end"]} for e in t["Exon"]), key=lambda e: e["start"])
    return {"id": f"{t['id']}.{t['version']}", "base": t["id"], "name": t.get("display_name") or t["id"],
            "chrom": t["seq_region_name"], "strand": 1 if t["strand"] == 1 else -1,
            "biotype": t.get("biotype"), "exons": exons,
            "gstart": min(e["start"] for e in exons), "gend": max(e["end"] for e in exons)}


def genomic_seq(chrom, start, end):
    d = _get(f"/sequence/region/human/{chrom}:{start}..{end}:1")
    return d["seq"].upper()


def build_gene(cfg, tab1_ids, ascot, score5, score3):
    gene = cfg["gene"]
    ids = list(cfg.get("_ids", []))
    tms = []
    for enst in ids:
        try:
            tms.append(transcript_model(enst))
        except Exception as exc:  # a transcript may be retired in the current release
            print(f"    ! {gene} {enst}: {exc}")
    if not tms:
        return None
    strand = tms[0]["strand"]
    chrom = tms[0]["chrom"]
    tms = [t for t in tms if t["strand"] == strand and t["chrom"] == chrom]
    gstart = min(t["gstart"] for t in tms)
    gend = max(t["gend"] for t in tms)

    # genomic sequence of the whole locus, + 60 nt of slack for the score windows
    pad = 60
    seq = genomic_seq(chrom, gstart - pad, gend + pad)          # forward strand, 1-based gstart-pad
    seq_off = gstart - pad                                      # genomic 1-based coord of seq[0]

    def fwd(a, b):
        """forward-strand sequence for 1-based inclusive [a,b]"""
        i, j = a - seq_off, b - seq_off + 1
        if i < 0 or j > len(seq):
            return None
        return seq[i:j]

    # ---- transcript-orientation axis: 0 = 5'-most base of the locus
    def ax(g):
        return (g - gstart) if strand == 1 else (gend - g)

    def exon_axis(e):
        """(start,end) half-open on the transcript-orientation axis"""
        return (ax(e["start"]), ax(e["end"]) + 1) if strand == 1 else (ax(e["end"]), ax(e["start"]) + 1)

    # ---- exon table (union over transcripts), keyed by axis coords
    exon_tx = {}
    for t in tms:
        for e in t["exons"]:
            exon_tx.setdefault(exon_axis(e), []).append(t["id"])
    keys = sorted(exon_tx)
    eid = {k: f"E{i + 1}" for i, k in enumerate(keys)}

    # transcript exon chains in transcript order
    chains = {}
    for t in tms:
        chains[t["id"]] = [eid[exon_axis(e)] for e in sorted(
            t["exons"], key=lambda e: e["start"] if strand == 1 else -e["start"])]

    firsts = {c[0] for c in chains.values()}
    lasts = {c[-1] for c in chains.values()}
    starts = {}
    ends = {}
    for k in keys:
        starts.setdefault(k[0], set()).add(k[1])
        ends.setdefault(k[1], set()).add(k[0])

    exons_out = []
    for k in keys:
        s, e = k
        name = eid[k]
        txs = sorted(set(exon_tx[k]))
        if name in firsts:
            typ = "first"
        elif name in lasts:
            typ = "last"
        elif len(starts[s]) > 1:
            typ = "alt5"
        elif len(ends[e]) > 1:
            typ = "alt3"
        else:
            spanning = [t["id"] for t in tms if ax_span(t, ax, strand)[0] <= s and ax_span(t, ax, strand)[1] >= e]
            typ = "constitutive" if set(spanning) <= set(txs) else "cassette"
        exons_out.append({"id": name, "start": s, "end": e, "type": typ, "transcripts": txs,
                          "len": e - s, "g": f"chr{chrom}:{min(k2g(s, gstart, gend, strand), k2g(e - 1, gstart, gend, strand))}-"
                                              f"{max(k2g(s, gstart, gend, strand), k2g(e - 1, gstart, gend, strand))}"})

    # ---- introns (union of adjacent-exon pairs over transcripts)
    seen = {}
    for t in tms:
        ch = chains[t["id"]]
        pos = {x["id"]: x for x in exons_out}
        for a, b in zip(ch, ch[1:]):
            key = (pos[a]["end"], pos[b]["start"])
            seen.setdefault(key, []).append(t["id"])
    introns_out = []
    for i, key in enumerate(sorted(seen)):
        dpos, apos = key
        if strand == 1:
            dg, ag = gstart + dpos, gstart + apos - 1          # 1-based: first intron base / last intron base
            dseq = fwd(dg - 3, dg + 5)
            aseq = fwd(ag - 19, ag + 3)
        else:
            dg, ag = gend - dpos, gend - apos + 1
            dseq = fwd(dg - 5, dg + 3)
            dseq = revcomp(dseq) if dseq else None
            aseq = fwd(ag - 3, ag + 19)
            aseq = revcomp(aseq) if aseq else None
        introns_out.append({
            "id": f"I{i + 1}", "start": dpos, "end": apos, "transcripts": sorted(set(seen[key])),
            "donor": {"pos": dpos, "seq": dseq, "score": ss(dseq, score5, 9), "scale": "MaxEntScan 5'"},
            "acceptor": {"pos": apos, "seq": aseq, "score": ss(aseq, score3, 23), "scale": "MaxEntScan 3'"},
            "bp": {"pos": None, "score": None}})

    # ---- TSSs (distinct first-exon starts)
    tss_pos = {}
    for t in tms:
        p = exon_axis(sorted(t["exons"], key=lambda e: e["start"] if strand == 1 else -e["start"])[0])[0]
        tss_pos.setdefault(p, []).append(t["id"])
    tss_out = []
    tss_of = {}
    for i, p in enumerate(sorted(tss_pos)):
        tid = f"T{i + 1}"
        for t in tss_pos[p]:
            tss_of[t] = tid
        tss_out.append({"id": tid, "pos": p,
                        "support": "Ensembl-annotated transcript start of " + ", ".join(sorted(tss_pos[p])),
                        "cite": "ensembl"})

    # ---- isoforms
    names = cfg.get("iso_names", {})
    isos = []
    for t in tms:
        tab1 = tab1_ids.get(t["base"])
        isos.append({"id": tab1["id"] if tab1 else t["id"],
                     "name": names.get(t["base"]) or (tab1["name"] if tab1 else t["name"]),
                     "ensembl_name": t["name"], "tss_id": tss_of[t["id"]],
                     "exon_ids": chains[t["id"]], "utr5_id": tab1["id"] if tab1 else None})

    # documented exon names (coordinate-keyed, so independent of the shipped transcript set)
    for e in exons_out:
        a = EXON_ALIASES.get(gene, {}).get(e["g"])
        if a:
            e["alias"] = a["alias"]
            e["alias_cite"] = a["cite"]
            e["alias_note"] = a["note"]

    g = {"gene": gene, "chrom": chrom, "strand": "+" if strand == 1 else "-",
         "locus": f"{chrom}:{gstart}-{gend}", "span": gend - gstart + 1,
         "tss": tss_out, "exons": exons_out, "introns": introns_out, "isoforms": isos,
         "measured": {"psi": measured_psi(gene, exons_out, ascot), "isoform_fraction": []}}
    if cfg.get("event"):
        g["event"] = cfg["event"]
        g["cite"] = cfg.get("cite", [])
    return g


def ax_span(t, ax, strand):
    a, b = ax(t["gstart"]), ax(t["gend"])
    return (min(a, b), max(a, b) + 1)


def k2g(axis_pos, gstart, gend, strand):
    return gstart + axis_pos if strand == 1 else gend - axis_pos


def ss(seq, fn, n):
    if not seq or len(seq) != n or any(c not in "ACGT" for c in seq):
        return None
    try:
        return round(float(fn(seq)), 3)
    except Exception:
        return None


# --------------------------------------------------------------------------- measured PSI
def load_ascot():
    """ASCOT GTEx cassette-exon PSI -> {(gene_symbol, 'chrN:start-end'): {tissue: psi_percent}}"""
    if not os.path.exists(ASCOT_PSI):
        print(f"  ! ASCOT table missing ({ASCOT_PSI}); measured.psi will be empty. See pipeline/splice/README.md")
        return {}
    out = {}
    with gzip.open(ASCOT_PSI, "rt") as fh:
        hdr = fh.readline().lstrip("﻿").rstrip("\n").split(",")
        tissues = hdr[12:]
        isym, iloc = hdr.index("gene_symbol"), hdr.index("exon_location")
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) != len(hdr):
                continue
            vals = {}
            for t, v in zip(tissues, f[12:]):
                try:
                    x = float(v)
                except ValueError:
                    continue
                if x >= 0:
                    vals[t] = x
            if vals:
                out[(f[isym], f[iloc])] = vals
    print(f"  ASCOT GTEx PSI table: {len(out)} exons x {len(tissues)} tissues")
    return out


def measured_psi(gene, exons, ascot):
    psi = []
    for e in exons:
        vals = ascot.get((gene, e["g"]))
        if not vals:
            continue
        label = f"{e['id']} ({gene} {e['alias']})" if e.get("alias") else e["id"]
        for tissue, pct in sorted(vals.items()):
            psi.append({"event": f"{label} inclusion", "exon_id": e["id"], "value": round(pct / 100.0, 4),
                        "tissue": f"GTEx {tissue}", "n": None, "cite": "ling2020"})
    return psi


# --------------------------------------------------------------------------- hand-curated extras
def vastdb_psi():
    """VastDB whole-compendium PSI averages read off the hg38 gene page (see README). Gene -> entries."""
    return {"MAPT": [
        {"event": "MAPT exon 2 inclusion (VastDB HsaEX0037852)", "value": 0.4359,
         "tissue": "VastDB human compendium (average)", "n": None, "cite": "tapial2017",
         "note": "PSI Average 43.59 on https://vastdb.crg.eu/gene/ENSG00000186868@hg38 (PSI Range 100)"},
        {"event": "MAPT exon 3 inclusion (VastDB HsaEX0037853)", "value": 0.0620,
         "tissue": "VastDB human compendium (average)", "n": None, "cite": "tapial2017",
         "note": "PSI Average 6.20 on https://vastdb.crg.eu/gene/ENSG00000186868@hg38 (PSI Range 34.69)"}]}


def isoform_fractions():
    """Measured isoform fractions from primary papers. gene -> [entry]; ids filled in by match_label()."""
    return {"ERBB2": [
        {"id": None, "label": "ERBB2 exon 16 skipped (Delta16HER2)", "value": 0.09,
         "tissue": "human breast carcinoma (HER2 3+ to 0 by HercepTest)", "n": 46, "cite": "castiglioni2006",
         "note": "\"the DeltaHER2 splice variant represents about 9% of the HER2 mRNA\" (Castiglioni et al. 2006); "
                 "the complementary ~91% is full-length HER2. Not mapped to a shipped transcript id: the "
                 "exon-16-skipped form is not among the Ensembl transcripts in this file."}]}


def splice_findings(genes):
    """Literature 'which form is used more' pairs. higher/lower = shipped isoform id, or an event label."""
    by_gene = {g["gene"]: g for g in genes}

    def iso(gene, base):
        g = by_gene.get(gene)
        if not g:
            return None
        for i in g["isoforms"]:
            if i["id"].split(".")[0] == base:
                return i["id"]
        return None

    def ex(gene, alias):
        g = by_gene.get(gene)
        for e in (g["exons"] if g else []):
            if e.get("alias") == alias:
                return e["id"]
        return None

    def psi_median(gene, exon_id):
        g = by_gene.get(gene)
        v = sorted(p["value"] for p in (g["measured"]["psi"] if g else []) if p.get("exon_id") == exon_id)
        return round(v[len(v) // 2], 3) if v else None

    out = []
    tp63_dn, tp63_ta = iso("TP63", "ENST00000354600"), iso("TP63", "ENST00000264731")
    if tp63_dn and tp63_ta:
        out.append({"id": "tp63_dn_predominant", "gene": "TP63", "cite": "yang1998",
                    "title": "deltaNp63 is the predominant p63 form in many epithelial tissues",
                    "higher": tp63_dn, "lower": tp63_ta,
                    "evidence": "Yang et al. 1998: \"the predominant p63 isotypes in many epithelial tissues lack "
                                "an acidic N terminus corresponding to the transactivation domain of p53\". "
                                "deltaNp63alpha = ENST00000354600 (MANE Plus Clinical), TAp63alpha = "
                                "ENST00000264731 (MANE Select); no percentage is given in the paper.",
                    "value": None})
    if "MAPT" in by_gene:
        out.append({"id": "mapt_ex2_gt_ex3", "gene": "MAPT", "cite": "arikan2002",
                    "title": "MAPT exon 2 is included much more often than exon 3",
                    "higher": "MAPT exon 2 inclusion", "lower": "MAPT exon 3 inclusion",
                    "evidence": "Arikan et al. 2002: the two exons are alternatively spliced cassettes \"in which "
                                "exon 3 never appears independently of exon 2\", and \"exon 3 is intrinsically "
                                "suboptimal and its primary regulator is a weak branch point\". Consistent with "
                                "VastDB compendium PSI averages 43.59 (exon 2) vs 6.20 (exon 3) (tapial2017).",
                    "value": None,
                    "higher_exon_id": ex("MAPT", "exon 2"), "lower_exon_id": ex("MAPT", "exon 3")})
        out.append({"id": "mapt_ex10_h1_h2", "gene": "MAPT", "cite": "caffrey2006",
                    "title": "MAPT exon 10 inclusion is higher from the H1 than the H2 haplotype",
                    "higher": "MAPT exon 10 included (H1 haplotype)",
                    "lower": "MAPT exon 10 included (H2 haplotype)",
                    "evidence": "Caffrey et al. 2006, allele-specific expression in 14 H1/H2 heterozygous "
                                "post-mortem control brains: exon 10+ MAPT H1:H2 transcript ratio 1.425 (SD 0.205, "
                                "P < 0.0001) in globus pallidus and 1.291 (SD 0.315, P = 0.006) in frontal cortex; "
                                "total MAPT did not differ. This is a haplotype contrast, not two transcript ids.",
                    "value": 1.425, "higher_exon_id": ex("MAPT", "exon 10"),
                    "lower_exon_id": ex("MAPT", "exon 10")})
    if "ERBB2" in by_gene:
        out.append({"id": "erbb2_fl_gt_d16", "gene": "ERBB2", "cite": "castiglioni2006",
                    "title": "Full-length HER2 mRNA greatly exceeds the exon-16-skipped Delta16HER2 form",
                    "higher": iso("ERBB2", "ENST00000269571") or "ERBB2 full-length (exon 16 included)",
                    "lower": "ERBB2 exon 16 skipped (Delta16HER2)",
                    "evidence": "Castiglioni et al. 2006, 46 human breast carcinomas: \"the DeltaHER2 splice "
                                "variant represents about 9% of the HER2 mRNA\" (first described by Kwong & Hung "
                                "1998, kwong1998, as an in-frame deletion of a 16-codon extracellular-domain "
                                "exon). The exon-16-skipped transcript is not annotated in the Ensembl set shipped "
                                "here, so 'lower' is an event label; the exon itself is "
                                + str(ex("ERBB2", "exon 16")) + ", the 16th exon of ENST00000269571.",
                    "value": 0.09, "higher_exon_id": ex("ERBB2", "exon 16"),
                    "lower_exon_id": ex("ERBB2", "exon 16"),
                    "crosscheck": ("Independent of that paper: the ASCOT GTEx median PSI for this exon across 56 "
                                   "tissues is " + str(psi_median("ERBB2", ex("ERBB2", "exon 16"))) + " (breast "
                                   "0.97), i.e. a few percent skipping in normal tissue - the same direction as "
                                   "the ~9% Delta16 fraction reported in carcinomas."),
                    "crosscheck_cite": "ling2020"})
    if "MDM2" in by_gene:
        out.append({"id": "mdm2_p2_preferred", "gene": "MDM2", "cite": "landers1997",
                    "title": "Tumour cells preferentially use the internal (P2) MDM2 promoter, giving a distinct "
                             "5' UTR",
                    "higher": "MDM2 internal-promoter (P2) transcript",
                    "lower": "MDM2 upstream-promoter (P1) transcript",
                    "evidence": "Landers et al. 1997: in human tumour cells with stabilized wild-type p53 there is "
                                "\"a preferential increase in mdm2 transcription that is initiated from the "
                                "internal p53-responsive promoter region of this gene\", and those transcripts "
                                "\"contain a distinct 5' untranslated region and exhibit a significantly enhanced "
                                "translational efficiency\". No percentage is given; P1/P2 are not mapped to "
                                "single Ensembl ids here, so both sides are event labels.",
                    "value": None})
    return out


# --------------------------------------------------------------------------- validation
def validate(doc):
    errs = []
    src_keys = {s["key"] for s in doc["meta"]["sources"]}

    def cites(obj, where):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k == "cite":
                    for c in ([v] if isinstance(v, str) else (v or [])):
                        if c not in src_keys:
                            errs.append(f"{where}: cite '{c}' not in meta.sources")
                else:
                    cites(v, f"{where}.{k}")
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                cites(v, f"{where}[{i}]")

    cites(doc, "doc")
    for g in doc["genes"]:
        w = g["gene"]
        ex = {e["id"]: e for e in g["exons"]}
        if g["strand"] not in ("+", "-"):
            errs.append(f"{w}: bad strand")
        # exons ordered, non-overlapping-start, inside the locus
        prev = -1
        for e in g["exons"]:
            if not (0 <= e["start"] < e["end"] <= g["span"]):
                errs.append(f"{w} {e['id']}: positions {e['start']}..{e['end']} outside 0..{g['span']}")
            if e["start"] < prev:
                errs.append(f"{w} {e['id']}: exons not in transcript order")
            prev = e["start"]
            if e["type"] not in ("constitutive", "cassette", "alt5", "alt3", "first", "last"):
                errs.append(f"{w} {e['id']}: bad type {e['type']}")
            if not e["transcripts"]:
                errs.append(f"{w} {e['id']}: belongs to no transcript")
        for t in g["tss"]:
            if not (0 <= t["pos"] <= g["span"]):
                errs.append(f"{w} {t['id']}: tss pos out of range")
        tids = {t["id"] for t in g["tss"]}
        # introns: adjacency consistent with some transcript's exon chain
        pairs = set()
        for i in g["isoforms"]:
            if i["tss_id"] not in tids:
                errs.append(f"{w} {i['id']}: tss_id {i['tss_id']} unknown")
            chain = i["exon_ids"]
            if not chain:
                errs.append(f"{w} {i['id']}: empty exon chain")
                continue
            for a, b in zip(chain, chain[1:]):
                if a not in ex or b not in ex:
                    errs.append(f"{w} {i['id']}: unknown exon in chain")
                    continue
                if ex[a]["end"] > ex[b]["start"]:
                    errs.append(f"{w} {i['id']}: exons {a},{b} overlap / out of order")
                pairs.add((ex[a]["end"], ex[b]["start"]))
            for eid_ in chain:
                if eid_ in ex and i["id"] not in ex[eid_]["transcripts"]:
                    errs.append(f"{w} {i['id']}: exon {eid_} does not list this transcript")
        for intr in g["introns"]:
            if (intr["start"], intr["end"]) not in pairs:
                errs.append(f"{w} {intr['id']}: not adjacent in any isoform exon chain")
            if intr["donor"]["pos"] != intr["start"] or intr["acceptor"]["pos"] != intr["end"]:
                errs.append(f"{w} {intr['id']}: donor/acceptor pos != intron bounds")
            if not (0 <= intr["start"] < intr["end"] <= g["span"]):
                errs.append(f"{w} {intr['id']}: out of range")
            ds, as_ = intr["donor"]["seq"], intr["acceptor"]["seq"]
            if ds is not None:
                if len(ds) != 9:
                    errs.append(f"{w} {intr['id']}: donor seq not 9 nt")
                elif ds[3:5] != "GT":
                    print(f"    note {w} {intr['id']}: non-GT donor {ds} (minor/U12 or annotation edge)")
            if as_ is not None:
                if len(as_) != 23:
                    errs.append(f"{w} {intr['id']}: acceptor seq not 23 nt")
                elif as_[18:20] != "AG":
                    print(f"    note {w} {intr['id']}: non-AG acceptor {as_}")
        for p in g["measured"]["psi"]:
            if p["value"] is not None and not (0.0 <= p["value"] <= 1.0):
                errs.append(f"{w}: psi value {p['value']} outside 0-1")
    # every documented exon alias must land on an existing exon of the expected class
    built = {g["gene"]: g for g in doc["genes"]}
    for gene, al in EXON_ALIASES.items():
        if gene not in built:
            continue
        coords = {e["g"]: e for e in built[gene]["exons"]}
        for coord, spec in al.items():
            e = coords.get(coord)
            if not e:
                errs.append(f"{gene}: alias '{spec['alias']}' coordinate {coord} is not an exon in the built model")
            elif e["type"] != spec["expect"]:
                errs.append(f"{gene} {e['id']} ('{spec['alias']}'): type {e['type']} != expected {spec['expect']}")
    for g in doc["genes"]:
        for e in g["exons"]:
            if "alias_cite" in e and e["alias_cite"] not in src_keys:
                errs.append(f"{g['gene']} {e['id']}: alias_cite not in sources")
    for f in doc.get("splice_findings", []):
        for k in ("higher_exon_id", "lower_exon_id"):
            if f.get(k) is None:
                continue
            if f[k] not in {e["id"] for e in built[f["gene"]]["exons"]}:
                errs.append(f"finding {f['id']}: {k} {f[k]} is not an exon of {f['gene']}")

    iso_ids = {i["id"] for g in doc["genes"] for i in g["isoforms"]}
    for f in doc.get("splice_findings", []):
        for side in ("higher", "lower"):
            v = f[side]
            if v in iso_ids or isinstance(v, str):
                continue
            errs.append(f"finding {f['id']}: {side} is not an isoform id or label")
    return errs


# --------------------------------------------------------------------------- main
def main():
    global USE_CACHE
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--quick", action="store_true", help="first 2 genes only (fast schema-correct stub)")
    ap.add_argument("genes", nargs="*")
    args = ap.parse_args()
    USE_CACHE = not args.no_cache

    print("splice_data build")
    score5, score3 = validate_maxent()

    with open(UTR_DATA) as fh:
        utr = json.load(fh)
    tab1 = {}       # unversioned ENST -> {id, name}
    tab1_by_gene = {}
    for g in utr["genes"]:
        tab1_by_gene[g["gene"]] = [i["id"] for i in g["isoforms"]]
        for i in g["isoforms"]:
            tab1[i["id"].split(".")[0]] = {"id": i["id"], "name": i["name"], "gene": g["gene"]}

    cfgs = json.load(open(GENES_CFG))["genes"]
    if args.genes:
        cfgs = [c for c in cfgs if c["gene"] in args.genes]
    if args.quick:
        cfgs = cfgs[:2]
    ascot = load_ascot()
    vast = vastdb_psi()
    fracs = isoform_fractions()

    genes = []
    for cfg in cfgs:
        ids = []
        if cfg.get("from_tab1"):
            ids += [x.split(".")[0] for x in tab1_by_gene.get(cfg["gene"], [])]
        for x in cfg.get("extra", []):
            if x not in ids:
                ids.append(x)
        cfg = dict(cfg, _ids=ids)
        print(f"  {cfg['gene']}: {len(ids)} transcripts")
        g = build_gene(cfg, tab1, ascot, score5, score3)
        if not g:
            continue
        g["measured"]["psi"] += vast.get(g["gene"], [])
        g["measured"]["isoform_fraction"] += fracs.get(g["gene"], [])
        linked = sum(1 for i in g["isoforms"] if i["utr5_id"])
        scored = [x for i in g["introns"] for x in (i["donor"]["score"], i["acceptor"]["score"]) if x is not None]
        types = {}
        for e in g["exons"]:
            types[e["type"]] = types.get(e["type"], 0) + 1
        print(f"    {len(g['exons'])} exons {types}, {len(g['introns'])} introns, "
              f"{len(scored)} scored sites, {len(g['tss'])} TSS, {linked}/{len(g['isoforms'])} linked to tab 1, "
              f"{len(g['measured']['psi'])} measured PSI")
        genes.append(g)

    doc = {"meta": {"generated": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
                    "release": ensembl_release(),
                    "sources": SOURCES,
                    "methods": {"ss_score": SS_METHOD, "psi": PSI_METHOD, "exon_type": EXON_TYPE_METHOD,
                                "coords": "All positions are 0-based half-open on a per-gene pre-mRNA axis in "
                                          "transcript orientation: 0 = the 5'-most annotated base of the locus "
                                          "(gene.locus / gene.span give the GRCh38 span). exon.g is the exon's "
                                          "GRCh38 coordinate (1-based inclusive). intron.start = first intron "
                                          "base = donor.pos; intron.end = first base of the next exon = "
                                          "acceptor.pos.",
                                "join": "isoform.utr5_id is the tab-1 (data/utr_data.json) transcript id when the "
                                        "same versioned ENST is present there, else null; isoform.id equals that "
                                        "id whenever it exists."}},
           "genes": genes,
           "splice_findings": []}
    doc["splice_findings"] = splice_findings(genes)

    errs = validate(doc)
    if errs:
        print("\nVALIDATION FAILED:")
        for e in errs[:40]:
            print("  -", e)
        sys.exit(1)
    with open(OUT_DATA, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    shutil.copyfile(OUT_DATA, OUT_APP)
    sz = os.path.getsize(OUT_DATA)
    print(f"\nvalidate OK. wrote {OUT_DATA} ({sz / 1024:.0f} KB) and {OUT_APP}")
    print(f"  {len(genes)} genes, {sum(len(g['isoforms']) for g in genes)} isoforms, "
          f"{sum(len(g['introns']) for g in genes)} introns, "
          f"{sum(len(g['measured']['psi']) for g in genes)} measured PSI values, "
          f"{len(doc['splice_findings'])} splice_findings")


if __name__ == "__main__":
    main()
