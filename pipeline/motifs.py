"""Sequence motifs and interpretable features of a 5' UTR (pure functions, no I/O).

API (see data/SPEC_v2.md):
    find_motifs(utr, cds_head) -> list[dict]
    features(utr, cds_head, fold=None) -> dict[str, float]

Input may be RNA or DNA (T -> U). Positions are 0-based, half-open, within the 5' UTR
(cap = 0). cds_head = first >= 30 nt of the CDS including the start codon (longer is fine; uORF
stops are searched through it, so pass the whole CDS+3'UTR if you have it).

Every motif type has a cite key; the citations live in SOURCES (merged into meta.sources).
"""
import math
import re

STOPS = ("UAA", "UAG", "UGA")
NEARCOG = ("CUG", "GUG", "ACG", "UUG")
PYR = set("CU")

# cite key -> {citation, url}; verified against the abstract (see SOURCE_NOTES)
SOURCES = {
    "kozak1986": {
        "citation": "Kozak M. Point mutations define a sequence flanking the AUG initiator codon that modulates "
                    "translation by eukaryotic ribosomes. Cell 44:283-292 (1986)",
        "url": "https://doi.org/10.1016/0092-8674(86)90762-2"},
    "calvo2009": {
        "citation": "Calvo SE, Pagliarini DJ, Mootha VK. Upstream open reading frames cause widespread reduction of "
                    "protein expression and are polymorphic among humans. PNAS 106:7507-7512 (2009)",
        "url": "https://doi.org/10.1073/pnas.0810916106"},
    "kearse2017": {
        "citation": "Kearse MG, Wilusz JE. Non-AUG translation: a new start for protein synthesis in eukaryotes. "
                    "Genes Dev 31:1717-1731 (2017)",
        "url": "https://doi.org/10.1101/gad.305250.117"},
    "thoreen2012": {
        "citation": "Thoreen CC et al. A unifying model for mTORC1-mediated regulation of mRNA translation. "
                    "Nature 485:109-113 (2012)",
        "url": "https://doi.org/10.1038/nature11083"},
    "elfakess2008": {
        "citation": "Elfakess R, Dikstein R. A translation initiation element specific to mRNAs with very short "
                    "5'UTR that also regulates transcription. PLoS One 3:e3094 (2008)",
        "url": "https://doi.org/10.1371/journal.pone.0003094"},
    "huppert2005": {
        "citation": "Huppert JL, Balasubramanian S. Prevalence of quadruplexes in the human genome. "
                    "Nucleic Acids Res 33:2908-2916 (2005)",
        "url": "https://doi.org/10.1093/nar/gki609"},
    "kumari2007": {
        "citation": "Kumari S, Bugaut A, Huppert JL, Balasubramanian S. An RNA G-quadruplex in the 5' UTR of the "
                    "NRAS proto-oncogene modulates translation. Nat Chem Biol 3:218-221 (2007)",
        "url": "https://doi.org/10.1038/nchembio864"},
}

# feature name -> human label (ridge agent uses these labels)
FEATURES = {
    "log_len": "log(1 + 5' UTR length)",
    "gc": "GC fraction of the 5' UTR",
    "gc_first50": "GC fraction of the first 50 nt (cap-proximal)",
    "n_uAUG": "upstream AUGs (any context)",
    "n_uAUG_strong": "upstream AUG, strong Kozak context",
    "n_uAUG_adequate": "upstream AUG, adequate Kozak context",
    "n_uAUG_weak": "upstream AUG, weak Kozak context",
    "n_uORF": "uORFs ending inside the 5' UTR",
    "n_uORF_overlap": "uORFs overlapping the main start (out of frame)",
    "n_uORF_inframe_ext": "upstream AUGs in frame with the CDS, no stop (N-terminal extension)",
    "max_uORF_aa": "longest uORF (codons, excl. stop), 0 if none",
    "dist_last_uAUG_to_start": "log(1 + nt from last upstream AUG to main start); log(1+len) if none",
    "n_nearcog_strong": "upstream CUG/GUG/ACG/UUG in strong Kozak context",
    "has_TOP": "5' TOP / TOP-like motif (>= 5 pyrimidines within 4 nt of the 5' end)",
    "n_G4": "putative G-quadruplexes (G3+N1-7 x4)",
    "has_TISU": "TISU element at the main start (UTR <= 30 nt)",
    "main_kozak_strong": "main start: strong Kozak context",
    "main_kozak_adequate": "main start: adequate Kozak context (weak = baseline)",
    "start_minus3_purine": "main start: purine at -3",
    "start_plus4_G": "main start: G at +4",
    "mfe_per_nt": "MFE per nt, local fold (max bp span 100) if given (kcal/mol/nt)",
    "cap_dG": "MFE of the first 50 nt (kcal/mol)",
    "acc_start": "mean unpaired probability, 30 nt around the main start",
}

G4_RE = re.compile(r"(?=((?:G{3,}[ACGU]{1,7}){3,}G{3,}))")
TISU_RE = re.compile(r"^AA[CG]AUGGCGGC$")  # AAS-AUG-GCGGC (Elfakess 2008): -3..+8, A of AUG = +1


def _rna(s):
    return s.upper().replace("T", "U")


def kozak_class(full, a):
    """full = UTR+CDS; a = index of first base of the start codon. -3 purine / +4 G rule."""
    m3 = full[a - 3] if a - 3 >= 0 else "N"
    p4 = full[a + 3] if a + 3 < len(full) else "N"
    return ["weak", "adequate", "strong"][(m3 in "AG") + (p4 == "G")]


def _orf_end(full, i):
    for j in range(i + 3, len(full) - 2, 3):
        if full[j:j + 3] in STOPS:
            return j
    return None


def _uaugs(utr, cds):
    """Each upstream AUG with stop position and class (same rules as build_data.find_uorfs)."""
    full, L = utr + cds, len(utr)
    out = []
    for i in range(L - 2):
        if utr[i:i + 3] != "AUG":
            continue
        stop = _orf_end(full, i)
        inframe = (L - i) % 3 == 0
        if stop is None:
            typ, aa = ("in-frame-extension" if inframe else "overlapping"), None
        else:
            aa = (stop - i) // 3
            if inframe and stop >= L:
                typ = "in-frame-extension"
            elif stop + 3 <= L:
                typ = "uORF"
            else:
                typ = "overlapping"
        out.append({"start": i, "stop": stop, "aa": aa, "type": typ, "kozak": kozak_class(full, i)})
    return out


def find_motifs(utr, cds_head):
    utr, cds = _rna(utr), _rna(cds_head)
    full, L = utr + cds, len(utr)
    m = []
    for u in _uaugs(utr, cds):
        i = u["start"]
        m.append({"type": "uAUG", "start": i, "end": i + 3, "seq": utr[i:i + 3],
                  "label": f"uAUG ({u['kozak']}) -> {u['type']}", "strength": u["kozak"], "cite": "calvo2009",
                  "note": f"{u['type']}" + (f", {u['aa']} codons" if u["aa"] is not None else ", no stop in sequence")})
        if u["type"] == "uORF":
            s = u["stop"]
            m.append({"type": "uORF_stop", "start": s, "end": s + 3, "seq": utr[s:s + 3],
                      "label": f"uORF stop (uAUG at {i})", "strength": None, "cite": "calvo2009",
                      "note": f"ends uORF {i}-{s + 3} ({u['aa']} codons)"})
    for i in range(L - 2):
        c = utr[i:i + 3]
        if c in NEARCOG and kozak_class(full, i) == "strong":
            m.append({"type": "nearcog", "start": i, "end": i + 3, "seq": c,
                      "label": f"near-cognate {c} (strong context)", "strength": "strong", "cite": "kearse2017",
                      "note": "non-AUG start candidate; only strong -3R/+4G context reported"})
    # TOP / TOP-like (Thoreen 2012): a stretch of >= 5 pyrimidines starting within 4 nt of the TSS
    for s0 in range(min(5, L)):
        k = s0
        while k < L and utr[k] in PYR:
            k += 1
        if k - s0 >= 5:
            m.append({"type": "TOP", "start": s0, "end": k, "seq": utr[s0:k], "label": "5' TOP / TOP-like motif",
                      "strength": None, "cite": "thoreen2012",
                      "note": ">= 5 pyrimidines within 4 nt of the 5' end (Thoreen 2012 TOP-like rule); depends on "
                              "the exact TSS, Ensembl 5' ends are approximate"})
            break
    for mm in G4_RE.finditer(utr):
        s, seq = mm.start(), mm.group(1)
        if any(x["type"] == "G4" and x["start"] <= s < x["end"] for x in m):
            continue  # nested match inside a reported G4
        m.append({"type": "G4", "start": s, "end": s + len(seq), "seq": seq, "label": "putative G-quadruplex",
                  "strength": None, "cite": "huppert2005", "note": "pattern (G3+N1-7){3,}G3+; not structure-tested"})
    if 3 <= L <= 30:
        ctx = full[L - 3:L + 8]
        if TISU_RE.match(ctx):
            m.append({"type": "TISU", "start": L - 3, "end": L, "seq": ctx, "label": "TISU at main start",
                      "strength": None, "cite": "elfakess2008",
                      "note": "AAS-AUG-GCGGC around the main AUG, UTR <= 30 nt (end clipped to the UTR; seq spans into CDS)"})
    ctx = ("N" * max(0, 6 - L) + utr[-6:]) + cds[:4]
    kc = kozak_class(full, L)
    m.append({"type": "kozak_main", "start": max(0, L - 6), "end": L, "seq": ctx, "label": f"main start Kozak: {kc}",
              "strength": kc, "cite": "kozak1986", "note": "-6..+4 context; seq includes the start codon and +4"})
    m.sort(key=lambda x: (x["start"], x["type"]))
    return m


def _gc(s):
    return sum(s.count(x) for x in "GC") / len(s) if s else 0.0


def features(utr, cds_head, fold=None):
    """Flat numeric features. fold (optional) = utr_data.json fold dict: mfe_local|mfe, cap_dG, acc_start|punp."""
    utr, cds = _rna(utr), _rna(cds_head)
    full, L = utr + cds, len(utr)
    ua = _uaugs(utr, cds)
    mot = find_motifs(utr, cds)
    kc = kozak_class(full, L)
    f = {
        "log_len": math.log1p(L),
        "gc": _gc(utr),
        "gc_first50": _gc(utr[:50]),
        "n_uAUG": float(len(ua)),
        "n_uAUG_strong": float(sum(u["kozak"] == "strong" for u in ua)),
        "n_uAUG_adequate": float(sum(u["kozak"] == "adequate" for u in ua)),
        "n_uAUG_weak": float(sum(u["kozak"] == "weak" for u in ua)),
        "n_uORF": float(sum(u["type"] == "uORF" for u in ua)),
        "n_uORF_overlap": float(sum(u["type"] == "overlapping" for u in ua)),
        "n_uORF_inframe_ext": float(sum(u["type"] == "in-frame-extension" for u in ua)),
        "max_uORF_aa": float(max([u["aa"] for u in ua if u["type"] == "uORF"] or [0])),
        "dist_last_uAUG_to_start": math.log1p(L - ua[-1]["start"] if ua else L),
        "n_nearcog_strong": float(sum(x["type"] == "nearcog" for x in mot)),
        "has_TOP": float(any(x["type"] == "TOP" for x in mot)),
        "n_G4": float(sum(x["type"] == "G4" for x in mot)),
        "has_TISU": float(any(x["type"] == "TISU" for x in mot)),
        "main_kozak_strong": float(kc == "strong"),
        "main_kozak_adequate": float(kc == "adequate"),
        "start_minus3_purine": float(L >= 3 and full[L - 3] in "AG"),
        "start_plus4_G": float(len(full) > L + 3 and full[L + 3] == "G"),
    }
    if fold:
        mfe = fold.get("mfe_local", fold.get("mfe"))  # utr_data.json fold: mfe_local = span<=100 structure
        if mfe is not None and L:
            f["mfe_per_nt"] = mfe / L
        if fold.get("cap_dG") is not None:
            f["cap_dG"] = float(fold["cap_dG"])
        if fold.get("acc_start") is not None:  # precomputed (covers +-15 nt incl. CDS side)
            f["acc_start"] = float(fold["acc_start"])
        elif fold.get("punp"):
            p = fold["punp"]
            w = p[max(0, L - 15):L]  # UTR side; CDS side only if punp covers it
            if len(p) > L:
                w = p[max(0, L - 15):L + 15]
            if w:
                f["acc_start"] = sum(w) / len(w)
    return f
