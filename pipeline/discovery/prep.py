"""Discovery step 1-2: same-protein, different-leader isoform groups from GENCODE v47, folded.

Inputs (read-only):
  pipeline/ridge/downloads/gencode.v47.pc_transcripts.fa.gz   (UTR5/CDS ranges in headers)
  pipeline/discovery/dl/gencode.v47.basic.annotation.gtf.gz     (tags: basic, MANE_Select, TSL, *_NF)
Output: pipeline/discovery/work/isoforms.json  (list of groups; each isoform in the utr_data.json
isoform shape that app/src/analytic.js needs: len, tags, main_kozak, uorfs, hairpins)

Filters (documented in data/discovery.json meta):
  protein_coding transcript, GENCODE 'basic' tag, no mRNA_start_NF / cds_start_NF / readthrough tag,
  not a PAR_Y copy, header has UTR5 starting at 1 and a CDS that starts with ATG,
  5' UTR length 10..MAXLEN nt.  Group = (gene, identical CDS nucleotide sequence). Identical UTRs are
  collapsed (best tag kept). Groups with >= 2 distinct UTRs are kept.
Folding (two methods, calls must agree across both):
  'global' = ViennaRNA MFE of the whole UTR (same as build_data.py v1, used for the app data the
             feasible parameter region was computed on);
  'local'  = ViennaRNA MFE with max_bp_span = 100 (local pairs only; cheap stand-in for SPEC_v2 plfold).
  Hairpins from the dot-bracket with the same rule as build_data.find_hairpins (copied, not imported,
  so another agent's edits can't change this run).
"""
import gzip
import json
import os
import re
import sys
from multiprocessing import Pool

import RNA

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FA = os.path.join(ROOT, "pipeline/ridge/downloads/gencode.v47.pc_transcripts.fa.gz")
GTF = os.path.join(ROOT, "pipeline/discovery/dl/gencode.v47.basic.annotation.gtf.gz")
OUT = os.path.join(ROOT, "pipeline/discovery/work/isoforms.json")
MINLEN, MAXLEN = 10, 1500
NPROC = 8


def parse_gtf():
    info = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 9)
            if c[2] != "transcript":
                continue
            a = c[8]
            tid = re.search(r'transcript_id "([^"]+)"', a).group(1)
            tags = re.findall(r'tag "([^"]+)"', a)
            tsl = re.search(r'transcript_support_level "([^"]+)"', a)
            info[tid] = {
                "type": re.search(r'transcript_type "([^"]+)"', a).group(1),
                "gene": re.search(r'gene_name "([^"]+)"', a).group(1),
                "tags": tags,
                "tsl": tsl.group(1) if tsl else "NA",
                "chrom": c[0], "strand": c[6],
            }
    return info


def parse_fa():
    recs = []
    with gzip.open(FA, "rt") as f:
        h, seq = None, []
        for line in f:
            if line.startswith(">"):
                if h:
                    recs.append((h, "".join(seq)))
                h, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
        if h:
            recs.append((h, "".join(seq)))
    return recs


def kozak_class(full, a):
    m3 = full[a - 3] if a - 3 >= 0 else "N"
    p4 = full[a + 3] if a + 3 < len(full) else "N"
    n = (m3 in "AG") + (p4 == "G")
    return ["weak", "adequate", "strong"][n]


def find_uorfs(utr, cds):  # identical rule to build_data.find_uorfs
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


def find_hairpins(db, seq, min_bp=3):  # identical rule to build_data.find_hairpins
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
            continue
        i, j = k, l
        while True:
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


def fold_both(utr):
    db_g, mfe_g = RNA.fold_compound(utr).mfe()
    md = RNA.md()
    md.max_bp_span = 100
    db_l, mfe_l = RNA.fold_compound(utr, md).mfe()
    return {"global": {"mfe": round(mfe_g, 2), "hairpins": find_hairpins(db_g, utr)},
            "local": {"mfe": round(mfe_l, 2), "hairpins": find_hairpins(db_l, utr)}}


def tag_rank(t):
    tg = t["tags"]
    return (("MANE_Select" in tg) * 4 + ("Ensembl_canonical" in tg) * 2 +
            (t["tsl"] in ("1", "2")) * 1, -int(t["tsl"]) if t["tsl"].isdigit() else -9)


def main():
    info = parse_gtf()
    print("gtf transcripts", len(info), file=sys.stderr)
    stats = {"fa_records": 0, "not_in_basic_gtf": 0, "not_pc": 0, "bad_tags": 0, "no_utr5": 0,
             "non_atg": 0, "len_filter": 0, "kept": 0}
    groups = {}
    for h, seq in parse_fa():
        stats["fa_records"] += 1
        f = h.split("|")
        tid = f[0]
        if "PAR_Y" in tid:
            stats["bad_tags"] += 1
            continue
        t = info.get(tid)
        if t is None:
            stats["not_in_basic_gtf"] += 1
            continue
        if t["type"] != "protein_coding":
            stats["not_pc"] += 1
            continue
        if "basic" not in t["tags"] or any(x in t["tags"] for x in ("mRNA_start_NF", "cds_start_NF", "readthrough_transcript")):
            stats["bad_tags"] += 1
            continue
        rng = {k: tuple(map(int, v.split("-"))) for k, v in (x.split(":") for x in f[7:] if ":" in x)}
        if "UTR5" not in rng or "CDS" not in rng or rng["UTR5"][0] != 1:
            stats["no_utr5"] += 1
            continue
        s = seq.upper().replace("T", "U")
        utr = s[:rng["CDS"][0] - 1]
        cds = s[rng["CDS"][0] - 1:rng["CDS"][1]]
        rest = s[rng["CDS"][0] - 1:]
        if not cds.startswith("AUG") or len(utr) != rng["UTR5"][1]:
            stats["non_atg"] += 1
            continue
        if not (MINLEN <= len(utr) <= MAXLEN) or set(utr) - set("ACGU"):
            stats["len_filter"] += 1
            continue
        stats["kept"] += 1
        key = (t["gene"], cds)
        groups.setdefault(key, []).append({"id": tid, "name": f[4], "gene": t["gene"], "tags": t["tags"],
                                           "tsl": t["tsl"], "utr5": utr, "rest": rest})
    out = []
    for (gene, cds), txs in groups.items():
        by_utr = {}
        for t in txs:
            if t["utr5"] not in by_utr or tag_rank(t) > tag_rank(by_utr[t["utr5"]]):
                by_utr[t["utr5"]] = t
        if len(by_utr) < 2:
            continue
        out.append({"gene": gene, "cds_len": len(cds), "isoforms": sorted(by_utr.values(), key=tag_rank, reverse=True)})
    stats["groups"] = len(out)
    stats["isoforms_in_groups"] = sum(len(g["isoforms"]) for g in out)
    stats["genes"] = len({g["gene"] for g in out})
    print(json.dumps(stats), file=sys.stderr)
    utrs = sorted({i["utr5"] for g in out for i in g["isoforms"]}, key=len, reverse=True)
    print("folding", len(utrs), "unique UTRs, total nt", sum(map(len, utrs)), file=sys.stderr)
    with Pool(NPROC) as pool:
        folds = dict(zip(utrs, pool.map(fold_both, utrs, chunksize=4)))
    for g in out:
        for i in g["isoforms"]:
            u, rest = i["utr5"], i.pop("rest")
            L = len(u)
            i["len"] = L
            i["main_kozak"] = {"class": kozak_class(u + rest, L)}
            i["cds_context"] = ("N" * max(0, 6 - L) + u[-6:]) + rest[:4]
            i["uorfs"] = find_uorfs(u, rest)
            i["fold"] = folds[u]
            i["gc"] = round(sum(u.count(x) for x in "GC") / L, 4)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump({"stats": stats, "filters": {"minlen": MINLEN, "maxlen": MAXLEN}, "groups": out}, open(OUT, "w"))
    print("wrote", OUT, file=sys.stderr)


if __name__ == "__main__":
    main()
