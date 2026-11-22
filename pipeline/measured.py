"""MEASURED layer: transcript-level translation efficiency (TE) from Cao et al. 2021 (Nat Commun 12:4106),
tables from https://github.com/zzz2010/5UTR_Optimizer (MIT), data/df_counts_and_len.TE_sorted.*.with_annot.txt,
downloaded into pipeline/external/.

TE = rpkm_riboseq / rpkm_rnaseq (column `te`). Matching: Ensembl transcript ID without version.
Their IDs come from an older Ensembl release than ours, so a transcript may be missing or its model may
have changed between releases even if the stable ID matches.

Ambiguity flag: the source quantification does not always resolve reads between sibling isoforms; some
sibling transcripts of the same gene carry bit-identical TE values. If any other transcript of the same
gene (external_gene_id) in the same source table has a TE equal to ours (to 10 significant digits),
the value is flagged ambiguous for that cell type.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(HERE, "external")
TABLES = {  # cell label -> file
    "HEK293T": "df_counts_and_len.TE_sorted.HEK_Andrev2015.with_annot.txt",
    "PC3": "df_counts_and_len.TE_sorted.pc3.with_annot.txt",
    "muscle": "df_counts_and_len.TE_sorted.Muscle.with_annot.txt",
}
METRIC = "TE (Ribo-seq/RNA-seq RPKM)"
NOTE = ("Cao et al. 2021 transcript-level TE (rpkm_riboseq/rpkm_rnaseq) from the 5UTR_Optimizer repo; "
        "matched by Ensembl transcript ID ignoring version (source uses an older Ensembl release than "
        "this file). ambiguous=true: a sibling isoform of the same gene has an identical TE in that cell "
        "type, i.e. reads were not resolved between isoforms; treat as gene-level, not isoform-level.")

_DATA = None


def _key(v):
    return f"{v:.10g}"


def load():
    global _DATA
    if _DATA is not None:
        return _DATA
    data = {}
    for cell, fn in TABLES.items():
        path = os.path.join(EXT, fn)
        with open(path) as fh:
            header = fh.readline().split()
            rows = {}
            by_gene = {}
            for line in fh:
                p = line.split()
                # first field is the row name (versioned ID), then the header columns
                rec = dict(zip(header, p[1:]))
                tid = rec["ensembl_tx_id"].split(".")[0]
                te = rec["te"]
                try:
                    te = float(te)
                    if te != te:  # NaN
                        te = None
                except ValueError:
                    te = None
                g = rec.get("external_gene_id")
                rows[tid] = (te, g, p[0])
                by_gene.setdefault(g, []).append((tid, te))
        data[cell] = (rows, by_gene)
    _DATA = data
    return data


def lookup(tx_id):
    """Return the `measured` dict for a (versioned) Ensembl transcript id, or None if absent in all tables."""
    data = load()
    tid = tx_id.split(".")[0]
    values, amb_cells, src_ids = {}, [], {}
    found = False
    for cell, (rows, by_gene) in data.items():
        r = rows.get(tid)
        if r is None:
            values[cell] = None
            continue
        found = True
        te, g, src = r
        src_ids[cell] = src
        values[cell] = None if te is None else round(te, 6)
        if te is not None:
            k = _key(te)
            if any(o != tid and v is not None and _key(v) == k for o, v in by_gene.get(g, [])):
                amb_cells.append(cell)
    if not found or all(v is None for v in values.values()):
        return None
    return {"metric": METRIC, "values": values, "cite": "cao2021", "ambiguous": bool(amb_cells),
            "ambiguous_cells": amb_cells, "source_ids": src_ids, "note": NOTE}
