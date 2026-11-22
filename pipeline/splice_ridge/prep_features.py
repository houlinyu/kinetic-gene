"""Interpretable per-transcript features for predicting a transcript's share of its gene.

Sources (all read-only, under /data/mdl_data/hyu/projects/public_isoform/isopedia/table_wide/):
  isoforms.parquet          isoform_uid, exons (chrom:(strand)[[s,e],...]), locus_id, transcript_id
  catalog_features.parquet  seq_uid(=isoform_uid) -> sequence/UTR/ORF/TSS features
  splice/shards/*.parquet   seq_uid -> Pangolin splice-site usage probabilities per isoform
                            (meta.json: {"model":"pangolin"}) -- NOT MaxEntScan; see README.
Exon structure, alt-event flags and TSS distances are computed here from `exons`.
Output: work/features.tsv
"""
import os, sys, glob, json, ast, collections
import pyarrow.parquet as pq

BUILD = os.environ.get("ISOPEDIA_BUILD", "table_wide2")   # table_wide2 = the 18,839-gene panel
B = f"/data/mdl_data/hyu/projects/public_isoform/isopedia/{BUILD}/"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "work", "features.tsv")

def parse_exons(s):
    chrom, rest = s.split(":", 1)
    strand = rest[1]
    ex = ast.literal_eval(rest[rest.index("["):])
    ex = sorted([(int(a), int(b)) for a, b in ex])
    return chrom, strand, ex

iso = pq.read_table(B + "isoforms.parquet").to_pydict()
struct = {}
for u, e, lid, tid in zip(iso["isoform_uid"], iso["exons"], iso["locus_id"], iso["transcript_id"]):
    try: chrom, strand, ex = parse_exons(e)
    except Exception: continue
    struct[u] = dict(locus=lid, tid=tid, chrom=chrom, strand=strand, ex=ex)
print("structures", len(struct), file=sys.stderr)

# catalog features
keep = ["gc", "log_len", "utr5_len", "utr3_len_nt", "orf_len_nt", "n_uatg", "nmd_ptc",
        "tss_gc", "tss_cpg_island", "tss_tata", "tes_pas_canonical", "exons_vs_gene_max",
        "n_siblings", "frac_share_tss", "frac_share_tes", "intron_retention", "likely_5p_fragment"]
cf = pq.read_table(B + "catalog_features.parquet", columns=["seq_uid"] + keep).to_pydict()
cat = {u: {k: cf[k][i] for k in keep} for i, u in enumerate(cf["seq_uid"])}
print("catalog features", len(cat), file=sys.stderr)

# Pangolin splice-site scores
pang = {}
_single = B + "splice/isoform_splice_scores.parquet"
_files = [_single] if os.path.exists(_single) else glob.glob(B + "splice/shards/*.parquet")
for f in _files:
    d = pq.read_table(f).to_pydict()
    for i, u in enumerate(d["seq_uid"]):
        pang[u] = dict(ss_p_mean=d["site_p_mean"][i], ss_p_min=d["site_p_min"][i],
                       ss_use_mean=d["site_use_mean"][i], ss_use_min=d["site_use_min"][i],
                       ss_logit_sum=d["logit_p_sum"][i], ss_njunc=d["n_junctions"][i])
print("pangolin isoforms", len(pang), file=sys.stderr)

# group by locus, pick structural reference = longest exonic length (no usage used -> no leakage)
by_locus = collections.defaultdict(list)
for u, s in struct.items(): by_locus[s["locus"]].append(u)

def five_end(s):  return s["ex"][0][0] if s["strand"] == "+" else s["ex"][-1][1]
def exlen(s):     return sum(b - a + 1 for a, b in s["ex"])
def donors(s):
    # genomic coords of the exon-end side of each junction, strand-aware
    j = [(s["ex"][i][1], s["ex"][i + 1][0]) for i in range(len(s["ex"]) - 1)]
    return ([a for a, b in j], [b for a, b in j]) if s["strand"] == "+" else ([b for a, b in j], [a for a, b in j])

rows = []
for lid, members in by_locus.items():
    if len(members) < 2: continue
    ref = max(members, key=lambda u: (exlen(struct[u]), struct[u]["tid"]))
    R = struct[ref]
    rd, ra = donors(R)
    r_introns = [(R["ex"][i][1] + 1, R["ex"][i + 1][0] - 1) for i in range(len(R["ex"]) - 1)]
    r_internal = R["ex"][1:-1]
    r5 = five_end(R)
    # dominant TSS among members = the 5' end shared by the most transcripts of the locus
    tss_count = collections.Counter(five_end(struct[u]) for u in members)
    dom_tss = tss_count.most_common(1)[0][0]
    for u in members:
        S = struct[u]
        if u not in cat: continue
        ex = S["ex"]; n_ex = len(ex)
        intr = [ex[i + 1][0] - ex[i][1] - 1 for i in range(n_ex - 1)]
        d, a = donors(S)
        my_introns = [(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(n_ex - 1)]
        # cassette: a reference internal exon swallowed by one of my introns
        cassette = int(any(s <= xa and xb <= e for (s, e) in my_introns for (xa, xb) in r_internal))
        alt5 = int(any(x not in set(rd) for x in d))
        alt3 = int(any(x not in set(ra) for x in a))
        altfirst = int(abs(five_end(S) - r5) > 10)
        r = dict(isoform_uid=u, locus_id=lid, transcript_id=S["tid"], is_ref=int(u == ref),
                 n_exons=n_ex,
                 log_tx_len=round(__import__("math").log10(max(exlen(S), 1)), 6),
                 log_mean_intron=round(__import__("math").log10(max(sum(intr) / len(intr), 1)), 6) if intr else 0.0,
                 log_min_intron=round(__import__("math").log10(max(min(intr), 1)), 6) if intr else 0.0,
                 log_max_intron=round(__import__("math").log10(max(max(intr), 1)), 6) if intr else 0.0,
                 log_mean_exon=round(__import__("math").log10(max(exlen(S) / n_ex, 1)), 6),
                 log_dist_dom_tss=round(__import__("math").log10(1 + abs(five_end(S) - dom_tss)), 6),
                 is_dom_tss=int(five_end(S) == dom_tss),
                 cassette_skip=cassette, alt_5ss=alt5, alt_3ss=alt3, alt_first_exon=altfirst,
                 is_longest=int(u == ref), n_iso_locus=len(members))
        r.update({k: cat[u][k] for k in keep})
        p = pang.get(u)
        for k in ("ss_p_mean", "ss_p_min", "ss_use_mean", "ss_use_min", "ss_logit_sum", "ss_njunc"):
            r[k] = (p[k] if p else "")
        r["has_pangolin"] = int(p is not None)
        rows.append(r)

hdr = list(rows[0].keys())
with open(OUT, "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows: fh.write("\t".join(str(r[k]) for k in hdr) + "\n")
print("wrote", OUT, len(rows), "rows,", len({r['locus_id'] for r in rows}), "multi-isoform loci",
      "pangolin-covered:", sum(r['has_pangolin'] for r in rows), file=sys.stderr)
