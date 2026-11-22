"""python pipeline/test_motifs.py  (plain asserts, no pytest needed)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motifs import find_motifs, features, FEATURES, SOURCES

CDS = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA"


def by(m, t):
    return [x for x in m if x["type"] == t]


def test_uorf_upstream():
    utr = "CCCAAAAUGGCAUAACCCCCC"  # uAUG at 6 (A at -3 = A, +4 G -> strong), stop UAA at 12
    m = find_motifs(utr, CDS)
    u = by(m, "uAUG")
    assert [x["start"] for x in u] == [6] and u[0]["strength"] == "strong", u
    s = by(m, "uORF_stop")
    assert [(x["start"], x["end"]) for x in s] == [(12, 15)], s
    f = features(utr, CDS)
    assert f["n_uORF"] == 1 and f["n_uAUG_strong"] == 1 and f["max_uORF_aa"] == 2


def test_overlapping_and_inframe():
    utr = "CCAUGCC"  # uAUG at 2, L-2=5 not multiple of 3 -> overlapping (stop in CDS, out of frame)
    m = find_motifs(utr, CDS)
    assert "overlapping" in by(m, "uAUG")[0]["note"] and not by(m, "uORF_stop")
    f = features(utr, CDS)
    assert f["n_uORF_overlap"] == 1 and f["n_uORF"] == 0
    utr2 = "CCAUGCCC"  # 8-2 = 6 -> in frame, no stop before main AUG
    f2 = features(utr2, CDS)
    assert f2["n_uORF_inframe_ext"] == 1, f2


def test_top_g4_kozak():
    utr = "CUUCCUCAGGGAGGGAUGGGAGGGACAAAAAAAAAAGCCACC"
    m = find_motifs(utr, CDS)
    t = by(m, "TOP")
    assert t and t[0]["start"] == 0 and t[0]["end"] == 7, t
    g = by(m, "G4")
    assert len(g) == 1 and g[0]["start"] == 8 and g[0]["seq"] == "GGGAGGGAUGGGAGGG", g
    k = by(m, "kozak_main")
    assert k[0]["strength"] == "strong" and k[0]["seq"] == "GCCACCAUGG"
    t2 = by(find_motifs("AGCUUCCUCA", CDS), "TOP")  # starts at 2 -> TOP-like
    assert t2 and t2[0]["start"] == 2, t2
    assert not by(find_motifs("AAAAACUUCCUC", CDS), "TOP")  # starts at 5 -> no
    assert not by(find_motifs("CUUCAG", CDS), "TOP")  # only 4 pyrimidines


def test_nearcog_tisu():
    utr = "AAAACUGGAAAAAA"  # CUG at 4: -3 = A, +4 = G -> strong
    m = find_motifs(utr, CDS)
    assert [x["start"] for x in by(m, "nearcog")] == [4]
    tisu_utr = "CCCCCCCCCCGAAC"
    m = find_motifs(tisu_utr, "AUGGCGGCCCCCCCCCCCCCCCCCCCCCCCCCC")
    assert by(m, "TISU"), m


def test_features_fold_and_labels():
    utr = "CCCAAAAUGGCAUAACCCCCC"
    f = features(utr, CDS, fold={"mfe": -2.1, "cap_dG": -1.0, "punp": [0.5] * len(utr)})
    assert abs(f["mfe_per_nt"] - (-0.1)) < 1e-9 and f["acc_start"] == 0.5 and f["cap_dG"] == -1.0
    assert set(f) <= set(FEATURES), set(f) - set(FEATURES)
    for x in find_motifs("CUUCCUCAGGGAGGGAUGGGAGGGAAUGAAAUAA", CDS):
        assert x["cite"] in SOURCES, x
    assert features("ACGT" * 5, CDS.replace("U", "T"))["gc"] == 0.5  # DNA input ok


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("ok", k)
    print("all motif tests passed")
