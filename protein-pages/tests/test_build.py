"""Checks of build.py's parsers on small hand-written files.

The test set has no hhsearch or InterProScan output, so their parsers are checked here on
files written by hand in each tool's documented layout. Run with: python -m pytest tests/
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build  # noqa: E402

HHR = """Query         FUN000001_FUN000001
Match_columns 120
No_of_seqs    1 out of 1
Neff          1
Searched_HMMs 1000
Date          Tue Oct 27 10:00:00 2026
Command       hhsearch -i q.a3m -d db

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 sp|P10415|BCL2_HUMAN Apoptosis  99.5 1.2E-25 3E-30  180.2   0.0   40   11-50     91-130 (239)
  2 sp|Q07817|B2CL1_HUMAN Bcl-2-lik  45.0   0.002 5E-08   30.1   0.0   20   60-79      5-24 (233)

No 1
>sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2
Probab=99.50  E-value=1.2e-25  Score=180.20  Aligned_cols=40  Identities=50%  Similarity=0.9

Q ss_pred             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Q FUN000001_FUN   11 MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPP   50 (120)
Q Consensus      11 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   50 (120)
                    |.|.|..|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|
T Consensus      91 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  130 (239)
T sp|P10415|BCL   91 MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPP  130 (239)
T ss_dssp           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

No 2
>sp|Q07817|B2CL1_HUMAN Bcl-2-like protein 1
Probab=45.00  E-value=0.002

Q FUN000001_FUN   60 MSQSNRELVV-DFLSYKLSQK   79 (120)
T sp|Q07817|B2C    5 MSQSNRELVVDFLSYKLSQKG   24 (233)
"""


def test_hhr(tmp_path):
    p = tmp_path / "FUN000001_FUN000001.hhr"
    p.write_text(HHR)
    q, hits = build.parse_hhr(p)
    assert q == "FUN000001_FUN000001"
    assert [h["target"] for h in hits] == ["P10415", "Q07817"]
    a, b = hits
    assert (a["q0"], a["q1"], a["t0"], a["t1"]) == (11, 50, 91, 130)
    assert a["passes"] and not b["passes"]  # b: Prob 45 is below 50
    assert a["aln"]["q"] == "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPP"
    assert a["aln"]["qs"] == 11 and a["aln"]["ts"] == 91
    assert b["aln"]["q"] == "MSQSNRELVV-DFLSYKLSQK"


def row(target, query, full_e, i_e, hmm=(1, 50), ali=(1, 50)):
    return [target, "-", "200", query, "-", "100", str(full_e), "50", "0", "1", "1", "1e-9",
            str(i_e), "50", "0", str(hmm[0]), str(hmm[1]), str(ali[0]), str(ali[1]),
            str(ali[0]), str(ali[1]), "0.9", "desc"]


def test_last_round_keeps_final_round():
    rounds = [
        [row("A", "q", 1e-30, 1e-30), row("B", "q", 1e-5, 1e-5)],
        [row("A", "q", 1e-40, 1e-40), row("B", "q", 1e-20, 1e-20), row("C", "q", 1e-4, 1e-4)],
        [row("A", "q", 1e-45, 1e-45), row("C", "q", 1e-25, 1e-25), row("B", "q", 1e-21, 1e-21)],
    ]
    rows = [r for rd in rounds for r in rd]
    assert [r[0] for r in build.last_round(rows)] == ["A", "C", "B"]


def test_domtbl_seq_uses_hmm_columns_for_query(tmp_path):
    p = tmp_path / "phmmer.domtblout"
    p.write_text("# header\n" + " ".join(row("sp|P1|X_HUMAN", "q", 1e-9, 1e-9, hmm=(5, 40), ali=(100, 135))) + "\n")
    h = build.materialize(build.parse_domtbl_seq(p, "phmmer")["q"][0])
    assert (h["q0"], h["q1"], h["t0"], h["t1"], h["target"]) == (5, 40, 100, 135, "P1")


def test_interproscan(tmp_path):
    p = tmp_path / "ips.tsv"
    p.write_text(
        "q\tmd5\t300\tPfam\tPF00452\tApoptosis regulator proteins, Bcl-2 family\t90\t190\t1.1E-20\tT\t27-10-2026\tIPR046371\tBcl-2 BH1-BH3\t-\t-\n"
        "q\tmd5\t300\tMobiDBLite\tmobidb-lite\tconsensus disorder prediction\t1\t30\t-\tT\t27-10-2026\t-\t-\t-\t-\n"
        "P10415\tmd5\t239\tPfam\tPF00452\tBcl-2\t93\t191\t1E-30\tT\t27-10-2026\tIPR046371\tBcl-2 BH1-BH3\t-\t-\n")
    out, target_ann, names = build.parse_interproscan(p, {"q": {}})
    a, b = out["q"]
    assert a["passes"] and a["names"] == {"PF00452", "IPR046371"}
    assert not b["passes"]  # a disorder prediction is not a family
    assert target_ann["P10415"][0][:2] == (93, 191)


def test_agreement_needs_same_family_and_overlap():
    h1 = build.hit("phmmer", 10, 100, "P1", 1, 90, 1e-9, [], True)
    h2 = build.hit("hmmscan", 20, 110, "Pkinase", None, None, 1e-9, [], True)
    h2["names"] = {"PF00069"}
    ann = {"P1": [(1, 80, {"PF00069"})]}
    for h in (h1, h2):
        h["_names"] = build.target_names(h, ann)
    progs, status, name = build.agreement([h1, h2], ann)
    assert status == "agree" and name == "PF00069"
    h2["q0"], h2["q1"] = 200, 300  # no longer overlapping
    assert build.agreement([h1, h2], ann)[1] == "apart"


def _two(names_a, names_b, a=(10, 100), b=(20, 110)):
    """A phmmer hit on P1 and an hmmscan hit, with the given families on each."""
    h1 = build.hit("phmmer", a[0], a[1], "P1", 1, 90, 1e-9, [], True)
    h2 = build.hit("hmmscan", b[0], b[1], "Pf", None, None, 1e-9, [], True)
    h2["names"] = names_b
    ann = {"P1": [(1, 80, names_a)]} if names_a is not None else {}
    for h in (h1, h2):
        h["_names"] = build.target_names(h, ann)
    return [h1, h2], ann


def test_status_agree():
    hits, ann = _two({"PF00069"}, {"PF00069"})
    assert build.agreement(hits, ann)[1:] == ("agree", "PF00069")


def test_status_disagree_lists_families_with_programs():
    hits, ann = _two({"PF00069"}, {"PF07714"})
    assert build.agreement(hits, ann)[1] == "disagree"
    assert build.DISAGREE_NAMES[0] == [("PF00069", ["phmmer"]), ("PF07714", ["hmmscan"])]
    s = build.header_sentence(2, "disagree", None, {"PF00069": "Pkinase", "PF07714": "PK_Tyr_Ser-Thr"})
    assert s.endswith("name different Pfam families: Pkinase PF00069 (phmmer), PK_Tyr_Ser-Thr PF07714 (hmmscan with Pfam).")


def test_status_overlap_unknown():
    hits, ann = _two(None, {"PF00069"})  # P1 has no Pfam output
    assert build.agreement(hits, ann)[1] == "overlap_unknown"


def test_status_apart():
    hits, ann = _two({"PF00069"}, {"PF07714"}, b=(200, 300))
    assert build.agreement(hits, ann)[1] == "apart"


def test_status_none_and_one():
    hits, ann = _two({"PF00069"}, {"PF00069"})
    assert build.agreement([], ann)[1] == "none"
    assert build.agreement(hits[:1], ann)[1] == "one"


def test_weak_hits_do_not_count():
    hits, ann = _two({"PF00069"}, {"PF00069"})
    hits[1]["evalue"] = 0.5  # drawn (at most 10) but above the counting threshold 0.001
    assert build.agreement(hits, ann)[1] == "one"



def test_coverage_counts_programs_once():
    a = build.hit("mmseqs2", 1, 10, "P1", 1, 10, 1e-9, [], True)
    b = build.hit("mmseqs2_3rounds", 5, 15, "P2", 1, 10, 1e-9, [], True)
    c = build.hit("phmmer", 8, 12, "P3", 1, 10, 1e-9, [], True)
    cov = build.coverage([a, b, c], 20)
    assert cov[0] == 1 and cov[7] == 2 and cov[14] == 1 and cov[15] == 0


def test_alphabets_from_kmerseek_source():
    src = Path("/Users/olga/code/kmerseek/src/rust/alphabets.rs")
    if not src.exists():
        return
    al = build.parse_alphabets(src)
    assert al["hp_lehninger2"]["G"] == "h" and al["hp_lehninger2"]["S"] == "p"
    assert len(set(al["gbmr4"].values())) == 4


def test_agreement_with_tied_hits_and_unknown_targets():
    # Two programs with a hit at the same range and the same E-value, targets with no known
    # family: the pair search must not fail on the tie (it did on the full Botryllus run).
    a = build.hit("phmmer", 10, 100, "P1", 1, 90, 1e-9, [], True)
    b = build.hit("jackhmmer", 10, 100, "P2", 1, 90, 1e-9, [], True)
    c = build.hit("jackhmmer", 200, 300, "P3", 1, 90, 1e-5, [], True)
    for h in (a, b, c):
        h["_names"] = build.target_names(h, {})
    assert build.agreement([a, b, c], {})[1] == "overlap_unknown"
