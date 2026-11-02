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
    # template length and description, for templates not in swissprot.fasta (PDB70 chains)
    assert (a["tlen"], b["tlen"]) == (239, 233)
    assert a["desc"] == "Apoptosis regulator Bcl-2" and b["desc"] == "Bcl-2-like protein 1"


def test_hhr_pdb70_template(tmp_path):
    """A PDB70 template keeps its chain name, description and length from the .hhr."""
    p = tmp_path / "q.hhr"
    p.write_text(HHR.replace("sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2", "6XYZ_A Apoptosis regulator Bcl-2; BH3 domain")
                 .replace("sp|P10415|BCL2_HUMAN Apoptosis ", "6XYZ_A Apoptosis regulator "))
    q, hits = build.parse_hhr(p)
    assert hits[0]["target"] == "6XYZ_A" and hits[0]["tlen"] == 239
    assert hits[0]["desc"] == "Apoptosis regulator Bcl-2; BH3 domain"


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
    # a disorder prediction is drawn, labelled "disorder", and never counts as a hit
    assert b["feature"] and b["target"] == "disorder" and not build.counts(b)
    b["_names"] = None
    assert not build.agreement([b], target_ann)[0]  # alone, the protein has no program with a hit
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


def test_weak_hits_do_not_count(monkeypatch):
    monkeypatch.setattr(build, "COUNT_E", 1e-3)
    hits, ann = _two({"PF00069"}, {"PF00069"})
    hits[1]["evalue"] = 0.5  # drawn (at most 10) but above a counting threshold of 0.001
    assert build.agreement(hits, ann)[1] == "one"


def test_sentence_follows_the_slider():
    hits, ann = _two({"PF00069"}, {"PF00069"})
    hits[1]["evalue"] = 0.5
    ctx = dict(target_ann=ann, names={"PF00069": "Pkinase"}, n_run=5)
    steps = build.sentence_steps(hits, ctx, -6)
    at = lambda v: [x for x in steps if x[0] <= v][-1][1]
    assert at(-3).startswith("1 of the 5 programs run reports a hit")
    assert at(0).startswith("2 of the 5 programs run report a hit. Their hits along the protein agree on Pkinase")
    assert build.COUNT_E == 10.0  # put back afterwards



def test_coverage_counts_programs_once():
    a = build.hit("mmseqs2", 1, 10, "P1", 1, 10, 1e-9, [], True)
    b = build.hit("mmseqs2_3rounds", 5, 15, "P2", 1, 10, 1e-9, [], True)
    c = build.hit("phmmer", 8, 12, "P3", 1, 10, 1e-9, [], True)
    cov = build.coverage([a, b, c], 20)
    assert cov[0] == 1 and cov[7] == 2 and cov[14] == 1 and cov[15] == 0


DATA = Path(__file__).resolve().parents[1] / "test_inputs" / "data"


def test_alphabets_from_kmerseek_source():
    al = build.parse_alphabets(DATA / "kmerseek_alphabets.rs")
    assert al["hp_lehninger2"]["G"] == "h" and al["hp_lehninger2"]["S"] == "p"
    assert len(set(al["gbmr4"].values())) == 4


def test_site_from_test_data(tmp_path):
    """Build the 51-protein site from test_inputs/data and run check_site.py on it."""
    import json
    import subprocess
    here = Path(build.__file__).resolve().parent
    out = tmp_path / "site"
    r = subprocess.run([sys.executable, str(here / "build.py"), str(DATA), "--out", str(out), "--jobs", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    summary = json.loads((out / "summary.json").read_text())
    assert summary["total"] == 51 and summary["programs_run"] == 8  # Foldseek only as ProstT5 + Foldseek
    r = subprocess.run([sys.executable, str(here / "tests" / "check_site.py"), str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    # The BHF page has the hit search box and the folding note, and its page data carry what the
    # search looks through: 5 ProstT5 + Foldseek hits are on RRP36 proteins (entry name RRP36_...).
    import re
    html = (out / "proteins" / "FUN008084_FUN008084.html").read_text()
    assert '<input id="hq" type="search"' in html and 'id="hqn"' in html
    assert "A program with more than 100 hits starts folded into one row, shaded by how many hits cover each residue. Click its name or the folded row to see each hit; click a name to fold any row." in html
    assert '<div id="ht"></div>' in html
    page = json.loads(re.search(r"const PAGE=(.*?);</script>", html).group(1).replace("<\\/", "</"))
    rrp36 = [h for h in page["hits"] if "rrp36" in " ".join(str(h.get(f) or "") for f in "ltngo").lower()]
    assert len(rrp36) == 5 and {h["k"] for h in rrp36} == {"prostt5"}
    # Biomni's programs with no position are rows of the drawing and count in the one sentence:
    # BHF's only one is FANTASIA, at distance 0.22, above the 0.15 that counts, so 5 of 12.
    for k in ("diamond", "jackhmmer", "fantasia", "unifire"):
        assert 'data-k="w_%s"' % k in html
    assert "not close enough to count" in html and "P32845" in html
    assert '<p class="summary" id="sent">5 of 12 programs report a hit.' in html
    assert "What the whole-proteome annotation found" not in html
    # KIF11 has a Swiss-Prot hit, so FANTASIA was not run on it: 11 programs run
    kif = (out / "proteins" / "FUN016555_KIF11.html").read_text()
    assert '<p class="summary" id="sent">10 of the 11 programs run report a hit.' in kif
    assert "P52732" in kif and "not run: FANTASIA was run only on proteins" in kif
    # The BHF sentence and the note under it link each Pfam accession to its InterPro page, as
    # built and in the slider's copies the page swaps in (PAGE.sent).
    vpr = ('VPR (<a class="mono" href="https://www.ebi.ac.uk/interpro/entry/pfam/PF00522/" '
           'target="_blank" rel="noopener">PF00522</a>)')
    assert vpr in re.search(r'<p class="summary" id="sent">(.*?)</p>', html).group(1)
    note = re.search(r'<div id="sentnote">(.*?)</div>', html).group(1)
    assert 'GrpE (<a class="mono" href="https://www.ebi.ac.uk/interpro/entry/pfam/PF01025/"' in note
    assert any(vpr in s[1] for s in page["sent"])
    assert not [s for s in page["sent"] if re.search(r"\(PF\d{5}\)", s[1] + s[2])]


def test_acc_html_links_and_escapes():
    s = build.acc_html('A<b> & "x" PF00069.32, IPR000719; PF1234 xPF00069')
    assert s.startswith("A&lt;b&gt; &amp; &quot;x&quot; ")
    assert '<a class="mono" href="https://www.ebi.ac.uk/interpro/entry/pfam/PF00069/" target="_blank" rel="noopener">PF00069.32</a>' in s
    assert 'href="https://www.ebi.ac.uk/interpro/entry/InterPro/IPR000719/"' in s
    assert s.count("<a ") == 2  # PF1234 is too short; xPF00069 is part of a longer word


def test_accession_keeps_isoform_number():
    """An AlphaFold model of an isoform names the isoform: C3VPR6-2, not the canonical C3VPR6."""
    assert build.af_accession("AF-P12345-F1-model_v6") == "P12345"
    assert build.af_accession("AF-C3VPR6-2-F1-model_v6") == "C3VPR6-2"
    assert build.af_accession("AF-A0A024R1R8-F1-model_v6") == "A0A024R1R8"
    assert build.accession("AF-C3VPR6-2-F1-model_v6") == "C3VPR6-2"
    assert build.accession("sp|C3VPR6-2|NLRC5_MOUSE") == "C3VPR6-2"
    # the committed ProstT5 + Foldseek table: its last column is the accession of the target name
    import gzip
    rows = [line.rstrip("\n").split("\t") for line in gzip.open(DATA / "prostt5_foldseek.m8.gz", "rt")]
    assert all(f[-1] == build.af_accession(f[1]) for f in rows)
    assert sum("-" in f[-1] for f in rows) > 0


def test_agreement_with_tied_hits_and_unknown_targets():
    # Two programs with a hit at the same range and the same E-value, targets with no known
    # family: the pair search must not fail on the tie (it did on the full Botryllus run).
    a = build.hit("phmmer", 10, 100, "P1", 1, 90, 1e-9, [], True)
    b = build.hit("jackhmmer", 10, 100, "P2", 1, 90, 1e-9, [], True)
    c = build.hit("jackhmmer", 200, 300, "P3", 1, 90, 1e-5, [], True)
    for h in (a, b, c):
        h["_names"] = build.target_names(h, {})
    assert build.agreement([a, b, c], {})[1] == "overlap_unknown"


def test_sentence_counts_programs_run():
    assert build.header_sentence(2, "apart", None, {}, 5).startswith("2 of the 5 programs run report a hit.")
    assert build.header_sentence(1, "one", None, {}, 8) == "1 of 8 programs reports a hit."
    assert build.header_sentence(0, "none", None, {}, 5) == "None of the 5 programs run reports a hit."


def test_out_folder_that_is_not_a_site_is_kept(tmp_path):
    import subprocess
    out = tmp_path / "mine"
    out.mkdir()
    (out / "notes.txt").write_text("keep me")
    r = subprocess.run([sys.executable, str(Path(build.__file__)), str(tmp_path), "--out", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 2 and (out / "notes.txt").exists()


def test_seed_runs_finds_k_columns_in_one_class():
    """BHF (FUN008084) against PG03_MYCTO (P9WIG2), kmerseek gbmr4 k = 19: residues 63-119 of
    BHF against 882-938 of P9WIG2. The first 25 columns are in the same gbmr4 class in both, so
    seven 19-column windows start there; nowhere else is there a run of 19."""
    al = build.parse_alphabets(DATA / "kmerseek_alphabets.rs")
    prots = {h.split()[0]: s for h, s in build.read_fasta(DATA / "proteins.fasta").items()}
    sp = {build.accession(h.split()[0]): s for h, s in build.read_fasta(DATA / "swissprot.fasta").items()}
    q, t = prots["FUN008084_FUN008084"][62:119], sp["P9WIG2"][881:938]
    assert build.seed_runs(q, t, al["gbmr4"], 19) == ([(0, 24)], 7)
    assert build.seed_runs(q, t, al["gbmr4"], 26) == ([], 0)


def test_seed_runs_small_cases():
    tab = {"A": "a", "D": "a", "L": "l", "I": "l", "G": "g"}
    # one break in the middle: two runs of 3, each one 3-mer
    assert build.seed_runs("AALGAAL", "DDIADDI", tab, 3) == ([(0, 2), (4, 6)], 2)
    # a residue outside the alphabet matches nothing, not even itself
    assert build.seed_runs("AXA", "AXA", tab, 2) == ([], 0)
    assert build.seed_runs("ALA", "DID", tab, 3) == ([(0, 2)], 1)


def test_dayhoff6_is_sourmash_groups():
    d = build.parse_alphabets(DATA / "kmerseek_alphabets.rs")["dayhoff6"]
    assert len(set(d.values())) == 6 and len(d) == 20
    assert d["A"] == d["G"] == d["P"] == d["S"] == d["T"] and d["H"] == d["K"] == d["R"] and d["C"] != d["A"]


def test_read_target_pfam(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("# made by hand\naccession\tpfam\tpfam_name\tstart\tend\tsource\n"
                 "P9WIG2\tPF00934\tPE family\t4\t94\tInterPro API\n"
                 "Q00001\t-\t\t\t\tInterPro API: no Pfam match\n")
    doms, names = build.read_target_pfam(p)
    assert doms == {"P9WIG2": [(4, 94, {"PF00934"})], "Q00001": []}
    assert names == {"PF00934": "PE family"}
    # asked, and InterPro has no Pfam match: known to have none (an empty set), not unknown (None)
    assert build.target_names(dict(target="Q00001", t0=1, t1=50), doms) == set()
    assert build.target_names(dict(target="Q99999", t0=1, t1=50), doms) is None
    assert build.target_names(dict(target="P9WIG2", t0=880, t1=940), doms) == set()
    assert build.target_names(dict(target="P9WIG2", t0=50, t1=60), doms) == {"PF00934"}


def test_kmerseek_hits_carry_alphabet_ksize_and_seed(tmp_path):
    """Each kmerseek hit keeps its own alphabet and k-mer size, a page lists the classes of
    exactly the alphabets its kmerseek hits use, and the BHF / PG03_MYCTO hit has its seed."""
    import json
    import re
    import subprocess
    here = Path(build.__file__).resolve().parent
    out = tmp_path / "site"
    r = subprocess.run([sys.executable, str(here / "build.py"), str(DATA), "--out", str(out), "--jobs", "2"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    page = (out / "proteins" / "FUN008084_FUN008084.html").read_text()
    P = json.loads(re.search(r"const PAGE=(\{.*?\});</script>", page).group(1))
    km = [h for h in P["hits"] if h["k"] == "kmerseek"]
    assert set(P["alpha"]) == {h["m"] for h in km}
    pg = [h for h in km if h.get("t") == "P9WIG2" and h["m"] == "gbmr4" and h["kk"] == 19]
    assert pg and pg[0]["sd"] == [[63, 87]] and pg[0]["sn"] == 7
    # P9WIG2's Pfam domains come from target_pfam_interpro.tsv: PE family at residues 4-94
    assert P["td"]["P9WIG2"] == [[4, 94, "PE family (PF00934)"]] and P["ds"]["P9WIG2"] == "interpro"
