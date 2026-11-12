"""Checks of build_full.py and of what it assumes about the Biomni tables, on the committed
tables in full_inputs/. Run with: python -m pytest tests/
"""
import base64
import csv
import gzip
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import build  # noqa: E402
import build_full  # noqa: E402

sys.path.insert(0, str(HERE / "tests"))
import check_hover  # noqa: E402  (needs Playwright: pip install -r requirements-test.txt)

FULL = HERE / "full_inputs"


def tsv(name):
    return list(csv.DictReader(build_full.opentext(build_full.inp(FULL, name)), delimiter="\t"))


def test_name_status_is_a_comparison_of_the_two_names():
    """The index's wording for name_status rests on this rule holding on every row."""
    for r in tsv("botryllus_protein_annotations.tsv"):
        e, p = r["existing_name"], r["proposed_name"]
        want = ("agree" if e.lower() == p.lower() else "conflict") if e and p else "unconfirmed" if e else "new" if p else ""
        assert r["name_status"] == want, r["protein_id"]


def test_isoform_ids_follow_file_order():
    """proteome_lengths.tsv names the k-th copy of a FASTA header <gene>.k; DIAMOND's query
    length in the allorecognition table agrees on all 832 rows, so the tables use the same names."""
    L = {r["protein_id"]: int(r["length"]) for r in tsv("proteome_lengths.tsv")}
    rows = list(csv.reader(open(FULL / "botryllus_allorecognition_hits.tsv"), delimiter="\t"))
    assert len(rows) == 832
    assert all(int(f[6]) == L[f[0]] for f in rows)


def test_shard_hash_is_the_same_in_python_and_javascript():
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")
    ids = ["FUN000001", "FUN019272.2", "FUN999999novelgene12345b6aa78", "FUN008084"]
    js = open(HERE / "assets" / "viewer.js").read()
    fn = js[js.index("function shardOf"):js.index("function esc")]
    out = subprocess.run([node, "-e", "var N_SHARDS=%d;%s;console.log(JSON.stringify(%s.map(shardOf)))"
                          % (build_full.N_SHARDS, fn, json.dumps(ids))], capture_output=True, text=True, check=True)
    assert json.loads(out.stdout) == [build_full.shard_of(i) for i in ids]


@pytest.fixture(scope="module")
def site(tmp_path_factory):
    out = tmp_path_factory.mktemp("full")
    summ = build_full.build_site(FULL, out)
    return out, summ


def read_js(path, prefix):
    """A data file's JSON: gzip and base64 between its first and last quote (build_full.packed)."""
    s = path.read_text()
    return json.loads(gzip.decompress(base64.b64decode(s[s.index('="') + 2:s.rindex('"')])))


def test_every_protein_in_its_one_data_file(site):
    out, summ = site
    F = read_js(out / "full" / "rows.js", "window.FULL=")
    ids = [r[0] for r in F["rows"]]
    assert len(ids) == len(set(ids)) == 45339 == summ["total"]
    where = Counter()
    for s in range(build_full.N_SHARDS):
        shard = read_js(out / "full" / ("%03d.js" % s), "window.FULL_SHARD=")
        for pid in shard:
            assert build_full.shard_of(pid) == s
            where[pid] += 1
    assert set(where) == set(ids) and set(where.values()) == {1}


def test_counts_match_the_tables(site):
    _, s = site
    ann = tsv("botryllus_protein_annotations.tsv")
    assert s["swissprot_hit"] == sum(1 for r in ann if r["sp_acc"])
    assert s["pfam_domain"] == sum(1 for r in ann if r["pfam_domains"])
    # The category set by a Fester or BHF hit is dropped (a search within Botryllus); its 53 proteins
    # are in no other category, so they are no longer candidates
    kept = [[c for c in build_full.split_list(r["categories"]) if c != "allorecognition"] for r in tsv("botryllus_immune_genes.tsv")]
    assert s["immune_candidates"] == sum(1 for c in kept if c) == 2835 - 53
    assert "allorecognition" not in s["immune_categories"] and "allorecognition_hit" not in s
    assert s["clusters"] == len({f[0] for f in csv.reader(open(FULL / "botryllus_mmseqs_clusters.tsv"), delimiter="\t")})
    assert sum(s["name_check"].values()) == s["total"]
    want = Counter(c for cats in kept for c in cats)
    assert s["immune_categories"] == dict(want)


def test_drawn_proteins_matched_by_sequence(tmp_path):
    """Build the 51-protein site, then link it: every drawn page is reached, and the TTN page goes
    to FUN019272.2, the isoform with its 2,160 residues, not to FUN019272 (1,991)."""
    import build
    drawn = tmp_path / "site"
    subprocess.run([sys.executable, str(HERE / "build.py"), str(HERE / "test_inputs" / "data"), "--out", str(drawn), "--jobs", "1"],
                   check=True, capture_output=True)
    out = tmp_path / "full"
    out.mkdir()
    s = build_full.build_site(FULL, out, drawn, HERE / "test_inputs" / "data" / "proteins.fasta")
    assert s["drawn_pages"] == 51
    F = read_js(out / "full" / "rows.js", "window.FULL=")
    page = {r[0]: r[14] for r in F["rows"]}
    assert page["FUN019272.2"] == "proteins/FUN019272_TTN.1.html" and page["FUN019272"] == ""
    assert page["FUN008084"] == "proteins/%s.html" % build.safe("FUN008084_FUN008084")
    assert (out / "drawn.html").exists() and (out / "proteins" / "FUN019272_TTN.1.html").exists()
    # in a browser: pointing at a hit on a test protein's page moves nothing, and no prediction label is misplaced
    res = check_hover.check(out, ["proteins/FUN036660_MYO10.2.html", "proteins/FUN019272_TTN.1.html"], max_hits=60)
    assert len(res) == 2 * len(check_hover.WIDTHS) and all(r["tried"] >= 50 for _, _, r in res)
    assert not check_hover.failures(res), check_hover.failures(res)


def test_out_guard(tmp_path):
    d = tmp_path / "x"
    d.mkdir()
    (d / "keep.txt").write_text("not a site")
    r = subprocess.run([sys.executable, str(HERE / "build_full.py"), str(FULL), "--out", str(d)], capture_output=True)
    assert r.returncode == 2 and (d / "keep.txt").exists()


def test_protein_missing_from_clusters_is_a_clear_error(tmp_path):
    d = tmp_path / "in"
    shutil.copytree(FULL, d)
    f = d / "botryllus_mmseqs_clusters.tsv"
    lines = f.read_text().splitlines(keepends=True)
    f.write_text("".join(lines[1:]))
    with pytest.raises(SystemExit) as e:
        build_full.load(d)
    assert "no row in botryllus_mmseqs_clusters.tsv" in str(e.value)


def test_unifire_on_the_pages(site):
    """UniFIRE predictions reach the right protein: a later copy of a repeated header
    (FUN000051.2 is UniFIRE's FUN000051_GALNT6.2), and proteins over 5,000 aa read "not run"."""
    out, summ = site
    uf = {}
    for s in range(build_full.N_SHARDS):
        for pid, d in read_js(out / "full" / ("%03d.js" % s), "window.FULL_SHARD=").items():
            uf[pid] = d["uf"]
    assert sum(1 for v in uf.values() if v and v != "not run") == 19001
    assert sum(1 for v in uf.values() if v == "not run") == 50
    raw = [r for r in csv.DictReader(gzip.open(FULL / "unifire_predictions.tsv.gz", "rt"), delimiter="\t")
           if r["protein_id"] == "FUN000051.2"]
    assert raw and len(uf["FUN000051.2"]) == len(raw)


def test_foldseek_on_the_pages(site):
    """The best ProstT5 + Foldseek hit reaches its protein at any E-value in the table (42,035 proteins);
    25,804 are at FOLDSEEK_MAX_E or below, as in the Biomni report (v3), and only those are counted.
    BHF's best is at E-value 0.06, so the pages show it and the counts leave it out."""
    out, summ = site
    fs = {}
    for s in range(build_full.N_SHARDS):
        for pid, d in read_js(out / "full" / ("%03d.js" % s), "window.FULL_SHARD=").items():
            fs[pid] = d["fs"]
    assert sum(1 for v in fs.values() if v) == 42035
    assert summ["foldseek_hit"] == sum(1 for v in fs.values() if v and float(v[4]) <= build_full.FOLDSEEK_MAX_E) == 25804
    assert float(fs["FUN008084"][4]) > build_full.FOLDSEEK_MAX_E
    raw = {r["protein_id"]: r for r in csv.DictReader(gzip.open(FULL / "foldseek_best_per_protein.tsv.gz", "rt"), delimiter="\t")}
    pid = next(p for p, v in fs.items() if v)
    assert build.af_accession(raw[pid]["target"]) == fs[pid][0] and int(raw[pid]["qstart"]) == fs[pid][2]


def test_foldseek_on_structures_on_the_pages(site):
    """The best hit of Foldseek on Biomni's predicted structures: 13,946 proteins, 3,729 of them at
    FOLDSEEK_MAX_E or below, all 200 aa or less; every longer protein says "not run", BHF (252 aa) among them."""
    out, summ = site
    fst, length = {}, {}
    for s in range(build_full.N_SHARDS):
        for pid, d in read_js(out / "full" / ("%03d.js" % s), "window.FULL_SHARD=").items():
            fst[pid] = d["fst"]
    for r in csv.DictReader(gzip.open(FULL / "proteome_lengths.tsv.gz", "rt"), delimiter="\t"):
        length[r["protein_id"]] = int(r["length"])
    assert sum(1 for v in fst.values() if v and v != "not run") == 13946
    assert summ["foldseek_struct_hit"] == sum(1 for v in fst.values() if v and v != "not run" and float(v[4]) <= build_full.FOLDSEEK_MAX_E) == 3729
    assert {p for p, v in fst.items() if v == "not run"} == {p for p, n in length.items() if n > build_full.STRUCT_MAX_AA}
    assert fst["FUN008084"] == "not run"
    raw = {r["protein_id"]: r for r in csv.DictReader(gzip.open(FULL / "foldseek_struct_best_per_protein.tsv.gz", "rt"), delimiter="\t")}
    pid = next(p for p, v in fst.items() if v and v != "not run")
    assert build.af_accession(raw[pid]["target"]) == fst[pid][0] and int(raw[pid]["qstart"]) == fst[pid][2]
    assert fst[pid][8:10] == [float(raw[pid]["alntmscore"]), float(raw[pid]["lddt"])]
    assert fst[pid][10:] == [int(raw[pid]["tstart"]), int(raw[pid]["tend"])]  # the range on the other protein


def test_known_proteins_link_only_when_drawn():
    page = {"FUN008084": "proteins/FUN008084_FUN008084.html"}
    html = build_full.known_html(page)
    assert 'href="proteins/FUN008084_FUN008084.html"' in html and "BHF" in html
    assert build_full.known_html({}) == ""


def test_index_position_track(site):
    """The first page's track: each row carries InterProScan's Pfam matches and MMseqs2's best hit,
    read from the same tables the protein pages use, and the headline counts match the summary."""
    out, summ = site
    data = read_js(out / "full" / "rows.js", "window.FULL=")
    rows = {r[0]: r for r in data["rows"]}
    ips = {}
    for f in csv.reader(gzip.open(FULL / "interproscan.tsv.gz", "rt"), delimiter="\t"):
        if f[3] == "Pfam":
            ips.setdefault(f[0], []).append((int(f[6]), int(f[7]), f[4]))
    assert sum(1 for r in rows.values() if r[15]) == len(ips)
    for pid in list(ips)[:200]:
        got = sorted((s, e, data["pfam"][i][0]) for s, e, i in rows[pid][15])
        assert got == sorted(ips[pid])
    mm = {r["protein_id"]: r for r in csv.DictReader(gzip.open(FULL / "mmseqs2_best_per_protein.tsv.gz", "rt"), delimiter="\t")}
    assert sum(1 for r in rows.values() if r[16]) == len(mm)
    pid = next(iter(mm))
    assert rows[pid][16] == [int(mm[pid]["qstart"]), int(mm[pid]["qend"]), round(100 * float(mm[pid]["fident"]), 1), mm[pid]["target"]]
    html = (out / "index.html").read_text()
    n = summ["total"]
    assert ("%s of %s proteins (%d%%) match a Swiss-Prot protein" % ("{:,}".format(summ["swissprot_hit"]), "{:,}".format(n),
                                                                     round(100 * summ["swissprot_hit"] / n))) in html
    assert 'href="assets/dense-bhf.css' in html and (out / "assets" / "index.css").exists()
    # 17: the other programs with a hit, as bits in the order of data["other"]; bit 0 is MMseqs2
    assert [o[0] for o in data["other"]][0] == "MMseqs2"
    assert all(bool(r[17] & 1) == (r[16] is not None) for r in data["rows"])


def test_protein_html_is_the_dense_page(site):
    """protein.html, the page of every protein without a drawn page, is the dense page: the same
    header, grid and Swiss-Prot table as a drawn page, with the constants protein-dense.js builds
    the page data from (the tracks in build.py's order, the InterProScan predictions, the cutoffs)."""
    out, _ = site
    import re
    html = (out / "protein.html").read_text()
    for el in ('<header class="top">', 'id="ladder"', '<div id="grid"></div>', 'id="biomni" data-base="" data-undrawn="1"',
               "assets/protein-dense.js", '<section class="below" id="swissprot">', '<aside class="help" id="help">'):
        assert el in html, el
    const = json.loads(re.search(r"window\.DENSE_CONST=(.*?);</script>", html).group(1))
    assert list(const["tracks"]) == [t["key"] for t in build.TRACKS] + ["diamond", "unifire"]
    assert const["struct_max_aa"] == build_full.STRUCT_MAX_AA and const["fantasia_max"] == build.FANTASIA_MAX_DIST
    assert const["ips_drawn"]["MobiDBLite"] == {"*": "disorder"}
    assert (out / "assets" / "protein-dense.js").exists() and (out / "assets" / "dense-bhf.js").exists()


def test_read_proteome_names_isoforms_and_checks_md5(tmp_path):
    """--proteome-fasta: the second copy of a gene's header is FUN1.2, as make_proteome_lengths.py names
    it, and a sequence whose length or MD5 differs from proteome_lengths.tsv stops the build."""
    import hashlib
    fa = tmp_path / "p.fasta"
    fa.write_text(">FUN1_A\nMKV\n>FUN1_A\nMKVL*\n>FUN2_B\nMA\n")
    md5 = lambda s: hashlib.md5(s.encode()).hexdigest()
    lengths = {"FUN1": (3, md5("MKV")), "FUN1.2": (4, md5("MKVL")), "FUN2": (2, md5("MA"))}
    assert build_full.read_proteome(fa, lengths) == {"FUN1": "MKV", "FUN1.2": "MKVL", "FUN2": "MA"}
    with pytest.raises(SystemExit):
        build_full.read_proteome(fa, dict(lengths, FUN2=(2, md5("MG"))))


def test_names_conflict():
    assert not build_full.names_conflict("NLRC5", "Nlrc5-like")
    assert not build_full.names_conflict("LET70", "let-70")
    assert build_full.names_conflict("ISM2", "adt-2")
    assert build_full.names_conflict("LRP2", "Cd320")
    assert not build_full.names_conflict("", "Cd320")


def test_unifire_rule_points_at_its_signatures():
    rules = {"R1": ({"IPR1": [(0, 0)], "PF1": [(0, 1)], "IPR9": [(1, 0)]}, [2, 1])}
    m = [["Pfam", "PF1", "", 5, 9, "", "", ""], ["CDD", "cd1", "", 1, 3, "", "IPR1", ""], ["CDD", "cd2", "", 7, 8, "", "IPR2", ""]]
    # set 0 needs both conditions and has them; set 1 does not hold; a row with a position keeps it
    assert build_full.rule_places([["ARBA", "R1", "keyword", "x", None, None], ["UniRule", "R1", "site", "y", 4, 4]], m, rules) == [
        ["ARBA", "R1", "keyword", "x", None, None, [0, 1]], ["UniRule", "R1", "site", "y", 4, 4]]
    assert build_full.rule_places([["ARBA", "R1", "keyword", "x", None, None]], m[2:], rules)[0][6] == []
    # no set holds in full: the nearest one still says where the rule points
    assert build_full.rule_places([["ARBA", "R1", "keyword", "x", None, None]], m[:1], rules)[0][6] == [0]


def test_unifire_hits_one_per_rule_and_place():
    m = [["Pfam", "PF1", "", 5, 9, "", "IPR1", ""]]
    p = [["ARBA", "R1", "keyword", "x", None, None, [0]], ["ARBA", "R1", "go", "y", None, None, [0]],
         ["UniRule", "R2", "site", "z", 4, None], ["ARBA", "R3", "keyword", "w", None, None, []]]
    h = {tuple(x["q"]): x for x in build_full.unifire_hits(p, m, 50)}
    assert set(h) == {(5, 9), (4, 4), (1, 50)}
    assert h[(5, 9)]["l"] == "x +1" and h[(5, 9)]["n"] == "x; y" and h[(5, 9)]["s"][3] == "PF1 / IPR1"
    assert h[(1, 50)]["s"][3].startswith("whole protein")


def test_gene_groups_match_the_main_name_only():
    r = dict(sp_acc="P1", sp_evalue="1e-30", sp_name="Calnexin")
    uni, aging = {"P1": dict(gene="Canx"), "P3": dict(gene="CD74")}, {"canx"}
    mm = {"x": ["P2", 0.3, 1, 9, "2e-5", "50", "Isoform 2 of Calnexin (Major histocompatibility complex class I antigen-binding protein p88)", ""]}
    km = {"x": ["hp2, k = 20", 0.3, "Q9Z301", "Period circadian protein homolog 2", "Rattus norvegicus", 35]}
    fs = {"x": ["P3", 0.2, 1, 9, "4.0", "30", "Isoform p41 of HLA class II histocompatibility antigen gamma chain (HLA-DR antigens-associated invariant chain)", ""]}
    got = build_full.gene_groups(r, "x", uni, aging, mm, fs, None, km)
    names = [g for g, _ in build_full.GENE_GROUPS]
    # Calnexin's other name does not make it an MHC protein; a weak hit counts; best E-value first
    assert [(names[g], prog, e) for g, prog, e, _, _ in got] == [
        ("aging (GenAge human genes)", "DIAMOND", 1e-30), ("Period", "kmerseek", 0.3),
        ("HLA and other MHC", "ProstT5 + Foldseek", 4.0), ("MHC class II", "ProstT5 + Foldseek", 4.0)]
    assert got[0][4] == "Calnexin, gene Canx"  # the symbol matches whatever its capitals
    assert build_full.gene_groups(r, "y", uni, set(), None, None, None, None) == []
    class_one = build_full.gene_groups(dict(r, sp_name="HLA class I histocompatibility antigen, alpha chain E"), "y", uni, set(), None, None, None, None)
    assert [names[g] for g, *_ in class_one] == ["HLA and other MHC", "MHC class I"]
    capital = build_full.gene_groups(dict(r, sp_name="Class I histocompatibility antigen, Gogo-C*0201 alpha chain"), "y", uni, set(), None, None, None, None)
    assert [names[g] for g, *_ in capital] == ["HLA and other MHC", "MHC class I"]


def test_genage_table_is_the_307_human_genes():
    assert len(build_full.read_genage(FULL)) == 307 and {"tp53", "ghr", "sirt1"} <= build_full.read_genage(FULL)


def test_gene_group_counts_are_in_the_summary_and_the_page(site):
    out, s = site
    assert list(s["gene_groups"]) == [g for g, _ in build_full.GENE_GROUPS]
    html = (out / "index.html").read_text()
    assert 'id="grp"' in html and all(build_full.build.esc(g) in html for g in s["gene_groups"])


def test_allorecognition_group_is_the_best_hit_of_each_of_the_53():
    allo = build_full.read_allorecognition(FULL)
    assert len(allo) == 53 == len({f[0] for f in csv.reader(open(FULL / "botryllus_allorecognition_hits.tsv"), delimiter="\t")})
    assert allo["FUN006619"] == [2.27e-232, "XTL10616.1", allo["FUN006619"][2]] and "Fester" in allo["FUN006619"][2]
    names = [g for g, _ in build_full.GENE_GROUPS]
    r = dict(sp_acc="", sp_evalue="", sp_name="")
    got = build_full.gene_groups(r, "FUN006619", {}, set(), None, None, None, None, allo)
    assert [(names[g], prog, acc) for g, prog, _, acc, _ in got] == [("allorecognition (Fester, Fester co-receptor, BHF)", "DIAMOND against Botryllus proteins", "XTL10616.1")]
    assert build_full.gene_groups(r, "FUN000001", {}, set(), None, None, None, None, allo) == []


def test_isoform_counts_on_the_index(site):
    out, s = site
    ids = [r["protein_id"] for r in tsv("proteome_lengths.tsv")]
    per_gene = Counter(i.split(".")[0] for i in ids)
    assert s["genes"] == len(per_gene) and s["genes_with_isoforms"] == sum(1 for k in per_gene.values() if k > 1)
    assert s["most_isoforms"] == max(per_gene.values())
    # every isoform after the first is named gene.N, so the counts and the IDs agree
    assert sum(k - 1 for k in per_gene.values()) == sum(1 for i in ids if "." in i)
    assert "have 2 to %d isoforms" % s["most_isoforms"] in (out / "index.html").read_text()
    assert s["genes_by_isoforms"] == dict(sorted(Counter(per_gene.values()).items()))
    html = (out / "index.html").read_text()
    assert html.count('data-iso="') == len(s["genes_by_isoforms"]) and "Genes by number of proteins (isoforms)" in html


def test_every_citation_is_on_the_index_with_its_doi(site):
    out, _ = site
    html = (out / "index.html").read_text()
    dois = [c[7] for c in build_full.CITATIONS]
    assert len(set(dois)) == len(dois) == 20 and all(d.startswith("10.") for d in dois)
    assert all('href="https://doi.org/%s"' % d in html for d in dois)
    assert all(build_full.build.esc(what) in html for what, _ in build_full.NO_PAPER)


def test_text_that_changes_on_hover_has_a_fixed_place_in_the_header():
    """Pointing at a hit fills in the overlap note and, for a kmerseek hit, changes the alphabet's
    classes in the key. Both sit in the fixed header, so each must keep one line whatever its text:
    hidden-while-empty or free to wrap, it made the header taller and moved every row (2026-11-10).
    tests/check_hover.py measures the same thing in a browser."""
    css = (Path(build_full.__file__).parent / "assets" / "dense-bhf.css").read_text()
    assert "#overlap:empty" not in css
    rule = lambda sel: next(l for l in css.split("\n") if l.startswith(sel + "{"))
    note = rule(".filters #overlap,.sentence #overlap")
    assert all(w in note for w in ("height:17px", "white-space:nowrap", "overflow:hidden"))
    key = rule(".key #alphakey")
    assert all(w in key for w in ("flex:0 0 100%", "height:16px", "flex-wrap:nowrap", "overflow:hidden"))


def test_pointing_at_a_hit_moves_nothing_on_protein_html(site):
    """In a browser, at three window widths: on protein.html, the page of every protein without a
    drawn page, pointing at any hit box leaves the fixed header and every row where they were, and
    the labels of the dashed prediction boxes lie over no box, no other label and not off the edge.
    FUN004057 has a narrow disorder box and a coiled coil at its last residues; FUN005955.2 has
    five predictions and kmerseek hits; FUN000001 has no prediction."""
    out, _ = site
    pages = ["protein.html#FUN004057", "protein.html#FUN005955.2", "protein.html#FUN000001", "protein.html#FUN001586"]
    res = check_hover.check(out, pages)
    assert len(res) == len(pages) * len(check_hover.WIDTHS)
    by = {(pg, w): r for pg, w, r in res}
    assert all(r["tried"] == r["hits"] > 0 for r in by.values())  # under MAX_HITS, so every box was tried
    assert by[("protein.html#FUN004057", 1400)]["labels"] == 2 and by[("protein.html#FUN005955.2", 1400)]["labels"] == 5
    assert not check_hover.failures(res), check_hover.failures(res)


def test_a_click_on_a_bar_or_a_domain_brings_up_its_hits(site):
    """In a browser, on protein.html: a click on a bar of a program's histogram keeps only the hits
    that reach the bar's residues and unfolds that program; a click on a box of the domain lane does
    the same for the box and unfolds every program with such a hit. Clearing folds them again. The
    link back to all proteins is the first thing in the title bar, at the top left."""
    r = check_hover.clicks(site[0], "protein.html#FUN005955.2")
    assert r["before"] == 0  # every program starts folded
    assert r["back"]["text"].endswith("All proteins") and r["back"]["href"] == "index.html"
    assert r["back"]["first"] and r["back"]["left"] < 40 and r["back"]["top"] < 40
    for key in ("bar", "dom"):
        c = r[key]
        assert c["rows"] > 0 and len(c["ranges"]) == c["rows"] and all(c["ranges"])
        assert all(lo <= c["b"] and hi >= c["a"] for lo, hi in c["ranges"]), key  # every hit shown reaches the residues clicked
        assert "%d\u2013%d" % (c["a"], c["b"]) in c["count"]
    assert r["bar"]["open"] == [r["bar"]["k"]]
    assert len(r["dom"]["open"]) >= 1 and r["dom"]["on"] and r["dom"]["l"] in r["dom"]["count"]
    for key in ("barCleared", "domCleared"):
        assert r[key]["rows"] == 0 and r[key]["open"] == [] and "\u2013" not in r[key]["count"]


def test_a_protein_only_kmerseek_found_shows_kmerseek_e_value(site):
    """The first page's E-value column is DIAMOND's. A protein that only kmerseek found has no DIAMOND
    match, so the page writes kmerseek's E-value there with the program's name under it."""
    js = (HERE / "assets" / "full.js").read_text()
    assert "!r[9] && r[19] ? '<div class=\"num ekm\"" in js and '<span class="s">kmerseek</span>' in js
    F = read_js(site[0] / "full" / "rows.js", "window.FULL=")
    only = [r for r in F["rows"] if r[19]]
    assert len(only) == site[1]["only_kmerseek_3"] > 0
    assert all(not r[9] and 0 <= r[19][0] <= build_full.ONLY_KMERSEEK_E for r in only)  # none has a DIAMOND E-value to hide


def test_s4pred_table_has_every_protein_and_no_residues(site):
    """s4pred_states.tsv.gz gives each of the 45,339 proteins its helix, strand and coil as runs that
    add up to its length, and nothing else: no residue letter is in the file. The shares of all
    residues are those of the run's own summary (olgabot/botryllus-mhc PR 13: 31.5, 11.4, 57.1)."""
    import re
    csv.field_size_limit(10 ** 9)
    L = {r["protein_id"]: int(r["length"]) for r in tsv("proteome_lengths.tsv")}
    rows = tsv("s4pred_states.tsv")
    assert [r["protein_id"] for r in rows] == list(L) and len(rows) == 45339
    for r in rows:
        runs = re.findall(r"([HEC])(\d+)", r["states"])
        assert "".join(c + n for c, n in runs) == r["states"], r["protein_id"]  # runs only, so no sequence
        assert sum(int(n) for _, n in runs) == L[r["protein_id"]], r["protein_id"]
        assert all(a[0] != b[0] for a, b in zip(runs, runs[1:])), r["protein_id"]  # neighbouring runs differ
    s = site[1]["s4pred"]
    assert (s["proteins"], s["helix_pct"], s["strand_pct"], s["coil_pct"]) == (45339, 31.5, 11.4, 57.1)
    assert s["residues"] == sum(L.values())
    # each protein's page data carries its own runs
    by = {r["protein_id"]: r["states"] for r in rows}
    for pid in ("FUN000001", "FUN005955.2", "FUN019272.2"):
        d = read_js(site[0] / "full" / ("%03d.js" % build_full.shard_of(pid)), "")
        assert d[pid]["ss"] == by[pid]


def test_s4pred_strip_is_drawn_under_the_ruler(site):
    """In a browser: protein.html has a row right under the ruler, named "helix and strand, predicted"
    with s4pred and this protein's shares at the left, with one bar per helix and strand run; the
    key says what the two bars are."""
    import re
    r = check_hover.s4pred(site[0], "protein.html#FUN005955.2")
    rows = {x["protein_id"]: x["states"] for x in tsv("s4pred_states.tsv")}
    runs = re.findall(r"([HEC])(\d+)", rows["FUN005955.2"])
    assert r["helix"] == sum(1 for c, _ in runs if c == "H") > 0 and r["strand"] == sum(1 for c, _ in runs if c == "E") > 0
    assert r["inside"] and r["afterRuler"] and "helix" in r["key"] and "strand" in r["key"] and "s4pred" in r["key"]
    L = sum(int(n) for _, n in runs)
    assert "helix and strand, predicted" in r["name"] and "s4pred" in r["name"]
    assert "helix %d%%" % round(100 * sum(int(n) for c, n in runs if c == "H") / L) in r["name"]


def test_a_row_with_no_diamond_match_still_says_what_else_was_found(site):
    """The first page's Swiss-Prot column is DIAMOND's, so its empty state says "no DIAMOND match",
    not "no match", and the row carries the Foldseek best hits weaker than E-value 0.001 that the
    protein's own page shows (FUN000001: ProstT5 + Foldseek 0.69 and Foldseek on structure 5.1)."""
    js = (HERE / "assets" / "full.js").read_text()
    assert "no DIAMOND match</span>" in js and ">no match<" not in js
    F = read_js(site[0] / "full" / "rows.js", "window.FULL=")
    names = [o[0] for o in F["other"]]
    r = next(r for r in F["rows"] if r[0] == "FUN000001")
    assert not r[6] and r[17] >> names.index("kmerseek") & 1  # no DIAMOND match, kmerseek has a hit
    assert [(names[i], e) for i, e in r[21]] == [("ProstT5 + Foldseek", 0.6941), ("Foldseek on structure", 5.123)]
    # field 21 holds only hits weaker than the cutoff, and never a program that field 17 already counts
    for r in F["rows"]:
        for i, e in r[21] or []:
            assert e > build_full.FOLDSEEK_MAX_E and not r[17] >> i & 1, r[0]



def test_immune_programs_panel(site):
    site, _ = site
    """The panel's data has every protein of immune_candidates_by_program.tsv.gz, each new-proteome ID has a page,
    and the counts drawn on the index are the table's."""
    raw = (site / "full" / "immune_programs.js").read_text()
    d = json.loads(gzip.decompress(base64.b64decode(raw.split('"')[1])))
    keys = [k for k, _, _ in build_full.IMMUNE_PROGRAMS]
    assert [p[0] for p in d["programs"]] == keys
    table = tsv("immune_candidates_by_program.tsv.gz")
    assert len(d["rows"]) == len(table) == 10948
    lengths = {r["protein_id"] for r in tsv("proteome_lengths.tsv")}
    assert all(r[0] in lengths for r in d["rows"] if r[1])
    assert not any(r[0] in lengths for r in d["rows"] if not r[1])
    html = (site / "index.html").read_text()
    for i, k in enumerate(keys):
        n = sum(1 for r in table if r[k])
        assert sum(1 for r in d["rows"] if r[2] >> i & 1) == n
        assert '<span class="ip-n mono">%s</span>' % "{:,}".format(n) in html, k
    for r in d["rows"]:  # one piece of evidence per program that calls it immune
        assert len(r[3]) == bin(r[2]).count("1")
    assert "1,105 by three or more" in html
