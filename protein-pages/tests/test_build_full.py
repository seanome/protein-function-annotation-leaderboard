"""Checks of build_full.py and of what it assumes about the Biomni tables, on the committed
tables in full_inputs/. Run with: python -m pytest tests/
"""
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
    """A data file's JSON; a shard (full/NNN.js) with its strings put back from full/str.js, as FULL_X does."""
    s = path.read_text()
    assert s.startswith(prefix)
    if prefix == "window.FULL=":
        return json.loads(s[len(prefix):-1])
    d = json.loads(s[s.index("FULL_X(") + len("FULL_X("):-3])
    st = path.parent / "str.js"
    S = json.loads(st.read_text().split(";window.FULL_X=")[0][len("window.FULL_STR="):])
    x = lambda i: S[i] if isinstance(i, int) else i
    for e in d.values():
        e["go"] = [[x(a), x(b)] for a, b in e["go"]]
        e["kw"] = [x(k) for k in e["kw"]]
        if isinstance(e.get("ips"), list):
            e["ips"] = [[x(m[0]), x(m[1]), x(m[2]), m[3], m[4], m[5], x(m[6]), x(m[7])] for m in e["ips"]]
        if isinstance(e.get("uf"), list):
            e["uf"] = [[x(u[0]), x(u[1]), x(u[2]), x(u[3]), u[4], u[5]] for u in e["uf"]]
        if e.get("fa"):
            e["fa"][6] = [[x(a), x(b)] for a, b in e["fa"][6]]
    return d


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
    """The best ProstT5 + Foldseek hit reaches its protein, kept at FOLDSEEK_MAX_E or below:
    25,804 proteins, as in the Biomni report (v3); BHF (FUN008084, best E-value 0.06) has none."""
    out, summ = site
    fs = {}
    for s in range(build_full.N_SHARDS):
        for pid, d in read_js(out / "full" / ("%03d.js" % s), "window.FULL_SHARD=").items():
            fs[pid] = d["fs"]
    assert summ["foldseek_hit"] == sum(1 for v in fs.values() if v) == 25804
    assert fs["FUN008084"] == 0
    raw = {r["protein_id"]: r for r in csv.DictReader(gzip.open(FULL / "foldseek_best_per_protein.tsv.gz", "rt"), delimiter="\t")}
    pid = next(p for p, v in fs.items() if v)
    assert build.af_accession(raw[pid]["target"]) == fs[pid][0] and int(raw[pid]["qstart"]) == fs[pid][2]


def test_foldseek_on_structures_on_the_pages(site):
    """The best hit of Foldseek on Biomni's predicted structures: 3,729 proteins at FOLDSEEK_MAX_E
    or below, all 200 aa or less; every longer protein says "not run", BHF (252 aa) among them."""
    out, summ = site
    fst, length = {}, {}
    for s in range(build_full.N_SHARDS):
        for pid, d in read_js(out / "full" / ("%03d.js" % s), "window.FULL_SHARD=").items():
            fst[pid] = d["fst"]
    for r in csv.DictReader(gzip.open(FULL / "proteome_lengths.tsv.gz", "rt"), delimiter="\t"):
        length[r["protein_id"]] = int(r["length"])
    assert summ["foldseek_struct_hit"] == sum(1 for v in fst.values() if v and v != "not run") == 3729
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
    assert list(const["tracks"]) == [t["key"] for t in build.TRACKS]
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
