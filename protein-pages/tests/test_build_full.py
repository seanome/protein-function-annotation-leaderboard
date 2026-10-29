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
    s = path.read_text()
    assert s.startswith(prefix)
    return json.loads(s[len(prefix):-1]) if prefix == "window.FULL=" else json.loads(s[s.index(",", s.index("||{}")) + 1:-2])


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
    assert s["immune_candidates"] == len(tsv("botryllus_immune_genes.tsv"))
    assert s["allorecognition_hit"] == len({f[0] for f in csv.reader(open(FULL / "botryllus_allorecognition_hits.tsv"), delimiter="\t")})
    assert s["clusters"] == len({f[0] for f in csv.reader(open(FULL / "botryllus_mmseqs_clusters.tsv"), delimiter="\t")})
    assert sum(s["name_check"].values()) == s["total"]
    want = Counter(c for r in tsv("botryllus_immune_genes.tsv") for c in build_full.split_list(r["categories"]))
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
