#!/usr/bin/env python3
"""Pick boxes at random from a built site and check them against the raw input lines.

    python tests/check_boxes.py SITE_DIR INPUT_DIR [N] [SEED]

For each box: find the line of the tool's own file that it came from, and print the query
range, target range and score from that line next to what the page holds. For a box with an
alignment, count identical residues again from the two aligned rows, and check that the
aligned residues, with gaps removed, are the residues of proteins.fasta and swissprot.fasta at
the stated coordinates.
"""
import json
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build  # noqa: E402


def main():
    site, inp = Path(sys.argv[1]), Path(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    prots = {h.split()[0]: s for h, s in build.read_fasta(inp / "proteins.fasta").items()}
    sp = {build.accession(h.split()[0]): s for h, s in build.read_fasta(inp / "swissprot.fasta").items()}
    boxes = []
    for page in sorted((site / "proteins").glob("*.html")):
        m = re.search(r"const PAGE=(\{.*?\});</script>", page.read_text())
        data = json.loads(m.group(1).replace("<\\/", "</"))
        aln_file = page.with_name(page.stem + ".aln.js")
        alns = json.loads(aln_file.read_text()[len("window.ALN="):-1]) if aln_file.exists() else {}
        for i, h in enumerate(data["hits"]):
            if h.get("a"):
                h = dict(h, a=alns.get(str(i)))
            if h["k"] in ("mmseqs2_3rounds", "phmmer", "jackhmmer", "kmerseek") and i < 200:
                boxes.append((page.stem, h, data))
    fails = []
    for pid, h, data in rng.sample(boxes, n):
        print("=" * 70)
        print(pid, h["k"], h["l"], "page: query %s target %s %s scores %s" % (h["q"], h["t"], h["r"], h["s"]))
        f = {"mmseqs2_3rounds": "mmseqs2_3rounds.m8", "phmmer": "phmmer.domtblout",
             "jackhmmer": "jackhmmer.domtblout", "kmerseek": "kmerseek.tsv"}[h["k"]]
        for line in build.opentext(inp / f):
            if pid not in line or h["t"] not in line:
                continue
            c = line.split("\t") if f.endswith((".m8", ".tsv")) else line.split()
            if h["k"] == "mmseqs2_3rounds":
                got = ([int(c[4]), int(c[5])], [int(c[6]), int(c[7])], c[8])
            elif h["k"] == "kmerseek":
                continue
            else:
                got = ([int(c[15]), int(c[16])], [int(c[17]), int(c[18])], c[12])
            if got[0] == h["q"] and got[1] == h["r"]:
                print("input line: query %s target %s E %s  -> matches the page" % got)
                break
        else:
            if h["k"] != "kmerseek":
                print("input line: NOT FOUND")
                fails.append((pid, h["l"], "no input line"))
        a = h.get("a")
        if a:
            ident = sum(1 for x, y in zip(a["q"], a["t"]) if x != "-" and x.upper() == y.upper())
            qres, tres = a["q"].replace("-", "").upper(), a["t"].replace("-", "").upper()
            q_ok = prots[pid][a["qs"] - 1:a["qs"] - 1 + len(qres)].upper() == qres
            t_ok = sp[h["t"]][a["ts"] - 1:a["ts"] - 1 + len(tres)].upper() == tres
            inside = a["qs"] >= h["q"][0] and a["qs"] + len(qres) - 1 <= h["q"][1]
            print("alignment (%s): %d of %d identical (recounted); query %d-%d matches proteins.fasta: %s; "
                  "target %d-%d matches swissprot.fasta: %s; inside the box: %s" % (
                      a.get("src", "")[:40], ident, len(a["q"]), a["qs"], a["qs"] + len(qres) - 1, q_ok,
                      a["ts"], a["ts"] + len(tres) - 1, t_ok, inside))
            if not (q_ok and t_ok and inside):
                fails.append((pid, h["l"], "alignment"))
        elif h["k"] == "kmerseek":
            q = prots[pid][h["q"][0] - 1:h["q"][1]]
            t = sp[h["t"]][h["r"][0] - 1:h["r"][1]]
            for line in open(inp / f):
                if line.startswith(pid + "\t") and "\t" + h["t"] + "\t" in line:
                    head = open(inp / f).readline().rstrip("\n").split("\t")
                    d = dict(zip(head, line.rstrip("\n").split("\t")))
                    print("input line: region %s-%s target %s-%s (0-based, end excluded); region_subseq "
                          "equals the page's query residues: %s; target_subseq equals: %s; identical %d of %d" % (
                              d["region_start"], d["region_end"], d["target_start"], d["target_end"],
                              d["region_subseq"] == q, d["target_subseq"] == t,
                              sum(x == y for x, y in zip(q, t)), len(q)))
                    if not (d["region_subseq"] == q and d["target_subseq"] == t):
                        fails.append((pid, h["l"], "kmerseek residues"))
        else:
            print("no alignment stored (beyond the best 5; the page aligns it in the browser)")
    print("=" * 70)
    print("%d of %d boxes checked failed%s" % (len(fails), n, (": %s" % fails) if fails else ""))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
