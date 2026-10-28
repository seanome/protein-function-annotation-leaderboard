#!/usr/bin/env python3
"""Checks on a built site.

    python tests/check_site.py SITE_DIR

1. Every colour drawn on a protein page is one the legend explains (program colours, the protein
   line, the agreement line, pLDDT and expression colours) or white / text ink.
2. The six class numbers on the index add up to the total: no hit, exactly one program, agree,
   disagree, overlap but no family known, apart. summary.json holds the same numbers.
3. Each protein page's sentence agrees with its row in the index table, including "N of the M
   programs run" with M the programs that have an input file.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build  # noqa: E402

LEGEND = {c.lower() for c in [f["color"] for f in build.PROGRAMS.values()]
          + [build.PROTEIN_LINE, build.AGREEMENT, build.EXPRESSION] + [p[2] for p in build.PLDDT]}
NEUTRAL = {"#fff", "#ffffff", "#1a1a1a", "none", "#c3c2b7"}  # #c3c2b7: frame of the residue panel's map


STATUS_WORDS = {"none": "", "one": "", "agree": "They agree on",
                "disagree": "name different Pfam families",
                "overlap_unknown": "no Pfam family is known",
                "apart": "do not overlap"}


def SENT_LEAD(n, m):
    """How the sentence starts for n programs with a hit, of m run (build.header_sentence)."""
    if n == 0:
        return "No tool reports a hit." if m == 8 else "None of the %d programs run" % m
    return "%d %s" % (n, build.of_programs(m))


def main():
    site = Path(sys.argv[1])
    bad = {}
    for page in sorted((site / "proteins").glob("*.html")):
        s = page.read_text()
        used = {c.lower() for c in re.findall(r'(?:fill|stroke)="([^"]+)"', s)}
        js = (site / "assets" / "protein.js").read_text()
        used |= {c.lower() for c in re.findall(r'(?:fill|stroke)="(#[0-9A-Fa-f]{3,6})"', js)}
        extra = {c for c in used - LEGEND - NEUTRAL if not c.startswith("url(#hatch-")}  # hatch: a kind's colour
        if extra:
            bad[page.name] = extra
    print("1. colours outside the legend:", bad or "none")

    idx = (site / "index.html").read_text()
    nums = [int(x.replace(",", "")) for x in re.findall(r'<span class="big">([\d,]+)</span>', idx)]
    total, classes = nums[0], dict(zip(build.STATUSES, nums[1:]))
    summary = json.loads((site / "summary.json").read_text())
    sums_ok = len(nums) == 1 + len(build.STATUSES) and sum(classes.values()) == total
    sums_ok = sums_ok and all(summary[k] == v for k, v in classes.items()) and summary["total"] == total
    print("2. index numbers %s: %s = %d, total %d, summary.json the same -> %s" % (
        nums, " + ".join("%d %s" % (v, k) for k, v in classes.items()), sum(classes.values()), total,
        "OK" if sums_ok else "FAIL"))

    rows = re.findall(r'<tr data-n="\d" data-s="(\w+)" data-f="([^"]*)"><td><a class="mono" href="proteins/([^"]+)\.html">'
                      r'.*?</td><td class="r" data-v="\d+">[^<]*</td><td class="r" data-v="(\d+)">', idx)
    if len(rows) != total:  # a parse that reads no rows must not pass
        sys.exit("3. read %d table rows, expected %d: the row pattern no longer fits index.html" % (len(rows), total))
    mism = []
    for status, progs, pid, n in rows:
        s = (site / "proteins" / (pid + ".html")).read_text()
        sent = re.search(r'<p class="summary">(.*?)</p>', s).group(1)
        n = int(n)
        lead = SENT_LEAD(n, summary["programs_run"])
        if n == 0:
            ok = sent.startswith(lead) and status == "none"
        else:
            ok = sent.startswith(lead) and (status == "one") == (n == 1)
            ok = ok and STATUS_WORDS[status] in sent
        if not ok:
            mism.append((pid, status, n, sent))
    print("3. pages whose sentence disagrees with the index (%d rows read):" % len(rows), mism or "none")
    sys.exit(1 if bad or mism or not sums_ok or len(rows) != total else 0)


if __name__ == "__main__":
    main()
