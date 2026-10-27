#!/usr/bin/env python3
"""Checks on a built site.

    python tests/check_site.py SITE_DIR

1. Every colour drawn on a protein page is one the legend explains (family colours, the protein
   line, the agreement line, pLDDT and expression colours) or white / text ink.
2. The five numbers on the index add up: no hit + exactly one + agree + disagree = total.
3. Each protein page's sentence agrees with its row in the index table.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build  # noqa: E402

LEGEND = {c.lower() for c in [f["color"] for f in build.FAMILIES.values()]
          + [build.PROTEIN_LINE, build.AGREEMENT, build.EXPRESSION] + [p[2] for p in build.PLDDT]}
NEUTRAL = {"#fff", "#ffffff", "#1a1a1a", "none"}


def main():
    site = Path(sys.argv[1])
    bad = {}
    for page in sorted((site / "proteins").glob("*.html")):
        s = page.read_text()
        used = {c.lower() for c in re.findall(r'(?:fill|stroke)="([^"]+)"', s)}
        js = (site / "assets" / "protein.js").read_text()
        used |= {c.lower() for c in re.findall(r'(?:fill|stroke)="(#[0-9A-Fa-f]{3,6})"', js)}
        extra = used - LEGEND - NEUTRAL
        if extra:
            bad[page.name] = extra
    print("1. colours outside the legend:", bad or "none")

    idx = (site / "index.html").read_text()
    nums = [int(x.replace(",", "")) for x in re.findall(r'<span class="big">([\d,]+)</span>', idx)]
    total, none, one, agree, disagree = nums
    print("2. index numbers %s: %d + %d + %d + %d = %d, total %d -> %s" % (
        nums, none, one, agree, disagree, none + one + agree + disagree, total,
        "OK" if none + one + agree + disagree == total else "FAIL"))

    rows = re.findall(r'<tr data-s="(\w+)" data-f="([^"]*)"><td class="mono"><a href="proteins/([^"]+)\.html">'
                      r'.*?</a></td><td class="r" data-v="\d+">[^<]*</td><td class="r" data-v="(\d+)">', idx)
    mism = []
    for status, fams, pid, n in rows:
        s = (site / "proteins" / (pid + ".html")).read_text()
        sent = re.search(r'<p class="summary">(.*?)</p>', s).group(1)
        n = int(n)
        if n == 0:
            ok = sent == "No tool reports a hit." and status == "none"
        else:
            ok = sent.startswith("%d of 8 tool families" % n)
            if status == "agree":
                ok = ok and "They agree on" in sent
            elif status == "disagree":
                ok = ok and ("do not agree" in sent or "no Pfam family is known" in sent)
            elif status == "one":
                ok = ok and n == 1
        if not ok:
            mism.append((pid, status, n, sent))
    print("3. pages whose sentence disagrees with the index (%d rows read):" % len(rows), mism or "none")
    sys.exit(1 if bad or mism or none + one + agree + disagree != total or len(rows) != total else 0)


if __name__ == "__main__":
    main()
