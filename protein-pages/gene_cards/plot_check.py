"""Figure for the 100-protein gene-card test, from out/check.tsv (check_cards.py).

a: on the conflicts, the first model's name verdict (rows) against the second's (columns); the
   number in each square is the count of proteins.
b: the first model's confidence for each kind of protein in the test.

  python gene_cards/plot_check.py gene_cards/out claude-haiku-5-5 claude-opus-5-5
  python gene_cards/plot_check.py gene_cards/out DeepSeek-R1-Distill-Qwen-32B claude-opus-5-5 deepseek_check
"""
import collections
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pubfig as pf  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

VERDICTS = ["existing", "proposed", "same_gene", "family_only", "neither"]
VLAB = ["header name", "proposed name", "same gene", "family only", "neither"]
KINDS = ["conflict_different_gene", "conflict_like_only", "agree", "unconfirmed", "new", "immune", "fantasia_only"]
KLAB = ["conflict, two genes", "conflict, \"-like\" only", "names agree", "header name only",
        "proposed name only", "immune, no name", "FANTASIA only"]
CONF = ["high", "medium", "low", "none"]
NAMES = {"claude-haiku-5-5": "Haiku 5.5", "claude-opus-5-5": "Opus 5.5",
         "DeepSeek-R1-Distill-Qwen-32B": "DeepSeek-R1 32B"}


def heat(ax, m, rows, cols, xlabel, ylabel):
    cmap = plt.get_cmap("mako_r").copy()
    cmap.set_bad("white")  # an empty square is white, not the lightest count colour
    ax.imshow(np.ma.masked_equal(m, 0), cmap=cmap, vmin=0, vmax=max(1, m.max()), aspect="equal")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if m[i, j]:
                ax.text(j, i, str(int(m[i, j])), ha="center", va="center",
                        color="white" if m[i, j] > m.max() / 2 else "black")
    ax.set_xticks(range(len(cols)), cols, rotation=35, ha="right")
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def main(out, a, b, name="pilot_check"):
    rows = list(csv.DictReader(open(Path(out) / "check.tsv"), delimiter="\t"))
    pf.use_style()
    fig, (ax1, ax2) = pf.figure(pf.TWO_COLUMN_MM, 70, ncols=2, gridspec_kw=dict(width_ratios=[5, 4]))
    m = np.zeros((len(VERDICTS), len(VERDICTS)))
    for r in rows:
        if r["kind"].startswith("conflict") and r.get(a + ":verdict") in VERDICTS and r.get(b + ":verdict") in VERDICTS:
            m[VERDICTS.index(r[a + ":verdict"]), VERDICTS.index(r[b + ":verdict"])] += 1
    heat(ax1, m, VLAB, VLAB, NAMES.get(b, b) + " verdict", NAMES.get(a, a) + " verdict")
    c = collections.Counter((r["kind"], r[a + ":confidence"]) for r in rows if r.get(a + ":confidence"))
    m2 = np.array([[c[k, x] for x in CONF] for k in KINDS], float)
    heat(ax2, m2, KLAB, CONF, NAMES.get(a, a) + "\nconfidence", "Kind of protein")
    pf.panel_label(ax1, "a", dx_pt=-62)
    pf.panel_label(ax2, "b", dx_pt=-82)
    pf.save(fig, str(Path(out) / name))


if __name__ == "__main__":
    main(*sys.argv[1:])
