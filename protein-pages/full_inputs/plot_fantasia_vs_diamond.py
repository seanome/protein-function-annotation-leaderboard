"""Is FANTASIA's distance a usable cutoff? Plots it against DIAMOND's E-value and agreement.

    python plot_fantasia_vs_diamond.py . ../figures/fantasia_vs_diamond

a: each protein with both a FANTASIA answer and a DIAMOND Swiss-Prot hit (E-value 0.001 or
   below), at its FANTASIA distance and DIAMOND E-value, filled when the two name a gene of the
   same family (gene name without its trailing number and letter, as fantasia_vs_diamond.py).
b: in windows of 0.025 in distance (each point at the median distance of its window), the share
   of proteins whose FANTASIA GO terms share at least one GO term with the DIAMOND hit's (any term,
   so broad ones such as "cytoplasm" count), and the share whose two answers name a gene of the
   same family. Windows with fewer than 20 proteins are left out.
c: how FANTASIA's distance is spread over all 26,676 proteins it was run on, split into those with
   and without a DIAMOND hit. BHF (FUN008084) is marked.
"""
import csv
import gzip
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "gene_cards"))
import pubfig as pf  # noqa: E402

inp, out = Path(sys.argv[1]), sys.argv[2]
gene = {}
for line in gzip.open(inp / "uniprot_names.tsv.gz", "rt"):
    acc, _n, g, _o = line.rstrip("\n").split("\t")
    if g:
        gene[acc] = g.upper()
family = lambda g: re.sub(r"\d+[A-Z]?$", "", g)

both, d_hit, d_nohit, bhf = [], [], [], None
for r in csv.DictReader(gzip.open(inp / "botryllus_protein_annotations_v2.tsv.gz", "rt"), delimiter="\t"):
    if not r["fantasia_distance"]:
        continue
    d = float(r["fantasia_distance"])
    if r["protein_id"].startswith("FUN008084"):
        bhf = d
    if r["sp_acc"]:
        d_hit.append(d)
        if r["sp_acc"] in gene and r["fantasia_uniprot"] in gene:
            a, b = gene[r["sp_acc"]], gene[r["fantasia_uniprot"]]
            ga = {x.strip() for x in r["go_terms"].split(";") if x.strip()}
            gb = {x.split("|")[0].strip() for x in r["fantasia_go_terms"].split(";") if x.strip()}
            go_ok = bool(ga & gb) if ga and gb else np.nan
            both.append((d, max(float(r["sp_evalue"]), 1e-180), a == b, family(a) == family(b), go_ok))
    else:
        d_nohit.append(d)
both = np.array(both, dtype=float)
dist, ev, same, fam, goo = both.T

pf.use_style()
fig, axs = pf.figure(pf.TWO_COLUMN_MM, 62, ncols=3)
blue, ver = pf.OKABE_ITO["blue"], pf.OKABE_ITO["vermillion"]

ax = axs[0]
ax.scatter(dist[fam == 0], ev[fam == 0], s=5, facecolor="none", edgecolor=pf.GREY, linewidth=0.5,
           label="different gene family")
ax.scatter(dist[fam == 1], ev[fam == 1], s=5, color=blue, linewidth=0, label="same gene family")
ax.set_yscale("log")
ax.set_yticks([1e-180, 1e-150, 1e-100, 1e-50, 1e-3])
ax.set_yticklabels(["0 (as reported)", "1e-150", "1e-100", "1e-50", "0.001"])
ax.set_ylim(3e-181, 1e-1)
ax.invert_yaxis()
ax.set_xlabel("FANTASIA distance (0 = identical)")
ax.set_ylabel("DIAMOND E-value (smaller = stronger)")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, handletextpad=0.2, ncol=2, markerscale=2)

ax = axs[1]
edges = np.arange(0, 0.45001, 0.025)
xs, g_pct, f_pct, go_pct, ns = [], [], [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (dist >= lo) & (dist < hi)
    if m.sum() >= 20:
        xs.append(np.median(dist[m])); g_pct.append(100 * same[m].mean()); f_pct.append(100 * fam[m].mean()); ns.append(m.sum())
        go_pct.append(100 * np.nanmean(goo[m]))
green = pf.OKABE_ITO["bluish_green"]
ax.plot(xs, go_pct, "-o", color=green, ms=3, label="share a GO term")
ax.plot(xs, f_pct, "-o", color=blue, ms=3, label="same gene family")
if bhf is not None:
    ax.axvline(bhf, color="#1a1a1a", linewidth=0.5, linestyle=":")
    ax.text(bhf + 0.005, 92, "BHF, %.2f" % bhf)
ax.set_ylim(0, 100)
ax.set_xlabel("FANTASIA distance (0 = identical)")
ax.set_ylabel("Proteins whose FANTASIA pick\nmatches DIAMOND's hit (%)")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, ncol=2)

ax = axs[2]
bins = np.arange(0, 0.62, 0.01)
ax.hist(d_nohit, bins=bins, color=pf.GREY, label="no DIAMOND hit (%d)" % len(d_nohit))
ax.hist(d_hit, bins=bins, color=blue, label="DIAMOND hit (%d)" % len(d_hit))
if bhf is not None:
    ax.axvline(bhf, color="#1a1a1a", linewidth=0.5, linestyle=":")
    ax.text(bhf + 0.008, ax.get_ylim()[1] * 0.92, "BHF")
ax.set_xlabel("FANTASIA distance (0 = identical)")
ax.set_ylabel("Proteins")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, ncol=1)

for a, letter in zip(axs, "abc"):
    pf.panel_label(a, letter)
pf.save(fig, out)
print("pairs %d; windows (median distance, n, same gene %%, same family %%):" % len(dist))
for x, n, g, f in zip(xs, ns, g_pct, f_pct):
    print("%.3f\t%d\t%.0f\t%.0f" % (x, n, g, f))
print("GO-term share per window:", ["%.0f" % v for v in go_pct])
from scipy.stats import spearmanr  # noqa: E402
print("Spearman distance vs log10 E-value: %.2f" % spearmanr(dist, np.log10(ev)).correlation)
