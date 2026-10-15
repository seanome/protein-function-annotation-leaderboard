"""Step 6: the same three labels on four proteins, drawn on their real residues.

For each protein: every InterPro match with its coordinates (InterPro API), every MMseqs2
hit to another Swiss-Prot entry above 30% identity with its query coordinates (the
protein's own entry is left out, since all four are Swiss-Prot entries), and the
AlphaFold DB per-residue pLDDT.

    python3 -I scripts/why/four.py --raw ~/data/why-numbers/2026_03 --sherlock ~/data/why-numbers/2026_03/sherlock
Writes results/why/four_proteins_calls.parquet and results/why/maps/four_<accession>.{svg,png,pdf}.
"""

import argparse
import json
import pathlib
import sys

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.patches import Patch, Rectangle

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import pubfig as pf  # noqa: E402
from labels import MIN_FIDENT, acc_of  # noqa: E402

FOUR = {
    "Q25460": "Mussel foot protein 1, Mytilus edulis",
    "P24856": "Antifreeze glycoprotein polyprotein, Notothenia neglecta",
    "P05484": "Omega-conotoxin MVIIA, Conus magus",
    "P42212": "Green fluorescent protein, Aequorea victoria",
}
PER_LINE = 100
# AlphaFold DB's own pLDDT bands and colours.
BANDS = [(90, 101, "#0053D6", "pLDDT > 90"), (70, 90, "#65CBF3", "70-90"), (50, 70, "#FFDB13", "50-70"), (0, 50, "#FF7D45", "< 50")]
IPR, HIT = "#009E73", "#CC79A7"  # Okabe-Ito green: InterPro match; reddish purple: Swiss-Prot hit


def interpro_calls(path, acc):
    d = json.loads(path.read_text())
    rows = []
    for r in d.get("results", []):
        m = r["metadata"]
        for p in r["proteins"]:
            for loc in p.get("entry_protein_locations") or []:
                for f in loc["fragments"]:
                    rows.append({
                        "accession": acc, "tool": "InterPro", "label": f"{m['accession']} {m['name']} ({m['source_database']})",
                        "start": f["start"], "end": f["end"], "score": loc.get("score"),
                    })
    return rows


def mmseqs_calls(path, acc):
    cols = ["query", "target", "fident", "alnlen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen"]
    hits = pl.read_csv(path, separator="\t", has_header=False, new_columns=cols).with_columns(
        query=acc_of("query"), target_acc=acc_of("target"),
    ).filter(pl.col("query") == acc)
    rows = []
    for h in hits.sort("bits", descending=True).iter_rows(named=True):
        rows.append({
            "accession": acc, "tool": "MMseqs2 vs Swiss-Prot",
            "label": f"{h['target']} {100 * h['fident']:.0f}% identity, E = {h['evalue']:.1e}" + (" (itself)" if h["target_acc"] == acc else ""),
            "start": h["qstart"], "end": h["qend"], "score": h["fident"], "self_hit": h["target_acc"] == acc,
            "above_30pct": h["fident"] > MIN_FIDENT,
        })
    return rows


def draw(acc, seq, calls, plddt, out):
    pf.use_style()
    ipr = [c for c in calls if c["tool"] == "InterPro"]
    hits = [c for c in calls if c["tool"].startswith("MMseqs2") and not c["self_hit"] and c["above_30pct"]]
    lines = [(i, seq[i : i + PER_LINE]) for i in range(0, len(seq), PER_LINE)]
    char_mm = 1.27  # Courier at 6 pt is 3.6 pt wide
    row_mm = 2.6
    tracks = 1 + max(1, len(ipr)) + max(1, len(hits))  # pLDDT, InterPro rows, hit rows
    block_mm = row_mm * (tracks + 1.6)
    height = 30 + block_mm * len(lines)
    width = pf.TWO_COLUMN_MM
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.set_axis_off()
    x0 = 30
    ax.text(4, 5, f"{acc}  {FOUR[acc]}, {len(seq)} residues", fontsize=7, weight="bold", va="center")
    legend = [Patch(facecolor=c, label=f"AlphaFold DB {lab}") for _, _, c, lab in BANDS]
    legend += [Patch(facecolor=IPR, label="InterPro match"), Patch(facecolor=HIT, label=f"MMseqs2 hit to another Swiss-Prot entry, identity > {100 * MIN_FIDENT:.0f}%")]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(4 / width, 1 - 8 / height), ncol=3, frameon=False)
    y = 26
    for start, chunk in lines:
        ax.text(x0 - 2, y, f"{start + 1}", fontsize=6, ha="right", va="center", family="monospace")
        for j, aa in enumerate(chunk):
            ax.text(x0 + j * char_mm + char_mm / 2, y, aa, fontsize=6, ha="center", va="center", family="monospace")
        end = start + len(chunk)
        ax.text(x0 + len(chunk) * char_mm + 2, y, f"{end}", fontsize=6, ha="left", va="center", family="monospace")
        ty = y + row_mm
        ax.text(x0 - 2, ty, "pLDDT", fontsize=5, ha="right", va="center")
        if plddt:
            for j in range(len(chunk)):
                v = plddt[start + j]
                col = next(c for lo, hi, c, _ in BANDS if lo <= v < hi)
                ax.add_patch(Rectangle((x0 + j * char_mm, ty - row_mm * 0.35), char_mm, row_mm * 0.7, facecolor=col, linewidth=0))
        else:
            ax.text(x0, ty, "no AlphaFold DB model", fontsize=5, va="center")
        for group, color, empty in ((ipr, IPR, "no InterPro match"), (hits, HIT, "no Swiss-Prot hit above 30% identity but itself")):
            for k in range(max(1, len(group))):
                ty += row_mm
                if not group:
                    if start == 0:
                        ax.text(x0, ty, empty, fontsize=5, va="center")
                    continue
                c = group[k]
                s, e = max(c["start"] - 1, start), min(c["end"], end)
                if s < e:
                    ax.add_patch(Rectangle((x0 + (s - start) * char_mm, ty - row_mm * 0.35), (e - s) * char_mm, row_mm * 0.7, facecolor=color, linewidth=0))
                if start == 0 or s < e:
                    ax.text(x0 + PER_LINE * char_mm + 10, ty, f"{c['label']}: {c['start']}-{c['end']}", fontsize=5, va="center")
        y += block_mm
    pf.save(fig, str(out), formats=("svg", "png", "pdf"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    ap.add_argument("--sherlock", required=True, type=pathlib.Path)
    args = ap.parse_args()
    all_calls = []
    for acc in FOUR:
        entry = json.loads((args.raw / f"uniprot_{acc}.json").read_text())
        seq = entry["sequence"]["value"]
        calls = interpro_calls(args.raw / f"interpro_matches_{acc}.json", acc)
        calls += mmseqs_calls(args.sherlock / "mmseqs_vs_sprot.tsv.gz", acc)
        conf = json.loads((args.raw / f"afdb_confidence_{acc}.json").read_text())
        plddt = conf["confidenceScore"] if conf else None
        if plddt and len(plddt) != len(seq):
            raise ValueError(f"{acc}: {len(plddt)} pLDDT values for {len(seq)} residues")
        if plddt:
            calls.append({"accession": acc, "tool": "AlphaFold DB", "label": f"mean pLDDT {sum(plddt) / len(plddt):.1f}", "start": 1, "end": len(seq), "score": sum(plddt) / len(plddt)})
        draw(acc, seq, calls, plddt, pathlib.Path("results/why/maps") / f"four_{acc}")
        all_calls += calls
    pl.DataFrame(all_calls, infer_schema_length=None).write_parquet("results/why/four_proteins_calls.parquet")
    print(pl.DataFrame(all_calls, infer_schema_length=None).select("accession", "tool", "label", "start", "end").to_pandas().to_string())


if __name__ == "__main__":
    main()
