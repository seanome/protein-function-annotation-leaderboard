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
from matplotlib.lines import Line2D
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
            "label": f"{h['target']} {100 * h['fident']:.1f}% identity, E = {h['evalue']:.1e}" + (" (itself)" if h["target_acc"] == acc else ""),
            "start": h["qstart"], "end": h["qend"], "score": h["fident"], "self_hit": h["target_acc"] == acc,
            "above_30pct": h["fident"] > MIN_FIDENT,
        })
    return rows


MAX_HITS = 5  # hits to other Swiss-Prot entries drawn per protein, strongest (bit score) first


def spans(text):
    return ", ".join(f"{a}-{b}" for a, b in text)


def draw(acc, seq, calls, plddt, out):
    """Residues 100 per line; under each line, one row per call, a bar over the residues it covers."""
    pf.use_style()
    # One row per InterPro entry or signature, with all its pieces on that row.
    ipr = {}
    for c in calls:
        if c["tool"] == "InterPro":
            ipr.setdefault(c["label"], []).append((c["start"], c["end"]))
    mm = [c for c in calls if c["tool"].startswith("MMseqs2")]
    self_hit = [c for c in mm if c["self_hit"]]
    others = [c for c in mm if not c["self_hit"] and c["above_30pct"]]  # already sorted by bit score
    rows = [("InterPro " + lab.split()[0], IPR, sp) for lab, sp in ipr.items()]
    rows += [("Swiss-Prot " + c["label"].split()[0], HIT, [(c["start"], c["end"])]) for c in others[:MAX_HITS]]

    # Header lines: what each row is, in full, with its coordinates.
    notes = [f"{acc}  {FOUR[acc]}, {len(seq)} residues"]
    notes.append("AlphaFold DB: " + (f"model, mean pLDDT {sum(plddt) / len(plddt):.1f}" if plddt else "no model"))
    notes.append("InterPro: " + ("no match" if not ipr else f"{len(ipr)} entries or signatures"))
    notes += [f"    {lab}: residues {spans(sp)}" for lab, sp in ipr.items()]
    mm_note = "MMseqs2 against Swiss-Prot 2026_03: "
    if not mm:
        mm_note += "no hit at all, not even to the protein's own Swiss-Prot entry"
    else:
        mm_note += ("hits its own entry (not drawn); " if self_hit else "") + f"{len(others)} other entries above {100 * MIN_FIDENT:.0f}% identity"
        if len(others) > MAX_HITS:
            mm_note += f"; the {MAX_HITS} with the highest bit score are drawn"
    notes.append(mm_note)
    notes += [f"    {c['label']}: residues {c['start']}-{c['end']}" for c in others[:MAX_HITS]]

    lines = [(i, seq[i : i + PER_LINE]) for i in range(0, len(seq), PER_LINE)]
    char_mm, row_mm, note_mm = 1.27, 2.4, 2.6  # Courier 6 pt is 3.6 pt = 1.27 mm wide
    head_mm = 4 + note_mm * len(notes) + 10  # notes, then the legend
    block_mm = row_mm * (2 + len(rows)) + 2
    width, height = pf.TWO_COLUMN_MM, head_mm + block_mm * len(lines) + 2
    # These are for the web page, where a long protein may run past Nature's 170 mm print
    # height; 875 residues at 100 per line need about 185 mm.
    if height > 260:
        raise ValueError(f"{acc}: figure would be {height:.0f} mm tall")
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.set_axis_off()
    for i, t in enumerate(notes):
        ax.text(3, 4 + i * note_mm, t, fontsize=7 if i == 0 else 6, weight="bold" if i == 0 else "normal", va="center")
    legend = [Patch(facecolor=c, label=f"AlphaFold DB {lab}") for _, _, c, lab in BANDS]
    legend += [Line2D([], [], color="#999999", linewidth=0.6, label="the protein (residues on the line above)"), Patch(facecolor=IPR, label="InterPro match (one bar per piece)"), Patch(facecolor=HIT, label=f"MMseqs2 hit to another Swiss-Prot entry, > {100 * MIN_FIDENT:.0f}% identity")]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(3 / width, 1 - (4 + note_mm * len(notes)) / height), ncol=4, frameon=False, borderaxespad=0)
    x0 = 34
    y = head_mm
    for start, chunk in lines:
        end = start + len(chunk)
        ax.text(x0 - 2, y, f"{start + 1}", fontsize=6, ha="right", va="center", family="monospace")
        for j, aa in enumerate(chunk):
            ax.text(x0 + j * char_mm + char_mm / 2, y, aa, fontsize=6, ha="center", va="center", family="monospace")
        ax.text(x0 + len(chunk) * char_mm + 1.5, y, f"{end}", fontsize=6, ha="left", va="center", family="monospace")
        ty = y + row_mm
        ax.text(x0 - 2, ty, "AlphaFold DB pLDDT", fontsize=5, ha="right", va="center")
        if plddt:
            for j in range(len(chunk)):
                v = plddt[start + j]
                col = next(c for lo, hi, c, _ in BANDS if lo <= v < hi)
                ax.add_patch(Rectangle((x0 + j * char_mm, ty - row_mm * 0.35), char_mm, row_mm * 0.7, facecolor=col, linewidth=0))
        for name, color, sp in rows:
            ty += row_mm
            ax.text(x0 - 2, ty, name, fontsize=5, ha="right", va="center")
            ax.plot([x0, x0 + len(chunk) * char_mm], [ty, ty], color="#999999", linewidth=0.3)  # the protein, as a line
            for a, b in sp:
                s0, e0 = max(a - 1, start), min(b, end)
                if s0 < e0:
                    ax.add_patch(Rectangle((x0 + (s0 - start) * char_mm, ty - row_mm * 0.35), (e0 - s0) * char_mm, row_mm * 0.7, facecolor=color, edgecolor="white", linewidth=0.4))
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
