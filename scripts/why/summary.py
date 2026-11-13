"""Write results/why/SUMMARY.md: every number of the "Why this project" page with its unit,
source file and release, read from results/why/numbers.parquet. Sentences take their
numbers from the table too ({id} below), so no number in the text is typed by hand.

    python3 -I scripts/why/summary.py --raw ~/data/why-numbers/2026_03
"""

import argparse
import gzip
import json
import pathlib
import re

import polars as pl

import sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from four import FILE  # noqa: E402

OUT = pathlib.Path("results/why")
STEPS = {
    1: "Swiss-Prot by organism",
    2: "Proteins with an experimental structure (PDB)",
    3: "Described species, and species with a Swiss-Prot entry",
    4: "Human and Pacific oyster, one row per protein",
    7: "What labelling the two proteomes cost",
}
FOUR = {
    "Q25460": "Mussel foot protein 1 (Mytilus edulis)",
    "P24856": "Antifreeze glycoprotein polyprotein (Notothenia neglecta)",
    "P05484": "Omega-conotoxin MVIIA (Conus magus)",
    "P42212": "Green fluorescent protein (Aequorea victoria)",
}


def fmt(v, unit):
    if unit.startswith("%") or unit == "pLDDT" or unit == "% points":
        return f"{v:,.2f}"
    if abs(v) >= 1e7:
        return f"{v:.3g}"
    return f"{v:,.0f}"


def sprot_names(raw, accs):
    """Protein name and organism from the FASTA headers of the same Swiss-Prot release."""
    out = {}
    with gzip.open(raw / "uniprot_sprot.fasta.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line.split("|")[1]
                if acc in accs:
                    m = re.match(r">\S+ (.*?) OS=(.*?) OX=", line)
                    out[acc] = f"{m.group(1)}, {m.group(2)}"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, type=pathlib.Path)
    args = ap.parse_args()
    nums = pl.read_parquet(OUT / "numbers.parquet")
    verdicts = OUT / "verification" / "verdicts.parquet"  # from the second agent's recompute.py
    vmap = dict(pl.read_parquet(verdicts).select("id", "verdict").iter_rows()) if verdicts.exists() else {}
    sources = json.loads((OUT / "sources.json").read_text())
    V = {r["id"]: fmt(r["value"], r["unit"]) for r in nums.iter_rows(named=True)}

    def t(s):
        return s.format(**V)

    md = []
    md.append("# Why this project: the numbers\n")
    md.append(t(
        "Every number below was computed from a file downloaded on 2026-10-15 or made from one, and sits next to "
        "that file and its release. The table holding them is `results/why/numbers.parquet`; the files, their URLs "
        "and sha256 sums are in `results/why/sources.json`. UniProt is release 2026_03 (2 September 2026) throughout. "
        "The last column is the verdict of a second, independent recomputation (HOLDS, WRONG or UNVERIFIABLE).\n"
    ))
    md.append("## What the numbers say\n")
    md.append(t(
        "- Swiss-Prot 2026_03 holds {sp_entries} reviewed entries from {sp_species} species. The five largest organisms "
        "hold {sp_top5_share_by_organisms}% of them.\n"
        "- Animal species that the World Register of Marine Species (WoRMS) lists as marine supply {sp_marine_animal_entries} "
        "entries, {sp_marine_animal_share}% of Swiss-Prot.\n"
        "- {kb_pdb_entries} UniProt entries from cellular life cite an experimental structure in the PDB. Counted by organism "
        "name, as the deck does, they come from {kb_pdb_organism_name} organisms, and the five largest hold "
        "{kb_pdb_top5_share_by_organism_name}%. The deck said 62,253 entries, 5,918 organisms and 26.4%.\n"
        "- Of {col_species_Animalia} described animal species (Catalogue of Life), {sp_species_Animalia} have a Swiss-Prot entry "
        "({sp_species_share_Animalia}%).\n"
        "- In the Pacific oyster, {oyster_no_interpro_no_sprot_pct}% of proteins have neither an InterPro match nor a "
        "Swiss-Prot hit above 30% identity; in human, {human_no_interpro_no_sprot_pct}%. Both count a hit to the protein's own "
        "Swiss-Prot entry; without it, human rises to {human_no_interpro_no_sprot_nonself_pct}% and oyster is "
        "{oyster_no_interpro_no_sprot_nonself_pct}%, since few oyster proteins are Swiss-Prot entries. An AlphaFold DB model for the protein's accession exists for "
        "{oyster_af_model_pct}% of oyster and {human_af_model_pct}% of human proteins, so with the model counted as a label "
        "only {oyster_no_label_pct}% and {human_no_label_pct}% have none.\n"
    ))

    for step, title in STEPS.items():
        md.append(f"## {step}. {title}\n")
        md.append("| Number | Value | Unit | Source file | Release | Check |\n|---|---|---|---|---|---|")
        for r in nums.filter(pl.col("step") == step).sort("id").iter_rows(named=True):
            desc = r["description"] + (f" ({r['note']})" if r["note"] else "")
            md.append(
                f"| {desc} | {fmt(r['value'], r['unit'])} | {r['unit']} | [`{r['source_file']}`]({r['source_url']}) "
                f"| {r['release']} | {vmap.get(r['id'], 'not yet checked')} |"
            )
        md.append("")
        if step == 4:
            md.append("![Proteome maps, function labels](maps/proteome_map_human_oyster_interpro_or_swissprot_uniprot_2026_03.png)\n")
            md.append("![Proteome maps, all three labels](maps/proteome_map_human_oyster_any_label_uniprot_2026_03.png)\n")
            md.append(
                "Each tile is one protein; both maps use 150 tiles per row and the same tile size. A tile is lit when the "
                "protein has the labels named in the legend. The counts printed above each map come from "
                "`results/why/proteome_labels.parquet`, the same table as the numbers above.\n"
            )
            if "cost_mmseqs2_reserved_core_h" in V:
                md.append(t(
                    "Labelling both proteomes this way took {cost_mmseqs2_reserved_core_h} reserved core-hours and "
                    "{cost_mmseqs2_peak_rss_gb} GB of memory on {cost_mmseqs2_cpus} CPUs for MMseqs2 (one search for both), "
                    "and {cost_busco_oyster_reserved_core_h} reserved core-hours for BUSCO on the oyster, Sherlock job 47961379. "
                    "Step 7 has both species and the CPU-hours used.\n"
                ))

    # Step 6
    calls = pl.read_parquet(OUT / "four_proteins_calls.parquet")
    mm = calls.filter(pl.col("tool").str.starts_with("MMseqs2"), ~pl.col("self_hit"), pl.col("above_30pct"))
    names = sprot_names(args.raw, set(mm["label"].str.split(" ").list.get(0).to_list()))
    md.append("## 6. Four proteins, the same three labels\n")
    md.append(
        "Each figure shows the protein's residues, 100 per line, with a bar under the residues each call covers. "
        "All four are Swiss-Prot entries, so the MMseqs2 search finds the protein's own entry; that hit is not drawn "
        "and not counted below.\n"
    )
    for acc, name in FOUR.items():
        c = calls.filter(pl.col("accession") == acc)
        af = c.filter(pl.col("tool") == "AlphaFold DB")
        ipr = c.filter(pl.col("tool") == "InterPro")
        hits = mm.filter(pl.col("accession") == acc)
        selfh = c.filter(pl.col("tool").str.starts_with("MMseqs2"), pl.col("self_hit"))
        seq_len = json.loads((args.raw / f"uniprot_{acc}.json").read_text())["sequence"]["length"]
        md.append(f"### {acc} {name}, {seq_len} residues\n")
        md.append(f"![{acc}](maps/{FILE[acc]}.png)\n")
        lines = []
        lines.append("AlphaFold DB: " + (af["label"][0] if af.height else "no model") + ".")
        if ipr.height:
            covered = set()
            for s, e in ipr.select("start", "end").iter_rows():
                covered.update(range(s, e + 1))
            sigs = ipr["label"].unique(maintain_order=True).to_list()
            lines.append(f"InterPro labels it: {len(sigs)} entries or signatures, covering {len(covered)} of {seq_len} residues ("
                         + "; ".join(f"{lab} at {', '.join(f'{s}-{e}' for s, e in ipr.filter(pl.col('label') == lab).select('start', 'end').iter_rows())}" for lab in sigs) + ").")
        else:
            lines.append("InterPro: no match anywhere on the protein.")
        if hits.height:
            top = hits.head(3)
            lines.append(f"MMseqs2 finds {hits.height} other Swiss-Prot entries above 30% identity; the strongest: "
                         + "; ".join(f"{r['label'].split(' ')[0]} ({names.get(r['label'].split(' ')[0], '?')}), "
                                     f"{100 * r['score']:.1f}% identity over residues {r['start']}-{r['end']}" for r in top.iter_rows(named=True)) + ".")
        elif selfh.height:
            lines.append("MMseqs2 finds no other Swiss-Prot entry above 30% identity, only the protein's own entry.")
        else:
            lines.append("MMseqs2 finds no hit at all, not even the protein's own Swiss-Prot entry.")
        md.append("\n".join(f"- {x}" for x in lines) + "\n")

    md.append("## Notes on the sources\n")
    md.append(
        "- The deck's 62,253 / 5,918 / 26.4% recompute from the table it names (`uniprot_combined_sprot_trembl_taxonomy_summary.parquet` "
        "in seanome/2025-biodata-wealth-inequality, UniProt as parsed on 2025-03-11) as " + V["deck_entries_recomputed"] + " entries, "
        + V["deck_organisms_recomputed"] + " organism names and " + V["deck_top5_pct_recomputed"] + "%; the slide rounds 26.45 down. "
        "Its 'species' are organism names, which count each named strain separately; counted as NCBI species there are " + V["kb_pdb_species_taxid"] + ".\n"
        "- The InterPro API timed out (HTTP 408) on every attempt to list the human proteome. The InterPro column of step 4 therefore "
        "uses UniProt's cross-references to InterPro and its member databases (Pfam, PANTHER, CDD, Gene3D, SUPERFAMILY, SMART, "
        "PROSITE, PRINTS, PIRSF, HAMAP, SFLD, NCBIfam), which UniProt takes from InterPro's own runs. The four proteins' match "
        "coordinates come from the InterPro API (release " + sources["interpro_matches_P42212"]["release"].replace("InterPro ", "") + ").\n"
        "- Catalogue of Life lists " + V["col_species_Protozoa"] + " described Protozoa species, against " + V["col_species_Chromista"] + " Chromista. "
        "Check COL's coverage of protists before quoting the Protozoa ratio in step 3.\n"
        "- Larsen et al. 2017 (https://doi.org/10.1086/693564), the source the deck uses, gives ESTIMATED total species, not "
        "described species. Those rows are marked as estimates and are read from the notebook that tables them.\n"
        "- WoRMS matched " + V["worms_matched_species"] + " of " + V["sp_metazoa_species"] + " animal species names. WoRMS lists "
        "marine species and many others, but not all land and freshwater ones, so an unmatched name was counted as not marine.\n"
        "- Larsen et al. 2017 estimate " + V["larsen2017_estimated_species_Plants"] + " plant species in all, fewer than the "
        + V["col_species_Plantae"] + " Plantae species the Catalogue of Life already lists as described.\n"
        "- The Archaea (" + V["sp_species_share_Archaea"] + "%) and Protozoa (" + V["sp_species_share_Protozoa"] + "%) ratios in step 3 divide "
        "NCBI species counts by Catalogue of Life counts; the two sources draw these groups differently, so the two ratios are not comparable to the others.\n"
        "- 'Above 30% identity' is MMseqs2's identity strictly greater than 0.30. Counting 0.30 itself adds " + V["oyster_sprot_hit_at_exactly_30pct"] + " oyster and " + V["human_sprot_hit_at_exactly_30pct"] + " human proteins.\n"
        "- The antifreeze glycoprotein P24856 is a Swiss-Prot entry, yet MMseqs2 returned no hit for it, not even to itself. "
        "Its sequence is mostly the repeat AAT; low-complexity masking in MMseqs2 is the likely cause, not tested here.\n"
        "- The InterPro API listing of the oyster proteome, meant as a check on UniProt's cross-references, timed out (HTTP 408) on "
        "three tries over two days, so that check was not done.\n"
        "- In `proteome_labels.parquet`, `sprot_best_target` is the hit with the highest bit score, not the highest identity.\n"
        "- MMseqs2 18.8cc5c (-s 7.5, E ≤ 0.001) and BUSCO 6.1.0 (metazoa_odb12) ran on Sherlock as job 47961379, from "
        "`scripts/why/sherlock_labels.sbatch`.\n"
    )
    (OUT / "SUMMARY.md").write_text("\n".join(md))
    print("wrote", OUT / "SUMMARY.md", len(md), "blocks;", len(vmap), "verdicts")


if __name__ == "__main__":
    main()
