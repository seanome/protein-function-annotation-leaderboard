"""Turn the tolerance atlas's species column into the atlas group, by PLAN.md's
"Atlas list decisions, 2026-10-19".

Writes, in --out-dir:
  atlas_species_list.csv  every species named in the atlas, its binomial, whether it is in
                          the atlas group, and the reason when it is not;
  atlas_species_input.csv the binomials in the atlas group (input to proteome_table.py);
  atlas_named_species.txt every binomial named anywhere in the atlas (left out of the
                          control pool).

    python3 -I atlas_species.py --atlas atlas.csv --out-dir results/nature-solved-it/unreadable
"""

import argparse
import csv
import pathlib

# PLAN.md decisions 1 and 2.
EXCLUDED = {
    "Deinococcus radiodurans": "bacterium; the question is about animals and the control pool is Eukaryota",
    "Flavobacterium frigoris": "bacterium; the question is about animals and the control pool is Eukaryota",
    "Danio rerio": "model organism; kept in the model-organism group only",
    "Drosophila melanogaster": "model organism; kept in the model-organism group only",
}
NO_NUMBER = ("no number", "not found", "none found")


def binomial(text):
    """Decision 4: drop text in brackets or after a semicolon; keep genus and species."""
    words = text.split("(")[0].split(";")[0].split()
    return " ".join(words[:2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, type=pathlib.Path)
    ap.add_argument("--species-column", default="species")
    ap.add_argument("--number-column", default="measured_ability")
    ap.add_argument("--out-dir", required=True, type=pathlib.Path)
    args = ap.parse_args()
    with open(args.atlas, newline="") as fh:
        rows = list(csv.DictReader(fh))

    listed = {}
    for r in rows:
        written = r[args.species_column].strip()
        if not written or written.lower().startswith("none found"):
            continue  # decision 5
        b = binomial(written)
        has_number = not r[args.number_column].strip().lower().startswith(NO_NUMBER)
        e = listed.setdefault(b, {"binomial": b, "as_written": written, "rows": 0, "rows_with_number": 0})
        e["rows"] += 1
        e["rows_with_number"] += has_number

    for b, e in listed.items():
        genus_only = len(b.split()) < 2 or b.split()[1].lower() in ("spp.", "sp.", "spp")
        if genus_only:
            e.update(in_atlas_group="no", reason="genus only, not a species (decision 5)")
        elif b in EXCLUDED:
            e.update(in_atlas_group="no", reason=EXCLUDED[b])
        elif not e["rows_with_number"]:
            e.update(in_atlas_group="no", reason="no row with a measured number (decision 3)")
        else:
            e.update(in_atlas_group="yes", reason="")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    cols = ["binomial", "as_written", "rows", "rows_with_number", "in_atlas_group", "reason"]
    entries = sorted(listed.values(), key=lambda e: e["binomial"])
    with open(out / "atlas_species_list.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, lineterminator="\n")
        w.writeheader()
        w.writerows(entries)
    group = [e for e in entries if e["in_atlas_group"] == "yes"]
    with open(out / "atlas_species_input.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, ["species", "as_written"], lineterminator="\n")
        w.writeheader()
        w.writerows({"species": e["binomial"], "as_written": e["as_written"]} for e in group)
    named = [e["binomial"] for e in entries]
    (out / "atlas_named_species.txt").write_text("\n".join(named) + "\n")
    with_number = sum(1 for e in entries if e["rows_with_number"] and len(e["binomial"].split()) == 2 and e["binomial"].split()[1] != "spp.")
    print(f"{len(entries)} names in the atlas; {with_number} binomials with a measured number; {len(group)} in the atlas group")
    for e in entries:
        if e["in_atlas_group"] == "no":
            print(f"  left out: {e['binomial']}: {e['reason']}")


if __name__ == "__main__":
    main()
