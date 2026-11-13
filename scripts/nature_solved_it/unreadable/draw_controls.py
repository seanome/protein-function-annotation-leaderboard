"""Draw the random non-model controls, matched to the atlas species on divergence time from
human in 100-million-year bins, and write the run list proteomes.tsv.

Rules (PLAN.md, group 3 and "Details fixed before the draw"):
  - candidates: Eukaryota rows of the UniProt reference proteomes README, minus every atlas
    species and the seven model organisms, matched by binomial (first two words of the
    README's species name);
  - shuffled once with random.Random(seed).shuffle, then read in that order;
  - a candidate joins if TimeTree gives it a time (adjusted_age, from human) and its bin
    (floor(time / 100)) still needs a species; a candidate TimeTree cannot place is skipped
    and counted; a second proteome of an already drawn binomial is skipped;
  - each atlas species with a reference proteome and a TimeTree time asks for one control
    in its bin; a bin short of candidates is reported, never filled from another bin.

Also writes the "leaderboard reference" row for the Pacific oyster (UP000005408), unless the
oyster is already in a group, so that the run reproduces the leaderboard's 19.41% and the
figure can label it. That row is in no group median.

    python3 -I scripts/nature_solved_it/unreadable/draw_controls.py \
        --atlas results/nature-solved-it/unreadable/atlas_proteomes.csv \
        --models results/nature-solved-it/unreadable/model_organisms.csv \
        --out-dir results/nature-solved-it/unreadable

Outputs in --out-dir: controls.csv, draw_record.json, proteomes.tsv. TimeTree answers are
cached in timetree_cache.json there; a redo with the same cache and seed gives the same draw.
"""

import argparse
import csv
import gzip
import json
import pathlib
import random
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import proteome_table  # noqa: E402
import today  # noqa: E402

SEED = 20261019
BIN_MY = 100
OYSTER = ("Magallana gigas", "UP000005408")


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def readme_rows(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("Proteome_ID\t"))
    return [r for r in csv.DictReader(lines[start:], delimiter="\t") if r.get("Tax_ID")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, type=pathlib.Path)
    ap.add_argument("--models", required=True, type=pathlib.Path)
    ap.add_argument("--out-dir", required=True, type=pathlib.Path)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--readme", type=pathlib.Path, help="saved README (.gz); default <out-dir>/reference_proteomes_README.txt.gz")
    ap.add_argument("--exclude-species", type=pathlib.Path, help="file of binomials, one a line, also left out of the pool (every species named in the atlas)")
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    atlas = read_csv(args.atlas)
    models = read_csv(args.models)
    readme_path = args.readme or out / "reference_proteomes_README.txt.gz"
    candidates_all = readme_rows(readme_path)
    readme = {r["Proteome_ID"]: r for r in candidates_all}
    tt = today.TimeTree(out / "timetree_cache.json")

    # What each bin needs: one control per atlas species that has a proteome and a time.
    need, atlas_unbinned = {}, []
    for r in atlas:
        if not r["upid"]:
            continue
        if r["divergence_my"] == "":
            atlas_unbinned.append(r["species"])
            continue
        b = int(float(r["divergence_my"]) // BIN_MY)
        need[b] = need.get(b, 0) + 1
    need_start = dict(need)

    excluded = {r["species"].lower() for r in atlas} | {r["species"].lower() for r in models}
    if args.exclude_species:
        excluded |= {line.strip().lower() for line in args.exclude_species.read_text().splitlines() if line.strip()}
    pool = [
        r for r in candidates_all
        if r["SUPERREGNUM"] == "eukaryota" and today.binomial_of(r["Species Name"]).lower() not in excluded
    ]
    pool.sort(key=lambda r: r["Proteome_ID"])  # fixed order before the shuffle
    rng = random.Random(args.seed)
    rng.shuffle(pool)

    drawn, drawn_binomials, read, no_time, same_species, bin_full = [], set(), 0, 0, 0, 0
    for i, r in enumerate(pool):
        if not any(need.values()):
            break
        if i % 40 == 0:  # read TimeTree 40 candidates ahead, 4 at a time
            tt.prefetch(int(x["Tax_ID"]) for x in pool[i:i + 40])
        read += 1
        name = today.binomial_of(r["Species Name"])
        if name.lower() in drawn_binomials:
            same_species += 1
            continue
        t = tt.lookup(int(r["Tax_ID"]))
        if t["adjusted_age"] is None:
            no_time += 1
            continue
        b = int(t["adjusted_age"] // BIN_MY)
        if need.get(b, 0) <= 0:
            bin_full += 1
            continue
        need[b] -= 1
        drawn.append((name, r["Proteome_ID"], b))
        drawn_binomials.add(name.lower())

    # Full rows for the drawn controls (UniProt lineage for the BUSCO dataset, FASTA URL).
    controls = []
    for name, upid, b in drawn:
        prot, release, url = today.reference_proteomes(f"upid:{upid}")
        controls.append(proteome_table.row_for("random control", name, prot, release, url, readme, tt))
    with open(out / "controls.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=proteome_table.COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(controls)

    record = {
        "seed": args.seed,
        "bin_width_my": BIN_MY,
        "readme": str(readme_path.name),
        "readme_release_line": next((l for l in gzip.open(readme_path, "rt").read().splitlines() if l.startswith("Release ")), ""),
        "excluded_binomials": sorted(excluded),
        "candidates_eukaryota_after_exclusion": len(pool),
        "candidates_read": read,
        "skipped_no_timetree_time": no_time,
        "skipped_same_binomial_already_drawn": same_species,
        "skipped_bin_already_full": bin_full,
        "atlas_species_with_proteome": sum(1 for r in atlas if r["upid"]),
        "atlas_species_without_timetree_time": atlas_unbinned,
        "needed_per_bin": {str(k): v for k, v in sorted(need_start.items())},
        "short_per_bin": {str(k): v for k, v in sorted(need.items()) if v > 0},
        "drawn": [{"species": n, "upid": u, "bin": b} for n, u, b in drawn],
    }
    (out / "draw_record.json").write_text(json.dumps(record, indent=1) + "\n")

    # The run list: every row with a proteome, in group order.
    run = [r for r in models if r["upid"]] + [r for r in atlas if r["upid"]] + [r for r in controls if r["upid"]]
    if OYSTER[1] not in {r["upid"] for r in run}:
        prot, release, url = today.reference_proteomes(f"upid:{OYSTER[1]}")
        run.append(proteome_table.row_for("leaderboard reference", OYSTER[0], prot, release, url, readme, tt))
    dup = {r["upid"] for r in run if [x["upid"] for x in run].count(r["upid"]) > 1}
    if dup:
        raise SystemExit(f"a proteome is in two groups: {sorted(dup)}")
    with open(out / "proteomes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=proteome_table.COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(run)
    print(json.dumps({k: record[k] for k in ("seed", "needed_per_bin", "short_per_bin", "candidates_read", "skipped_no_timetree_time")}))
    print(f"{out / 'proteomes.tsv'}: {len(run)} proteomes")


if __name__ == "__main__":
    main()
