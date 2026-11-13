"""One table of the ESMFold2 proteome run, and what is left to fold.

    python3 esmfold2_proteome_collect.py esmfold2_proteome_chunks.tsv $SCRATCH/protein-pages-esmfold2/proteome

Reads chunk_NNN/esmfold2_summary.tsv under the run folder and writes
esmfold2_proteome_summary.tsv there: the summary columns plus chunk and status (folded, out of
memory, or not run). A protein rerun keeps its last line. Prints the counts, the chunks to rerun
as an --array value, and the proteins whose median distance between neighbouring CA atoms is
far from the 3.8 angstrom of a real chain. Standard library only, for the login node's Python 3.6.
"""
import csv
import os
import sys

COLUMNS = ["id", "length", "mean_plddt", "ptm", "seconds", "ca_ca_median", "gpu_gb"]
CA_LOW, CA_HIGH = 3.6, 4.0


def main(chunks_tsv, run):
    plan = list(csv.DictReader(open(chunks_tsv), delimiter="\t"))
    seen = {}
    for c in sorted({int(r["chunk"]) for r in plan}):
        path = os.path.join(run, "chunk_%03d" % c, "esmfold2_summary.tsv")
        if os.path.exists(path):
            for r in csv.DictReader(open(path), delimiter="\t"):
                seen[r["id"]] = r
    counts = {"folded": 0, "out of memory": 0, "not run": 0}
    rerun, broken, seconds = set(), [], 0.0
    with open(os.path.join(run, "esmfold2_proteome_summary.tsv"), "w") as fh:
        fh.write("\t".join(COLUMNS + ["chunk", "status"]) + "\n")
        for p in plan:
            r = seen.get(p["protein_id"])
            if r is None:
                status = "not run"
                r = {"id": p["protein_id"], "length": p["length"]}
                rerun.add(int(p["chunk"]))
            elif r["ca_ca_median"] == "out of memory":
                status = "out of memory"
                seconds += float(r["seconds"])
            else:
                status = "folded"
                seconds += float(r["seconds"])
                if not CA_LOW <= float(r["ca_ca_median"]) <= CA_HIGH:
                    broken.append(r["id"])
            counts[status] += 1
            fh.write("\t".join([r.get(k, "") for k in COLUMNS] + [p["chunk"], status]) + "\n")
    print("%d proteins: %d folded, %d out of memory, %d not run" % (len(plan), counts["folded"], counts["out of memory"], counts["not run"]))
    print("%.0f GPU hours of folding so far, %.0f predicted for the whole run" % (seconds / 3600, sum(float(p["predicted_seconds"]) for p in plan) / 3600))
    print("CA spacing outside %.1f to %.1f angstrom: %d %s" % (CA_LOW, CA_HIGH, len(broken), " ".join(broken[:20])))
    print("chunks to rerun: %s" % (",".join(str(c) for c in sorted(rerun)) or "none"))


if __name__ == "__main__":
    main(*sys.argv[1:3])
