"""The Biomni whole-proteome results that have no position on the protein, for the drawn proteins.

    python test_inputs/make_biomni_best.py full_inputs test_inputs/data/proteins.fasta > test_inputs/data/biomni_best.json

build.py draws these as rows that name the best hit and say "whole protein, no position". The
drawn proteins are matched to the proteome by sequence MD5, as build_full.py does. One entry per
drawn protein ID, each key 0 when the program was run and found nothing:
    diamond   [accession, protein name, gene, organism, identity %, coverage %, E-value]
    jackhmmer [seed family, seed protein, round, E-value], the best of the hits from immune seeds
    fantasia  [accession, protein name, gene, organism, distance], or "not run" (FANTASIA was run only
              on proteins with no Swiss-Prot match or in the immune table)
    unifire   [what it predicts (a name, family or keywords), number of predictions], or "not run"
No sequence is written.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import build_full  # noqa: E402

folder, fasta = Path(sys.argv[1]), Path(sys.argv[2])
ann, lengths, immune, jh, rep_of, members, uni, go, uf, km, ips, mm, fs = build_full.load(folder)
by_md5 = defaultdict(list)
for pid, (_, md5) in lengths.items():
    by_md5[md5].append(pid)
row = {r["protein_id"]: r for r in ann}


def unifire(pid):
    if uf is None:
        return None
    u = uf.get(pid, [])
    if not u:
        return "not run" if lengths[pid][0] > build_full.UNIFIRE_MAX_AA else 0
    pick = lambda t: [x[3] for x in u if x[2] == t]
    what = (pick("protein.recommendedName.fullName") or pick("comment.similarity") or
            [", ".join(list(dict.fromkeys(pick("keyword")))[:5])])[0]
    return [what, len(u)]


out = {}
for did, md5 in sorted(build_full.read_fasta_md5(fasta).items()):
    pids = sorted(by_md5.get(md5, []))
    if not pids:
        continue
    pid = pids[0]  # proteins with the same sequence have the same results
    r = row[pid]
    u = uni.get(r["sp_acc"], {})
    j = sorted(jh.get(pid, []), key=lambda x: float(x[3]))
    fu = uni.get(r.get("fantasia_uniprot", ""), {})
    out[did] = dict(
        proteome_id=pid,
        diamond=[r["sp_acc"], r["sp_name"], u.get("gene", ""), u.get("organism", ""), build_full.num(r["sp_pident"]),
                 build_full.num(r["sp_qcov"]), r["sp_evalue"]] if r["sp_acc"] else 0,
        jackhmmer=j[0][:4] if j else 0,
        fantasia=[r["fantasia_uniprot"], fu.get("protein_name", ""), fu.get("gene", ""), fu.get("organism", ""),
                  round(float(r["fantasia_distance"]), 3)] if r.get("fantasia_uniprot") else "not run",
        unifire=unifire(pid))
json.dump(out, sys.stdout, indent=0, sort_keys=True)
sys.stdout.write("\n")
