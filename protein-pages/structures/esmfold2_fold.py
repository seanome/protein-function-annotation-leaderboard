"""ESMFold2 structures for proteins of a FASTA file, one protein at a time.

    python esmfold2_fold.py [--kernels fused] <proteins.fasta> <out folder> [protein ID ...]

With no IDs, every protein of the FASTA is folded. For each protein it writes
<out folder>/<ID>.cif (the model's mmCIF) and <ID>.pdb (the same atoms; pLDDT, 0 to 100, in the
B-factor column, which build.py reads from structures/<ID>.pdb), and adds a line to
<out folder>/esmfold2_summary.tsv: ID, length, mean pLDDT (0 to 100), pTM, seconds, the
median distance between neighbouring CA atoms (about 3.8 angstrom in any real chain, so a value
far from it means the coordinates are broken), and the most GPU memory the fold used (GB). A
protein that does not fit in GPU memory gets a line with "out of memory", and a rerun tries it again. A protein whose .pdb exists is skipped, so a rerun
folds only what is missing. Files are written under a .part name and then renamed.

Settings are the package's own defaults for a fold (esm 3.4.2, ESMFold2InputBuilder.fold):
20 loops, 200 diffusion steps, one sample, language-model dropout 0.3, the ESMC backbone in fp8
(fp8 instructions need compute capability 8.9 or newer, so H100 or L40S; other GPUs run it slower
in another format); seed 0, so a rerun gives the same structure.

--kernels fused switches on the package's Triton kernels (EsmFold2Model.set_kernel_backend), which
its documentation says give the same numbers as the default code. Without it the default code runs.
"""
import argparse
import math
import statistics
import sys
import time
from pathlib import Path

import torch
from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model, ProteinInput, StructurePredictionInput


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs.setdefault(name, [])  # an ID repeated for isoforms keeps its first sequence only
            name = name if not seqs[name] else None
        elif name:
            seqs[name].append(line)
    return {k: "".join(v).rstrip("*") for k, v in seqs.items()}


def ca_ca_median(pdb):
    ca = [tuple(float(l[i:i + 8]) for i in (30, 38, 46)) for l in pdb.splitlines()
          if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    d = [math.dist(a, b) for a, b in zip(ca, ca[1:])]
    return statistics.median(d) if d else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernels", choices=["fused"])
    ap.add_argument("fasta", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("ids", nargs="*")
    a = ap.parse_args()
    fasta, out, ids = a.fasta, a.out, a.ids
    seqs = read_fasta(fasta)
    missing = [i for i in ids if i not in seqs]
    if missing:
        sys.exit("not in %s: %s" % (fasta, ", ".join(missing)))
    todo = [i for i in (ids or list(seqs)) if not (out / (i + ".pdb")).exists()]
    out.mkdir(parents=True, exist_ok=True)
    print("%d to fold, %d already done" % (len(todo), len(ids or seqs) - len(todo)), flush=True)
    if not todo:
        return
    t0 = time.time()
    model = EsmFold2Model.from_pretrained("biohub/ESMFold2", device="cuda").eval()
    if a.kernels:
        model.set_kernel_backend(a.kernels)
    print("model loaded in %.0f s; kernels: %s" % (time.time() - t0, a.kernels or "default"), flush=True)
    builder = ESMFold2InputBuilder()
    summary = out / "esmfold2_summary.tsv"
    if not summary.exists():
        summary.write_text("id\tlength\tmean_plddt\tptm\tseconds\tca_ca_median\tgpu_gb\n")
    for pid in todo:
        seq, t0 = seqs[pid], time.time()
        spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=seq)])
        torch.cuda.reset_peak_memory_stats()
        try:
            with torch.inference_mode():
                res = builder.fold(model, spi, num_loops=20, num_sampling_steps=200, num_diffusion_samples=1, seed=0)
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            line = "%s\t%d\t\t\t%.1f\tout of memory\t%.1f\n" % (pid, len(seq), time.time() - t0, torch.cuda.max_memory_allocated() / 1e9)
            with open(summary, "a") as fh:
                fh.write(line)
            print(line, end="", flush=True)
            continue
        res = res[0] if isinstance(res, list) else res
        cif = res.complex.to_mmcif()
        pdb = res.complex.to_protein_complex().to_pdb_string()
        for ext, text in ((".cif", cif), (".pdb", pdb)):  # the .pdb last: it marks the protein done
            part = out / (pid + ext + ".part")
            part.write_text(text)
            part.rename(out / (pid + ext))
        secs = time.time() - t0
        line = "%s\t%d\t%.1f\t%s\t%.1f\t%.2f\t%.1f\n" % (pid, len(seq), 100 * float(res.plddt.float().mean()),
                                                           "%.3f" % res.ptm if res.ptm is not None else "", secs, ca_ca_median(pdb),
                                                           torch.cuda.max_memory_allocated() / 1e9)
        with open(summary, "a") as fh:
            fh.write(line)
        print(line, end="", flush=True)


if __name__ == "__main__":
    main()
