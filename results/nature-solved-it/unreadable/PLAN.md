# Are the animals that solved it the ones the tools cannot read? Plan, fixed before any run

Written 2026-10-19, before any proteome was downloaded or searched. This file is not changed
after the first run starts. If a rule here turns out wrong, the result is still reported
under it, and the problem is written in the README next to the result.

## Question

Do species that tolerate a human insult (the tolerance atlas, Prompt 1) have a larger share
of proteins with no function label than model organisms, and than random species at the
same distance from human?

## Three groups, chosen before any run

1. **Atlas species.** Every species in column 2 of the tolerance atlas CSV that has a UniProt
   reference proteome, read from UniProt on the day the list arrives, with the proteome
   identifier (UP...) and the UniProt release.
2. **Model organisms.** *Homo sapiens*, *Mus musculus*, *Danio rerio*, *Drosophila
   melanogaster*, *Caenorhabditis elegans*, *Saccharomyces cerevisiae* (strain S288C),
   *Arabidopsis thaliana*. Proteome identifiers and release are in `model_organisms.csv`,
   read from UniProt on 2026-10-19.
3. **Random controls.** As many species as the atlas group, drawn at random by
   `scripts/nature_solved_it/draw_controls.py` with a recorded seed, from UniProt reference
   proteomes of species not in the atlas. Each atlas species takes one control from the same
   100-million-year bin of divergence time from human (bin 0 is 0 to under 100 My, bin 1 is
   100 to under 200 My, and so on). Divergence times come from TimeTree, read on the day of
   the draw. If a bin has fewer unused candidates than atlas species, the draw takes every
   candidate in that bin and the shortfall is reported; it is not filled from another bin.

## Measurement, identical for every proteome

The measurement is the one that gave 19.41% for the Pacific oyster and 0.88% for human, in
seanome/protein-function-annotation-leaderboard, branch `olgabot/why-numbers` (PR 12),
`scripts/why/labels.py` with the search in `scripts/why/sherlock_labels.sbatch`. There is no
notebook; those two scripts are the pipeline.

**Unlabelled share, the leaderboard's definition verbatim:** "share of proteins with neither
an InterPro match nor a Swiss-Prot hit above 30% identity (a hit to its own entry counts)".

- The proteins are the UniProt reference proteome's one-protein-per-gene FASTA
  (`<UP>_<taxid>.fasta.gz` from the release's `reference_proteomes` FTP folder).
- "InterPro match": UniProt cross-references the protein to any InterPro entry or
  member-database signature (`xref_interpro, xref_pfam, xref_panther, xref_cdd, xref_gene3d,
  xref_supfam, xref_smart, xref_prosite, xref_prints, xref_pirsf, xref_hamap, xref_sfld,
  xref_ncbifam`).
- "Swiss-Prot hit above 30% identity": MMseqs2 `fident > 0.30` against Swiss-Prot of the same
  UniProt release, a hit to the protein's own entry included.

The search command, copied verbatim from `scripts/why/sherlock_labels.sbatch` (MMseqs2
18.8cc5c, the image pinned by sha256 in that script):

```
    run mmseqs easy-search queries.fasta.gz "$SPROT" mmseqs_vs_sprot.tsv "$tmp" \
        -s 7.5 -e 1e-3 --max-seqs 1000 --threads "$SLURM_CPUS_PER_TASK" \
        --format-output query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen
```

The BUSCO command, copied verbatim from the same script (BUSCO 6.1.0):

```
    run busco -i "$fa.fasta" -m proteins -l "busco_downloads_${tag}/lineages/$LINEAGE" --offline \
        -c "$SLURM_CPUS_PER_TASK" -o "busco_${tag}" -f --opt-out-run-stats
```

### How this differs from the prompt's wording

Prompt 3 says "no InterPro entry, no GO term and no named function". The leaderboard measured
something narrower: no InterPro match and no Swiss-Prot hit above 30% identity. It did not
look at GO terms or at the protein name. This plan keeps the leaderboard's definition, so that
the oyster and human numbers are comparable with every new proteome. The share named
"Uncharacterized protein" in UniProt is recorded as a separate column; it is not part of the
unlabelled share.

### Also recorded per proteome

- number of proteins in the one-protein-per-gene FASTA;
- BUSCO completeness (Complete single + Complete duplicated, as % of the lineage's BUSCOs) and
  the lineage dataset name and date;
- share of proteins whose UniProt protein name is exactly "Uncharacterized protein";
- divergence time from human (TimeTree, with the URL and the date read);
- CPU-hours (Slurm TotalCPU) and peak memory (MaxRSS) of its MMseqs2 and BUSCO tasks.

## BUSCO lineage, chosen by kingdom before any run

One broad lineage per kingdom, so that completeness is counted against the same gene set
within each kingdom, and animals use the lineage that gave the oyster and human values:

| UniProt taxonomy lineage contains | BUSCO dataset (version 2026-05-22) |
|---|---|
| Metazoa | `metazoa_odb12` |
| Fungi | `fungi_odb12` |
| Viridiplantae | `viridiplantae_odb12` |
| none of these | `eukaryota_odb12` |

Every proteome also gets `eukaryota_odb12`, so that all three groups can be compared on one
gene set. The 90% cut in the rule uses the kingdom dataset in the table. Dataset names and
dates were read from https://busco-data.ezlab.org/v5/data/file_versions.tsv on 2026-10-19.

## Figure

x axis: divergence time from human (million years). y axis: unlabelled share (% of
proteins). One mark per proteome. Shape is the group (model organism, atlas species, random
control). Fill is assembly quality: filled if BUSCO completeness is at least 90%, hollow
otherwise. Legend above the plot. Oyster and human labelled with leader lines.

## Rule, copied verbatim from Prompt 3

Three sentences are possible and the figure decides which one is true:

* "The solvers are unreadable": the median unlabelled share of atlas species is at least 5 percentage points above the median of the random controls, among proteomes with BUSCO of at least 90%.
* "Non-model animals are unreadable": atlas species and random controls are within 5 points of each other, and both are at least 5 points above the model organisms. This is still a measured result about the annotation gap, and it is the honest sentence if it is what the data show; it is not the sentence about the thesis.
* "The gap is assembly quality": the difference between groups disappears when restricted to BUSCO of at least 90%, or divergence time explains it (report the correlation of unlabelled share with divergence time within each group). Report which sentence holds. Do not report the first unless the random controls were run.

## Compute budget

Two compute days on Sherlock, no GPUs. The resource table, measured from Slurm's record of
job 47961379, is in `resources.md`.
