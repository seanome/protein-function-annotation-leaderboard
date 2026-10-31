# Inputs of build_full.py: annotation tables for all 45,339 Botryllus proteins

These tables were made by a Biomni agent on its machine `worker-0`, not in this repository, and
downloaded by Olga on 2026-10-28. `biomni_PLAN.md` is the plan the agent followed and
`biomni_worker-0.ipynb` is the notebook from that machine, both as downloaded. The notebook
holds only tool lookups and an ESMCFold test run; the commands that made the tables below are
not in it, so the settings given here come from the plan and are not checked against a log.

The proteome is `botryllus_new_genome_protein.fasta` (45,339 proteins; on the laptop at
`~/data/botryllus/`). It is not committed: the genome is unpublished.

| File | What it holds | How it was made, per the plan |
|---|---|---|
| `botryllus_protein_annotations.tsv.gz` | One row per protein: the name in the FASTA header, a proposed name, how the two compare (`name_status`), the best Swiss-Prot hit (accession, name, % identity, % coverage of your protein, E-value), that entry's GO terms and UniProt keywords, the Pfam domains, the best Fester or BHF hit, and `mmseqs_cluster` | DIAMOND blastp `--sensitive` against Swiss-Prot, top 25 hits, E-value at most 0.001; Pfam-A by pyhmmer hmmsearch with the gathering thresholds (`--cut_ga`); the proposed name from the best Swiss-Prot hit or the Pfam domains |
| `botryllus_immune_genes.tsv` | 2,835 immune-gene candidates, their categories and the evidence for each | A category is set by a Pfam domain, the gene name, an immune keyword or GO term of the best Swiss-Prot hit, or a Fester or BHF hit |
| `botryllus_jackhmmer_immune.tsv` | 7,881 jackhmmer hits from 20 immune seed proteins (human C1QA, NLRC5 and others), on 3,589 proteins | Not in the plan; the number of rounds and the E-value cut are not recorded. The rounds in the table are 2, 3 and 5; the weakest E-value is 9.1e-06 |
| `botryllus_allorecognition_hits.tsv` | 832 DIAMOND hits of 53 proteins to 96 reference proteins: 39 Fester family, 56 Fester co-receptor and 1 BHF (partial). No header; the columns are query, reference, % identity, aligned columns, E-value, bit score, query length, reference length | DIAMOND against a small database of Fester and BHF proteins from NCBI |
| `botryllus_mmseqs_clusters.tsv` | 27,483 MMseqs2 clusters: representative, member (45,339 rows) | Not in the plan; the clustering settings are not recorded |
| `allorecognition_references.tsv` | The NCBI title of each of the 96 reference proteins | NCBI E-utilities `esummary`, db=protein, fetched 2026-10-28 |
| `proteome_lengths.tsv.gz` | Each protein's ID as the tables name it, its length and the MD5 of its sequence | `make_proteome_lengths.py` on the proteome FASTA |
| `botryllus_protein_annotations_v2.tsv.gz`, `botryllus_immune_genes_v2.tsv.gz` | The same tables with four FANTASIA columns added: the nearest experimentally annotated UniProt protein by ProtT5 embedding, its cosine distance, a reliability index and the GO terms copied from it. FANTASIA was run on the 26,676 proteins with no Swiss-Prot match or in the immune table. Every other column is the same as in the first version | Biomni, downloaded 2026-10-30; `biomni_report_v2.md` is its report |
| `uniprot_names.tsv.gz` | Protein name, gene and organism of every Swiss-Prot match, FANTASIA reference protein and best MMseqs2 hit (21,473 of 21,610; the other 137 are TrEMBL entries no longer in UniProt) | `../gene_cards/fetch_names.py`, UniProt REST, 2026-10-30 |
| `go_names.tsv.gz` | Name and namespace of every GO term in the tables, UniFIRE's included (17,241) | `../gene_cards/fetch_names.py`, go-basic.obo, 2026-10-31 |
| `interproscan.tsv.gz` | InterProScan 5.59 matches, all member databases, under the Biomni protein IDs: 391,212 rows on 33,389 proteins; the 50 proteins over 5,000 aa were not run | `make_interproscan_table.py` from the XML of olgabot/botryllus-mhc PR 13's run, via `../test_inputs/interproscan_xml_to_tsv.py` |
| `mmseqs2_best_per_protein.tsv.gz` | Best MMseqs2 hit against Swiss-Prot per protein (lowest E-value, then highest bit score): 22,043 proteins, all at E-value 0.001 or below | `botryllus_mmseqs2_sp.m8` from Olga, 2026-10-31 (3,566,205 rows; settings not recorded with it) |
| `kmerseek_best_per_protein.tsv.gz` | Best kmerseek hit against all of Swiss-Prot 2026_03 per protein, over the 35 settings whose E-value fit worked: 45,317 proteins have one | `make_kmerseek_table.py` from `$SCRATCH/botryllus-kmerseek-swissprot/results/summary/kmerseek_botryllus_swissprot.per_protein.tsv.gz` on Sherlock |
| `unifire_predictions.tsv.gz` | UniFIRE 2026.3 predictions (ARBA, UniRule and PIRSR rules on InterProScan 5.59 matches) under the Biomni protein IDs: 231,595 rows on 19,001 proteins (ARBA 16,372, UniRule 11,575, PIRSR 1,069). The 50 proteins over 5,000 aa were not run. | `make_unifire_table.py`, from olgabot/botryllus-mhc PR 13, commit 257832d, `nextflow_runs/timing-pilot/unifire_whole_proteome` |

Two things about the tables that build_full.py works around:

- **Isoform IDs.** The FASTA repeats a header for each isoform of a gene. The tables call the
  first copy `FUN019272` and the second `FUN019272.2`, in file order. The 51-protein test set
  of `build.py` was named differently: its `FUN019272_TTN.1` is the tables' `FUN019272.2`. So
  build_full.py links a protein to its drawn page by sequence (MD5), not by ID.
- **The `mmseqs_cluster` column is wrong.** It is filled only on each cluster's representative,
  and for the 5,808 clusters of two or more proteins it names one other member instead of the
  representative. build_full.py takes clusters from `botryllus_mmseqs_clusters.tsv` only.

None of the tables gives where on the protein a hit or domain lies, so the pages list them and
draw nothing along the protein. The plan says the Pfam step kept domain positions; asking for
that file, and for DIAMOND's `qstart qend`, would let the pages draw them. The request is
drafted in `biomni_request_positions.md`.
