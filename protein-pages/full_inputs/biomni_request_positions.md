# Request to Biomni: where on each protein the hits and domains lie

Drafted 2026-10-29 for Olga to paste into the Biomni chat that made the tables in this folder.
Why: the tables say which Swiss-Prot entry, Pfam domain and Fester or BHF protein each
Botryllus protein matches, but not where on the protein. Without positions, build_full.py can
list the matches but cannot draw them along the protein. Step 3 of `biomni_PLAN.md` says the
Pfam step collected positions, so they may already be on worker-0 under `/mnt/shared-workspace/`.

The text below is the message to send.

---

Thank you for the five tables for the 45,339 Botryllus proteins (annotations, immune genes,
jackhmmer hits, Fester/BHF hits, MMseqs2 clusters). I need where on each protein every hit and
domain lies, so I can draw them. Please send these as tab-separated files with a header row.
Use the same protein IDs as in `botryllus_protein_annotations.tsv` (the first copy of a repeated
FASTA header is `FUN019272`, the second `FUN019272.2`). All positions 1-based and inclusive, as
HMMER and DIAMOND print them.

1. **`botryllus_pfam_domains.tsv`**: one row per Pfam domain hit that passed the gathering
   threshold (`--cut_ga`), from the pyhmmer hmmsearch of step 3 of your plan. If the hits from
   that run are saved (checkpoints in `/mnt/shared-workspace/`), please export those rather than
   running it again. Columns:
   `protein_id, pfam_acc, pfam_name, ali_from, ali_to, env_from, env_to, hmm_from, hmm_to, model_length, score_bits, i_evalue, c_evalue`

2. **`botryllus_diamond_swissprot.tsv`**: the DIAMOND blastp `--sensitive` hits against
   Swiss-Prot from step 2, all hits kept (top 25 per protein, E-value at most 0.001), not only
   the best. Columns as DIAMOND `--outfmt 6` with these fields, in this order:
   `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen`

3. **`botryllus_diamond_allorecognition.tsv`**: the same 14 columns for the DIAMOND search
   against the Fester and BHF references (the 832 rows of `botryllus_allorecognition_hits.tsv`,
   with `qstart qend sstart send` added).

4. **The commands that made each table**: the exact command line or pyhmmer call, the tool
   versions, and the Swiss-Prot and Pfam release used. Also, for the jackhmmer table: the
   number of rounds, the E-value cuts (`-E`, `--incE`, `--domE`, `--incdomE`) and the database
   searched; and for the MMseqs2 clusters: the command and its `--min-seq-id`, `-c` and
   `--cov-mode`. Your notebook on worker-0 holds tool lookups and an ESMCFold test, but not
   these commands.

5. If it is quick: **`botryllus_jackhmmer_domains.tsv`**, the per-domain table of the jackhmmer
   hits (`--domtblout`), with `seed, iteration, target, ali_from, ali_to, env_from, env_to, hmm_from, hmm_to, score_bits, i_evalue`.

Two things I found in the tables, in case they matter elsewhere:

- The `mmseqs_cluster` column of `botryllus_protein_annotations.tsv` is filled only on
  cluster representatives, and for the 5,808 clusters of two or more proteins it names one
  other member instead of the representative. `botryllus_mmseqs_clusters.tsv` looks right, so I
  use that.
- jackhmmer hits do not set an immune category: 2,183 of the 3,589 proteins with a jackhmmer
  hit are not in `botryllus_immune_genes.tsv`. Was that on purpose?

---

## What build_full.py will do with the files (not written yet)

- Pfam: a box per domain from `env_from` to `env_to`, on a line the length of the protein.
- DIAMOND: a bar per hit from `qstart` to `qend` under the protein line, for Swiss-Prot and
  for Fester or BHF.
- Checks before drawing: every `protein_id` is in `proteome_lengths.tsv.gz`, every `qlen`
  equals that length, and every end position is at most the length. A mismatch would mean the
  IDs are numbered differently, so the build should stop.
