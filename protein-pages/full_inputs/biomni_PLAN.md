# Plan: Functional annotation of the Botryllus proteome (immune-focused)

## Summary
Annotate all 45,339 proteins in `botryllus_new_genome_protein.fasta` (~16.9 M aa; 15,829 already carry prior BLAST-based header names, 29,510 unnamed `FUN` IDs) using two evidence layers — DIAMOND similarity vs UniProt Swiss-Prot and Pfam domain scanning (HMMER) — then curate a prioritized immune-gene table covering the conserved deuterostome immune repertoire **plus** Botryllus allorecognition genes (fester, BHF). All prior header names are re-derived from evidence; disagreements with existing names are flagged, not silently overwritten.

## Workflow
1. **Setup (sandbox, worker-0)**: install `diamond` (static binary) and `pyhmmer` via `uv pip`; download UniProt Swiss-Prot (`uniprot_sprot.fasta.gz` for search + `uniprot_sprot.dat.gz` for GO terms/keywords/DE) and Pfam-A.hmm (+ descriptions). Download fester and BHF (Botryllus histocompatibility factor) reference protein sequences from NCBI.
2. **Similarity annotation**: `diamond makedb` on Swiss-Prot; `diamond blastp --sensitive` of all 45,339 proteins, top 25 hits, e-value ≤ 1e-3, report identity, coverage, e-value, bit score. Also run a small second DIAMOND db (fester + BHF + key allorecognition references) to catch allorecognition candidates directly.
3. **Domain annotation**: `pyhmmer` hmmsearch of all proteins against Pfam-A with trusted cutoffs (`--cut_ga`), 16 threads, run in background with checkpointing; collect per-protein domain architecture (accession, name, score, e-value, positions).
4. **Name assignment & conflict flagging**: for each protein, propose a gene name from best Swiss-Prot hit (preferably a named DE product) and/or diagnostic Pfam architecture; compare to the existing header name (case-insensitive, exact match) and flag conflicts in a dedicated column. Existing names are never overwritten.
5. **Immune curation**: merge evidence into an immune-gene table using (a) Pfam immune-domain list (TIR, NACHT/LRR/NLR, death domains, IgSF, C-type lectin, C1q, fibrinogen/ficolin, defensins/TIL, helicases, JAK/STAT, TNF family, caspases, complement components, etc.), (b) Swiss-Prot hits to immune-annotated proteins (immune keyword/GO), and (c) allorecognition hits (fester/BHF). Assign each immune gene a category (e.g. complement, PRR/TLR, NLR, lectin, antiviral, JAK/STAT, NF-κB pathway, IgSF/adaptive-like, allorecognition, effector/antimicrobial).

## Deliverables (in `/mnt/results/`)
- `botryllus_protein_annotations.tsv` — one row per protein: ID, existing header name, proposed name, best Swiss-Prot hit (accession, name, e-value, %identity, coverage), GO terms, Pfam domains, conflict flag.
- `botryllus_immune_genes.tsv` — curated immune subset with category, evidence (similarity/domain/both), and key domains.
- `report_botryllus_immune_annotation.md` — concise report: methods, immune-repertoire summary by family, allorecognition candidates, notable findings, caveats.
- Annotated FASTA (header names kept, annotations in tables only) if useful.

## Compute / execution target
- Default machine worker-0 (16 CPU / 16 GB RAM), all local — no HPC needed (Swiss-Prot is small; DIAMOND ~15–30 min).
- Pfam-A.hmm download ~1.9 GB; hmmsearch of 45k proteins ≈ 1–3 h → run as tracked background job with intermediate checkpointing to `/mnt/shared-workspace/`.
- Total wall time ≈ 3–4 h, dominated by Pfam scanning. Disk: ~5 GB scratch.

## Test cases / acceptance criteria
- Annotation table has exactly 45,339 rows, all input IDs present exactly once.
- Spot-check known immune genes: NLRC5-family proteins get NACHT/LRR + effector domains; TLRs get TIR + LRR; complement C3 gets the C3 domain architecture; fester/BHF hits recovered by the allorecognition db.
- Conflict flags populated for all proteins where proposed ≠ existing name; agreement rate reported.
- Immune table categories all non-empty and traceable to evidence columns.

## Assumptions
- Swiss-Prot (curated) is the preferred similarity target over full TrEMBL/nr for interpretability; divergent Botryllus genes with no good Swiss-Prot hit are covered by Pfam domains.
- Existing header names are treated as prior hypotheses to verify, not ground truth.
- Report is exploratory-grade (not publication-ready figures) unless requested.
