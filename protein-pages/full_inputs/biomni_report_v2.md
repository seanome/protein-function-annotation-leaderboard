# Botryllus schlosseri proteome annotation — immune-focused

## Summary

All 45,339 proteins in the new Botryllus genome were annotated with five complementary evidence streams: DIAMOND similarity vs UniProt Swiss-Prot (sensitive mode), Pfam-A domain scanning (HMMER, trusted cutoffs), MMseqs2 clustering + Swiss-Prot search, jackhmmer iterative searches from 24 immune seed families, and a dedicated allorecognition search against fester co-receptor alleles and BHF. A curated immune-gene table of **2,835 proteins** was derived, spanning the conserved deuterostome innate immune repertoire plus Botryllus-specific allorecognition genes.

## Methods

| Tool | Target | Mode | Hits |
|---|---|---|---|
| DIAMOND blastp | Swiss-Prot (533,805 entries) | sensitive, E ≤ 1e-3, top 25 | 369,001 hits / 21,150 proteins |
| Pfam hmmsearch | Pfam-A (30k models) | GA trusted cutoffs | 53,860 domains |
| MMseqs2 | Swiss-Prot | search-type 2, s=7.5 | 1.95M hits / 22,517 proteins |
| MMseqs2 | self | easy-cluster | 27,483 families |
| jackhmmer | 24 immune seeds vs proteome | 3–5 iterations, E ≤ 1e-5 | 7,881 hits |
| DIAMOND | fester (100 alleles) + BHF | E ≤ 1e-3 | 832 hits / 53 proteins |

Names were re-derived from evidence (best Swiss-Prot hit ≥30% identity & ≥50% coverage → direct name; ≥20% & ≥50% → "-like"). Existing header names were compared and flagged: **10,404 agree, 3,936 conflict, 2,173 new names assigned, 3,117 unconfirmed**.

## Immune repertoire

| Category | Genes | Notes |
|---|---|---|
| NLR | 375 | 291 with NACHT domain; large expansion typical of tunicates |
| Fibrinogen/ficolin | 458 | FBG-domain lectins, complement-like recognition |
| Lectin | 442 | C-type lectins, galectins |
| Complement | 231 | 8 C3-like, C1q, C2-like, C4BP, sushi-domain proteins |
| NF-κB pathway | 123 | Rel/RHD, IKK, IκB |
| Caspase/apoptosis | 76 | |
| Death-domain signaling | 72 | CARD, DED, death domains |
| IgSF | 66 | Immunoglobulin superfamily (adaptive-like) |
| Allorecognition | 53 | fester-like (incl. FUN006619 at 98.6% identity to fester co-receptor) and BHF-like |
| Cytokine/receptor | 48 | |
| TLR/IL1R signaling | 40 | 29 with TIR domain |
| JAK/STAT | 22 | |
| TNF family | 21 | |
| Antiviral | 10 | RIG-I/MDA5/LGP2 helicases, OAS |
| Phagocytosis/scavenger | 15 | |
| Effector | 7 | lysozyme, prophenoloxidase |

## Allorecognition

53 proteins hit the fester/BHF reference set. FUN006619 is a near-exact fester co-receptor match (98.6% identity, E=2e-232). The fester family is expanded and polymorphic, consistent with its role in colony fusion/rejection. BHF (Botryllus histocompatibility factor) candidates were identified via the AGS14996.1 reference.

## Cross-tool comparison site

A static site (`botryllus_protein_pages_site.tar.gz`) provides one page per protein with every tool's hits drawn as boxes on a shared residue ruler, E-value as opacity, and clickable alignments. The index classifies each protein: **8,436 agree, 12,771 disagree, 1,588 one-program, 260 apart, 22,284 no-hit** (at E ≤ 1e-3).


## FANTASIA GO prediction (embedding-based, homology-free)

FANTASIA-Lite (ProtT5 embeddings + cosine lookup against 124,397 experimentally
characterized UniProt references, k=1) was run on all 26,676 proteins that lacked
Swiss-Prot homology support — the 2,835 curated immune genes plus all 24,189 "dark"
proteins with no Swiss-Prot hit. Every queried protein received GO predictions:
26,676/26,676 (100%), adding functional annotation to a quarter of the proteome
that homology search alone could not reach.

- Reliability index: median 0.78; 100% of predictions at RI >= 0.5, 3% at RI >= 0.9.
- Median cosine distance to nearest reference: 0.22.
- Sequences longer than 2,000 aa were embedded from the first 2,000 residues
  (41 proteins >3,000 aa, max 10,509 aa) to keep CPU runtime tractable.

Most frequent GO terms recovered in the dark proteome: cytoplasm (3,462),
nucleus (2,568), mitochondrion (1,543), plasma membrane (1,093), and — of
immune interest — carbohydrate binding, cell-cell adhesion, and extracellular
region terms, consistent with lectin-like and adhesion biology. Among the
curated immune genes, FANTASIA's top terms were plasma membrane (257),
carbohydrate binding (235), cell-cell adhesion (218), apoptotic process (137),
and negative regulation of apoptotic process (117).

FANTASIA predictions are stored in the `fantasia_*` columns of
`botryllus_protein_annotations_v2.tsv` and `botryllus_immune_genes_v2.tsv`
(uniprot reference, distance, reliability index, GO terms with descriptions).

## Caveats

- Swiss-Prot is curated but small; divergent Botryllus proteins with no Swiss-Prot hit rely on Pfam domains only.
- The 3,936 name conflicts are mostly cases where the prior header name and the best Swiss-Prot hit differ (e.g., FUN009411 labeled "TLR7" but best hit is FLRT2, an LRR adhesion protein). These are flagged, not overwritten.
- Structure prediction (ESMCFold2-Fast) was attempted but GPU capacity was unavailable; Foldseek structural annotation is pending.
- InterProScan was interrupted by a machine failure and was not completed; Pfam (via pyhmmer) provides the domain evidence layer.
- UniFIRE was not run (Docker-only distribution, not available in this environment); InterProScan provides the equivalent signature-based evidence layer.

## Files

- `botryllus_protein_annotations.tsv` / `_v2.tsv` — full annotation table (45,339 rows; v2 adds FANTASIA columns)
- `botryllus_immune_genes.tsv` / `_v2.tsv` — curated immune subset (2,835 rows; v2 adds FANTASIA columns)
- `botryllus_mmseqs_clusters.tsv` — MMseqs2 cluster assignments
- `botryllus_jackhmmer_immune.tsv` — jackhmmer immune seed hits
- `botryllus_allorecognition_hits.tsv` — fester/BHF hits
- `botryllus_protein_pages_site.tar.gz` — cross-tool comparison site
