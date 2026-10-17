"""Build data/known_cases/cases.json: fourteen protein pairs whose shared function was proven
by experiment, the input each tool would be given, the decoy, and whatever result each tool has.

The case list below is fixed. A case that cannot be run keeps its row and says why. Sequences
come from UniProt's REST API and are cached in data/known_cases/sequences.json, so a rebuild
without network access uses the cached copy.

Run: .venv/bin/python scripts/known_cases.py --out data/known_cases
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UNIPROT = "https://rest.uniprot.org/uniprotkb/"
# The release every search in this page's run searches. The leaderboard's own run searched the
# Quest for Orthologs proteomes instead; Methods says so.
TARGET_DB = "UniProtKB/Swiss-Prot release 2026_03"
SHUFFLE_SEED = 1

# Evidence codes from the Evidence and Conclusion Ontology (https://evidenceontology.org),
# each checked against the ontology on 2026-10-17.
ECO = {
    "ECO:0000012": "functional complementation evidence",
    "ECO:0005670": "x-ray crystallography evidence used in manual assertion",
    "ECO:0000314": "direct assay evidence used in manual assertion",
    "ECO:0000315": "mutant phenotype evidence used in manual assertion",
    "ECO:0000353": "physical interaction evidence used in manual assertion",
    "ECO:0005801": "enzymatic activity assay evidence used in manual assertion",
    "ECO:0000269": "experimental evidence used in manual assertion",
}

BOXES = {
    1: {
        "title": "Shared ancestor, same job",
        "right": "the target's family and its function",
        "words": "The two proteins descend from one ancestral protein and still do the same job. "
                 "A call is right when it names the target's family and transfers its function.",
    },
    2: {
        "title": "Shared ancestor, new job",
        "right": "the family only; transferring the function is a wrong label",
        "words": "The two proteins descend from one ancestral protein, but the query protein has "
                 "taken on a different job. A call is right when it names the family, and wrong "
                 "when it carries the target's function across.",
    },
    3: {
        "title": "No shared ancestor, whole protein",
        "right": "the function only; a family call is wrong by ancestry and right by function",
        "words": "The two proteins reached the same job from different ancestors, so there is no "
                 "family truth. A call that names the target's family is wrong by ancestry and "
                 "right by function, and the card shows both marks.",
    },
    4: {
        "title": "No shared ancestor, short motif",
        "right": "a call covering the right few residues",
        "words": "A short stretch of the query protein copies a short stretch of a human protein. "
                 "A call is right when it lands on those residues.",
    },
}

REFERENCES = [
    {"n": 1, "text": 'Hengartner, M. O., and H. R. Horvitz. "C. elegans Cell Survival Gene ced-9 '
                     'Encodes a Functional Homolog of the Mammalian Proto-oncogene bcl-2." '
                     '<i>Cell</i> 76, no. 4 (1994): 665–76.',
     "doi": "10.1016/0092-8674(94)90506-1"},
    {"n": 2, "text": 'Vaux, D. L., I. L. Weissman, and S. K. Kim. "Prevention of Programmed Cell '
                     'Death in Caenorhabditis elegans by Human bcl-2." <i>Science</i> 258, '
                     'no. 5090 (1992): 1955–57.',
     "doi": "10.1126/science.1470921"},
    {"n": 3, "text": 'Kvansakul, M., M. F. van Delft, E. F. Lee, J. M. Gulbis, W. D. Fairlie, '
                     'D. C. S. Huang, and P. M. Colman. "A Structural Viral Mimic of Prosurvival '
                     'Bcl-2: A Pivotal Role for Sequestering Proapoptotic Bax and Bak." '
                     '<i>Molecular Cell</i> 25, no. 6 (2007): 933–42.',
     "doi": "10.1016/j.molcel.2007.02.004"},
    {"n": 4, "text": 'Graham, S. C., M. W. Bahar, S. Cooray, R. A.-J. Chen, D. M. Whalen, '
                     'N. G. A. Abrescia, D. Alderton, et al. "Vaccinia Virus Proteins A52 and B14 '
                     'Share a Bcl-2-like Fold but Have Evolved to Inhibit NF-κB Rather Than '
                     'Apoptosis." <i>PLoS Pathogens</i> 4, no. 8 (2008): e1000128.',
     "doi": "10.1371/journal.ppat.1000128"},
    {"n": 5, "text": 'Löwe, J., and L. A. Amos. "Crystal Structure of the Bacterial Cell-Division '
                     'Protein FtsZ." <i>Nature</i> 391, no. 6663 (1998): 203–6.',
     "doi": "10.1038/34472"},
    {"n": 6, "text": 'van den Ent, F., L. A. Amos, and J. Löwe. "Prokaryotic Origin of the Actin '
                     'Cytoskeleton." <i>Nature</i> 413, no. 6851 (2001): 39–44.',
     "doi": "10.1038/35092500"},
    {"n": 7, "text": 'Bork, P., C. Sander, and A. Valencia. "An ATPase Domain Common to '
                     'Prokaryotic Cell Cycle Proteins, Sugar Kinases, Actin, and hsp70 Heat Shock '
                     'Proteins." <i>Proceedings of the National Academy of Sciences of the United '
                     'States of America</i> 89, no. 16 (1992): 7290–94.',
     "doi": "10.1073/pnas.89.16.7290"},
    {"n": 8, "text": 'Schumacher, M. A., R. F. Pearson, T. Møller, P. Valentin-Hansen, and '
                     'R. G. Brennan. "Structures of the Pleiotropic Translational Regulator Hfq '
                     'and an Hfq–RNA Complex: A Bacterial Sm-like Protein." <i>The EMBO Journal</i> '
                     '21, no. 13 (2002): 3546–56.',
     "doi": "10.1093/emboj/cdf322"},
    {"n": 9, "text": 'Essuman, K., D. W. Summers, Y. Sasaki, X. Mao, A. K. Y. Yim, A. DiAntonio, '
                     'and J. Milbrandt. "TIR Domain Proteins Are an Ancient Family of '
                     'NAD-Consuming Enzymes." <i>Current Biology</i> 28, no. 3 (2018): 421–430.e4.',
     "doi": "10.1016/j.cub.2017.12.024"},
    {"n": 10, "text": 'Buchwald, G., A. Friebel, J. E. Galán, W.-D. Hardt, A. Wittinghofer, and '
                      'K. Scheffzek. "Structural Basis for the Reversible Activation of a Rho '
                      'Protein by the Bacterial Toxin SopE." <i>The EMBO Journal</i> 21, no. 13 '
                      '(2002): 3286–95.',
     "doi": "10.1093/emboj/cdf329"},
    {"n": 11, "text": 'Stebbins, C. E., and J. E. Galán. "Modulation of Host Signaling by a '
                      'Bacterial Mimic: Structure of the Salmonella Effector SptP Bound to Rac1." '
                      '<i>Molecular Cell</i> 6, no. 6 (2000): 1449–60.',
     "doi": "10.1016/s1097-2765(00)00141-6"},
    {"n": 12, "text": 'Cameron, C. M., J. W. Barrett, M. Mann, A. Lucas, and G. McFadden. "Myxoma '
                      'Virus M128L Is Expressed as a Cell Surface CD47-like Virulence Factor That '
                      'Contributes to the Downregulation of Macrophage Activation in Vivo." '
                      '<i>Virology</i> 337, no. 1 (2005): 55–67.',
     "doi": "10.1016/j.virol.2005.03.037"},
    {"n": 13, "text": 'Tal, M. C., P. S. Hansen, H. A. Ogasawara, Q. Feng, R. F. Volk, B. Lee, '
                      'S. E. Casebeer, et al. "P66 Is a Bacterial Mimic of CD47 That Binds the '
                      'Anti-phagocytic Receptor SIRPα and Facilitates Macrophage Evasion by '
                      'Borrelia burgdorferi." bioRxiv preprint, April 30, 2024.',
     "doi": "10.1101/2024.04.29.591704", "note": "Preprint, not peer reviewed."},
    {"n": 14, "text": 'Marazzi, I., J. S. Y. Ho, J. Kim, B. Manicassamy, S. Dewell, R. A. Albrecht, '
                      'C. W. Seibert, et al. "Suppression of the Antiviral Response by an '
                      'Influenza Histone Mimic." <i>Nature</i> 483, no. 7390 (2012): 428–33.',
     "doi": "10.1038/nature10892"},
    {"n": 15, "text": 'Chen, L., A. L. DeVries, and C. H. Cheng. "Convergent Evolution of '
                      'Antifreeze Glycoproteins in Antarctic Notothenioid Fish and Arctic Cod." '
                      '<i>Proceedings of the National Academy of Sciences of the United States of '
                      'America</i> 94, no. 8 (1997): 3817–22.',
     "doi": "10.1073/pnas.94.8.3817"},
    {"n": 16, "text": 'Kumar, M., S. Michael, J. Alvarado-Valverde, B. Mészáros, '
                      'H. Sámano-Sánchez, A. Zeke, L. Dobson, et al. "The Eukaryotic Linear Motif '
                      'Resource: 2022 Release." <i>Nucleic Acids Research</i> 50, no. D1 (2022): '
                      'D497–D508.',
     "doi": "10.1093/nar/gkab975"},
    {"n": 17, "text": 'Cheng, H., B.-H. Kim, and N. V. Grishin. "MALISAM: A Database of '
                      'Structurally Analogous Motifs in Proteins." <i>Nucleic Acids Research</i> '
                      '36, Database issue (2008): D211–17.',
     "doi": "10.1093/nar/gkm698"},
    {"n": 18, "text": 'Davey, N. E., G. Travé, and T. J. Gibson. "How Viruses Hijack Cell '
                      'Regulation." <i>Trends in Biochemical Sciences</i> 36, no. 3 (2011): '
                      '159–69.',
     "doi": "10.1016/j.tibs.2010.10.002"},
]

# ---------------------------------------------------------------- the fixed case list
# Each case: the query protein searched, the human (or Arctic cod) protein whose label is at
# stake, the box, the experiment that proved the shared function, and the input each tool gets.
# "structure" is the experimental structure the cited paper reports, or None for an AlphaFold
# Database model. Coordinates are 1-based and inclusive, as UniProt writes them.
CASES = [
    {
        "id": 1, "box": 1,
        "name": "CED-9 against BCL-2",
        "query": {"acc": "P41958", "label": "CED-9", "species": "Caenorhabditis elegans"},
        "target": {"acc": "P10415", "label": "BCL-2", "species": "human"},
        "family": {"id": "PF00452", "name": "Bcl-2"},
        "experiment": "Human bcl-2 put into C. elegans stopped the programmed cell death that "
                      "loss of ced-9 causes, and ced-9 itself carries the same activity.",
        "takeaway": "Two proteins that do the same job can sit near the edge of what a sequence "
                    "search reports. A weak hit to a protein of the right family is worth "
                    "reading, not throwing away.",
        "earlier_runs": "An earlier search of CED-9 against the human GENCODE v49 canonical "
                        "proteome is stored in seanome/2024-kmerseek-analysis, branch "
                        "olgabot/case-gallery-reference-pairs, commit 3e3767c1, files "
                        "tables/case_gallery/other_tools/*.tsv. phmmer puts CED-9 residues 111 to "
                        "232 on BCL-2 residues 68 to 235, E = 7.9e-05. jackhmmer puts 100 to 229 "
                        "on 73 to 208, E = 7.7e-41. MMseqs2 with three rounds puts 101 to 236 on "
                        "82 to 211, E = 8.1e-20, and with one round reports no BCL-2 row at all. "
                        "Those runs had no decoy and searched a different database, so they are "
                        "not in the table above.",
        "eco": "ECO:0000012", "refs": [1, 2],
        "structure": None,
    },
    {
        "id": 2, "box": 1,
        "name": "Myxoma virus M11L against BCL-2",
        "query": {"acc": "Q85295", "label": "M11L", "species": "myxoma virus"},
        "target": {"acc": "P10415", "label": "BCL-2", "species": "human"},
        "family": {"id": "PF00452", "name": "Bcl-2"},
        "experiment": "The crystal structure of M11L bound to a Bak BH3 peptide shows the Bcl-2 "
                      "fold and the same groove that BCL-2 uses to hold Bax and Bak.",
        "takeaway": "A virus protein can hold a human fold with almost no sequence left in "
                    "common. A tool that reads only letters may report nothing here, even though "
                    "the label is right.",
        "eco": "ECO:0005670", "refs": [3],
        "structure": {"pdb": "2JBX", "chain": "A", "what": "M11L, 1.7 Å crystal structure"},
    },
    {
        "id": 3, "box": 1,
        "name": "Myxoma virus M128L against CD47",
        "query": {"acc": "Q9Q8H3", "label": "M128L", "species": "myxoma virus"},
        "target": {"acc": "Q08722", "label": "CD47", "species": "human"},
        "family": {"id": "PF07686", "name": "V-set immunoglobulin domain"},
        "experiment": "M128L sits on the surface of infected cells as a CD47-like protein, and "
                      "rabbits infected with virus lacking it had milder disease and more "
                      "activated macrophages.",
        "takeaway": "When a virus carries a copy of a host surface protein, the right label is "
                    "the host protein's own family. A tool that misses it leaves a virulence "
                    "factor unlabelled.",
        "eco": "ECO:0000315", "refs": [12],
        "structure": None,
    },
    {
        "id": 4, "box": 1,
        "name": "FtsZ against alpha- and beta-tubulin",
        "query": {"acc": "Q57816", "label": "FtsZ", "species": "Methanocaldococcus jannaschii"},
        "target": {"acc": "P07437", "label": "beta-tubulin", "species": "human"},
        "extra_targets": [{"acc": "Q71U36", "label": "alpha-tubulin", "species": "human"}],
        "family": {"id": "PF00091", "name": "Tubulin/FtsZ GTPase domain"},
        "experiment": "The crystal structure of FtsZ has the tubulin fold and the same "
                      "GTP-binding site. A bacterial cell-division protein and a human "
                      "cytoskeleton protein are one family.",
        "takeaway": "One ancient family can keep its fold and its binding site while the sequence "
                    "drifts far apart. Check how much of the call sits on the real domain, not "
                    "only whether the family is named.",
        "eco": "ECO:0005670", "refs": [5],
        "structure": {"pdb": "1FSZ", "chain": "A", "what": "FtsZ, 2.8 Å crystal structure"},
    },
    {
        "id": 5, "box": 1,
        "name": "MreB against actin",
        "query": {"acc": "Q9WZ57", "label": "MreB", "species": "Thermotoga maritima"},
        "target": {"acc": "P60709", "label": "beta-actin", "species": "human"},
        "family": {"id": "PF06723", "name": "MreB/Mbl, the actin-like ATPase fold"},
        "experiment": "The crystal structure of MreB has the actin fold, and MreB forms protein "
                      "filaments with the same subunit spacing as an actin filament.",
        "takeaway": "The same fold carries actin's job in a bacterium, so a family call here is "
                    "both correct and far from obvious from the sequence alone.",
        "eco": "ECO:0005670", "refs": [6, 7],
        "structure": {"pdb": "1JCE", "chain": "A", "what": "MreB, 1.9 Å crystal structure"},
    },
    {
        "id": 6, "box": 1,
        "name": "Hfq against the Sm proteins",
        "query": {"acc": "A0A0H3JV59", "label": "Hfq", "species": "Staphylococcus aureus"},
        "target": {"acc": "P62314", "label": "Sm D1", "species": "human"},
        "family": {"id": "PF01423", "name": "LSM domain"},
        "experiment": "The crystal structures of Hfq alone and bound to RNA show the Sm fold. Its "
                      "subunits form a ring that holds single-stranded RNA, as the human Sm "
                      "proteins do.",
        "takeaway": "Both proteins here are 77 residues long, which leaves a search little to "
                    "work with. A miss says more about the length than about the pair.",
        "eco": "ECO:0005670", "refs": [8],
        "structure": {"pdb": "1KQ1", "chain": "A", "what": "Hfq, 1.55 Å crystal structure"},
    },
    {
        "id": 7, "box": 2,
        "name": "Vaccinia virus A52 and B14 against BCL-2",
        "query": {"acc": "Q01220", "label": "A52", "species": "vaccinia virus"},
        "extra_queries": [{"acc": "P24772", "label": "B14", "species": "vaccinia virus",
                           "structure": {"pdb": "2VVY", "chain": "A",
                                         "what": "B14, 2.7 Å crystal structure"}}],
        "target": {"acc": "P10415", "label": "BCL-2", "species": "human"},
        "family": {"id": "PF00452", "name": "Bcl-2"},
        "experiment": "The crystal structures of A52 and B14 have the Bcl-2 fold but no groove "
                      "for a BH3 peptide. Neither protein changes programmed cell death; both "
                      "block NF-κB signalling.",
        "takeaway": "A52 and B14 have the Bcl-2 fold and none of the apoptosis job. Reading the "
                    "family off the fold is safe here; reading the job off the family is not.",
        "eco": "ECO:0005670", "refs": [4],
        "structure": {"pdb": "2VVW", "chain": "A", "what": "A52, 1.9 Å crystal structure"},
        "wrong_label": "stops programmed cell death (the job of BCL-2)",
    },
    {
        "id": 8, "box": 2,
        "name": "A bacterial TIR-domain protein against a human TIR domain",
        "query": {"acc": "P0DTS9", "label": "TirS", "species": "Staphylococcus aureus"},
        "target": {"acc": "Q6SZW1", "label": "SARM1", "species": "human"},
        "family": {"id": "PF01582", "name": "TIR domain"},
        "experiment": "Purified TirS cleaves NAD, and the full-length protein lowers NAD in human "
                      "cells, which the catalytically dead E216A protein does not. The human TIR "
                      "domains of the Toll-like receptors hold protein complexes together "
                      "instead.",
        "takeaway": "A domain shared between a bacterium and a human can carry a different job in "
                    "each. The safe transfer is the domain name, not the human protein's role.",
        "eco": "ECO:0005801", "refs": [9],
        "structure": None,
        "wrong_label": "holds a Toll-like receptor signalling complex together",
    },
    {
        "id": 9, "box": 3,
        "name": "P66 against CD47",
        "query": {"acc": "H7C7N8", "label": "P66", "species": "Borrelia burgdorferi"},
        "target": {"acc": "Q08722", "label": "CD47", "species": "human"},
        "family": None,
        "experiment": "P66 binds the human receptor SIRPα, the receptor CD47 binds, and Borrelia "
                      "lacking P66 is taken up by macrophages more readily.",
        "takeaway": "Two proteins can do the same job with no shared ancestry, so a search that "
                    "finds nothing here is not wrong about the family. It says nothing about the "
                    "job.",
        "earlier_runs": "An earlier search of P66 against the human GENCODE v49 canonical "
                        "proteome, same branch and commit, returns no CD47 row from phmmer, "
                        "jackhmmer or MMseqs2. Foldseek on the two AlphaFold models gives one "
                        "hit, E = 3.047, bit score 3, over 22 aligned residues "
                        "(seanome/2024-kmerseek-analysis notebook 276, commit 7043bd98). Those "
                        "runs had no decoy and searched a different database.",
        "eco": "ECO:0000353", "refs": [13],
        "structure": None,
        "preprint": True,
        "wrong_label": None,
    },
    {
        "id": 10, "box": 3,
        "name": "SopE against a Dbl-family exchange factor",
        "query": {"acc": "O52623", "label": "SopE", "species": "Salmonella typhimurium"},
        "target": {"acc": "Q13009", "label": "TIAM1", "species": "human"},
        "family": None,
        "experiment": "The crystal structure of SopE bound to Cdc42 shows SopE swapping the "
                      "nucleotide of a Rho protein. That is the job of the human Dbl-family "
                      "exchange factors, done from a fold unlike theirs.",
        "takeaway": "A bacterial toxin can take over a human enzyme's job with an unrelated fold. "
                    "No search here is expected to find it, and no hit is not evidence that the "
                    "protein is harmless.",
        "eco": "ECO:0005670", "refs": [10],
        "structure": {"pdb": "1GZS", "chain": "B", "what": "SopE with Cdc42, 2.3 Å crystal structure"},
    },
    {
        "id": 11, "box": 3,
        "name": "SptP against a RhoGAP",
        "query": {"acc": "P74873", "label": "SptP", "species": "Salmonella typhimurium"},
        "target": {"acc": "Q07960", "label": "p50-RhoGAP", "species": "human"},
        "family": None,
        "experiment": "The crystal structure of SptP bound to Rac1 shows SptP switching a Rho "
                      "protein off with an arginine finger. That is the job of the human RhoGAP "
                      "proteins, done from a fold unlike theirs.",
        "takeaway": "The same job with a different fold again. A tool that names a human RhoGAP "
                    "family here is right about what the protein does and wrong about where it "
                    "came from.",
        "eco": "ECO:0005670", "refs": [11],
        "structure": {"pdb": "1G4U", "chain": "S", "what": "SptP with Rac1, 2.3 Å crystal structure"},
    },
    {
        "id": 12, "box": 3,
        "name": "Antifreeze glycoprotein, Antarctic notothenioid against Arctic cod",
        "query": {"acc": "P24856", "label": "antifreeze glycoprotein polyprotein",
                  "species": "Notothenia neglecta (Antarctic)"},
        "target": {"acc": "A0A481T0N3", "label": "antifreeze glycoprotein polypeptide 9",
                   "species": "Boreogadus saida (Arctic cod)"},
        "family": None,
        "experiment": "Both fish make the same ice-binding glycoprotein from a repeat of the same "
                      "three residues. The notothenioid gene came from a trypsinogen and the cod "
                      "gene did not, and the two genes differ in their introns, spacers and codon "
                      "use.",
        "takeaway": "Two fish lineages arrived at the same repeated glycoprotein from different "
                    "genes, so a strong sequence match is not by itself evidence of common "
                    "ancestry.",
        "eco": "ECO:0000269", "refs": [15],
        "structure": None,
        "target_rule": "Of the Boreogadus saida antifreeze glycoprotein polypeptides in UniProt, "
                       "the longest, so that its length is closest to the notothenioid "
                       "polyprotein's.",
    },
    {
        "id": 13, "box": 4,
        "name": "Influenza A NS1 against the histone H3 tail",
        "query": {"acc": "P03495", "label": "NS1", "species": "influenza A H3N2, A/Udorn/307/1972"},
        "query_region": [226, 229],
        "target": {"acc": "P68431", "label": "histone H3.1", "species": "human"},
        "target_region": [2, 5],
        "family": None,
        "experiment": "The last residues of NS1 copy the start of the histone H3 tail. NS1 binds "
                      "the human PAF1 transcription complex through them, and changing them ends "
                      "the binding.",
        "takeaway": "Four residues are shorter than the shortest word any of these tools indexes. "
                    "Every search here is expected to be silent, and a motif scanner is the right "
                    "tool.",
        "eco": "ECO:0000315", "refs": [14],
        "structure": {"pdb": "4O42", "chain": "B",
                      "what": "the NS1 C-terminal tail, 15 residues, bound to the CHD1 "
                              "chromodomains"},
        "structure_note": "The only experimental structure of this stretch holds 15 residues of "
                          "NS1, not the whole protein, and it is the tail of a different H3N2 "
                          "strain (UniProt Q38SQ2).",
        "motif_note": "UniProt numbers the histone from its initial methionine, so the ARTK of "
                      "the H3 tail is residues 2 to 5 there and residues 1 to 4 in the "
                      "chromatin literature.",
    },
    {
        "id": 14, "box": 4,
        "name": "Ten pathogen motifs from ELM against their human partners",
        "query": None,
        "target": None,
        "family": None,
        "experiment": "Each instance is a short motif in a virus or bacterial protein. The "
                      "Eukaryotic Linear Motif resource records the experiments that showed it "
                      "binds a human protein.",
        "takeaway": "Short pathogen motifs are the part of annotation these nine tools are least "
                    "suited to. That is why the leaderboard reports region length as its own "
                    "axis.",
        "eco": "ECO:0000353", "refs": [16, 18],
        "structure": None,
        "elm": True,
        "pick_rule": "The ten ELM instances whose protein is from a virus or a bacterium "
                     "(UniProt lineage), whose recorded partner is a human protein, with the "
                     "lowest ELM instance identifiers. Chosen before any search was run.",
    },
]

# The nine tools and the setting Table 1 calls best for each, read from the leaderboard tables by
# build_site.card_numbers on 2026-10-17 and written here so this script does not need the parquet
# files. tests/check_known_cases.py checks each one against the tables.
TOOL_SETTINGS = {
    "kmerseek": "gbmr4_k21_lcFalse",
    "hmmer3_phmmer": "default",
    "hmmer3_jackhmmer": "n3",
    "mmseqs2_seqseq": "s7",
    "mmseqs2_iterative": "s7",
    "hhblits": "single_seq",
    "foldseek": "3di_aa",
    "prostt5": "3di_from_seq",
    "hmmscan": "direct",
}

# What each tool is given for a case: its sequence, or the structure named on the card.
STRUCTURE_TOOLS = {"foldseek"}

RESULT_WORDS = {
    "found": "a call on the right target, overlapping the truth region with an IoU of 0.5 or more",
    "wrong place": "a call on the right target, overlapping the truth region with an IoU below 0.5",
    "wrong label": "the call carries the human protein's function to a protein that does not have it",
    "missed": "no call on the target above the tool's own cutoff",
    "not run": "the search has not been run yet; the card says why",
}


# ---------------------------------------------------------------- sequences and decoys

def fetch_fasta(acc: str) -> str:
    url = f"{UNIPROT}{urllib.parse.quote(acc)}.fasta"
    txt = urllib.request.urlopen(url, timeout=120).read().decode()
    return "".join(l.strip() for l in txt.splitlines()[1:])


def fetch_meta(acc: str) -> dict:
    url = (f"{UNIPROT}search?" + urllib.parse.urlencode(
        {"query": f"accession:{acc}", "fields": "accession,id,protein_name,organism_name,length,reviewed",
         "format": "tsv", "size": 1}))
    txt = urllib.request.urlopen(url, timeout=120).read().decode().splitlines()
    if len(txt) < 2:
        return {}
    f = txt[1].split("\t")
    return {"entry_name": f[1], "protein_name": f[2], "organism": f[3], "length": int(f[4]),
            "reviewed": f[5] == "reviewed"}


def dipeptide_shuffle(seq: str, rng: random.Random) -> str:
    """Shuffle a sequence keeping every dipeptide count, by walking a random Eulerian path over
    the graph whose nodes are residues and whose edges are the dipeptides of the sequence
    (Altschul and Erickson 1985). The first and last residue stay where they are."""
    if len(seq) < 3:
        return seq
    edges: dict[str, list[str]] = {}
    for x, y in zip(seq, seq[1:]):
        edges.setdefault(x, []).append(y)
    last = seq[-1]

    def reaches_last(nxt: dict[str, str]) -> bool:
        """Following the chosen last edge from every residue must end at the final residue."""
        for v in nxt:
            seen, u = {v}, v
            while u != last:
                u = nxt.get(u)
                if u is None or u in seen:
                    return False
                seen.add(u)
        return True

    for _ in range(500):
        # One outgoing edge per residue, kept for the end of that residue's edge list; the walk
        # can then never leave a residue it still has edges to reach (Altschul and Erickson 1985).
        nxt = {v: rng.choice(e) for v, e in edges.items() if v != last}
        if reaches_last(nxt):
            break
    else:
        return seq

    pool = {v: list(e) for v, e in edges.items()}
    for v, u in nxt.items():
        pool[v].remove(u)
    for v in pool:
        rng.shuffle(pool[v])
    for v, u in nxt.items():
        pool[v].append(u)

    out = [seq[0]]
    v = seq[0]
    for _ in range(len(seq) - 1):
        if not pool.get(v):
            break
        u = pool[v].pop(0)
        out.append(u)
        v = u
    if len(out) != len(seq):  # the walk did not use every dipeptide; keep the real order
        return seq
    return "".join(out)


# ---------------------------------------------------------------- alignment

def blosum62(path: Path) -> tuple[dict, list[str]]:
    rows = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    cols = rows[0].split()
    m = {}
    for line in rows[1:]:
        f = line.split()
        for c, v in zip(cols, f[1:]):
            m[(f[0], c)] = int(v)
    return m, cols


def smith_waterman(a: str, b: str, mat: dict, gap_open: int = 11, gap_extend: int = 1):
    """Best local alignment of two sequences, with BLOSUM62 and affine gaps. Returns the two
    aligned strings, the 1-based start of each, and the score. Written out here because the page
    must show the residues of pairs no tool in the run aligned."""
    n, m = len(a), len(b)
    neg = -10**6
    M = np.full((n + 1, m + 1), neg, dtype=np.int32)
    Ix = np.full((n + 1, m + 1), neg, dtype=np.int32)
    Iy = np.full((n + 1, m + 1), neg, dtype=np.int32)
    M[0, :] = 0
    M[:, 0] = 0
    ptr = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 stop, 1 diagonal, 2 up (gap in b), 3 left
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            s = mat.get((ai, b[j - 1]), -4)
            diag = max(M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1], 0) + s
            Ix[i, j] = max(M[i - 1, j] - gap_open, Ix[i - 1, j] - gap_extend)
            Iy[i, j] = max(M[i, j - 1] - gap_open, Iy[i, j - 1] - gap_extend)
            cell = max(diag, Ix[i, j], Iy[i, j], 0)
            M[i, j] = cell
            ptr[i, j] = 0 if cell == 0 else (1 if cell == diag else (2 if cell == Ix[i, j] else 3))
            if cell > best:
                best, bi, bj = int(cell), i, j
    ra, rb = [], []
    i, j = bi, bj
    while i > 0 and j > 0 and ptr[i, j] != 0:
        if ptr[i, j] == 1:
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif ptr[i, j] == 2:
            ra.append(a[i - 1]); rb.append("-"); i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1
    return "".join(reversed(ra)), "".join(reversed(rb)), i + 1, j + 1, best


def match_line(a: str, b: str) -> tuple[str, int]:
    line = "".join("|" if x == y and x != "-" else " " for x, y in zip(a, b))
    return line, line.count("|")


# ---------------------------------------------------------------- build

def wrap60(s: str) -> str:
    return "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60))


def alphafold_model(acc: str) -> dict:
    """The AlphaFold Database file for an accession, or an empty answer when it has none. The
    model version moves (v4, v6), so the URL is read from the database rather than written here.
    The server answers 403 to a user agent naming curl, hence the browser-style string."""
    req = urllib.request.Request(
        f"https://alphafold.ebi.ac.uk/api/prediction/{urllib.parse.quote(acc)}",
        headers={"User-Agent": "Mozilla/5.0 (compatible; known-cases; "
                               "+https://github.com/seanome/"
                               "protein-function-annotation-leaderboard)"})
    try:
        d = json.load(urllib.request.urlopen(req, timeout=60))
    except Exception:
        return {}
    if not d:
        return {}
    return {"url": d[0].get("cifUrl"), "version": d[0].get("latestVersion")}


def tool_input(case: dict, tool: str, query: dict) -> dict:
    if tool in STRUCTURE_TOOLS:
        st = query.get("structure") or case.get("structure")
        if st:
            return {"kind": "experimental structure",
                    "what": f'PDB {st["pdb"]} chain {st["chain"]}: {st["what"]}'}
        af = query.get("alphafold") or {}
        if af.get("url"):
            return {"kind": "AlphaFold Database model",
                    "what": f'the AlphaFold Database model of {query["acc"]}, version '
                            f'{af.get("version")}', "url": af["url"]}
        return {"kind": "no structure",
                "what": f'no experimental structure in the papers cited here, and the AlphaFold '
                        f'Database has no model of {query["acc"]}'}
    if tool == "prostt5":
        return {"kind": "sequence",
                "what": "the sequence; ProstT5 predicts the 3-D shape letters from it"}
    return {"kind": "sequence", "what": "the sequence"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/known_cases")
    ap.add_argument("--elm", default="data/known_cases/elm_picked.json")
    ap.add_argument("--matrix", default="data/known_cases/BLOSUM62")
    ap.add_argument("--offline", action="store_true", help="use the cached sequences only")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    seq_f = out / "sequences.json"
    seqs = json.loads(seq_f.read_text()) if seq_f.exists() else {}

    def seq_of(acc: str) -> dict:
        if acc not in seqs:
            if a.offline:
                raise SystemExit(f"{acc} is not in {seq_f} and --offline was given")
            s = fetch_fasta(acc)
            seqs[acc] = {"sequence": s, "sha256": hashlib.sha256(s.encode()).hexdigest(),
                         **fetch_meta(acc)}
        return seqs[acc]

    af_f = out / "alphafold.json"
    af_cache = json.loads(af_f.read_text()) if af_f.exists() else {}
    elm_picked = json.loads(Path(a.elm).read_text())
    mat, _ = blosum62(Path(a.matrix))
    rng = random.Random(SHUFFLE_SEED)

    cases = []
    for c in CASES:
        row = {k: v for k, v in c.items() if k not in ("extra_queries", "extra_targets")}
        queries = []
        if c.get("elm"):
            for p in elm_picked:
                queries.append({"acc": p["Primary_Acc"], "label": p["ProteinName"],
                                "species": p["Organism"], "elm_instance": p["Accession"],
                                "elm_class": p["ELMIdentifier"],
                                "region": [int(p["Start"]), int(p["End"])],
                                "partner": {"acc": p["partner"]["interactorDomain"],
                                            "pfam": p["partner"]["Domain"],
                                            "pmid": p["partner"]["PMID"]}})
        else:
            queries.append(dict(c["query"]))
            queries += [dict(q) for q in c.get("extra_queries", [])]
        targets = []
        if c.get("elm"):
            for p in elm_picked:
                targets.append({"acc": p["partner"]["interactorDomain"], "species": "human"})
        else:
            targets.append(dict(c["target"]))
            targets += [dict(t) for t in c.get("extra_targets", [])]

        for q in queries:
            s = seq_of(q["acc"])
            q["length"] = s["length"] if s.get("length") else len(s["sequence"])
            q["protein_name"] = s.get("protein_name", "")
            q["reviewed"] = s.get("reviewed")
            q["decoy"] = {"method": "dipeptide-preserving shuffle (Altschul and Erickson 1985)",
                          "seed": SHUFFLE_SEED,
                          "sequence": dipeptide_shuffle(s["sequence"], rng)}
            if not (q.get("structure") or c.get("structure")):
                if q["acc"] not in af_cache:
                    if a.offline:
                        raise SystemExit(f"no AlphaFold answer cached for {q['acc']} and "
                                         "--offline was given")
                    af_cache[q["acc"]] = alphafold_model(q["acc"])
                q["alphafold"] = af_cache[q["acc"]]
            q["input_per_tool"] = {t: tool_input(c, t, q) for t in TOOL_SETTINGS}
        for t in targets:
            s = seq_of(t["acc"])
            t["length"] = s["length"] if s.get("length") else len(s["sequence"])
            t["protein_name"] = t.get("protein_name") or s.get("protein_name", "")
            t.setdefault("label", s.get("entry_name", t["acc"]).split("_")[0])

        # The residues, aligned. Where no tool in this page's run has an alignment, the block is
        # the best local alignment of the two sequences with BLOSUM62, computed here.
        q0, t0 = queries[0], targets[0]
        qs, ts = seqs[q0["acc"]]["sequence"], seqs[t0["acc"]]["sequence"]
        if c.get("query_region") and c.get("target_region"):
            qa, qb = c["query_region"]
            ta, tb = c["target_region"]
            # Set the two motifs under each other, with as much flank as both proteins have.
            flank = 6
            lf = min(flank, qa - 1, ta - 1)
            rf = min(flank, len(qs) - qb, len(ts) - tb, )
            n_motif = min(qb - qa, tb - ta) + 1
            qlo, tlo = qa - lf, ta - lf
            asub = qs[qlo - 1:qa - 1 + n_motif + rf]
            bsub = ts[tlo - 1:ta - 1 + n_motif + rf]
            line, ident = match_line(asub, bsub)
            aln = {"source": "the motif and six residues on each side, set side by side",
                   "query_start": qlo, "target_start": tlo, "query": asub, "target": bsub,
                   "match": line, "identical": ident, "columns": len(asub),
                   "motif_query": [qa, qb], "motif_target": [ta, tb],
                   "motif_columns": [lf + 1, lf + n_motif],
                   "identical_in_motif": line[lf:lf + n_motif].count("|")}
        else:
            qa_, tb_, qstart, tstart, score = smith_waterman(qs, ts, mat)
            line, ident = match_line(qa_, tb_)
            aln = {"source": "best local alignment of the two sequences, BLOSUM62, "
                             "gap open 11, gap extend 1, computed by scripts/known_cases.py",
                   "query_start": qstart, "target_start": tstart, "query": qa_, "target": tb_,
                   "match": line, "identical": ident, "columns": len(qa_), "score": score}
        row["alignment"] = aln
        row["queries"] = queries
        row["targets"] = targets
        not_yet = f"the search of this query against {TARGET_DB} has not been run yet"
        no_model = [q["acc"] for q in queries
                    if q["input_per_tool"]["foldseek"]["kind"] == "no structure"]
        row["tools"] = {t: {"setting": TOOL_SETTINGS[t], "result": "not run", "reason": not_yet,
                            "score": None, "decoy_score": None, "call": None, "iou": None}
                        for t in TOOL_SETTINGS}
        if no_model:
            row["tools"]["foldseek"]["reason"] = (
                f"{not_yet}; and {len(no_model)} of {len(queries)} query proteins have neither "
                f"an experimental structure in the papers cited here nor an AlphaFold Database "
                f"model ({', '.join(no_model)}), so Foldseek cannot search them at all")
        cases.append(row)

    doc = {
        "built_by": "scripts/known_cases.py",
        "target_db": TARGET_DB,
        "direction": "each case's query protein is searched against Swiss-Prot. Table 1 searches "
                     "the other way: human proteins against one relative's proteome.",
        "shuffle_seed": SHUFFLE_SEED,
        "boxes": BOXES,
        "result_words": RESULT_WORDS,
        "eco": ECO,
        "references": REFERENCES,
        "tool_settings": TOOL_SETTINGS,
        "cases": cases,
    }
    (out / "cases.json").write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
    seq_f.write_text(json.dumps(seqs, indent=1) + "\n")
    af_f.write_text(json.dumps(af_cache, indent=1) + "\n")

    # The queries the run searches: every query protein and its decoy, in case order. The name
    # before the first space is what the tools report back, and scripts/known_cases_results.py
    # reads the case number and the decoy mark off it.
    fa = []
    for c in cases:
        for q in c["queries"]:
            s = seqs[q["acc"]]["sequence"]
            fa.append(f'>case{c["id"]}|{q["acc"]}|real {q["label"]}\n' + wrap60(s))
            fa.append(f'>case{c["id"]}|{q["acc"]}|decoy {q["label"]}, dipeptide-preserving '
                      f'shuffle, seed {SHUFFLE_SEED}\n' + wrap60(q["decoy"]["sequence"]))
    (out / "queries.fasta").write_text("".join(fa))
    n_q = sum(len(c["queries"]) for c in cases)
    print(f"wrote {out/'cases.json'}: {len(cases)} cases, {n_q} query proteins, "
          f"{len(seqs)} sequences cached")


if __name__ == "__main__":
    sys.exit(main())
