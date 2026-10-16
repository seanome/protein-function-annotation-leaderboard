"""Facts about each tool that do not come from the tables.

Every entry was read from the run itself on 2026-10-15: the command line from the task's
.command.sh on Sherlock (the zebrafish search), the container from the Nextflow trace, and
the container's version tag from quay.io's API by the digest the trace names. Numbers
(scores, cost, strong and weak subsets) are NOT here; build_site.py reads those from the
tables so a card cannot disagree with the leaderboard.
"""

QFO = "Quest for Orthologs reference proteomes, release 2020_04"

# Kind of method -> one colour across the whole site. Tools of one kind share the colour
# and differ by marker shape.
KINDS = {
    "letters": {"label": "aligns amino-acid letters", "color": "#0072B2"},
    "profiles": {"label": "builds a profile of a protein family", "color": "#E69F00"},
    "shapes": {"label": "compares 3-D shape letters (3Di)", "color": "#009E73"},
    "kmers": {"label": "counts shared k-mers in a reduced alphabet", "color": "#CC3311"},
    "models": {"label": "scans curated Pfam family models", "color": "#6B4C9A"},
}

TOOLS = {
    "kmerseek": {
        "name": "kmerseek",
        "kind": "kmers",
        "shape": "circle",
        "compares": (
            "Short words of k letters (k-mers) shared between two proteins. First, each "
            "amino acid is replaced by a letter from a smaller alphabet, for example H for "
            "hydrophobic and P for polar. Two proteins can then match where their residues "
            "differ but their pattern of residue types is the same."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": (
            f"An index of the relative's proteome ({QFO}), built once per alphabet and "
            "k-mer length."
        ),
        "score": (
            "region_enrichment: the number of k-mers the region shares with the target "
            "divided by the number expected by chance (fold enrichment). Higher ranks "
            "first. A region was kept when either test passed: its own score was at least "
            "1.3, or the whole query protein's p-value was at most 0.05."
        ),
        "caveats": [
            "A feature shorter than k residues contains no complete k-mer, so it cannot be "
            "found. Figure 2 above shows how Fmax changes with region "
            "length.",
            "The task script records a kmerseek search that exits with an error as a search "
            "with no hits, so a crash would look like a miss. How often that happened in "
            "this run: not measured.",
            "Setting hp_thomas_dill2 k33 was scored, but no trace row holds its search task, "
            "so its cost is not measured, and the scoring task for it was not found either.",
        ],
        "version": "kmerseek 0.3.1 (printed in the search log)",
        "container": "docker.io/olgabot/kmerseek:2026-08-24-reduced-alphabets",
        "command": """kmerseek search \\
    --alphabet hp_thomas_dill2 \\
    --ksize    24 \\
    --query    UP000005640_9606.fasta \\
    --target   zebrafish.hp_thomas_dill2.k24.lcfalse.kmerseek.rocksdb \\
    --threshold         0.0 \\
    --min-shared-kmers  2 \\
    --max-query-pvalue  0.05 \\
    --min-region-score  1.3""",
        "command_note": "Shown for setting hp_thomas_dill2 k24. Other settings change "
        "--alphabet and --ksize only.",
        "licence": 'MIT',
        "commercial": 'Yes.',
        "licence_tag": 'permissive licence',
        "licence_url": 'https://github.com/seanome/kmerseek',
        "requirement_question": None,
    },
    "hmmer3_phmmer": {
        "name": "phmmer",
        "kind": "letters",
        "shape": "circle",
        "compares": (
            "Amino-acid letters, one query protein against each protein of the relative, "
            "scored with a substitution matrix the way BLAST does but with HMMER's model "
            "of insertions and deletions."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": f"The relative's proteome ({QFO}), as a FASTA file.",
        "score": "Domain bit score: higher ranks first. The E-value cutoff was 10.",
        "caveats": ["None recorded for this run."],
        "version": "HMMER 3.4",
        "container": "quay.io/biocontainers/hmmer:3.4--hb6cb901_4 "
        "(sha256:7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c)",
        "command": """phmmer \\
    --domtblout /dev/stdout \\
    --tblout /dev/null \\
    -o /dev/stderr \\
    --noali \\
    -E 10.0 \\
    --cpu 8 \\
    UP000005640_9606.fasta UP000000437_7955.fasta""",
        "command_note": "Query: 998 human proteins. Target: the zebrafish proteome.",
        "licence": 'BSD 3-clause (HMMER)',
        "commercial": 'Yes.',
        "licence_tag": 'permissive licence',
        "licence_url": 'https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE',
        "requirement_question": None,
    },
    "hmmer3_jackhmmer": {
        "name": "jackhmmer",
        "kind": "profiles",
        "shape": "circle",
        "compares": (
            "A profile of a protein family: it searches with one protein, builds a profile "
            "from the hits, and searches again with that profile, three rounds in all."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": f"The relative's proteome ({QFO}), as a FASTA file.",
        "score": "Domain bit score: higher ranks first. The E-value cutoff was 10.",
        "caveats": [
            "Each round can pull in proteins from a different family, which then pull in "
            "more. How often that happened here: not measured."
        ],
        "version": "HMMER 3.4",
        "container": "quay.io/biocontainers/hmmer:3.4--hb6cb901_4 "
        "(sha256:7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c)",
        "command": """jackhmmer \\
    -N 3 \\
    --domtblout /dev/stdout \\
    --tblout /dev/null \\
    -o /dev/stderr \\
    --noali \\
    -E 10.0 \\
    --cpu 8 \\
    UP000005640_9606.fasta UP000000437_7955.fasta""",
        "command_note": "Query: 998 human proteins. Target: the zebrafish proteome.",
        "licence": 'BSD 3-clause (HMMER)',
        "commercial": 'Yes.',
        "licence_tag": 'permissive licence',
        "licence_url": 'https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE',
        "requirement_question": None,
    },
    "mmseqs2_seqseq": {
        "name": "MMseqs2",
        "kind": "letters",
        "shape": "square",
        "compares": (
            "Amino-acid letters: it first keeps target proteins that share similar short "
            "words with the query, then aligns only those, which makes it much faster than "
            "aligning every pair."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": f"The relative's proteome ({QFO}), as an MMseqs2 database.",
        "score": "Bit score: higher ranks first. The E-value cutoff was 10.",
        "caveats": [
            "At most 1,000 target proteins are kept per query (--max-seqs 1000), so a query "
            "from a very large family can lose hits. How often: not measured."
        ],
        "version": "MMseqs2 18.8cc5c",
        "container": "quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0 "
        "(sha256:3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef)",
        "command": """mmseqs search \\
    human-38fdd87d_mmdb/db zebrafish_mmdb/db result mmseqs_tmp \\
    --threads 16 \\
    -s 7 \\
    --max-seqs 1000 \\
    -e 10.0""",
        "command_note": "Followed by mmseqs convertalis to write the table of hits.",
        "licence": 'MIT',
        "commercial": 'Yes.',
        "licence_tag": 'permissive licence',
        "licence_url": 'https://github.com/soedinglab/MMseqs2',
        "requirement_question": None,
    },
    "mmseqs2_iterative": {
        "name": "MMseqs2, three rounds",
        "kind": "profiles",
        "shape": "square",
        "compares": (
            "A profile of a protein family, built by MMseqs2 from its own hits over three "
            "rounds of search, as jackhmmer does."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": f"The relative's proteome ({QFO}), as an MMseqs2 database.",
        "score": "Bit score: higher ranks first. The E-value cutoff was 10.",
        "caveats": [
            "At most 1,000 target proteins are kept per query (--max-seqs 1000). How often "
            "that cut hits: not measured.",
            "Each round can pull in proteins from a different family. How often: not "
            "measured.",
        ],
        "version": "MMseqs2 18.8cc5c",
        "container": "quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0 "
        "(sha256:3503bfe576d560e550df2872af86a1ad1bcc1c06cfb7caadd3e7a95649f5f0ef)",
        "command": """mmseqs search \\
    human-38fdd87d_mmdb/db zebrafish_mmdb/db result mmseqs_tmp \\
    --threads 16 \\
    -s 7 \\
    --num-iterations 3 \\
    --max-seqs 1000 \\
    -e 10.0""",
        "command_note": "Followed by mmseqs convertalis to write the table of hits.",
        "licence": 'MIT',
        "commercial": 'Yes.',
        "licence_tag": 'permissive licence',
        "licence_url": 'https://github.com/soedinglab/MMseqs2',
        "requirement_question": None,
    },
    "hhblits": {
        "name": "HH-suite (hhsearch)",
        "kind": "profiles",
        "shape": "triangle",
        "compares": (
            "Two profiles against each other. Here each profile was built from a single "
            "protein, not from a family alignment, because no UniRef30 database was given."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": (
            f"Single-protein profiles of the relative's proteome ({QFO}), as an HH-suite "
            "database."
        ),
        "score": "Bit score from hhsearch: higher ranks first.",
        "caveats": [
            "The run names this tool hhblits (its tool id), but it ran hhsearch on profiles "
            "of single proteins. HH-suite's usual use, with family alignments from UniRef30, "
            "was not run, so these numbers say nothing about it.",
        ],
        "version": "HH-suite 3.3.0",
        "container": "quay.io/biocontainers/hhsuite:3.3.0--h503566f_15 "
        "(sha256:4bf9bb5229de18f522a94f4443c19fdcbb0f0cb0e6ea92f5390aa170bcb0a24f)",
        "command": """ffindex_apply \\
    human-38fdd87d_hhdb/db_a3m.ffdata human-38fdd87d_hhdb/db_a3m.ffindex \\
    -d results.ffdata -i results.ffindex \\
    -- hhsearch -i stdin -d zebrafish_hhdb/db -blasttab /dev/stdout -cpu 16 -v 0""",
        "command_note": "One hhsearch per query protein, run through ffindex_apply.",
        "licence": 'GPL-3.0',
        "commercial": 'Yes. If you distribute a changed version, you must release its source under GPL-3.0 too.',
        "licence_tag": 'GPL-3.0',
        "licence_url": 'https://github.com/soedinglab/hh-suite',
        "requirement_question": None,
    },
    "foldseek": {
        "name": "Foldseek",
        "kind": "shapes",
        "shape": "circle",
        "compares": (
            "Predicted 3-D shapes. Each residue of an AlphaFold model gets one of 20 shape "
            "letters (3Di) from its local fold. The shape letters are aligned together with "
            "the amino acids."
        ),
        "needs": {"sequence only": False, "3D structure prediction": True, "GPU": False},
        "reference": (
            f"AlphaFold models of the relative's proteome ({QFO} accessions), as a Foldseek "
            "database."
        ),
        "score": "Bit score: higher ranks first. The E-value cutoff was 10.",
        "caveats": [
            "A protein with no AlphaFold model cannot be searched. Five relatives have no "
            "Foldseek rows at all (chimpanzee, cow, dog, opossum, western clawed frog): "
            "not run.",
            "AlphaFold models of long proteins come in overlapping pieces (F1, F2, ...). "
            "The task script moves piece n's coordinates 200 x (n - 1) residues along the "
            "protein; if a piece started elsewhere, its calls would land in the wrong "
            "place. Checked against the real piece starts: not measured.",
        ],
        "version": "Foldseek 10.941cd33",
        "container": "quay.io/biocontainers/foldseek:10.941cd33--h5021889_1 "
        "(sha256:1156a052f31b2afb85257c02e83a962f559c9752273fe1064ab735f90ac29d1a)",
        "command": """foldseek search \\
    human-38fdd87d_fsdb/db zebrafish_fsdb/db result foldseek_tmp \\
    --threads 16 \\
    -e 10.0 \\
    --max-seqs 1000""",
        "command_note": "Followed by foldseek convertalis to write the table of hits.",
        "licence": 'GPL-3.0',
        "commercial": 'Yes. If you distribute a changed version, you must release its source under GPL-3.0 too.',
        "licence_tag": 'GPL-3.0',
        "licence_url": 'https://github.com/steineggerlab/foldseek',
        "requirement_question": (
            "Do you have a 3D structure prediction for every protein you want to annotate?",
            "If not, Foldseek cannot search the proteins without one.",
        ),
    },
    "prostt5": {
        "name": "ProstT5 + Foldseek",
        "kind": "shapes",
        "shape": "square",
        "compares": (
            "Shape letters (3Di) predicted from the amino-acid sequence by the ProstT5 "
            "protein language model, then searched with Foldseek, so no AlphaFold model "
            "is needed."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": True},
        "reference": (
            f"3Di letters predicted by ProstT5 for the relative's proteome ({QFO}), as a "
            "Foldseek database."
        ),
        "score": "Bit score: higher ranks first. The E-value cutoff was 10.",
        "caveats": [
            "Proteins longer than 6,000 residues are skipped before the 3Di letters are "
            "predicted, and are listed in a skipped-proteins file.",
            "Predicting the 3Di letters ran on one GPU (foldseek createdb --gpu 1). For the "
            "998 queries that step is counted in the GPU column and the cost; for each "
            "relative's proteome it is a one-off database build and is not counted.",
        ],
        "version": "Foldseek 10.941cd33 with ProstT5 weights",
        "container": "quay.io/biocontainers/foldseek:10.941cd33--h5021889_1 "
        "(sha256:1156a052f31b2afb85257c02e83a962f559c9752273fe1064ab735f90ac29d1a)",
        "command": """foldseek createdb filtered.fasta zebrafish_prostt5/db \\
    --prostt5-model weights --gpu 1 --threads 8

foldseek search \\
    human-38fdd87d_prostt5/db zebrafish_prostt5/db result tmp \\
    --threads 16 \\
    -e 10.0 \\
    --max-seqs 1000""",
        "command_note": "The first command (one per proteome, on a GPU) predicts the "
        "3Di letters. The second is the search the cost columns measure.",
        "licence": 'GPL-3.0 (Foldseek); MIT (ProstT5 model weights)',
        "commercial": "Yes. Foldseek's GPL-3.0 applies to changed versions you distribute.",
        "licence_tag": 'GPL-3.0',
        "licence_url": 'https://huggingface.co/Rostlab/ProstT5',
        "requirement_question": (
            "Do you have a GPU to predict the 3Di letters for your proteins?",
            "The run used one GPU for that step; its cost is not counted here.",
        ),
    },
    "hmmscan": {
        "name": "hmmscan with Pfam",
        "kind": "models",
        "shape": "circle",
        "compares": (
            "Each query protein against the curated profile models of every Pfam family, "
            "without searching any relative."
        ),
        "needs": {"sequence only": True, "3D structure prediction": False, "GPU": False},
        "reference": "Pfam-A 38.2 family models (30,134 models).",
        "score": "Domain bit score: higher ranks first. E-value cutoffs 10 (sequence and domain).",
        "caveats": [
            "Pfam domains are defined by these same models, so scoring hmmscan against the "
            "Pfam truth set grades the models against their own output. Read its Pfam row "
            "as the best a search could do, not as a competitor.",
            "It was scored against the Pfam truth set only. Swiss-Prot features and "
            "Pfam-N: not run.",
        ],
        "version": "HMMER 3.4",
        "container": "quay.io/biocontainers/hmmer:3.4--hb6cb901_4 "
        "(sha256:7a2b317b8d2fd3650b4924a8482cddeb940d4a0746c6a1501ff03ac1b7439e0c)",
        "command": """hmmscan \\
    --domtblout /dev/stdout \\
    -o /dev/null \\
    --noali \\
    -E 10.0 --domE 10.0 \\
    --cpu 4 \\
    Pfam-A.hmm UP000005640_9606.fasta""",
        "command_note": "Query: the same 998 human proteins.",
        "licence": 'BSD 3-clause (HMMER); CC0 1.0 (Pfam-A models, per the Pfam 38.2 release notes)',
        "commercial": 'Yes.',
        "licence_tag": 'permissive licence',
        "licence_url": 'https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE',
        "requirement_question": None,
    },
}

# Common names for the relatives, in the order of time since they split from human.
SPECIES = {
    "ptroglodytes": "chimpanzee",
    "rnorvegicus": "rat",
    "clupus": "dog",
    "btaurus": "cow",
    "mouse": "mouse",
    "mdomestica": "opossum",
    "chicken": "chicken",
    "xtropicalis": "western clawed frog",
    "zebrafish": "zebrafish",
    "ciona": "sea squirt (Ciona)",
    "fly": "fruit fly",
    "worm": "nematode (C. elegans)",
    "yeast": "baker's yeast",
    "arabidopsis": "thale cress (Arabidopsis)",
    "ecoli": "E. coli",
    "all": "none: Pfam models (hmmscan only)",
}
