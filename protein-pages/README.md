# Protein pages: annotation tools side by side

`build.py` turns the output files of up to ten annotation tools into a small website with one
page per protein of your proteome. On each page, every tool's hits are drawn as boxes on one
ruler of residue positions, so you can see which tools found something, where they agree,
where they disagree, and where only one tool says anything. Click a box to see the aligned
residues behind it.

The tools are the ones on the
[protein function annotation leaderboard](https://seanome.github.io/protein-function-annotation-leaderboard/),
plus InterProScan.

## What a protein page shows

From top to bottom:

1. The protein ID, its length, and one sentence made by a fixed rule:
   "N of 8 programs report a hit", then whether they name the same Pfam family.
   There are 10 rows but 8 programs: MMseqs2 run once and run for three rounds is one
   program, and Foldseek on predicted structures and ProstT5 + Foldseek is another. Two runs
   of the same program are not two independent opinions, so they count once. (On these pages
   "family" always means a Pfam protein family.)
   Two programs **agree** when a hit from each covers the same part of the protein (they
   share at least half of the shorter hit) and both name the same Pfam family or InterPro
   entry. For a search hit, the family named is the Pfam domain on the matched part of the
   Swiss-Prot protein, read from InterProScan's output on that protein, or hmmscan's if
   InterProScan has none. When hits overlap but no Pfam family is known for the Swiss-Prot
   proteins they match, the sentence says that instead of "They do not name the same Pfam
   family". With one program there is nothing to agree on, and the sentence stops after the count.
2. The drawing of what each tool found, the first thing under the sentence, with the legend
   in a narrow column beside it (below it on a phone):
   - a ruler of residue positions;
   - one row per tool, starting with the tool's name and its mark. A row says "not run" when
     that tool's file is not in the input folder. The best 5 hits per tool are drawn; "Show all
     hits" draws the rest. Only hits at the threshold are drawn: E-value at most 10, the
     `-E 10` the searches were run with (hhsearch:
     also Prob at least 50);
   - if `plddt.tsv` is given, a strip coloured by AlphaFold's confidence scale;
   - "Programs with a hit": at each residue, how many of the 8 programs have a hit over it.

   Hovering over or tapping a box shows the matched Swiss-Prot protein's entry name, full name,
   gene and species (from the `swissprot.fasta` header), the Pfam families on the matched part
   of it, and the tool's scores. For hmmscan and InterProScan boxes, it shows the family's
   name, accession and description. Clicking a box shows the aligned residues.
   The legend gives, for each tool, its mark, the score that decides a hit, the setting, the
   database searched and which columns give the box, then the other marks on the page. Every
   colour and shape on the page is in it.
3. If `expression.tsv` is given, one bar per sample.
4. A footer listing every input file read, with its size and date.

The index page has five counts (all proteins, no hit, exactly one program, programs agree on
a Pfam family, programs disagree), filters, a search box, and a sortable table of all proteins.
A narrow column beside the table counts, for each program, the proteins on which only it has a
hit; clicking a program lists those proteins.

The search box looks through every hit counted on the protein's page: the names, genes and
species of the Swiss-Prot proteins the tools matched, the Pfam families on the matched part of
them, and the Pfam families hmmscan or InterProScan found on the protein itself. So "HLA",
"MHC" or "immune" finds the proteins whose hits carry those words. Every word typed must start
a word in that text ("ring" finds "RING finger", not "during"). A last column shows what
matched and which programs' hits carry it. The filter and the search words are part of the
page address (`index.html#only-kmerseek&q=MHC`), so a filtered table can be shared as a link.
The search text is in `site/search.js`, one copy of each distinct text for the whole site.

## Inputs

Put all files in one folder. Only the first one is required. A missing tool file is drawn as
"not run", never guessed.

| File | What it is |
|---|---|
| `proteins.fasta` | Your predicted proteome. The protein ID is the header up to the first space. |
| `swissprot.fasta` | The Swiss-Prot sequences the tools searched, with UniProt's whole headers (`sp\|P04439\|HLAA_HUMAN HLA class I histocompatibility antigen, A alpha chain OS=Homo sapiens OX=9606 GN=HLA-A ...`). Used to show the target's residues, and its name, gene and species in the tooltip and the search. |
| `mmseqs2.m8` | MMseqs2, `--format-output "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,theader"`, optionally followed by `qaln,taln`: then the residue panel shows MMseqs2's own alignment |
| `mmseqs2_3rounds.m8` | MMseqs2 with `--num-iterations 3`, same columns |
| `kmerseek.tsv` | kmerseek regions, tab-separated. Needs query (`query` or `query_name`), target (`target` or `target_name`), query start and end (`query_start`/`query_end` or `region_start`/`region_end`), `target_start`, `target_end`, `region_evalue`. A column named like `go_terms` or `pfam2go` is shown in the tooltip. Coordinates are read as kmerseek writes them: 0-based, end not included. |
| `foldseek.m8` | Foldseek on your predicted structures against AlphaFold Swiss-Prot, columns as MMseqs2 plus `alntmscore` before `theader` |
| `prostt5_foldseek.m8` | Foldseek with 3Di letters predicted by ProstT5, same columns |
| `phmmer.domtblout` | phmmer `--domtblout` against Swiss-Prot |
| `jackhmmer.domtblout` | jackhmmer `-N 3 --domtblout`. HMMER 3.4 writes only the last round to this table; if a file holds more than one, the last is drawn. |
| `hhsearch/` | One `.hhr` file per protein |
| `phmmer.out`, `jackhmmer.out` | Optional: HMMER's main output (`-o`, best with `--notextw`). The residue panel then shows the tool's own alignment. |
| `hmmscan.domtblout` | hmmscan `--domtblout` against Pfam. Rows for Swiss-Prot proteins (not in `proteins.fasta`) give the Pfam families of search targets: their domains with i-Evalue at most 0.001. |
| `hmmscan.out` | Optional: hmmscan's main output on your proteins, for the alignment of each Pfam model to your protein. |
| `interproscan.tsv` | InterProScan TSV (15 columns). Rows for Swiss-Prot proteins are used the same way. |
| `protein_names.tsv` | Optional: protein ID, then a name (for example the Botryllus histocompatibility factor). Shown on the page and searched. |
| `plddt.tsv` | Protein ID, then comma-separated pLDDT per residue |
| `expression.tsv` | Protein ID, then one column per sample |
| `leaderboard.tsv` | Columns `tool`, `fmax`, `cost_per_1000_queries`, `date_read`, copied from the leaderboard page by hand. Shown in the legend. |
| `run_settings.tsv` | Columns `track`, `setting`, `database`. Replaces the legend's setting line when your run differs from the leaderboard's. |
| `kmerseek_alphabets.rs` | kmerseek's `src/rust/alphabets.rs`. The residue classes of each alphabet are read from it, to mark same-class residues in the kmerseek panel. |

### Which columns are positions on your protein

For hmmscan and InterProScan, the box is the envelope (`env from`, `env to`) on your protein.
For phmmer and jackhmmer, your protein is turned into the search model, so `env` and `ali`
are positions on the Swiss-Prot protein and your protein's positions are the `hmm from` and
`hmm to` columns. The box uses `hmm from`-`hmm to`; the tooltip gives the target's
`env` range.

## How to build

Python 3.9 or newer, and Biopython 1.79 or newer:

```bash
pip install biopython
```

```bash
python build.py path/to/inputs --out site --jobs 8
```

Open `site/index.html` in a browser. The pages need no server.

## Speed, memory and size

Measured on Sherlock (8 CPUs of an Intel Xeon E5-2620 v4, Python 3.9, Biopython 1.79), on every
Botryllus FUN protein that has search output: 7,216 proteins and all 575,748 Swiss-Prot
sequences. hmmscan output covers only the test proteins and their targets in these runs.

| Hits kept | Hits | Read all files | Whole build | Main process | Whole job (Slurm) | Site |
|---|---|---|---|---|---|---|
| E-value at most 0.001 (job 49809991) | 16.8 million | 73 s | 3 min 7 s | 6.2 GB | 19.8 GB | 2.0 GB |
| E-value below 1 (job 49843507) | 20.3 million | 97 s | 5 min 12 s | 10.4 GB | 29.3 GB | 8.1 GB |

The build at E-value 10, the site's threshold now, has not been timed on the full set yet (see
the pull request). Scaled to 20,000 proteins with as many hits per protein, the E < 1 build
takes about 14 minutes and about 80 GB of memory: more than a laptop has. These 7,216 proteins
are a hit-rich set (notebook 281 picked them), so a whole proteome, a third of it with no
hit, needs less. The build was not run on a laptop: the input files (4.5 GB) stay on Sherlock.

Two changes made it this fast. The rule that decides whether programs agree stops at the first
overlapping pair instead of comparing every pair of hits (that took 117 of 119 seconds on the
50-protein set before). Python's cycle collector is turned off while reading, and each hit is
kept as its raw text line until a worker draws its page; the first full run, which kept every
hit as a record, was killed at 32 GB.

## How the residue panel is made

Clicking a box opens a panel with:

- the tool's scores (E-value, score or bits), how many aligned columns are identical and how
  many more are similar (BLOSUM62 score above 0), the number of gap columns, and the ranges on
  both proteins with how much of each they cover;
- the matched Swiss-Prot protein's full name, gene, species and Pfam families;
- a small map: your protein across, the matched protein down, every hit on that protein as a
  line in its tool's colour (the one shown drawn thick), and Pfam domains along both edges;
- the residues, 60 columns a line, coloured by the residue classes of the kmerseek alphabet in
  `kmerseek_alphabets.rs` (orange: hydrophobic, blue: polar in hp_lehninger2), identical
  residues in bold, and a match line (| the same residue, + similar). Pointing at a residue
  gives its position and class.

Which alignment is shown:

- MMseqs2 and Foldseek: the `qaln` and `taln` columns, when the `.m8` file has them.
- phmmer and jackhmmer: the alignment in `phmmer.out` / `jackhmmer.out`. The search model was
  built from your protein, so its row is written with your protein's residues at the model's
  positions. hmmscan: the alignment in `hmmscan.out`, with the Pfam model's consensus as the
  other row.
- hhsearch: the alignment printed in the `.hhr` file.
- kmerseek: the two segments side by side without gaps, as kmerseek matches them, with 10
  faded residues on each side and the class letters of both between them.
- When a tool's own alignment is not given: the two ranges aligned again (local, BLOSUM62, gap
  open -11, gap extend -1), with Biopython at build time for the best 5 hits per tool and in the
  browser for the rest. The panel says which.

Alignments are stored next to each page (`proteins/<ID>.aln.js`) and read on the first click.

## How to add a tool

1. Write a parser in `build.py` that returns `{protein_id: [hit, ...]}`, using `hit(...)` for
   each hit: the query range (1-based, inclusive), the target accession and range, the
   E-value, the scores for the tooltip as `(name, value)` pairs, and whether it passes the
   threshold.
2. Add a line to `TRACKS` with its file name, program, name, setting and database, and call the
   parser from `main()`.
3. If it is a new program (not a second run of one already there), add it to `PROGRAMS` with a
   colour no other mark on the page uses, and change "8 programs" in the text.

## Test set

`test_inputs/` holds the scripts that cut a 51-protein test set out of real *Botryllus
schlosseri* outputs on Sherlock: 50 drawn at random plus the Botryllus histocompatibility
factor (BHF, `FUN008084_FUN008084`). See `test_inputs/PROVENANCE.md`.

## Checks

```bash
python -m pytest tests/
```

```bash
python tests/check_site.py site
```

```bash
python tests/check_boxes.py site path/to/inputs 3
```

`check_site.py` checks that every colour on a page is one the legend explains, that the five
index numbers add up, and that each page's sentence matches its index row. `check_boxes.py`
picks boxes at random and checks them against the line of the tool's file they came from.
`tests/screenshots.py` takes full-page screenshots with Playwright.
