# Protein pages: annotation tools side by side

`build.py` turns the output files of up to ten annotation tools into a small website with one
page per protein of your proteome. On each page, every tool's hits are drawn as boxes on one
ruler of residue positions, so you can see which tools found something, where they agree,
where they disagree, and where only one tool says anything. Click a box to see the aligned
residues behind it.

The tools are the ones on the
[protein function annotation leaderboard](https://seanome.github.io/protein-function-annotation-leaderboard/),
plus InterProScan.

## Two E-value thresholds

- A hit is **drawn** when its E-value is at most 10, the `-E 10` the searches were run with.
- A hit **counts** when its E-value is at most 0.001 (`--count-at`, a power of 10). Counted
  hits make the sentence at the top of a protein page, the "programs with a hit" line, the
  index page's seven numbers, the table's program count, Pfam family, small line and gene
  column, and `summary.json`.
- Two things use every drawn hit instead: the table's "Best E-value, any tool" column, with
  the index slider that filters on it, and the search box. So a protein with no counted hit
  can still show a best E-value, such as BHF's 1.2e-03.
- hhsearch hits with Prob below 50 are not drawn at all. An InterProScan match from a member
  database whose score is not an E-value (PRINTS, PROSITE and others) is drawn and counts at
  any threshold, unless it is a coil, disorder, membrane or signal-peptide prediction.

Each protein page has an E-value slider that starts at 0.001. Moving it up draws weaker hits,
up to 10; moving it down draws fewer. The sentence and the index do not move with it.

## What a protein page shows

From top to bottom:

1. The protein ID, its length, and one sentence made by a fixed rule, for example
   "4 of the 5 programs run report a hit. They agree on Rep_fac_C (PF08542)."
   The number of programs run is the number with at least one input file; when all 8 were
   run the sentence says "of 8 programs". There are 10 rows but 8 programs: MMseqs2 run once
   and run for three rounds is one program, and Foldseek on predicted structures and ProstT5 +
   Foldseek is another. Two runs of one program are not two opinions, so they count once.
   ("Family" on these pages means a Pfam protein family, or an InterPro entry where
   InterProScan gives one.)
   Hits of two programs **overlap** when they share at least half of the shorter hit. For a
   search hit, the family is the Pfam domain on the matched part of the Swiss-Prot protein,
   read from InterProScan's output on that protein, or hmmscan's if InterProScan has none.
   With two or more programs, the sentence ends in one of four ways:
   - **agree**: "They agree on MHC_I (PF00129)." Overlapping hits name the same family.
   - **disagree**: "Their hits overlap but name different Pfam families: A (phmmer),
     B (hmmscan with Pfam)." Overlapping hits both name a family, and none is shared.
   - **overlap, but no family known**: "Their hits overlap, but for one of them no Pfam
     family is known on the matched part of the Swiss-Prot protein." The Swiss-Prot protein
     may have no Pfam output, or its domains may lie outside the matched part.
   - **apart**: "Their hits are on different parts of the protein: no two programs' hits
     share at least half of the shorter one." They may still touch.
   The first of these that holds is used, in this order: agree, disagree, overlap but no
   family known, apart. So an "agree" protein can also have overlapping hits that name other
   families; the page lists the other families that also meet the rule.
   With one program there is nothing to compare, and the sentence stops after the count.
2. The legend, in a column left of the drawing (above it on a phone), so it is read first.
   One line per tool, its mark and its name; a click opens what makes a hit, the setting, the
   database and which columns give the box.
3. The drawing. Each box is one hit, at the position the tool reported on your protein.
   - a ruler of residue positions, in aa;
   - one row per tool. A row says "not run" when that tool's file is not in the input
     folder. The best 5 hits per tool are drawn; "Show all hits" draws the rest;
   - if `plddt.tsv` is given, a strip coloured by AlphaFold's confidence scale;
   - "Programs with a hit": at each residue, how many of the programs run (0 to M) have a
     counted hit over it.

   Hovering over or tapping a box shows the matched Swiss-Prot protein's entry name, full
   name, gene and species, the Pfam families on the matched part of it, and the tool's
   scores. Clicking a box, or its label, shows the aligned residues.
4. If `expression.tsv` is given, one outlined bar per sample.
5. A footer listing every input file read, with its size and date.

### How the drawing encodes what

| What | How it is drawn |
|---|---|
| Position on your protein | the horizontal position, the only thing on that axis |
| Kind of program | colour: sequence search `#008B7A` (MMseqs2, phmmer, jackhmmer, hhsearch), structure `#6A3D9A` (Foldseek), k-mer `#C0392B` (kmerseek), Pfam-model scan dark grey `#3D3D3D` (hmmscan, InterProScan) |
| Program, inside a kind | box shape: filled, open, hatched or rounded |
| E-value | box opacity, on a log scale: 0.25 at E-value 10, full at 1e-30 or below |
| The protein | a thin light grey line `#C8C8C8` |

Three hues and one grey. In normal vision they differ by at least 7 in lightness (CIE L*).
Under a deuteranopia simulation (Machado et al. 2009) they differ by at least 24 in colour
(CIE76), but the sequence-search teal and the kmerseek red are then nearly the same lightness,
so in greyscale those two are told apart by their rows, not their colour. A label
goes inside its box when it fits; otherwise it goes beside the box, joined to it by a thin
line in the box's colour, and is never cut or made smaller. `assets/site.css` lists every
colour on the pages and its one meaning.

### The index page

Seven numbers: all proteins, then six classes that add up to all: no hit, exactly one
program, agree, disagree, overlap but no family known, apart. Each is also a filter. A
narrow column beside the table counts, for each program, the proteins on which only it has a
hit; clicking a program lists those proteins. The table's small line is each protein's
"programs with a hit" line at the same 90-pixel width and the same 0-to-M height for every
protein; hovering over it gives the position and the count, for example "position 120 aa:
3 of 5 programs".

The search box looks through every hit drawn on the protein pages (E-value at most 10): the
names, genes and species of the Swiss-Prot proteins the tools matched, the Pfam families on
the matched part of them, and the Pfam families hmmscan or InterProScan found on the protein
itself. So "HLA", "MHC" or "immune" finds the proteins whose hits carry those words. Every
word typed must start a word in that text ("ring" finds "RING finger", not "during"). A last
column shows what matched and which programs' hits carry it. The filter and the search words
are part of the page address (`index.html#only-kmerseek&q=MHC`), so a filtered table can be
shared as a link.

## Inputs

Put all files in one folder. Only the first one is required. A missing tool file is drawn as
"not run", never guessed. Any file can be gzipped: build.py reads `name.gz` when `name` is
not there.

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

Open `site/index.html` in a browser. The pages need no server. `--count-at 0.01` counts hits
at E-value 0.01 instead of 0.001. build.py replaces `--out` only when that folder is empty or
holds `summary.json` (a site it built); otherwise it stops with exit code 2.

To build the 51-protein test site, which is in this repository:

```bash
python build.py test_inputs/data --out site
```

The built site and the screenshots are not committed. The GitHub Actions workflow
`.github/workflows/protein-pages.yml` builds the test site on every push that changes
`protein-pages/`, runs the tests and both site checks, and keeps the site as a downloadable
artifact of the run (`protein-pages-site`).

## Speed, memory and size

Measured on Sherlock (8 CPUs of an Intel Xeon E5-2620 v4, Python 3.9, Biopython 1.79), on every
Botryllus FUN protein that has search output: 7,216 proteins and all 575,748 Swiss-Prot
sequences. hmmscan output covers only the test proteins and their targets in these runs.

| Hits kept | Hits | Read all files | Whole build | Main process | Whole job (Slurm) | Site |
|---|---|---|---|---|---|---|
| E-value at most 0.001 (job 49809991) | 16.8 million | 73 s | 3 min 7 s | 6.2 GB | 19.8 GB | 2.0 GB |
| E-value below 1 (job 49843507) | 20.3 million | 97 s | 5 min 12 s | 10.4 GB | 29.3 GB | 8.1 GB |

Both runs are older than the two thresholds: the build now reads every hit up to E-value 10,
and that has not been timed on the full set. Scaled to 20,000 proteins with as many hits per protein, the E < 1 build
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

- the tool's scores (E-value, score or bits) and the ranges on both proteins with how much of
  each they cover;
- the matched Swiss-Prot protein's full name, gene, species and Pfam families;
- a small map: your protein across, the matched protein down, both axes named with the
  protein and "aa"; every hit on that protein is a line in the colour of its kind of program
  (the one shown drawn thick), with Pfam domains along both edges;
- one summary line above the residues, counted from the alignment shown: "226 of 312 aligned
  columns identical (72%), 43 more similar (BLOSUM62 score above 0), 0 gap columns" (for
  kmerseek, "more in the same class"); and directly under it, where the alignment comes from;
- the residues, 60 columns a line, with both coordinates, coloured by the residue classes of
  the kmerseek alphabet in `kmerseek_alphabets.rs` (orange: hydrophobic, blue: polar in
  hp_lehninger2), identical residues in bold, and a match line (| the same residue,
  + similar). Pointing at a residue gives its position and class.

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
   E-value, the scores for the tooltip as `(name, value)` pairs, and whether it is drawn
   (E-value at most 10). Whether it counts is decided later, by `counts()`.
2. Add a line to `TRACKS` with its file name, program, name, setting and database, and call the
   parser from `main()`.
3. If it is a new program (not a second run of one already there), add it to `PROGRAMS` with
   its kind (`CATEGORIES`) and a shape no other program of that kind uses, and change "8
   programs" in the text.

## Test set

`test_inputs/data/` holds the 51-protein test set: 50 *Botryllus schlosseri* proteins drawn at
random plus the Botryllus histocompatibility factor (BHF, `FUN008084_FUN008084`), with the
tools' real outputs cut to these proteins. `test_inputs/` also holds the scripts that made it on
Sherlock, and `cut_test_data.py`, which cut the files small enough to commit. See
`test_inputs/PROVENANCE.md`.

## Checks

```bash
python -m pytest tests/
```

```bash
python tests/check_site.py site
```

```bash
python tests/check_boxes.py site test_inputs/data 9 1
```

The tests include building the test site from `test_inputs/data` and running
`check_site.py` on it. `check_site.py` checks that every colour on a page is one the legend
explains, that the six class numbers add up to the total on the index and in `summary.json`,
and that each page's sentence matches its index row ("N of the M programs run", and the
status). `check_boxes.py` picks boxes at random and checks them against the line of the
tool's file they came from, and the alignment's residues against the FASTA files; it exits 1
if any fails. `tests/screenshots.py` takes full-page screenshots with Playwright
(`PROTEIN_ID:kmerseek` opens that protein's first kmerseek box).
