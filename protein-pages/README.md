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
   "N of 8 tool families report a hit", then whether they agree on a Pfam family.
   A **tool family** is one program. MMseqs2 run once and run for three rounds is one family;
   Foldseek on predicted structures and ProstT5 + Foldseek is another. Two runs of the same
   program are not two independent opinions, so they count once.
   Two families **agree** when a hit from each covers the same part of the protein (they
   share at least half of the shorter hit) and both name the same Pfam family or InterPro
   entry. For a search hit, the family named is the Pfam domain on the matched part of the
   Swiss-Prot protein, read from InterProScan's output on that protein, or hmmscan's if
   InterProScan has none. When hits overlap but no Pfam family is known for the Swiss-Prot
   proteins they match, the sentence says that instead of "They do not agree on a family".
   With one family there is nothing to agree on, and the sentence stops after the count.
2. A legend with one row per tool: the mark, the setting, the database searched, the score
   shown and the threshold. Every colour and shape on the page is in it.
3. A ruler of residue positions.
4. One row per tool. A row says "not run" when that tool's file is not in the input folder.
   The best 5 hits per tool are drawn; "Show all hits" draws the rest. Only hits at the
   threshold are drawn: E-value at most 0.001 (hhsearch: also Prob at least 50).
5. If `plddt.tsv` is given, a strip coloured by AlphaFold's confidence scale.
6. "Families with a hit": at each residue, how many of the 8 families have a hit over it.
7. If `expression.tsv` is given, one bar per sample.
8. A footer listing every input file read, with its size and date.

The index page has a sortable table of all proteins, five counts (all proteins, no hit,
exactly one family, families agree, families disagree), and filters such as
"only proteins where kmerseek is the only family with a hit". The filter is part of the
page address (`index.html#only-kmerseek`), so a filtered table can be shared as a link.

## Inputs

Put all files in one folder. Only the first one is required. A missing tool file is drawn as
"not run", never guessed.

| File | What it is |
|---|---|
| `proteins.fasta` | Your predicted proteome. The protein ID is the header up to the first space. |
| `swissprot.fasta` | The Swiss-Prot sequences the tools searched (headers `sp\|P12345\|NAME_HUMAN ...`). Used only to show the target's residues. |
| `mmseqs2.m8` | MMseqs2, `--format-output "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,theader"` |
| `mmseqs2_3rounds.m8` | MMseqs2 with `--num-iterations 3`, same columns |
| `kmerseek.tsv` | kmerseek regions, tab-separated. Needs query (`query` or `query_name`), target (`target` or `target_name`), query start and end (`query_start`/`query_end` or `region_start`/`region_end`), `target_start`, `target_end`, `region_evalue`. A column named like `go_terms` or `pfam2go` is shown in the tooltip. Coordinates are read as kmerseek writes them: 0-based, end not included. |
| `foldseek.m8` | Foldseek on your predicted structures against AlphaFold Swiss-Prot, columns as MMseqs2 plus `alntmscore` before `theader` |
| `prostt5_foldseek.m8` | Foldseek with 3Di letters predicted by ProstT5, same columns |
| `phmmer.domtblout` | phmmer `--domtblout` against Swiss-Prot |
| `jackhmmer.domtblout` | jackhmmer `-N 3 --domtblout`. HMMER 3.4 writes only the last round to this table; if a file holds more than one, the last is drawn. |
| `hhsearch/` | One `.hhr` file per protein |
| `hmmscan.domtblout` | hmmscan `--domtblout` against Pfam. Rows for Swiss-Prot proteins (not in `proteins.fasta`) give the Pfam domains of search targets. |
| `interproscan.tsv` | InterProScan TSV (15 columns). Rows for Swiss-Prot proteins are used the same way. |
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

```bash
pip install biopython
```

```bash
python build.py path/to/inputs --out site --jobs 8
```

Open `site/index.html` in a browser. The pages need no server.

## How the residue panel is made

- MMseqs2, Foldseek, phmmer, jackhmmer: the two coordinate ranges are aligned again with
  Biopython's `PairwiseAligner` (local, BLOSUM62, gap open -11, gap extend -1). The panel says
  so; the tool's own alignment is not shown. For the best 5 hits per tool this is done when
  the site is built. For the rest, the page does the same alignment in the browser when the
  box is clicked, because doing it for every hit of 20,000 proteins takes too long.
- hhsearch: the alignment printed in the `.hhr` file.
- kmerseek: the two segments side by side without gaps, as kmerseek matches them, with a mark
  where two different residues fall in the same class of the alphabet, and the same segments
  written as class letters.

## How to add a tool

1. Write a parser in `build.py` that returns `{protein_id: [hit, ...]}`, using `hit(...)` for
   each hit: the query range (1-based, inclusive), the target accession and range, the
   E-value, the scores for the tooltip as `(name, value)` pairs, and whether it passes the
   threshold.
2. Add a line to `TRACKS` with its file name, family, name, setting and database, and call the
   parser from `main()`.
3. If it is a new program (not a second run of one already there), add it to `FAMILIES` with a
   colour no other mark on the page uses, and change "8 tool families" in the text.

## Test set

`test_inputs/` holds the scripts that cut a 50-protein test set out of real *Botryllus
schlosseri* outputs on Sherlock. See `test_inputs/PROVENANCE.md`.
