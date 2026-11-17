# Contributing

The leaderboard grows in two ways, and neither needs anyone here to run your tool or your
experiment:

- **Tool authors** send predictions as a pull request. A workflow grades them with the same
  scorer that made every number on the site and posts the rows on the pull request.
- **Labs** send a protein region they have tested, through the form at
  [/contribute/truth/](https://seanome.github.io/protein-function-annotation-leaderboard/contribute/truth/).
  A curator accepts or rejects it with a reason ([docs/curation.md](docs/curation.md)), and
  accepted rows become a test every tool is graded on.

No account on the site and nothing to install: a GitHub account is enough.

## Tool authors: five steps

### 1. Pick a query set and get its files

A query set is a fixed list of query proteins with known answers. Every submission on one
query set is graded on the same proteins against the same answers. They are in
[`data/query_sets/`](data/query_sets/). Today there is one:

| folder | queries | what they are |
|---|---|---|
| `human_midi_plus` | 998 human proteins | chromosome 6 plus partners of the MHC (the major histocompatibility complex) on other chromosomes: B2M, CD1A-E, MR1, KIR and LILR genes |

Each folder holds `queries.fasta` (the proteins to search with), `relatives.tsv` (the proteomes
a search tool may search, defined below), the answers (`truth/`), and `SOURCE.json` (where every
file came from, with its sha256). There are two sets of answers, called truth sets: the Pfam-A
domains on each query (Pfam is the database of protein domain families) and the Swiss-Prot features
(curated sites and regions, such as active sites and membrane-spanning segments, from the reviewed
part of UniProt).

A **relative** is another species whose annotated proteome your tool searches. When your tool
aligns part of a human protein to part of a relative's protein, the scorer copies the known
annotation of that part of the relative's protein (its Pfam domains, its Swiss-Prot features)
onto the human region. The call is right when the copied family or feature is there on
the human protein, in the right place. `relatives.tsv` lists 15 relatives, from chimpanzee
(6 million years ago, Mya, since its line split from the human line) to *Arabidopsis* (1,516 Mya),
and *E. coli*, for which TimeTree 5 gives no age (`mya_to_human` is empty). The ages are TimeTree 5's
pairwise divergence times from human. All are from the Quest for Orthologs (QfO) reference proteomes, release 2020_04. Get the proteome file the table names
and check its sha256; the answers were made for that file and no other:

```bash
curl -fL -o qfo.tar.gz https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2020_04_with_updated_UP000008143/QfO_release_2020_04_with_updated_UP000008143.tar.gz
tar -xzf qfo.tar.gz Eukaryota/UP000000437_7955.fasta
shasum -a 256 Eukaryota/UP000000437_7955.fasta
```

The tarball is 2.6 GB. The last command must print the sha256 in the `zebrafish` row of
`relatives.tsv` (`bb03b87a…`).

A tool that names the family itself, with no relative (hmmscan against the Pfam models, for
example), uses `none` as the relative and puts the family in `label`.

### 2. Write the predictions

`predictions.tsv`, or `predictions.tsv.gz` if it is large: tab-separated, this header, one row
per call. Positions count the first residue as 1 and include both ends, as BLAST and MMseqs2
report them. (The scorer measures overlap as end minus start for every tool alike; a one-residue
call, start equal to end, is widened to the next residue so it can match a one-residue answer.)

| column | what goes in it |
|---|---|
| `query_set` | the folder name in `data/query_sets/`, such as `human_midi_plus` |
| `relative` | a `relative` from `relatives.tsv`, or `none` |
| `query` | the query protein as `queries.fasta` names it (`sp\|O00624\|NPT3_HUMAN`) or its accession (`O00624`) |
| `call_start`, `call_end` | the called region on the query |
| `target` | the relative's protein the region aligned to, as its FASTA file names it (`tr\|E7F467\|E7F467_DANRE`) or its accession (`E7F467`); empty when `relative` is `none` |
| `target_start`, `target_end` | the aligned region on the target; empty when `relative` is `none` |
| `label` | empty for a search call; for `relative` = `none`, a Pfam accession (`PF00001`; a version suffix such as `PF00001.24` is dropped) or one of the 12 Swiss-Prot feature types the answers use: `ACT_SITE`, `BINDING`, `COILED`, `DNA_BIND`, `DOMAIN`, `INTRAMEM`, `MOTIF`, `REGION`, `REPEAT`, `SITE`, `TRANSMEM`, `ZN_FING` |
| `score` | your tool's score; higher means more confident. Only the order matters |
| `tool` | the same short name as `tool` in `tool.yaml` |
| `version` | the same version as `version` in `tool.yaml` |
| `setting` | a short name for the parameter set, such as `s7`. One file may hold several settings |

Keep every call your tool makes, down to its loosest cut-off. The scorer sweeps the score to
find the best balance of precision and recall (Fmax), so it needs the weak calls too.

### 3. Write the timing log

`timing.tsv`: one row per run, this header. These numbers are shown on the site as
**reported**, with your hardware named, never as measured by the leaderboard.

| column | what goes in it |
|---|---|
| `query_set`, `relative`, `setting` | which run |
| `step` | `search`; or `database_build` (once per relative) or `query_conversion` (such as predicting structures for the queries) |
| `wall_s` | wall-clock seconds |
| `cores` | CPU cores the run used |
| `peak_ram_gb` | peak memory in GB |
| `gpu_model` | the GPU, such as `NVIDIA A100 80GB`, or `none` |
| `hardware` | the machine: CPU model, and anything else that matters |

Every (query set, relative, setting) in the predictions needs a `search` row, `none` included.
The comment adds up the runs of one setting and step across relatives; `cores` shows a number
only when every run used the same count.

### 4. Write the tool card

`tool.yaml` becomes the tool's card on the site:

```yaml
tool: mytool            # lower case, letters, digits and _
name: MyTool
version: "1.2.0"       # in quotes: unquoted, 1.10 would be read as 1.1
kind: letters           # letters, profiles, shapes, kmers or models
compares: >-
  One or two sentences: what the tool compares between two proteins.
needs:
  sequence_only: true
  structure_prediction: false
  gpu: false
licence: MIT
licence_url: https://...
container: docker.io/me/mytool@sha256:...   # pinned by digest, the image's sha256, so the exact image can be pulled again
command: |
  the exact command that made predictions.tsv
```

`kind` is the group the site files the tool under: `letters` (aligns amino-acid letters, such
as MMseqs2), `profiles` (builds a family profile, such as jackhmmer), `shapes` (compares 3D
shape, such as Foldseek), `kmers` (counts shared short words, such as kmerseek), `models`
(a trained model, such as ProstT5).

### 5. Open a pull request

Add one folder, `submissions/<name>/` (the name in lower-case letters, digits, `_` and `-`), holding only `predictions.tsv` (or `.tsv.gz`),
`timing.tsv` and `tool.yaml`, and open a pull request to `main`. A file larger than 100 MB
cannot go into GitHub: split the work into one folder per relative.

The `score-submission` workflow then

1. checks every file and every row, and lists any problem in a comment;
2. grades the calls with `scorer/evaluate_domain_calls.py`, a byte-for-byte copy of the
   leaderboard's scorer (`scorer/SOURCE.json` names its commit), once per relative and truth
   set, with and without merging overlapping calls of one family;
3. writes the rows in the layout of the leaderboard's own tables, with every cost number
   marked `provenance = reported`;
4. posts the scores, the costs and the tool card as one comment on the pull request, and
   updates that comment on every push.

The comment reports Fmax (the best balance of precision and recall over all score cut-offs; 1 is
perfect), AUPRC (area under the precision-recall curve) and median IoU (for correct calls, the
overlap of call and answer divided by their union; 1 is an exact match). Raw keeps every call;
dedup merges overlapping calls that copy the same family onto the same protein. Cost numbers
carry `provenance = reported`: they are your numbers, not measured by the leaderboard.

The full rows are in the run's `scored-rows` artifact (a zip of files attached to the workflow run, on its GitHub Actions page). A maintainer merges once the comment is
there and the tool card reads right. Your folder stays in the repository and is scored again
at every release, so it is also graded on rows labs add later. The leaderboard table does not
show submitted rows yet; until it does, they are on the pull request.

The workflow runs the scorer from the branch your pull request goes into, never code from the pull request, so a pull request
from a fork is scored the same way.

## Worked example: MMseqs2 against zebrafish

This repeats the leaderboard's own MMseqs2 search, so its scores must match the site's.
It needs `mmseqs` (version 18.8cc5c, the one on the leaderboard), `python3`, `bc`,
`/usr/bin/time`, and the zebrafish proteome from step 1.

```bash
git clone https://github.com/seanome/protein-function-annotation-leaderboard.git
cd protein-function-annotation-leaderboard
examples/mmseqs2/run.sh zebrafish Eukaryota/UP000000437_7955.fasta submissions/mmseqs2-zebrafish 8
```

[`examples/mmseqs2/run.sh`](examples/mmseqs2/run.sh) runs

```bash
mmseqs createdb Eukaryota/UP000000437_7955.fasta target
mmseqs createdb data/query_sets/human_midi_plus/queries.fasta query
mmseqs search query target result tmp --threads 8 -s 7 --max-seqs 1000 -e 10.0
mmseqs convertalis query target result hits.tsv \
    --format-output "query,target,qstart,qend,tstart,tend,bits,evalue"
```

and turns each hit into one row: `qstart`/`qend` become the call, `target`/`tstart`/`tend`
the aligned target region, and the bit score the score. It times both steps and fills in
`tool.yaml` from [`examples/mmseqs2/tool.yaml`](examples/mmseqs2/tool.yaml). On an Apple M4
Max with 8 threads it wrote 94,835 calls; the search took 13.6 s and peaked at 1.1 GB.
The first lines of `predictions.tsv.gz`:

```text
query_set        relative   query       call_start  call_end  target  target_start  target_end  label  score  tool            version   setting
human_midi_plus  zebrafish  A0A1B0GTQ4  33          64        P0DP88  43            74                 29     mmseqs2_seqseq  18.8cc5c  s7
human_midi_plus  zebrafish  O00624      11          424       E7F467  35            435                294    mmseqs2_seqseq  18.8cc5c  s7
```

(The file is tab-separated; the spaces here are for reading.) If you want to see the scores
before opening the pull request, score it the way the workflow will. This step is optional
and needs Python with the packages in `requirements.txt` (`python3 -m venv .venv &&
.venv/bin/pip install -r requirements.txt`):

```bash
make score SUB=submissions/mmseqs2-zebrafish
```

`scored/mmseqs2-zebrafish/comment.md` is the comment the workflow will post. For this example
its first table reads, with the setting, region counts and leaderboard columns left out:

| relative | truth set | Fmax raw | Fmax dedup | AUPRC raw | median IoU raw |
|---|---|---:|---:|---:|---:|
| zebrafish | Pfam-A domain instances (2,435) | 0.327 | 0.455 | 0.023 | 0.741 |
| zebrafish | Swiss-Prot features (7,185) | 0.156 | 0.156 | 0.017 | 0.744 |

the same as the leaderboard's MMseqs2 row for zebrafish. `tests/check_scorer.py --scored
scored/mmseqs2-zebrafish` compares every row with the leaderboard's; on this example all 532
rows the two share agree to within 1e-9.

Then `git add submissions/mmseqs2-zebrafish`, commit, push to your fork, and open the pull
request.

## Labs

Use the form at
[/contribute/truth/](https://seanome.github.io/protein-function-annotation-leaderboard/contribute/truth/),
or start a "Add a region your lab has tested" issue on GitHub. One region per issue.
[docs/curation.md](docs/curation.md) has the rules a curator applies.
