# Curating wet-lab rows

Labs send protein regions they have tested through the form at `/contribute/truth/` or the
"Add a region your lab has tested" issue form. Each arrives as a GitHub issue labelled
`truth-submission`: that list of issues is the curation queue. This page is the rule book a
curator follows.

## What a row holds

One region of one protein: a UniProt accession (or the sequence and organism when there is no
UniProt entry), the first and last residue (counting from 1, both included), what the
experiment showed (binding, activity, localisation or structure), an evidence code from the
Evidence and Conclusion Ontology (ECO), a citation or a note on unpublished work, and the
submitter's ORCID iD.

## Accept a row when all of these hold

1. **The experiment is the submitter's own or is published, and it tested this region.** The
   citation or note names the experiment that shows the region does what the row says:
   a deletion or point mutant that removes the binding, an assay on the isolated region, a
   solved structure of it. A whole-protein result does not place the region.
2. **The ECO code is an experimental one** and fits the experiment described: a code under
   "experimental evidence" (ECO:0000006) in the ECO tree, such as ECO:0000269 (experimental
   evidence used in manual assertion). A code for sequence similarity, prediction or author
   statement without data (for example ECO:0000250, ECO:0000255, ECO:0000303) is not.
3. **The coordinates are on the sequence named.** For an accession, the region lies inside the
   current UniProt sequence; the build checks this against UniProt and refuses a row that
   runs past the end. If the experiment used another isoform or an older sequence, the
   submitter converts the positions first.
4. **What was shown matches the evidence.** "Binding" needs a binding experiment, "activity" an
   activity assay, "localisation" an experiment that moves or tracks the protein, "structure"
   a solved or tested structure.
5. **The ORCID iD is valid**: the form and the build both check its last character.
6. **The submitter ticked both boxes**: the region is released under CC0, and it comes from an
   experiment rather than only from a database.

## Reject a row when any of these holds

- **It only restates a database annotation.** If UniProt, Pfam, InterPro, PDB or another
  database already holds the same region (within a few residues) with the same meaning, and the
  submission brings no new experiment, reject it. The tools are already graded on those
  databases; a copy adds nothing and counts one result twice. A new experiment on a region a
  database already lists is accepted, and the reason says what is new.
- The region is a prediction (from a structure model, a sequence motif scan, a language model),
  however confident.
- The evidence does not place the region (see rule 1), or the ECO code is not experimental.
- The coordinates cannot be placed on the named sequence.
- One issue holds several regions: ask the submitter to send one issue per region.

## How a decision is recorded

Every decision carries a reason, written for the submitter. The curator

1. comments on the issue: "Accepted: …" or "Rejected: …", with the reason;
2. adds the label `truth-accepted` or `truth-rejected` and closes the issue;
3. records the decision in `data/curation/decisions.tsv`, in a pull request:

   ```bash
   scripts/truth_contributed.py decide --issue 12 --accept BINDING \
       --reason "Deletion mutant in Fig. 3 removes binding; not in UniProt" --curator 0000-0002-1825-0097
   scripts/truth_contributed.py decide --issue 13 --reject \
       --reason "Same region and meaning as UniProt BINDING 830-831; no new experiment" --curator 0000-0002-1825-0097
   ```

An accepted row is filed under one Swiss-Prot feature type (`--accept`): `BINDING`, `ACT_SITE`,
`SITE`, `MOTIF`, `REGION`, `DOMAIN`, `TRANSMEM`, `INTRAMEM`, `DNA_BIND`, `ZN_FING`, `REPEAT` or
`COILED`. A search tool is graded on the row by copying the relative's Swiss-Prot features of
that type, so the type must be the one a database would use for such a region. A guide:
binding of a small molecule or ion at a few residues is `BINDING`; a catalytic residue is
`ACT_SITE`; a binding or activity stretch with no better type is `REGION`; a localisation
signal is `MOTIF`; a solved structured domain is `DOMAIN`.

## From decision to table to rescoring

At each release, `make truth` (`scripts/truth_contributed.py build`) reads the queue and the
decisions, checks every accepted row again, and writes

- `data/truth_contributed.tsv`: one row per accepted region, with `provenance = wet_lab`, the
  submitter's ORCID iD, the curator, the date and the reason;
- `data/truth_contributed/<query set>.parquet`: the rows that fall on a query set's proteins,
  in the scorer's layout;
- `docs/data/contributed_counts.json`: accepted rows per organism, which the organism page shows.

`make rescore` then scores every submission folder again, adding a third truth set,
"Wet-lab rows contributed through the site", wherever a query set has accepted rows. The
leaderboard's own runs of the in-house tools (in seanome/2024-kmerseek-analysis) are scored on
it when that pipeline is next run with `data/truth_contributed/` as an extra truth set; that step
is not automatic yet.

A row on a protein outside every query set is kept in the table and counted for its organism,
but no tool is graded on it until a query set includes that protein.
