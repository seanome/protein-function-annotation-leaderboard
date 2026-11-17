"""The two ways in to the leaderboard: /contribute/ and the wet-lab form at /contribute/truth/.

Called from build_site.main with its layout() so the pages share the header, fonts and
footer. The form has no server behind it: it opens GitHub's new-issue page with the answers
filled in, and the person submits the issue there. docs/curation.md says what happens next.
"""

import json
from pathlib import Path

REPO = "seanome/protein-function-annotation-leaderboard"
ISSUE_LABEL = "truth-submission"
SHOWN = [
    ("binding", "Binding", "the region binds a molecule, ion or another protein"),
    ("activity", "Activity", "the region carries out a catalytic or other measured activity"),
    ("localisation", "Localisation", "the region sends the protein to a place in the cell"),
    ("structure", "Structure", "a 3D structure of the region was solved"),
]

# One step per box, read left to right (top to bottom on a phone). Used by Figure 1 on both pages.
FLOW_TOOL = [
    ("1", "Run your tool", "on a published query set, on your own hardware"),
    ("2", "Open a pull request", "predictions, timing log, tool card"),
    ("3", "Workflow scores it", "the leaderboard's own scorer, against the known regions of the query proteins"),
    ("4", "Rows posted", "as a comment on the pull request, before merge"),
    ("5", "Merge", "the folder stays in the repository and is scored again at every release"),
]
FLOW_LAB = [
    ("1", "Fill in the form", "one tested region of one protein"),
    ("2", "GitHub issue", "you submit it on GitHub; it joins the curation queue"),
    ("3", "Curator decides", "accepted or rejected, always with a reason"),
    ("4", "Table of tested regions", "an accepted region is added, with your ORCID iD"),
    ("5", "Tools graded on it", "at the next rebuild of the tables, if a query set holds the protein"),
]


def flow(steps, label) -> str:
    items = "".join(
        f'<li><span class="flow-n">{n}</span><b>{t}</b><span class="flow-d">{d}</span></li>'
        for n, t, d in steps
    )
    return f'<ol class="flow" aria-label="{label}">{items}</ol>'


def page_contribute(layout, src) -> str:
    gh = f"https://github.com/{REPO}"
    body = f"""<h1>Contribute</h1>
<p class="lede">The leaderboard grows in two ways without anyone here running your tool or
your experiment. Tool authors send predictions; labs send regions they have tested.</p>

<h2 id="tools">Tool authors: send predictions</h2>
<figure>
<figcaption><b>Figure 1.</b> What happens to a tool submission, in order. No account on this
site and nothing of ours to install: a GitHub account is enough.</figcaption>
{flow(FLOW_TOOL, "Steps of a tool submission")}
</figure>
<p>A query set is a fixed list of query proteins whose known regions (Pfam domains and Swiss-Prot
features) are the answers every tool is graded on. Run your tool on one of the query sets in
<a href="{gh}/tree/main/data/query_sets"><span class="mono">data/query_sets/</span></a>, then open a
pull request that adds one folder under <span class="mono">submissions/</span> with three files:
the predictions, a timing log (how long each run took, on what hardware), and a short tool card
(what the tool compares, its licence and command, shown on the site). A workflow grades the predictions with
the same scorer that made every number on this site and posts the rows on the pull request.
Cost numbers you send are shown as <em>reported</em>, with your hardware named, never as measured.
<a href="{gh}/blob/main/CONTRIBUTING.md">CONTRIBUTING.md</a> has the five steps and a worked example
with MMseqs2. The leaderboard table does not show submitted rows yet; until it does, they are on
the pull request and in the repository.</p>

<h2 id="labs">Labs: send a region you have tested</h2>
<figure>
<figcaption><b>Figure 2.</b> What happens to a lab's submission (one region of one protein), in order.</figcaption>
{flow(FLOW_LAB, "Steps of a lab submission")}
</figure>
<p>If your lab has shown by experiment that a stretch of a protein binds something, carries an
activity, sends the protein somewhere in the cell, or has a solved structure, that region is a
test every tool can be graded on, once a query set holds the protein. <a href="truth/">The form</a> takes one region at a time. A
curator checks it against the rules in
<a href="{gh}/blob/main/docs/curation.md">docs/curation.md</a>; a row that only repeats an
annotation already in a database is turned down, because the tools are already graded on those.</p>
"""
    return layout("Contribute", "contribute", body, src, depth=1,
                  description="How tool authors and labs add to the protein annotation leaderboard.")


def page_truth_form(layout, src) -> str:
    shown = "".join(
        f'<label class="tf-radio"><input type="radio" name="shown" value="{v}" required>'
        f'<span><b>{t}</b>: {d}</span></label>' for v, t, d in SHOWN
    )
    config = json.dumps({"repo": REPO, "label": ISSUE_LABEL})
    body = f"""<nav class="crumb" aria-label="Breadcrumb"><a href="../">Contribute</a> › Add a tested region</nav>
<h1>Add a region your lab has tested</h1>
<p class="lede">One protein region per form, with what your experiment showed about it. It goes to
a public curation queue. Once accepted, every tool on the leaderboard is graded on it at the next
rebuild of the tables, if the protein is in one of the query sets (the fixed lists of proteins the
tools are run on); otherwise it is kept and counted for its organism until a query set includes it.</p>

<figure>
<figcaption><b>Figure 1.</b> Where your answers go, in order. Nothing leaves this page until you
press <em>Submit new issue</em> on GitHub in step 2.</figcaption>
{flow(FLOW_LAB, "Steps of a lab submission")}
</figure>

<h2 id="how">How the form reaches GitHub</h2>
<p>This site is a set of static pages with no server, so it cannot store anything. The button at
the end opens GitHub's new-issue page for
<a href="https://github.com/{REPO}/issues?q=label%3A{ISSUE_LABEL}"><span class="mono">{REPO}</span></a>
with your answers already written in. You read them there and press <em>Submit new issue</em>.
That needs a free GitHub account, and the issue is public. The curator answers in the same issue.</p>

<form id="truth-form" class="tf" novalidate>
<fieldset><legend>Protein</legend>
<p class="tf-hint">Give the UniProt accession, or, only if there is none, the sequence and the organism.</p>
<label for="tf-acc">UniProt accession <span class="tf-hint">such as P01375</span></label>
<input id="tf-acc" name="accession" autocomplete="off" spellcheck="false" inputmode="latin">
<label for="tf-seq">Or the sequence, if the protein has no UniProt entry <span class="tf-hint">one-letter amino-acid code</span></label>
<textarea id="tf-seq" name="sequence" rows="4" spellcheck="false"></textarea>
<label for="tf-org">Organism <span class="tf-hint">scientific name or NCBI taxon ID; needed when you give a sequence</span></label>
<input id="tf-org" name="organism" autocomplete="off">
</fieldset>

<fieldset><legend>Region</legend>
<p class="tf-hint">Required. Residue positions on that sequence, counting the first residue as 1. Both ends are part of the region. If your experiment used another isoform or an older sequence, convert the positions to the UniProt sequence first.</p>
<div class="tf-pair"><span><label for="tf-start">First residue</label>
<input id="tf-start" name="start" type="number" min="1" step="1" required></span>
<span><label for="tf-end">Last residue</label>
<input id="tf-end" name="end" type="number" min="1" step="1" required></span></div>
</fieldset>

<fieldset><legend>What the experiment showed about this region (required)</legend>
{shown}
<label for="tf-eco">Evidence code from the ECO ontology <span class="tf-hint">such as ECO:0000269, "experimental evidence used in manual assertion". It must be an experimental code: codes for sequence similarity, prediction or author statement are turned down. Look one up at the
<a href="https://www.ebi.ac.uk/ols4/ontologies/eco">ECO browser</a></span></label>
<input id="tf-eco" name="eco" autocomplete="off" spellcheck="false" placeholder="ECO:0000000" required>
<label for="tf-cite">Citation or note <span class="tf-hint">a DOI or PubMed ID of the paper, or a note on unpublished work: method, controls, who did it. Say what the region binds, what activity was measured, or where in the cell it sends the protein</span></label>
<textarea id="tf-cite" name="citation" rows="4" required></textarea>
</fieldset>

<fieldset><legend>You</legend>
<label for="tf-orcid">Your ORCID iD <span class="tf-hint">such as 0000-0002-1825-0097; it goes into the table beside your row</span></label>
<input id="tf-orcid" name="orcid" autocomplete="off" spellcheck="false" placeholder="0000-0000-0000-0000" required>
<label class="tf-radio"><input type="checkbox" name="cc0" id="tf-cc0" required><span>I release this region under
<a href="https://creativecommons.org/publicdomain/zero/1.0/">CC0</a>, so anyone may reuse it.</span></label>
<label class="tf-radio"><input type="checkbox" name="newexp" id="tf-new" required><span>It comes from an
experiment, not only from an annotation already in UniProt, Pfam, InterPro or another database.</span></label>
</fieldset>

<p id="tf-problems" class="tf-problems" role="alert" hidden></p>
<p><button type="submit" class="tf-go">Open the filled-in issue on GitHub</button></p>
<details class="tf-preview"><summary>The issue text the button will fill in</summary>
<pre id="tf-body"></pre></details>
</form>
<script type="application/json" id="tf-config">{config}</script>
"""
    head = '<script defer src="../../assets/truth-form.js"></script>'
    return layout("Add a tested region", "contribute", body, src, depth=2, head=head,
                  description="Send one protein region your lab has tested; it becomes a test for every tool.")


def write(out: Path, src: dict, layout) -> None:
    (out / "contribute" / "truth").mkdir(parents=True, exist_ok=True)
    (out / "contribute" / "index.html").write_text(page_contribute(layout, src))
    (out / "contribute" / "truth" / "index.html").write_text(page_truth_form(layout, src))
