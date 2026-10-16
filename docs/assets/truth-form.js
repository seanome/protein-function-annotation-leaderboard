/* The wet-lab form at /contribute/truth/. Checks the answers, writes them as the issue text
   (the same "### heading" layout as .github/ISSUE_TEMPLATE/truth.yml produces, which
   scripts/truth_contributed.py reads back), and opens GitHub's new-issue page with that text
   filled in. Nothing is sent anywhere until the person presses Submit on GitHub. */
(function () {
  "use strict";
  const form = document.getElementById("truth-form");
  if (!form) return;
  const cfg = JSON.parse(document.getElementById("tf-config").textContent);
  const $ = (id) => document.getElementById(id);
  const AA = /^[ACDEFGHIKLMNPQRSTVWYUOBZXJ*]+$/;
  const ACC = /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$/;

  // ORCID iDs end in a check character (ISO 7064 MOD 11-2), so a typo is caught here.
  function orcidOk(s) {
    if (!/^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$/.test(s)) return false;
    const d = s.replace(/-/g, "");
    let total = 0;
    for (let i = 0; i < 15; i++) total = (total + Number(d[i])) * 2;
    const r = (12 - (total % 11)) % 11;
    return d[15] === (r === 10 ? "X" : String(r));
  }

  function values() {
    const shown = form.querySelector('input[name="shown"]:checked');
    return {
      accession: $("tf-acc").value.trim().toUpperCase(),
      sequence: $("tf-seq").value.replace(/\s+/g, "").toUpperCase(),
      organism: $("tf-org").value.trim(),
      start: $("tf-start").value.trim(),
      end: $("tf-end").value.trim(),
      shown: shown ? shown.value : "",
      eco: $("tf-eco").value.trim().toUpperCase(),
      citation: $("tf-cite").value.trim(),
      orcid: $("tf-orcid").value.trim().toUpperCase(),
      cc0: $("tf-cc0").checked,
    };
  }

  function problems(v) {
    const p = [];
    if (!v.accession && !v.sequence) p.push("Give a UniProt accession or the sequence.");
    if (v.accession && !ACC.test(v.accession)) p.push("The accession does not look like a UniProt accession (such as P01375).");
    if (v.sequence && !AA.test(v.sequence)) p.push("The sequence may hold only one-letter amino-acid codes.");
    if (!v.accession && !v.organism) p.push("Name the organism when you give a sequence.");
    const s = Number(v.start), e = Number(v.end);
    if (!Number.isInteger(s) || s < 1) p.push("The first residue must be a whole number of 1 or more.");
    if (!Number.isInteger(e) || e < s) p.push("The last residue must be a whole number, no smaller than the first.");
    if (v.sequence && Number.isInteger(e) && e > v.sequence.length) p.push(`The region ends at ${e}, past the end of the ${v.sequence.length}-residue sequence.`);
    if (!v.shown) p.push("Pick what the experiment showed.");
    if (!/^ECO:\d{7}$/.test(v.eco)) p.push("The evidence code must look like ECO:0000269.");
    if (!v.citation) p.push("Give a citation, or a note on unpublished work.");
    if (!orcidOk(v.orcid)) p.push("The ORCID iD is not valid; it looks like 0000-0002-1825-0097, and its last character is a check digit.");
    if (!v.cc0) p.push("Tick the box about CC0 and database copies.");
    return p;
  }

  function body(v) {
    const none = "_No response_";
    const sections = [
      ["UniProt accession", v.accession || none],
      ["Sequence", v.sequence || none],
      ["Organism", v.organism || none],
      ["First residue", v.start],
      ["Last residue", v.end],
      ["What was shown", v.shown],
      ["ECO evidence code", v.eco],
      ["Citation or note", v.citation],
      ["ORCID iD", v.orcid],
    ];
    return sections.map(([h, x]) => `### ${h}\n\n${x}`).join("\n\n") +
      "\n\n### Licence and database copies\n\n- [x] I release this row under CC0, and it is not only a copy of a database annotation\n";
  }

  function title(v) {
    const who = v.accession || v.organism;
    return `[truth] ${who} ${v.start}-${v.end} ${v.shown}`;
  }

  function preview() { $("tf-body").textContent = body(values()); }
  form.addEventListener("input", preview);
  preview();

  form.addEventListener("submit", (ev) => {
    ev.preventDefault();
    const v = values();
    const p = problems(v);
    const box = $("tf-problems");
    if (p.length) {
      box.hidden = false;
      box.innerHTML = "";
      const ul = document.createElement("ul");
      for (const x of p) { const li = document.createElement("li"); li.textContent = x; ul.appendChild(li); }
      box.appendChild(document.createTextNode("Please fix before opening the issue:"));
      box.appendChild(ul);
      box.focus && box.focus();
      return;
    }
    box.hidden = true;
    const url = `https://github.com/${cfg.repo}/issues/new?` + new URLSearchParams({
      title: title(v), body: body(v), labels: cfg.label,
    }).toString();
    window.open(url, "_blank", "noopener");
  });
})();
