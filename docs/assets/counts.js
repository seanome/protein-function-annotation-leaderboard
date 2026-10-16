/* Two small components for the organism page and the Why page. Include with
     <script defer src="<up>assets/counts.js"></script>
   on any page; each part does nothing unless its marker is on the page.

   1. Contributed rows per organism. Put
        <span data-contributed-count data-taxon="9606"></span>
      where the count should appear. It is filled from data/contributed_counts.json, which
      scripts/truth_contributed.py writes from the accepted wet-lab rows. 0 is a real count.

   2. Copy-command clicks. A button
        <button type="button" data-copy-command data-organism="9606" data-toolset="sequence-only"
                data-text="mmseqs search ...">Copy command</button>
      copies data-text to the clipboard. If data/counter.json names an endpoint, the click
      also sends one request to <endpoint>/hit/<organism>/<toolset>, with no cookies, no
      body and no identifier: the counter learns only which organism and tool set were
      copied. With endpoint null (the default) nothing is sent.

   3. The counts on the Why page. Put
        <div data-usage-counts></div>
      where the table should go. With no endpoint it says that copies are not counted yet,
      rather than showing 0, because 0 would be a count nobody measured. */
(function () {
  "use strict";
  const script = document.currentScript;
  const root = script ? script.src.replace(/assets\/counts\.js.*$/, "") : "";
  const getJSON = (path) => fetch(root + path, { credentials: "omit", cache: "no-cache" })
    .then((r) => (r.ok ? r.json() : null)).catch(() => null);

  function contributed() {
    const spots = document.querySelectorAll("[data-contributed-count]");
    if (!spots.length) return;
    getJSON("data/contributed_counts.json").then((c) => {
      for (const el of spots) {
        if (!c) { el.textContent = "count not available"; continue; }
        const byTaxon = c.by_taxon || {};
        const t = el.dataset.taxon;
        const n = t && byTaxon[t] ? byTaxon[t].n_rows : (c.by_name || {})[el.dataset.organism] || 0;
        el.textContent = `${n} contributed ${n === 1 ? "row" : "rows"}`;
      }
    });
  }

  let counter = null;
  const counterReady = getJSON("data/counter.json").then((c) => { counter = c; return c; });

  function copyButtons() {
    document.addEventListener("click", (ev) => {
      const b = ev.target.closest("[data-copy-command]");
      if (!b) return;
      const text = b.dataset.text || "";
      const done = () => { const old = b.textContent; b.textContent = "Copied"; setTimeout(() => { b.textContent = old; }, 1500); };
      (navigator.clipboard ? navigator.clipboard.writeText(text) : Promise.reject()).then(done, () => {});
      counterReady.then((c) => {
        if (!c || !c.endpoint) return;
        const org = encodeURIComponent(b.dataset.organism || "unknown");
        const set = encodeURIComponent(b.dataset.toolset || "unknown");
        fetch(`${c.endpoint}/hit/${org}/${set}`, { method: "POST", credentials: "omit", keepalive: true,
          referrerPolicy: "no-referrer", mode: "no-cors" }).catch(() => {});
      });
    });
  }

  function usageTable() {
    const spots = document.querySelectorAll("[data-usage-counts]");
    if (!spots.length) return;
    counterReady.then((c) => {
      const say = (el, msg) => { const p = document.createElement("p"); p.textContent = msg; el.replaceChildren(p); };
      for (const el of spots) {
        if (!c || !c.endpoint) { say(el, "Copies of the command are not counted yet: no counter service is set up, so there is no number to show."); continue; }
        fetch(`${c.endpoint}/counts`, { credentials: "omit", referrerPolicy: "no-referrer" })
          .then((r) => r.json())
          .then((d) => {
            const rows = [];
            for (const [org, sets] of Object.entries(d.by_organism || {}))
              for (const [set, n] of Object.entries(sets)) rows.push([org, set, n]);
            rows.sort((a, b) => b[2] - a[2]);
            const t = document.createElement("table");
            t.className = "lb";
            t.innerHTML = "<thead><tr><th>Organism</th><th>Tool set</th><th>Copies</th></tr></thead>";
            const tb = document.createElement("tbody");
            for (const [o, s, n] of rows) {
              const tr = document.createElement("tr");
              for (const [x, num] of [[o, false], [s, false], [n, true]]) {
                const td = document.createElement("td");
                td.textContent = String(x);
                if (num) td.className = "num";
                tr.appendChild(td);
              }
              tb.appendChild(tr);
            }
            t.appendChild(tb);
            el.replaceChildren(t);
          })
          .catch(() => say(el, "The counter service did not answer, so no counts are shown."));
      }
    });
  }

  const start = () => { contributed(); copyButtons(); usageTable(); };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start); else start();
})();
