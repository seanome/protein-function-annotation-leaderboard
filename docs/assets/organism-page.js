/* Organism pages in the browser: the "I work on" box with its suggestions, the page built here
   for a species with no pre-built page, and the Copy command buttons. The answer itself comes
   from organism.js, the same code the build runs. */
(function () {
  "use strict";
  const host = document.querySelector("[data-organism]");
  const page = host ? host.dataset.organism : null;
  const base = host ? host.dataset.base || "" : "";

  // Copy command buttons.
  document.addEventListener("click", async (e) => {
    const b = e.target.closest("button[data-copy]");
    if (!b) return;
    const el = document.getElementById(b.dataset.copy);
    try {
      await navigator.clipboard.writeText(el.textContent);
      b.textContent = "Copied";
    } catch (err) {
      const r = document.createRange(); r.selectNodeContents(el);
      const s = getSelection(); s.removeAllRanges(); s.addRange(r);
      b.textContent = "Selected: press Ctrl+C or ⌘C";
    }
    setTimeout(() => { b.textContent = "Copy command"; }, 2500);
  });

  // ---------------------------------------------------------------- the "I work on" box
  if (page === "search") {
    // A link made for the old home page (the leaderboard, with its filters after #) goes there.
    if (/(^|[#&])(truth|species|axis|calls)=/.test(location.hash)) location.replace("leaderboard.html" + location.hash);
    const input = document.getElementById("species");
    const list = document.getElementById("suggest");
    const status = document.getElementById("suggest-status");
    const shards = {};
    let active = -1, items = [], seq = 0;
    const key = (s) => { const w = s.toLowerCase().replace(/[^a-z]/g, "").slice(0, 2); return w.length === 2 ? w : "_"; };
    async function shard(k) {
      if (!(k in shards)) {
        shards[k] = fetch(`${base}data/organism/names/${k}.json`).then((r) => (r.ok ? r.json() : [])).catch(() => []);
      }
      return shards[k];
    }
    function hrefFor(it) { return it[2] ? `${base}organism/${it[2]}.html` : `${base}organism/live.html#taxid=${it[1]}`; }
    function draw(text) {
      list.innerHTML = "";
      items.forEach((it, i) => {
        const li = document.createElement("li");
        li.id = `opt-${i}`; li.setAttribute("role", "option");
        li.setAttribute("aria-selected", i === active ? "true" : "false");
        const a = document.createElement("a");
        a.href = it.href || hrefFor(it);
        a.textContent = it[0];
        const tag = document.createElement("span");
        tag.className = "tag";
        tag.textContent = it.ncbi ? "from NCBI, built when opened" : it[2] ? "page ready" : "built when opened";
        a.appendChild(tag);
        li.appendChild(a);
        list.appendChild(li);
      });
      list.hidden = items.length === 0;
      input.setAttribute("aria-expanded", items.length ? "true" : "false");
      status.textContent = text || "";
    }
    async function ncbi(text) {
      // NCBI's own taxonomy search, for a species missing from the site's index.
      const t = encodeURIComponent(`${text}*[All Names] AND species[Rank]`);
      const s = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmode=json&retmax=8&term=${t}`).then((r) => r.json());
      const ids = s.esearchresult.idlist;
      if (!ids.length) return [];
      const m = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&retmode=json&id=${ids.join(",")}`).then((r) => r.json());
      return ids.map((id) => {
        const r = m.result[id];
        const label = r.commonname ? `${r.commonname} (${r.scientificname})` : r.scientificname;
        const it = [label, Number(id), ""]; it.ncbi = true; return it;
      });
    }
    async function update() {
      const text = input.value.trim();
      const my = ++seq;
      if (text.length < 2) { items = []; draw(""); return; }
      const low = text.toLowerCase();
      const all = await shard(key(text.split(/\s+/)[0]));
      if (my !== seq) return;
      const seen = new Set();
      items = [];
      for (const it of all) {
        const l = it[0].toLowerCase();
        if ((l.startsWith(low) || l.split(/[\s(]+/).some((w) => w.startsWith(low))) && !seen.has(it[1])) {
          seen.add(it[1]); items.push(it);
          if (items.length >= 10) break;
        }
      }
      active = -1;
      draw(items.length ? `${items.length} suggestions` : "Searching NCBI Taxonomy…");
      if (items.length < 3 && text.length >= 3) {
        try {
          const more = (await ncbi(text)).filter((it) => !seen.has(it[1]));
          if (my !== seq) return;
          items = items.concat(more).slice(0, 10);
          draw(items.length ? `${items.length} suggestions` : "No species found");
        } catch (e) { if (!items.length) draw("No species found in the site's index, and NCBI did not answer"); }
      }
    }
    input.addEventListener("input", update);
    input.addEventListener("keydown", (e) => {
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault();
        if (!items.length) return;
        active = (active + (e.key === "ArrowDown" ? 1 : items.length - 1)) % items.length;
        draw(document.getElementById("suggest-status").textContent);
        input.setAttribute("aria-activedescendant", `opt-${active}`);
      } else if (e.key === "Enter") {
        e.preventDefault();
        const it = items[active >= 0 ? active : 0];
        if (it) location.href = it.href || hrefFor(it);
      } else if (e.key === "Escape") { items = []; draw(""); }
    });
  }

  // ---------------------------------------------------------------- a page built here
  if (page === "live") {
    const out = document.getElementById("answer");
    const prog = document.getElementById("progress");
    async function run() {
      const m = location.hash.match(/taxid=(\d+)/);
      if (!m) { out.innerHTML = `<p>No species chosen. <a href="${base}index.html">Choose one</a>.</p>`; return; }
      const taxid = Number(m[1]);
      try {
        const [S, reviewed, mrca, prebuilt] = await Promise.all(["strata", "reviewed", "mrca", "prebuilt"].map((n) =>
          fetch(`${base}data/organism/${n}.json`).then((r) => r.json())));
        // A pre-built page exists: go there.
        // (#taxid=…&live=1 builds it here anyway; tests/check_organism.py compares the two.)
        if (prebuilt[String(taxid)] && !/live=1/.test(location.hash)) { location.replace(`${base}organism/${prebuilt[String(taxid)]}.html`); return; }
        const L = await Organism.lookupSpecies(taxid, { reviewed, mrca }, (what) => { prog.textContent = `Looking up: ${what}…`; });
        const A = Organism.computeAnswer(L, S);
        A.built = "in this browser, " + new Date().toISOString().slice(0, 10);
        document.title = `${L.common_name ? L.common_name + " (" + L.name + ")" : L.name}: which tools to run`;
        out.innerHTML = Organism.renderAnswer(A, S);
        prog.textContent = "";
        const blob = new Blob([JSON.stringify({ page: "organism/live.html", built_in_browser: true, lookups: L, answer: A }, null, 1)], { type: "application/json" });
        const a = document.getElementById("json-link");
        a.href = URL.createObjectURL(blob); a.download = `${L.slug}.json`; a.hidden = false;
        window.__organismAnswer = { L, A };
      } catch (e) {
        prog.textContent = "";
        out.innerHTML = `<p class="notice">This page could not be built: ${String(e.message || e).replace(/</g, "&lt;")}. UniProt, AlphaFold DB or GoaT did not answer; try again later.</p>`;
        window.__organismAnswer = { error: String(e) };
      }
    }
    window.addEventListener("hashchange", run);
    run();
  }
})();
