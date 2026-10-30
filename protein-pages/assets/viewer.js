// protein.html: one page for any protein. The ID after "#" picks the protein; its data are in
// full/NNN.js, where NNN is the same hash of the ID that build_full.py's shard_of computes.
(function () {
  var F = window.FULL, CATS = F.cats, N_SHARDS = 96;
  var byId = {};
  F.rows.forEach(function (r) { byId[r[0]] = r; });
  function shardOf(id) {
    var h = 0;
    for (var i = 0; i < id.length; i++) h = (Math.imul(h, 31) + id.charCodeAt(i)) >>> 0;
    return h % N_SHARDS;
  }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function fmt(n) { return n.toLocaleString("en-US"); }
  var main = document.getElementById("main"), loaded = {};
  function load(s, then) {
    if (loaded[s]) return then();
    var el = document.createElement("script");
    el.src = "full/" + ("00" + s).slice(-3) + ".js";
    el.onload = function () { loaded[s] = true; then(); };
    el.onerror = function () { main.innerHTML = "<p>Could not load " + esc(el.src) + ".</p>"; };
    document.head.appendChild(el);
  }
  function table(head, rows) {
    return "<div class=\"tablewrap\"><table class=\"small\"><thead><tr>" + head.map(function (h) {
      return "<th" + (h[1] ? " class=\"r\"" : "") + ">" + esc(h[0]) + "</th>";
    }).join("") + "</tr></thead><tbody>" + rows.map(function (r) {
      return "<tr>" + r.map(function (c, i) { return "<td" + (head[i][1] ? " class=\"r\"" : "") + ">" + c + "</td>"; }).join("") + "</tr>";
    }).join("") + "</tbody></table></div>";
  }
  function link(id) { return "<a class=\"mono\" href=\"protein.html#" + encodeURIComponent(id) + "\">" + esc(id) + "</a>"; }
  function show() {
    var id = location.hash.slice(1), r;
    try { id = decodeURIComponent(id); } catch (e) { /* malformed escape: keep the raw text */ }
    r = byId[id];
    document.getElementById("pid").textContent = id || "No protein named";
    if (!r) { main.innerHTML = "<p>No protein " + esc(id) + " in this proteome. <a href=\"index.html\">Search all proteins</a>.</p>"; return; }
    document.title = id + (r[3] || r[2] ? " " + (r[3] || r[2]) : "") + " · Botryllus protein";
    load(shardOf(id), function () { render(id, r, window.FULL_SHARD[id]); window.scrollTo(0, 0); });
  }
  function uni(acc) { return "<a class=\"mono\" href=\"https://www.uniprot.org/uniprotkb/" + encodeURIComponent(acc) + "\">" + esc(acc) + "</a>"; }
  // The top of every protein page, in the same order as a drawn page: name line, one sentence on
  // which programs report a hit, then one row per program with its best hit.
  function top(id, r, d) {
    var cats = []; for (var i = 0; i < CATS.length; i++) if (r[11] >> i & 1) cats.push(CATS[i]);
    var progs = [], rows = [];
    function row(name, found, hit, score, where) {
      if (found) progs.push(name);
      rows.push(["<b>" + esc(name) + "</b>", found ? hit : "<i>no hit</i>", found ? score : "", found ? where : ""]);
    }
    var nowhere = "not in the table";
    row("DIAMOND against Swiss-Prot", !!r[6], r[6] ? uni(r[6]) + " " + esc(r[5]) + (d.sp_gene ? " (" + esc(d.sp_gene) + ")" : "") + (d.sp_org ? ", <i>" + esc(d.sp_org) + "</i>" : "") : "",
      r[6] ? (r[7] === null ? "" : r[7].toFixed(1) + "% identity over " + Math.round(r[8]) + "% of your protein, ") + "E-value " + esc(d.sp_e) : "", nowhere);
    var pf = {}; d.pfam.forEach(function (p) { pf[p[1]] = (pf[p[1]] || 0) + 1; });
    row("hmmsearch with Pfam", d.pfam.length > 0, Object.keys(pf).map(function (n) { return esc(n) + (pf[n] > 1 ? " ×" + pf[n] : ""); }).join(", "),
      "gathering threshold", nowhere);
    var j = d.jh[0];
    row("jackhmmer from immune seed proteins", d.jh.length > 0, j ? esc(j[0]) + " (seed " + esc(j[1]) + ")" + (d.jh.length > 1 ? " and " + (d.jh.length - 1) + " more" : "") : "",
      j ? "E-value " + esc(j[3]) + ", round " + j[2] : "", nowhere);
    var a = d.al[0];
    row("DIAMOND against Fester and BHF", d.al.length > 0, a ? esc(a[1] || a[0]) + " <a class=\"mono\" href=\"https://www.ncbi.nlm.nih.gov/protein/" + encodeURIComponent(a[0]) + "\">" + esc(a[0]) + "</a>" : "",
      a ? a[2].toFixed(1) + "% identity over " + fmt(a[3]) + " columns, E-value " + esc(a[4]) : "", nowhere);
    var f = d.fa;
    row("FANTASIA (nearest protein by embedding)", !!f, f ? uni(f[0]) + " " + esc(f[1] || "") + (f[2] ? " (" + esc(f[2]) + ")" : "") + (f[3] ? ", <i>" + esc(f[3]) + "</i>" : "") : "",
      f ? "cosine distance " + f[4].toFixed(2) + ", reliability " + f[5].toFixed(2) : "", f ? "whole protein, no alignment" : "");
    var name = (d.card && d.card.card_name) || r[3] || r[2] || "";
    var line = d.card ? d.card.one_line : (cats.length ? "immune-gene candidate: " + cats.join(", ") : "");
    document.getElementById("sub").textContent = (name ? name + " · " : "") + (line ? line + " · " : "") + fmt(r[1]) + " amino acids";
    var n = progs.length, sentence = n === 0 ? "None of the 5 programs reports a hit." :
      n + " of the 5 programs report" + (n === 1 ? "s" : "") + " a hit: " + progs.join(", ") + ".";
    var h = ["<p class=\"lead\">" + esc(sentence) + "</p>", "<h2>What each program found</h2>",
      "<p class=\"note\">The best hit of each program. The Biomni tables do not say where on the protein a hit lies, " +
      "so nothing is drawn along it" + (r[14] ? "; the drawing below comes from separate searches." : ".") + "</p>",
      table([["Program"], ["Best hit"], ["Score"], ["Where on your protein"]], rows)];
    if (r[14]) h.push("<h2>Every program's hits drawn along this protein</h2><p class=\"note\">From separate searches against Swiss-Prot, " +
      "matched to this protein by its sequence. <a href=\"" + esc(r[14]) + "\">Open the drawing on its own page</a>.</p>" +
      "<iframe class=\"drawn\" title=\"Hits drawn along " + esc(id) + "\" src=\"" + esc(r[14]) + "\" loading=\"lazy\"></iframe>");
    if (d.card) {
      var c = d.card, v = c.name_verdict;
      h.push("<h2>Gene card</h2><p class=\"note\">Written by a language model (" + esc(c.model) + ") from the evidence on this page only; " +
        "check it against the tables below.</p><p>" + esc(c.summary) + "</p><p><b>Confidence: " + esc(c.confidence) + ".</b> " + esc(c.confidence_reason) + "</p>");
      if (v && v.choice !== "not_applicable") h.push("<p><b>Header name " + esc(r[2]) + " or proposed name " + esc(r[3]) + "?</b> " +
        esc(F.verdicts[v.choice] || v.choice) + "; suggested name <b>" + esc(v.recommended_name) + "</b>. " + esc(v.reason) + "</p>");
    }
    return h;
  }
  function render(id, r, d) {
    var h = top(id, r, d);
    // Names
    h.push("<h2>Name</h2>" + table([["Existing name (FASTA header)"], ["Proposed name"], ["Comparison"]],
      [[esc(r[2]) || "none", esc(r[3]) || "none", esc(F.check[r[4]])]]));
    // Best Swiss-Prot match
    h.push("<h2>Best Swiss-Prot match (DIAMOND)</h2>");
    if (r[6]) {
      h.push(table([["Accession"], ["Protein"], ["Identity (%)", 1], ["Coverage of your protein (%)", 1], ["E-value", 1]],
        [["<a class=\"mono\" href=\"https://www.uniprot.org/uniprotkb/" + encodeURIComponent(r[6]) + "\">" + esc(r[6]) + "</a>",
          esc(r[5]), r[7] === null ? "" : r[7].toFixed(1), r[8] === null ? "" : Math.round(r[8]), esc(d.sp_e)]]));
      h.push("<p class=\"note\">Where on your protein the match lies is not in the table this page was built from.</p>");
      if (d.go.length) h.push("<h3>GO terms of the matched Swiss-Prot entry (not measured on your protein)</h3><p class=\"chips\">" + d.go.map(function (g) {
        return "<a href=\"https://www.ebi.ac.uk/QuickGO/term/" + encodeURIComponent(g[0]) + "\">" + esc(g[1] || g[0]) + "</a>"; }).join(" · ") + "</p>");
      if (d.kw.length) h.push("<h3>UniProt keywords of that entry</h3><p>" + d.kw.map(esc).join(", ") + "</p>");
    } else h.push("<p>No Swiss-Prot match at E-value 0.001 or below.</p>");
    // Pfam
    h.push("<h2>Pfam domains (hmmsearch, gathering thresholds)</h2>");
    if (d.pfam.length) {
      var seen = {}, rows = [];
      d.pfam.forEach(function (p) { seen[p[0]] = seen[p[0]] || { n: p[1], c: 0 }; seen[p[0]].c++; });
      Object.keys(seen).forEach(function (a) {
        rows.push(["<a class=\"mono\" href=\"https://www.ebi.ac.uk/interpro/entry/pfam/" + encodeURIComponent(a.split(".")[0]) + "\">" + esc(a) + "</a>",
          esc(seen[a].n), seen[a].c]);
      });
      h.push(table([["Accession"], ["Name"], ["Times listed", 1]], rows));
      h.push("<p class=\"note\">Listed in the table's order; where each domain lies is not in the table.</p>");
    } else h.push("<p>No Pfam domain.</p>");
    // Immune
    var cats = []; for (var i = 0; i < CATS.length; i++) if (r[11] >> i & 1) cats.push(CATS[i]);
    h.push("<h2>Immune-gene candidate</h2>");
    if (cats.length) h.push(table([["Category"], ["Evidence"]], d.ev.map(function (e) { return [esc(e[0]), esc(e[1])]; })));
    else h.push("<p>Not in the immune-gene table.</p>");
    // jackhmmer
    h.push("<h2>jackhmmer hits from immune seed proteins</h2>");
    if (d.jh.length) h.push(table([["Seed family"], ["Seed protein"], ["Round", 1], ["E-value", 1], ["Score (bits)", 1]],
      d.jh.map(function (j) { return [esc(j[0]), esc(j[1]), j[2], esc(j[3]), esc(j[4])]; })));
    else h.push("<p>None.</p>");
    // Allorecognition
    h.push("<h2>Matches to Fester and BHF (DIAMOND)</h2>");
    if (d.al.length) h.push(table([["Reference protein"], ["NCBI accession"], ["Identity (%)", 1], ["Aligned columns", 1],
      ["E-value", 1], ["Bit score", 1], ["Reference length (aa)", 1]],
      d.al.map(function (a) {
        return [esc(a[1] || "name not found at NCBI"), "<a class=\"mono\" href=\"https://www.ncbi.nlm.nih.gov/protein/" + encodeURIComponent(a[0]) + "\">" + esc(a[0]) + "</a>",
          a[2].toFixed(1), fmt(a[3]), esc(a[4]), a[5].toFixed(0), fmt(a[7])];
      })));
    else h.push("<p>None.</p>");
    // Cluster
    h.push("<h2>MMseqs2 cluster</h2>");
    if (r[13] > 1) {
      h.push("<p>" + fmt(r[13]) + " proteins, including this one; its representative, the protein MMseqs2 keeps to stand for the cluster, is " + link(d.rep) + ".</p>");
      h.push(table([["Other member"], ["Name"]], d.mem.map(function (m) { return [link(m[0]), esc(m[1])]; })));
      if (d.more) h.push("<p class=\"note\">and " + fmt(d.more) + " more.</p>");
    } else h.push("<p>Alone in its cluster.</p>");
    main.innerHTML = h.join("");
  }
  window.addEventListener("hashchange", show);
  show();
})();
