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
  function render(id, r, d) {
    var sub = fmt(r[1]) + " aa";
    document.getElementById("sub").textContent = sub;
    var h = [];
    // Names
    h.push("<h2>Name</h2>" + table([["Existing name (FASTA header)"], ["Proposed name"], ["Comparison"]],
      [[esc(r[2]) || "none", esc(r[3]) || "none", esc(F.check[r[4]])]]));
    if (r[14]) h.push("<p><a href=\"" + esc(r[14]) + "\">Every annotation program's hits drawn along this protein</a> " +
      "(a page made from separate searches, matched to this protein by its sequence).</p>");
    // Best Swiss-Prot match
    h.push("<h2>Best Swiss-Prot match (DIAMOND)</h2>");
    if (r[6]) {
      h.push(table([["Accession"], ["Protein"], ["Identity (%)", 1], ["Coverage of your protein (%)", 1], ["E-value", 1]],
        [["<a class=\"mono\" href=\"https://www.uniprot.org/uniprotkb/" + encodeURIComponent(r[6]) + "\">" + esc(r[6]) + "</a>",
          esc(r[5]), r[7] === null ? "" : r[7].toFixed(1), r[8] === null ? "" : Math.round(r[8]), esc(d.sp_e)]]));
      h.push("<p class=\"note\">Where on your protein the match lies is not in the table this page was built from.</p>");
      if (d.go.length) h.push("<h3>GO terms of the matched Swiss-Prot entry (not measured on your protein)</h3><p class=\"chips\">" + d.go.map(function (g) {
        return "<a class=\"mono\" href=\"https://www.ebi.ac.uk/QuickGO/term/" + encodeURIComponent(g) + "\">" + esc(g) + "</a>"; }).join(" ") + "</p>");
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
